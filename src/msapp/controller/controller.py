from matplotlib.figure import Figure
from typing import Any

from msapp.model.msa import MultiSeqAlignment
from msapp.view.msapp_gui import App
from msapp.view.visualization import (color_clusters, create_resized_mat_visualization, get_empty_plot, highlight_row,
                                      show_as_subimages, visualize_dendrogram, visualize_domains)


class Controller:
    """Manages the model (data) and view (GUI). It passes input and state changes between the two."""

    def __init__(self) -> None:
        self.msa = None
        self.gui = None

        self.split_mat_visualization = True
        self.hide_empty_cols = False
        self.reorder_rows = False
        self.colour_clusters = False
        self.highlight_selected_seq = False
        self.selected_seq = -1
        self.dendro_hcutoff = 0.75
        self.calc_domain_mode = ""

        # dirty flags for visualizations that need to be refreshed
        self.msa_changed = True
        self.dendro_changed = True
        self.domains_changed = True

    def start(self) -> None:
        """Starts the main loop of the GUI."""
        if self.gui is None:
            self.gui = App(controller=self)
            self.gui.mainloop()

    # --- configure/trigger state changes

    def initialize_from_file(self, filename: str) -> bool:
        """Creates and initializes a MultiSeqAlignment object from a given file name."""
        try:
            msa = MultiSeqAlignment(filename)
        except Exception as e:
            err_msg = str(e)
            if len(err_msg) > 0:
                self.gui.add_to_textbox(f"ERR: {err_msg}")
            return False

        fname_truncated = filename.split('/')[-1]
        if msa.initialized:
            self.msa = msa
            self.gui.add_to_textbox(f"Sucessfully loaded MSA from file '{fname_truncated}'.\n")

            # flag all plots as dirty
            self.msa_changed = True
            self.dendro_changed = True
            self.domains_changed = True
        else:
            self.gui.add_to_textbox(f"Could not load MSA from file '{fname_truncated}'.\n")
            return False
        return True

    def run_filtering_pipeline(self, filter_type: str = "standard") -> None:
        """Triggers running a matrix filtering pipeline on the MultiSeqAlignment object."""
        if not self.is_mat_initialized(): return
        filter_type = filter_type.lower()
        self.msa.run_filtering_pipeline(filter_type=filter_type)
        self.gui.add_to_textbox(f"-- OP: Running {filter_type} filtering pipeline.")

        self.msa_changed = True
        self.dendro_changed = True
        self.domains_changed = True

    def set_dendro_hcutoff(self, dendro_hcutoff: float) -> None:
        if self.dendro_hcutoff == dendro_hcutoff:
            return
        self.dendro_hcutoff = dendro_hcutoff
        self.dendro_changed = True
        self.domains_changed = True
        if self.colour_clusters:
            self.msa_changed = True

    # --- sequence selection

    def on_seq_selection(self, seq_id: str) -> None:
        seq_indexer = self.get_seq_indexer()
        if seq_indexer is None or seq_id == self.selected_seq:
            return
        info = seq_indexer.get_infos_for_seq_id(seq_id)
        if info is not None:
            self.selected_seq = seq_id
            self.gui.set_selected_seq_info(f"{seq_id}\n{info[0]}")
            idx = info[2] if self.reorder_rows else info[1]
            self.gui.set_slider_value(idx)
            if self.highlight_selected_seq:
                self.msa_changed = True

    def reset_selected_seq(self) -> None:
        self.highlight_selected_seq = False
        self.selected_seq = -1
        self.gui.set_slider_value(0)

    def select_ith_seqid(self, i) -> None:
        if self.reorder_rows:
            seq_id = self.msa.get_seq_indexer().get_ith_seqid_dendro(i)
        else:
            seq_id = self.msa.get_seq_indexer().get_ith_seqid_mat(i)
        self.on_seq_selection(seq_id)

    def get_selseq_idx(self) -> int:
        """Retrieves the index of the currently selected sequence in the MSA mat."""
        idx = -1
        if self.selected_seq != -1 and self.highlight_selected_seq:
            seq_indexer = self.msa.get_seq_indexer()
            if self.reorder_rows:
                idx = seq_indexer.get_dendroidx_from_seqid(self.selected_seq)
            else:
                idx = seq_indexer.get_matidx_from_seqid(self.selected_seq)
        return idx

    # --- toggle mat display options

    def toggle_split_mat_visualization(self, should_split: bool = False) -> None:
        self.split_mat_visualization = should_split
        self.msa_changed = True

    def toggle_hide_empty_cols(self, should_hide: bool = False) -> None:
        self.hide_empty_cols = should_hide
        self.msa_changed = True

    def toggle_reorder_mat_rows(self, should_reorder: bool = False) -> None:
        self.reorder_rows = should_reorder
        self.msa_changed = True
        if self.highlight_selected_seq and self.selected_seq != -1:
            self.gui.set_slider_value(self.get_selseq_idx())

    def toggle_colour_clusters(self, should_colour: bool = False) -> None:
        self.colour_clusters = should_colour
        self.msa_changed = True

    def toggle_highlight_selected_seq(self, should_highlight: bool = False) -> None:
        self.highlight_selected_seq = should_highlight
        if self.selected_seq != -1:
            self.msa_changed = True

    # --- display graphs

    def get_msa_figure(self, split_mat=True) -> Figure:
        mat = self.msa.get_mat(self.hide_empty_cols, self.reorder_rows)
        if mat is None:
            return None
        img = mat
        if self.colour_clusters:
            cluster_labels = self.msa.get_cluster_labels(self.dendro_hcutoff, self.reorder_rows)
            img = color_clusters(mat, cluster_labels)
        img = highlight_row(img, self.get_selseq_idx())
        if split_mat and self.split_mat_visualization:
            wh_ratio = self.gui.get_mat_frame_wh_ratio()
            img = create_resized_mat_visualization(img, wh_ratio, not self.colour_clusters)
        fig = show_as_subimages(img, "")
        return fig

    def on_show_msa_mat(self, force: bool = False) -> None:
        """Controller is ordered to fetch an MSA matrix visualization figure and forward it to the GUI."""
        if not self.is_mat_initialized(): return
        if not force and not self.msa_changed: return

        fig = self.get_msa_figure()
        if fig is not None:
            self.gui.show_matrix(fig)
            self.msa_changed = False

    def on_show_dendrogram(self, force: bool = False) -> None:
        """Controller is ordered to fetch a dendrogram visualization figure and forward it to the GUI."""
        if not self.is_mat_initialized(): return
        if not force and not self.dendro_changed: return

        fig = visualize_dendrogram(self.msa.get_linkage_mat(), self.dendro_hcutoff)
        self.gui.show_dendrogram(fig)
        self.dendro_changed = False

        cluster_labels = self.msa.get_cluster_labels(self.dendro_hcutoff)
        nclusters = len(set(cluster_labels))
        self.gui.set_cluster_count(nclusters)

    def clear_domains_view(self) -> None:
        empty_plot = get_empty_plot()
        self.gui.show_domains(empty_plot)

    def on_show_domains(self, calc_domain_mode: str) -> None:
        if not self.is_mat_initialized(): return
        if not self.domains_changed and self.calc_domain_mode == calc_domain_mode: return

        self.gui.add_to_textbox(f"Calculating domains with mode '{calc_domain_mode}'.")
        domains = self.msa.calculate_domains(self.dendro_hcutoff, calc_domain_mode)
        fig = visualize_domains(domains)
        self.gui.show_domains(fig)
        self.domains_changed = False
        self.calc_domain_mode = calc_domain_mode

    # --- query state

    def get_seq_indexer(self) -> Any:
        return self.msa.get_seq_indexer() if self.is_mat_initialized() else None

    def is_mat_initialized(self) -> bool:
        return self.msa is not None and self.msa.initialized

    # --- save

    def on_save(self, tight=True) -> None:
        msa_fig = self.get_msa_figure(split_mat=False)
        if tight:
            msa_fig.get_axes()[0].set_ylabel("")
            msa_fig.get_axes()[0].set_xlabel("")
        fname = "msa"
        if self.hide_empty_cols:
            fname += "_gapless"
        if self.reorder_rows:
            fname += "_reordered"
        if self.colour_clusters:
            fname += "_colored"
        if tight:
            msa_fig.savefig(f'out/{fname}.png', bbox_inches='tight', pad_inches=0, dpi='figure')
        else:
            msa_fig.savefig(f'out/{fname}.png', bbox_inches='tight', dpi='figure')
