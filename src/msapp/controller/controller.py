from msapp.model.msa import MultiSeqAlignment
from msapp.view.msapp_gui import App
from msapp.model.mat_manipulation import cross_convolve, gaussian_blur

from msapp.view.visualization import (color_clusters, create_resized_mat_visualization, highlight_row,
                                      show_as_subimages, visualize_dendrogram, visualize_domains)


class Controller:
    """Manages the model (data) and view (GUI). It passes input and state changes between the two."""

    def __init__(self):
        self.msa = None
        self.gui = None

        self.split_mat_visualization = True
        self.hide_empty_cols = False
        self.reorder_rows = False
        self.colour_clusters = False
        self.highlight_selected_seq = False
        self.selected_seq = -1
        self.dendro_hcutoff = 0.75

        # dirty flags for visualizations that need to be refreshed
        self.msa_changed = True
        self.dendro_changed = True
        self.domains_changed = True

    def start(self):
        """Starts the main loop of the GUI."""
        if self.gui is None:
            self.gui = App(controller=self)
            self.gui.mainloop()


    # --- configure/trigger state changes

    def initialize_from_file(self, filename: str):
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

    def run_filtering_pipeline(self, filter_type: str = ""):
        """Triggers running a matrix filtering pipeline on the MultiSeqAlignment object."""
        if not self.is_mat_initialized(): return

        self.gui.add_to_textbox("Running filtering pipeline.")
        self.msa.filter_by_length_statistic()
        self.gui.add_to_textbox("-- OP: Filtering by length statistic > 3 sigma.")
        self.msa.remove_isolated_connections()
        self.gui.add_to_textbox("-- OP: Removing isolated connections (likely errors).")

        if filter_type.lower() == "standard":
            self.msa.img_process(gaussian_blur, ksize=5)
            self.gui.add_to_textbox("-- OP: Cross-convolving (1x5)")
            self.msa.img_process(cross_convolve, col_size=5)
            self.gui.add_to_textbox("-- OP: Cross-convolving (1x5)")
            self.msa.img_process(cross_convolve, col_size=17, row_size=7)
            self.gui.add_to_textbox("-- OP: Cross-convolving (7x17)")

        if filter_type.lower() == "aggressive":
            vsize = self.msa.nrows // 10
            vsize = vsize if vsize % 2 == 1 else vsize + 1
            self.msa.img_process(cross_convolve, col_size=vsize, row_size=3)
            self.gui.add_to_textbox(f"-- OP: Cross-convolving ({3}x{vsize})")

        self.msa_changed = True
        self.dendro_changed = True
        self.domains_changed = True

    def set_dendro_hcutoff(self, dendro_hcutoff: float):
        if self.dendro_hcutoff == dendro_hcutoff:
            return
        self.dendro_hcutoff = dendro_hcutoff
        self.dendro_changed = True
        self.domains_changed = True
        if self.colour_clusters:
            self.msa_changed = True

    # --- sequence selection

    def on_seq_selection(self, seq_id: str):
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

    def reset_selected_seq(self):
        self.highlight_selected_seq = False
        self.selected_seq = -1
        self.gui.set_slider_value(0)

    def select_ith_seqid(self, i):
        if self.reorder_rows:
            seq_id = self.msa.get_seq_indexer().get_ith_seqid_dendro(i)
        else:
            seq_id = self.msa.get_seq_indexer().get_ith_seqid_mat(i)
        self.on_seq_selection(seq_id)

    def get_selseq_idx(self):
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

    def toggle_split_mat_visualization(self, should_split: bool = False):
        self.split_mat_visualization = should_split
        self.msa_changed = True

    def toggle_hide_empty_cols(self, should_hide: bool = False):
        self.hide_empty_cols = should_hide
        self.msa_changed = True

    def toggle_reorder_mat_rows(self, should_reorder: bool = False):
        self.reorder_rows = should_reorder
        self.msa_changed = True
        if self.highlight_selected_seq and self.selected_seq != -1:
            self.gui.set_slider_value(self.get_selseq_idx())

    def toggle_colour_clusters(self, should_colour: bool = False):
        self.colour_clusters = should_colour
        self.msa_changed = True

    def toggle_highlight_selected_seq(self, should_highlight: bool = False):
        self.highlight_selected_seq = should_highlight
        if self.selected_seq != -1:
            self.msa_changed = True


    # --- display graphs

    def on_show_msa_mat(self, force: bool = False):
        """Controller is ordered to fetch an MSA matrix visualization figure and forward it to the GUI."""
        if not self.is_mat_initialized(): return
        if not force and not self.msa_changed: return

        mat = self.msa.get_mat(self.hide_empty_cols, self.reorder_rows)
        if mat is not None:
            wh_ratio = self.gui.get_mat_frame_wh_ratio()
            img = mat
            if self.colour_clusters:
                # TODO: consider mat vs. dendro ordering!
                cluster_labels = self.msa.get_cluster_labels(self.dendro_hcutoff, self.reorder_rows)
                img = color_clusters(mat, cluster_labels)
            img = highlight_row(img, self.get_selseq_idx())
            if self.split_mat_visualization:
                img = create_resized_mat_visualization(img, wh_ratio, not self.colour_clusters)
            fig = show_as_subimages(img, "")
            self.gui.show_matrix(fig)
            self.msa_changed = False

    def on_show_dendrogram(self, force: bool = False):
        """Controller is ordered to fetch a dendrogram visualization figure and forward it to the GUI."""
        if not self.is_mat_initialized(): return
        if not force and not self.dendro_changed: return

        fig = visualize_dendrogram(self.msa.get_linkage_mat(), self.dendro_hcutoff)
        self.gui.show_dendrogram(fig)
        self.dendro_changed = False

        cluster_labels = self.msa.get_cluster_labels(self.dendro_hcutoff)
        nclusters = len(set(cluster_labels))
        self.gui.set_cluster_count(nclusters)


    def on_show_domains(self):
        if not self.is_mat_initialized(): return
        if not self.domains_changed: return

        self.gui.add_to_textbox("Calculating domains.")
        domains = self.msa.retrieve_domains(self.dendro_hcutoff)
        fig = visualize_domains(domains)
        self.gui.show_domains(fig)
        self.domains_changed = False

    # --- query state

    def get_seq_indexer(self):
        return self.msa.get_seq_indexer() if self.is_mat_initialized() else None


    def is_mat_initialized(self):
        return self.msa is not None and self.msa.initialized
