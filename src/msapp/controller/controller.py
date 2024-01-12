from msapp.model.msa import MultiSeqAlignment
from msapp.view.msapp_gui import App
from msapp.model.mat_manipulation import cross_convolve

from msapp.view.visualization import create_cluster_consensus_visualization, \
    create_dendrogram_height_cluster_count_plot, create_resized_mat_visualization, get_empty_plot, show_as_subimages, \
    visualize_dendrogram


class Controller:
    """Manages the model (data) and view (GUI). It passes input and state changes between the two."""

    def __init__(self):
        self.msa = None
        self.gui = None

        self.hide_empty_cols = False
        self.reorder_rows = False
        self.dendro_hcutoff = 0.75
        self.selected_seq = -1

        # dirty flags for visualizations that need to be refreshed
        self.msa_changed = True
        self.dendro_changed = True
        self.consensus_changed = True
        self.dclustercount_changed = True

    def start(self):
        """Starts the GUI."""
        if self.gui is None:
            self.gui = App(controller=self)
            self.gui.mainloop()

    def initialize_from_file(self, filename: str):
        """Creates and initializes a MultiSeqAlignment object from a given file name and if successful,
        shows the corresponding plots in the GUI."""
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
            self.consensus_changed = True
            self.dclustercount_changed = True
        else:
            self.gui.add_to_textbox(f"Could not load MSA from file '{fname_truncated}'.\n")
            return False
        return True

    def run_filtering_pipeline(self, aggressive: bool = False):
        """Triggers running a matrix filtering pipeline on the MultiSeqAlignment object."""
        if not self.is_mat_initialized(): return
        self.gui.add_to_textbox("Running filtering pipeline.")
        self.msa.filter_by_length_statistic()
        self.gui.add_to_textbox("-- OP: Filtering by length statistic > 3 sigma.")
        self.msa.img_process(cross_convolve, col_size=5)
        self.gui.add_to_textbox("-- OP: Cross-convolving (1x5)")
        self.msa.img_process(cross_convolve, col_size=17, row_size=7)
        self.gui.add_to_textbox("-- OP: Cross-convolving (7x17)")
        if aggressive:
            vsize = self.msa.nrows // 10
            vsize = vsize if vsize % 2 == 1 else vsize + 1
            self.msa.img_process(cross_convolve, col_size=vsize, row_size=3)
            self.gui.add_to_textbox(f"-- OP: Cross-convolving ({3}x{vsize})")
        self.msa_changed = True
        self.dendro_changed = True
        self.consensus_changed = True
        self.dclustercount_changed = True

    def is_mat_initialized(self):
        return self.msa is not None and self.msa.initialized

    def toggle_hide_empty_cols(self, should_hide: bool = False):
        self.hide_empty_cols = should_hide
        self.msa_changed = True

    def toggle_reorder_mat_rows(self, should_reorder: bool = False):
        self.reorder_rows = should_reorder
        self.msa_changed = True

    def set_dendro_hcutoff(self, dendro_hcutoff: float):
        if self.dendro_hcutoff == dendro_hcutoff:
            return
        self.dendro_hcutoff = dendro_hcutoff
        self.dendro_changed = True
        self.consensus_changed = True
        self.dclustercount_changed = True

    def set_selected_seq(self, selected_seq: str):
        if self.selected_seq == selected_seq:
            return
        self.selected_seq = selected_seq

    def get_selected_seq_mat_idx(self):
        idx = -1
        return idx

    def on_show_msa_mat(self, force: bool = False):
        """Controller is ordered to fetch an MSA matrix visualization figure and forward it to the GUI."""
        if not self.is_mat_initialized(): return
        if not force and not self.msa_changed: return

        mat = self.msa.get_mat(self.hide_empty_cols, self.reorder_rows)
        if mat is not None:
            wh_ratio = self.gui.get_mat_frame_wh_ratio()
            resized = create_resized_mat_visualization(mat, wh_ratio, self.get_selected_seq_mat_idx())
            fig = show_as_subimages(resized, "")
            self.gui.show_matrix(fig)
            self.msa_changed = False

    def on_show_dendrogram(self, force: bool = False):
        """Controller is ordered to fetch a dendrogram visualization figure and forward it to the GUI."""
        if not self.is_mat_initialized(): return
        if not force and not self.dendro_changed: return

        fig = visualize_dendrogram(self.msa.get_linkage_mat(), self.dendro_hcutoff)
        self.gui.show_dendrogram(fig)
        self.dendro_changed = False

    def on_show_consensus(self, force: bool = False):
        """Controller is ordered to fetch a consensus sequence visualization figure and forward it to the GUI."""
        if not self.is_mat_initialized(): return
        if not force and not self.consensus_changed: return

        consensus_clusters, cluster_sizes = self.msa.calc_consensus_clusters(perc_threshold=self.dendro_hcutoff)
        nclusters = len(consensus_clusters)
        self.gui.set_cluster_count(nclusters)
        fig = create_cluster_consensus_visualization(consensus_clusters, cluster_sizes)
        self.gui.show_consensus(fig)
        self.consensus_changed = False

    def on_show_dendro_clustercount(self, force: bool = False):
        """Controller is ordered to fetch a cluster per dendrogram height plot figure and forward it to the GUI."""
        if not self.is_mat_initialized(): return
        if not force and not self.dclustercount_changed: return

        fig = create_dendrogram_height_cluster_count_plot(self.msa.get_linkage_mat(), self.dendro_hcutoff)
        self.gui.show_dendro_clustercount(fig)
        self.consensus_changed = False

    def get_seq_indexer(self):
        return self.msa.get_seq_indexer() if self.is_mat_initialized() else None