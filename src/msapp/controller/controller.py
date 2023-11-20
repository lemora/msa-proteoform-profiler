import copy

from msapp.model.msa import MultiSeqAlignment
from msapp.view.msapp_gui import App
from msapp.model.mat_manipulation import cross_convolve

from msapp.view.visualization import create_dendogram, create_dendogram_height_cluster_count_plot, show_as_subimages, \
    create_resized_mat_visualization


class Controller:
    """Manages the model (data) and view (GUI). It passes input and state changes between the two."""

    def __init__(self):
        """Constructor."""
        self.msa = None
        self.gui = None
        self.msa_changed = True  # dirty flag for msa visualization
        self.dendro_changed = True  # dirty flag for dendrogram visualization
        self.hide_empty_cols = False
        self.reorder_rows = False

    def start(self):
        """Starts the GUI."""
        if self.gui is None:
            self.gui = App(controller=self)
            self.gui.mainloop()

    def initialize_from_file(self, filename: str):
        """Creates and initializes a msa object from a given file name and if successful, shows it in the GUI."""
        try:
            msa = MultiSeqAlignment(filename)
        except Exception as e:
            err_msg = str(e)
            if len(err_msg) > 0:
                self.gui.add_to_textbox(f"ERR: {err_msg}")
            return

        fname_truncated = filename.split('/')[-1]
        if msa.initialized:
            self.msa = msa
            self.gui.add_to_textbox(f"Sucessfully loaded MSA from file '{fname_truncated}'.\n")
            self.msa_changed = True
            self.dendro_changed = True
        else:
            self.gui.add_to_textbox(f"Could not load MSA from file '{fname_truncated}'.")

    def cross_convolve_mat(self):
        if not self.is_mat_initialized(): return
        col_size = 5
        self.msa.img_process(cross_convolve, col_size=col_size)
        self.gui.add_to_textbox(f"Convolving with {col_size}x{col_size} cross-kernel (1x{col_size}).")
        self.msa_changed = True
        self.dendro_changed = True

    def run_filtering_pipeline(self):
        if not self.is_mat_initialized(): return
        self.gui.add_to_textbox("Running filtering pipeline.")
        self.msa.filter_by_length_statistic()
        self.gui.add_to_textbox("-- OP: Filtering by length statistic > 3 sigma.")
        self.msa.img_process(cross_convolve, col_size=5)
        self.gui.add_to_textbox("-- OP: Cross-convolving (1x5)")
        self.msa.img_process(cross_convolve, col_size=17, row_size=7)
        self.gui.add_to_textbox("-- OP: Cross-convolving (7x17)")
        self.msa_changed = True
        self.dendro_changed = True

    def is_mat_initialized(self):
        return self.msa is not None and self.msa.initialized

    def toggle_hide_empty_cols(self, should_hide: bool = False):
        self.hide_empty_cols = should_hide
        self.msa_changed = True

    def toggle_reorder_mat_rows(self, should_reorder: bool = False):
        self.reorder_rows = should_reorder
        self.msa_changed = True

    def show_msa_mat_old(self, force: bool = False):
        if not self.is_mat_initialized(): return
        if not force and not self.msa_changed: return

        max_row_width = self.gui.get_mat_frame_width()
        fig = self.msa.get_mat_visualization(self.hide_empty_cols, self.reorder_rows, max_row_width)
        self.gui.show_matrix(fig)
        self.msa_changed = False

    def show_msa_mat(self, force: bool = False):
        if not self.is_mat_initialized(): return
        if not force and not self.msa_changed: return

        mat = self.msa.get_mat(self.hide_empty_cols, self.reorder_rows)
        if mat is not None:
            wh_ratio = self.gui.get_mat_frame_wh_ratio()
            resized = create_resized_mat_visualization(mat, wh_ratio)
            fig = show_as_subimages(resized, "")
            self.gui.show_matrix(fig)
            self.msa_changed = False

    def show_dendrogram(self, force: bool = False):
        if not self.is_mat_initialized(): return
        if not force and not self.dendro_changed: return

        fig = create_dendogram(self.msa.get_linkage_mat())
        # fig = create_dendogram_height_cluster_count_plot(self.msa.get_linkage_mat())
        self.gui.show_dendrogram(fig)
        self.dendro_changed = False
