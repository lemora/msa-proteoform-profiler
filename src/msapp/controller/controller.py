import copy

from msapp.model.msa import MultiSeqAlignment
from msapp.view.msapp_gui import App
from msapp.model.img_processing import cross_convolve

from msapp.view.visualization import create_dendogram

class Controller:
    """Manages the model (data) and view (GUI). It passes input and state changes between the two."""

    def __init__(self):
        """Constructor."""
        self.msa = None
        self.gui = None
        self.hide_empty_cols = False

    def start(self):
        """Starts the GUI."""
        if self.gui is None:
            self.gui = App(controller=self)
            self.gui.mainloop()

    def initialize_from_file(self, filename: str):
        """Creates and initializes a msa object from a given file name and if successful, shows it in the GUI."""
        msa = MultiSeqAlignment(filename)
        fname_truncated = filename.split('/')[-1]
        if msa.initialized:
            self.msa = msa
            self.show_msa_mat()
            self.gui.add_to_textbox(f"Sucessfully loaded MSA from file '{fname_truncated}'.")
            self.show_dendrogram()
        else:
            self.gui.add_to_textbox(f"Could not load MSA from file '{fname_truncated}'.")

    def cross_convolve_mat(self):
        if not self.is_mat_initialized(): return
        col_size = 5
        self.msa.img_process(cross_convolve, col_size=col_size)
        self.gui.add_to_textbox(f"Convolving with {col_size}x{col_size} cross-kernel (1x{col_size}).")

    def is_mat_initialized(self):
        return self.msa is not None and self.msa.initialized

    def toggle_hide_empty_cols(self, should_hide: bool = False):
        if not self.is_mat_initialized():
            self.hide_empty_cols = should_hide
            return

        if self.hide_empty_cols is not should_hide:
            # update matrix
            self.hide_empty_cols = should_hide
            self.show_msa_mat()

    def show_msa_mat(self):
        if not self.is_mat_initialized():
            return
        max_row_width = self.gui.get_mat_frame_width()
        fig = self.msa.get_mat_visualization(self.hide_empty_cols, max_row_width=max_row_width)
        self.gui.show_matrix(fig)

    def show_dendrogram(self):
        if not self.is_mat_initialized(): return
        fig = create_dendogram(self.msa.get_linkage_mat())
        self.gui.show_dendrogram(fig)


