from tkinter import filedialog, Listbox, END, ANCHOR
import customtkinter as ctk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from typing import Union

ctk.set_appearance_mode("System")
ctk.set_default_color_theme("blue")


class App(ctk.CTk):
    """The MSAPP GUI."""

    def __init__(self, controller):
        super().__init__()
        self.controller = controller

        # ------ configure window

        self.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.bind("<Configure>", self.on_resizing)

        self.title("MSA Proteoform Profiler")
        self.geometry(f"{1400}x{900}")
        self.minsize(width=1200, height=700)
        self.after(10, self._create_widgets)

    def _create_widgets(self):
        # configure grid layout (4x4)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure((2, 3), weight=0)
        self.grid_rowconfigure((0, 1, 2, 3), weight=1)

        # ------ left sidebar frame

        self.sidebar_frame = ctk.CTkFrame(self, width=140, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
        self.sidebar_frame.grid_rowconfigure(4, weight=1)
        self.logo_label = ctk.CTkLabel(self.sidebar_frame, text="MSAPP", font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_label.grid(row=0, column=0, padx=20, pady=(20, 10))
        self.sidebar_button_1 = ctk.CTkButton(self.sidebar_frame, text='Load FASTA', command=self.on_load_file)
        self.sidebar_button_1.grid(row=1, column=0, padx=20, pady=10)
        self.sidebar_button_save = ctk.CTkButton(self.sidebar_frame, text="Save", command=self.on_save)
        self.sidebar_button_save.grid(row=2, column=0, padx=20, pady=10)
        self.sidebar_button_close = ctk.CTkButton(self.sidebar_frame, text="Quit", command=self.on_closing)
        self.sidebar_button_close.grid(row=3, column=0, padx=20, pady=10)
        self.appearance_mode_label = ctk.CTkLabel(self.sidebar_frame, text="Appearance:", anchor="w")
        self.appearance_mode_label.grid(row=5, column=0, padx=20, pady=(10, 0))
        self.appearance_mode_optionemenu = ctk.CTkOptionMenu(self.sidebar_frame, values=["Light", "Dark", "System"],
                                                             command=self.on_change_appearance_mode_event)
        self.appearance_mode_optionemenu.grid(row=6, column=0, padx=20, pady=(10, 10))

        # ------ CENTRAL VIEW

        # ------ MAT DISPLAY FRAME (upper middle)

        self.tabview_mat = ctk.CTkTabview(self)
        self.tabview_mat.grid(row=0, rowspan=2, column=1, padx=(20, 0), pady=(20, 0), sticky="nsew")
        self.tabview_mat.grid_columnconfigure(0, weight=1)
        self.tabview_mat.add("MSA")
        self.tabview_mat.add("Consensus")
        self.tabview_mat.tab("MSA").grid_columnconfigure(0, weight=1)
        self.tabview_mat.tab("MSA").grid_rowconfigure(0, weight=1)
        self.tabview_mat.tab("Consensus").grid_columnconfigure(0, weight=1)
        self.tabview_mat.tab("Consensus").grid_rowconfigure(0, weight=1)

        # ------ mat view

        self.mat_display_frame = ctk.CTkFrame(self.tabview_mat.tab("MSA"), fg_color="transparent")
        self.mat_display_frame.grid(row=0, column=0, padx=(0, 0), pady=(0, 0), sticky="nsew")
        self.mat_display_frame.grid_columnconfigure(0, weight=1)
        self.mat_display_frame.grid_rowconfigure(0, weight=1)
        self.mat_label = ctk.CTkLabel(self.mat_display_frame, text="MSA Matrix", font=ctk.CTkFont(size=15))
        self.mat_label.grid(row=0, column=0, padx=20, pady=(0, 0), sticky="n")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fig.set_figheight(4)
        ax.axis("off")
        self.canvas_mat = FigureCanvasTkAgg(fig, self.mat_display_frame)
        self.canvas_mat.draw()
        self.canvas_mat.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

        # ------ consensus view

        self.consensus_frame = ctk.CTkFrame(self.tabview_mat.tab("Consensus"), fg_color="transparent")
        self.consensus_frame.grid(row=0, column=0, padx=(0, 0), pady=(0, 0), sticky="nsew")
        self.consensus_frame.grid_columnconfigure(0, weight=1)
        self.consensus_frame.grid_rowconfigure(0, weight=1)
        self.mat_label = ctk.CTkLabel(self.consensus_frame, text="MSA Consensus", font=ctk.CTkFont(size=15))
        self.mat_label.grid(row=0, column=0, padx=20, pady=(0, 0), sticky="n")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fig.set_figheight(4)
        ax.axis("off")
        self.canvas_consensus = FigureCanvasTkAgg(fig, self.consensus_frame)
        self.canvas_consensus.draw()
        self.canvas_consensus.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

        # ------ DENDROGRAM DISPLAY FRAME (lower middle)

        self.tabview_dendro = ctk.CTkTabview(self)
        self.tabview_dendro.grid(row=2, rowspan=2, column=1, padx=(20, 0), pady=(20, 20), sticky="nsew")
        self.tabview_dendro.grid_columnconfigure(0, weight=1)
        self.tabview_dendro.add("Dendrogram")
        self.tabview_dendro.add("Domains")
        self.tabview_dendro.tab("Dendrogram").grid_columnconfigure(0, weight=1)
        self.tabview_dendro.tab("Dendrogram").grid_rowconfigure(0, weight=1)
        self.tabview_dendro.tab("Domains").grid_columnconfigure(0, weight=1)
        self.tabview_dendro.tab("Domains").grid_rowconfigure(0, weight=1)

        # ------ dendrogram view

        self.dendrogram_display_frame = ctk.CTkFrame(self.tabview_dendro.tab("Dendrogram"), fg_color="transparent")
        self.dendrogram_display_frame.grid(row=0, column=0, padx=(0, 0), pady=(0, 0), sticky="nsew")
        self.dendrogram_display_frame.grid_columnconfigure(0, weight=1)
        self.dendrogram_display_frame.grid_rowconfigure(0, weight=1)
        self.dendrogram_label = ctk.CTkLabel(self.dendrogram_display_frame, text="Dendrogram",
                                             font=ctk.CTkFont(size=15))
        self.dendrogram_label.grid(row=0, column=0, padx=20, pady=(0, 0), sticky="n")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fig.set_figheight(4)
        ax.axis("off")
        self.canvas_dendro = FigureCanvasTkAgg(fig, self.dendrogram_display_frame)
        self.canvas_dendro.draw()
        self.canvas_dendro.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

        # ------ prediceted domains view

        self.domains_frame = ctk.CTkFrame(self.tabview_dendro.tab("Domains"),
                                          fg_color="transparent")
        self.domains_frame.grid(row=0, column=0, padx=(0, 0), pady=(0, 0), sticky="nsew")
        self.domains_frame.grid_columnconfigure(0, weight=1)
        self.domains_frame.grid_rowconfigure(0, weight=1)
        self.domains_label = ctk.CTkLabel(self.domains_frame, text="Predicted protein domains",
                                          font=ctk.CTkFont(size=15))
        self.domains_label.grid(row=0, column=0, padx=20, pady=(0, 0), sticky="n")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fig.set_figheight(4)
        ax.axis("off")
        self.canvas_domain_view = FigureCanvasTkAgg(fig, self.domains_frame)
        self.canvas_domain_view.draw()
        self.canvas_domain_view.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10),
                                                     sticky="nsew")

        # ------ options/operations frame

        self.tabview = ctk.CTkTabview(self, width=280)
        self.tabview.grid(row=0, rowspan=2, column=2, padx=(20, 0), pady=(20, 0), sticky="nsew")
        self.tabview.add("Global")
        self.tabview.add("SingleSeq")
        self.tabview.tab("Global").configure(width=280)
        self.tabview.tab("SingleSeq").configure(width=280)
        self.tabview.tab("Global").grid_columnconfigure(0, weight=1)
        self.tabview.tab("Global").grid_columnconfigure(1, weight=0)
        self.tabview.tab("SingleSeq").grid_columnconfigure(0, weight=1)
        self.tabview.tab("SingleSeq").grid_columnconfigure(1, weight=0)

        # Tab 1: Global

        self.filter_msa_selector = ctk.CTkOptionMenu(self.tabview.tab("Global"), width=110, dynamic_resizing=False,
                                                     values=["Standard", "Aggressive"])
        self.filter_msa_selector.grid(row=0, column=0, padx=(10, 5), pady=(15, 0), sticky="nw")

        self.button_filter_msa = ctk.CTkButton(self.tabview.tab("Global"), width=120, text="Filter MSA",
                                               command=self.on_filter_msa)
        self.button_filter_msa.grid(row=0, column=1, padx=(0, 10), pady=(15, 0))

        # dendro cutoff spinbox
        self.dendro_cutoff_spinbox = FloatSpinbox(self.tabview.tab("Global"), width=110, step_size=0.05, minval=0.25,
                                                  maxval=1.0)
        self.dendro_cutoff_spinbox.grid(row=1, column=0, padx=(10, 5), pady=(10, 0), sticky="nw")

        self.button_dendro_hcutoff = ctk.CTkButton(self.tabview.tab("Global"), width=120, text="Use dendro cutoff",
                                                   command=self.on_dendro_height_change)
        self.button_dendro_hcutoff.grid(row=1, column=1, padx=(0, 10), pady=(10, 0))

        # visualization checkboxes
        self.checkbox_var_hide_empty_cols = ctk.BooleanVar()
        self.checkbox_hide_empty_cols = ctk.CTkCheckBox(self.tabview.tab("Global"), text="Hide empty columns",
                                                        command=self.on_hide_empty_cols_switch,
                                                        variable=self.checkbox_var_hide_empty_cols, onvalue=True,
                                                        offvalue=False)
        self.checkbox_hide_empty_cols.grid(row=2, column=0, columnspan=2, pady=(20, 0), padx=(10, 10), sticky="nw")

        self.checkbox_var_reorder_mat_rows = ctk.BooleanVar()
        self.checkbox_reorder_mat_rows = ctk.CTkCheckBox(self.tabview.tab("Global"), text="Sort sequences",
                                                         command=self.on_reorder_mat_rows_switch,
                                                         variable=self.checkbox_var_reorder_mat_rows, onvalue=True,
                                                         offvalue=False)
        self.checkbox_reorder_mat_rows.grid(row=3, column=0, columnspan=2, pady=(10, 0), padx=(10, 10), sticky="nw")

        self.checkbox_var_colour_clusters = ctk.BooleanVar()
        self.checkbox_colour_clusters = ctk.CTkCheckBox(self.tabview.tab("Global"), text="Colour clusters",
                                                         command=self.on_colour_clusters_switch,
                                                         variable=self.checkbox_var_colour_clusters, onvalue=True,
                                                         offvalue=False)
        self.checkbox_colour_clusters.grid(row=4, column=0, columnspan=2, pady=(10, 0), padx=(10, 10), sticky="nw")

        # Tab 2: SingleSeq

        self.button_choose_seq = ctk.CTkButton(self.tabview.tab("SingleSeq"), width=120, text="Select Sequence",
                                               command=self.on_open_search_seq_dialog)
        self.button_choose_seq.grid(row=0, column=0, columnspan=2, padx=(0, 10), pady=(15, 0))

        # selected sequence
        self.seq_slider = ctk.CTkSlider(self.tabview.tab("SingleSeq"), from_=0, to=100, command=self.on_slider_event)
        self.seq_slider.grid(row=1, column=0, columnspan=2, padx=(5, 0), pady=(5, 0), sticky="ew")

        self.selected_seq_label = ctk.CTkLabel(self.tabview.tab("SingleSeq"), text="Selected sequence:",
                                               font=ctk.CTkFont(size=14), text_color="gray")
        self.selected_seq_label.grid(row=2, column=0, columnspan=2, padx=(10, 0), pady=(10, 0), sticky="nw")
        self.selected_seq_textfield = ctk.CTkTextbox(self.tabview.tab("SingleSeq"), font=ctk.CTkFont(size=14),
                                                     fg_color="transparent", height=100)
        self.selected_seq_textfield.grid(row=4, column=0, columnspan=2, padx=(5, 0), pady=(5, 0), sticky="ew")

        self.checkbox_var_highlight_selected_seq = ctk.BooleanVar()
        self.checkbox_highlight_selected_seq = ctk.CTkCheckBox(self.tabview.tab("SingleSeq"),
                                                               text="Show selected sequence",
                                                               command=self.on_highlight_selected_seq,
                                                               variable=self.checkbox_var_highlight_selected_seq,
                                                               onvalue=True, offvalue=False)
        self.checkbox_highlight_selected_seq.grid(row=5, column=0, columnspan=2, pady=(0, 0), padx=(10, 10),
                                                  sticky="nw")

        # ------ info frame

        self.scrollable_frame = ctk.CTkScrollableFrame(self, label_text="Info")
        self.scrollable_frame.grid(row=0, rowspan=2, column=3, padx=(20, 20), pady=(37, 0), sticky="nsew")
        self.scrollable_frame.grid_columnconfigure(0, weight=1)

        # general MSA info section
        self.info_msa_label = ctk.CTkLabel(self.scrollable_frame, text="General", width=120,
                                           font=ctk.CTkFont(size=14, weight="bold"), text_color="gray28",
                                           bg_color="gray90")
        self.info_msa_label.grid(row=0, column=0, columnspan=3, padx=(5, 5), pady=(5, 5), sticky="ew")

        self.seq_count_label = ctk.CTkLabel(self.scrollable_frame, text="Sequence count:", font=ctk.CTkFont(size=14),
                                            text_color="gray")
        self.seq_count_label.grid(row=1, column=0, padx=(10, 0), pady=(5, 0), sticky="nw")
        self.seq_count_value = ctk.CTkLabel(self.scrollable_frame, text="", font=ctk.CTkFont(size=14))
        self.seq_count_value.grid(row=1, column=1, padx=(0, 10), pady=(5, 0), sticky="nw")

        self.msa_length_label = ctk.CTkLabel(self.scrollable_frame, text="Alignment length:",
                                             font=ctk.CTkFont(size=14), text_color="gray")
        self.msa_length_label.grid(row=2, column=0, padx=(10, 0), pady=(5, 0), sticky="nw")
        self.msa_length_value = ctk.CTkLabel(self.scrollable_frame, text="", font=ctk.CTkFont(size=14))
        self.msa_length_value.grid(row=2, column=1, padx=(0, 10), pady=(5, 0), sticky="nw")

        # clustering info section
        self.info_clustering_label = ctk.CTkLabel(self.scrollable_frame, text="Clustering", width=120,
                                                  font=ctk.CTkFont(size=14, weight="bold"), text_color="gray28",
                                                  bg_color="gray90")
        self.info_clustering_label.grid(row=3, column=0, columnspan=3, padx=(5, 5), pady=(20, 5), sticky="ew")

        self.cluster_count_label = ctk.CTkLabel(self.scrollable_frame, text="Cluster count:",
                                                font=ctk.CTkFont(size=14), text_color="gray")
        self.cluster_count_label.grid(row=4, column=0, padx=(10, 0), pady=(5, 0), sticky="nw")
        self.cluster_count_value = ctk.CTkLabel(self.scrollable_frame, text="", font=ctk.CTkFont(size=14))
        self.cluster_count_value.grid(row=4, column=1, padx=(0, 10), pady=(5, 0), sticky="nw")

        # TODO: add suggested clickable clusterings

        # ------ textbox

        self.logging_frame = ctk.CTkFrame(self)
        self.logging_frame.grid(row=2, rowspan=2, column=2, columnspan=2, padx=(20, 20), pady=(20, 20), sticky="nsew")
        self.logging_frame.grid_columnconfigure((0, 1), weight=1)
        self.logging_frame.grid_rowconfigure(0, weight=1)
        self.textbox = ctk.CTkTextbox(self.logging_frame)
        self.textbox.grid(row=0, column=0, columnspan=2, padx=(10, 10), pady=(10, 10), sticky="nsew")

        # ------ set default values

        self.appearance_mode_optionemenu.set("System")
        self.sidebar_button_save.configure(state="disabled")

        self.button_filter_msa.configure(state="disabled")
        self.checkbox_hide_empty_cols.deselect()
        self.controller.hide_empty_cols = False

        self.dendro_cutoff_spinbox.set_highlighted(0.75)
        self.dendro_cutoff_spinbox.enable(False)
        self.button_dendro_hcutoff.configure(state="disabled")

        self.button_choose_seq.configure(state="disabled")
        self.selected_seq_textfield.configure(state="disabled")
        self.seq_slider.set(0)
        self.seq_slider.configure(state="disabled")

        self.textbox.configure(state="disabled")
        self.last_mat_frame_hwratio = 0

    # ------ general

    def on_closing(self):
        self.quit()

    def on_change_appearance_mode_event(self, new_appearance_mode: str):
        ctk.set_appearance_mode(new_appearance_mode)

    def on_resizing(self, event):
        if event.widget == self:
            if getattr(self, "_after_id", None):
                self.after_cancel(self._after_id)
            self._after_id = self.after(500, lambda: self.eval_significant_resize_event())

    def eval_significant_resize_event(self):
        curr_mfratio = self.get_mat_frame_wh_ratio()
        if abs(self.last_mat_frame_hwratio - curr_mfratio) > 0.2:
            self.controller.on_show_msa_mat(force=True),
            self.controller.on_show_dendrogram(force=True)
            self.last_mat_frame_hwratio = curr_mfratio

    # ------ callbacks triggered by user interactions

    def on_save(self):
        msg = "'Save' clicked, but not yet implemented."
        print(msg)
        self.add_to_textbox(msg)

    def on_load_file(self):
        filename = filedialog.askopenfilename(title="Select a File",
                                              filetypes=(("fasta files", "*.fa*"), ("all files", "*.*")))
        if len(filename) == 0: return
        success: bool = self.controller.initialize_from_file(filename)
        if not success: return
        self.dendro_cutoff_spinbox.set_highlighted(0.75)
        self.controller.set_dendro_hcutoff(0.75)
        self.controller.reset_selected_seq()

        self.controller.on_show_msa_mat()
        self.controller.on_show_dendrogram()
        # self.controller.on_show_consensus()
        # self.controller.on_show_domains()

        self.button_filter_msa.configure(state="normal")
        self.filter_msa_selector.configure(state="normal")
        self.button_dendro_hcutoff.configure(state="normal")
        self.dendro_cutoff_spinbox.enable(True)
        self.seq_count_value.configure(text=self.controller.msa.nrows)
        self.msa_length_value.configure(text=self.controller.msa.ncols)
        self.button_choose_seq.configure(state="normal")
        self.set_selected_seq_info(" ")
        self.seq_slider.configure(state="normal")
        self.seq_slider.configure(to=self.controller.msa.nrows-1)

    def on_filter_msa(self):
        msa_filter_type = self.filter_msa_selector.get()
        self.controller.run_filtering_pipeline(msa_filter_type == "Aggressive")
        self.controller.on_show_msa_mat()
        self.controller.on_show_dendrogram()
        # self.controller.on_show_consensus()
        # self.controller.on_show_domains()
        self.filter_msa_selector.configure(state="disabled")
        self.button_filter_msa.configure(state="disabled")

    def on_hide_empty_cols_switch(self):
        hide = self.checkbox_var_hide_empty_cols.get()
        self.controller.toggle_hide_empty_cols(hide)
        self.controller.on_show_msa_mat()

    def on_reorder_mat_rows_switch(self):
        reorder = self.checkbox_var_reorder_mat_rows.get()
        self.controller.toggle_reorder_mat_rows(reorder)
        self.controller.on_show_msa_mat()

    def on_colour_clusters_switch(self):
        colour = self.checkbox_var_colour_clusters.get()
        self.controller.toggle_colour_clusters(colour)
        self.controller.on_show_msa_mat()

    def on_highlight_selected_seq(self):
        should_highlight = self.checkbox_var_highlight_selected_seq.get()
        self.controller.toggle_highlight_selected_seq(should_highlight)
        self.controller.on_show_msa_mat()

    def on_dendro_height_change(self):
        val = self.dendro_cutoff_spinbox.get()
        highlight_val = self.dendro_cutoff_spinbox.highlight_val
        if highlight_val == val: return
        self.dendro_cutoff_spinbox.set_highlighted(val)
        self.add_to_textbox(f"Dendrogram height cutoff set to {val}")
        self.controller.set_dendro_hcutoff(val)
        self.controller.on_show_dendrogram()
        self.controller.on_show_msa_mat()
        # self.controller.on_show_consensus()
        # self.controller.on_show_domains()

    def on_open_search_seq_dialog(self):
        SeqSearchDialogue(self, title="Select a Sequence")

    def on_seq_selection(self, seq_id: str):
        self.controller.on_seq_selection(seq_id)
        self.controller.on_show_msa_mat()

    def set_slider_value(self, val):
        self.seq_slider.set(val)

    def on_slider_event(self, value):
        intval = int(value)
        self.controller.select_ith_seqid(intval)
        self.controller.on_show_msa_mat()

    # ------ called by controller

    def show_matrix(self, mat_fig: Figure):
        plt.close()
        ax = mat_fig.subplots()
        ax.axis("off")
        mat_fig.set_figheight(4)
        mat_fig.set_tight_layout(True)

        self.canvas_mat.get_tk_widget().destroy()
        self.canvas_mat = FigureCanvasTkAgg(mat_fig, self.mat_display_frame)
        self.canvas_mat.draw()
        self.canvas_mat.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

    def show_dendrogram(self, dendro_fig: Figure):
        plt.close()
        ax = dendro_fig.subplots()
        ax.axis("off")
        dendro_fig.set_figheight(4)
        dendro_fig.set_tight_layout(True)

        self.canvas_dendro.get_tk_widget().destroy()
        self.canvas_dendro = FigureCanvasTkAgg(dendro_fig, self.dendrogram_display_frame)
        self.canvas_dendro.draw()
        self.canvas_dendro.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="new")

    def show_consensus(self, consensus_fig: Figure):
        plt.close()
        ax = consensus_fig.subplots()
        ax.axis("off")
        consensus_fig.set_figheight(4)
        consensus_fig.set_tight_layout(True)

        self.canvas_consensus.get_tk_widget().destroy()
        self.canvas_consensus = FigureCanvasTkAgg(consensus_fig, self.consensus_frame)
        self.canvas_consensus.draw()
        self.canvas_consensus.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="new")

    def show_domains(self, dc_fig: Figure):
        plt.close()
        ax = dc_fig.subplots()
        ax.axis("off")
        dc_fig.set_figheight(4)
        dc_fig.set_tight_layout(True)

        self.canvas_domain_view.get_tk_widget().destroy()
        self.canvas_domain_view = FigureCanvasTkAgg(dc_fig, self.domains_frame)
        self.canvas_domain_view.draw()
        self.canvas_domain_view.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10),
                                                     sticky="new")

    def add_to_textbox(self, text_to_add: str):
        self.textbox.configure(state="normal")
        self.textbox.insert(0.0, f"{text_to_add}\n")
        self.textbox.configure(state="disabled")

    def set_cluster_count(self, cluster_count: int):
        self.cluster_count_value.configure(text=cluster_count)

    def set_selected_seq_info(self, info_text: str):
        if info_text is None or len(info_text) == 0: return
        self.selected_seq_textfield.configure(state="normal")
        self.selected_seq_textfield.delete(0.0, END)
        self.selected_seq_textfield.insert(0.0, info_text)
        self.selected_seq_textfield.configure(state="disabled")

    def get_mat_frame_wh_ratio(self):
        """Get the current matrix frame width to height ratio with reduced height (for plot title + x-axis label)."""
        w = self.mat_display_frame._current_width
        h = self.mat_display_frame._current_height - 150
        wh_ratio = w / h
        return wh_ratio


# ------------------------------------------------------------------------

class SeqSearchDialogue(ctk.CTkToplevel):
    """Dialogue window for searching/selecting a sequence."""

    def __init__(self, main_window, title: str):
        super().__init__()
        self.main_window = main_window
        self._title = title
        self.title(self._title)
        self.lift()
        self.attributes("-topmost", True)
        self.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.after(10, self._create_widgets)
        self.grab_set()

    def _create_widgets(self):
        self.geometry(f"{700}x{400}")
        self.minsize(width=400, height=400)
        self.grid_columnconfigure((0, 1), weight=1)
        self.rowconfigure(0, weight=1)

        self._entry_field = ctk.CTkEntry(master=self, width=230)
        self._entry_field.grid(row=0, column=0, columnspan=2, padx=10, pady=(0, 0), sticky="ew")

        self._list_box = Listbox(master=self)
        self._list_box.grid(row=1, rowspan=2, column=0, columnspan=2, padx=20, pady=20, sticky="ew")

        self._ok_button = ctk.CTkButton(master=self, width=100, border_width=0, text='Select',
                                        command=self.on_ok_event)
        self._ok_button.grid(row=3, column=0, columnspan=1, padx=(20, 10), pady=(0, 20), sticky="ew")

        self._cancel_button = ctk.CTkButton(master=self, width=100, border_width=0, text='Cancel',
                                            command=self.on_cancel_event)
        self._cancel_button.grid(row=3, column=1, columnspan=1, padx=(10, 20), pady=(0, 20), sticky="ew")

        self._entry_field.bind("<KeyRelease>", self.on_key_release_event)
        self._list_box.bind("<Double-Button-1>", self.on_click_select_event)
        self.after(150, lambda: self._entry_field.focus())

        self.seq_indexer = self.main_window.controller.get_seq_indexer()
        self.set_to_listbox(self.get_matching_seqs(""))

    def on_key_release_event(self, event) -> None:
        """Something has been typed. Finds and shows sequences that match."""
        typed = self._entry_field.get()
        if self.seq_indexer is None:
            self.set_to_listbox([])
            return
        matching_seqs = self.get_matching_seqs(typed)
        self.set_to_listbox(matching_seqs)

    def on_click_select_event(self, event) -> None:
        """A sequence has been selected by double-clicking on the list."""
        selected = self._list_box.get(ANCHOR)
        if len(selected) == 0:
            return

        self._entry_field.delete(0, END)
        self._entry_field.insert(0, selected)
        seq_id = selected.split()[0]
        matching_seqs = self.get_matching_seqs(seq_id)
        self.set_to_listbox(matching_seqs)

    def get_matching_seqs(self, the_str: str) -> list:
        """Returns a list of tuples (seq_id:str, text:str) where one or
        both entries contain the given string."""
        if self.seq_indexer is None:
            return []
        matching_seqs = self.seq_indexer.get_seq_infos_containing_string(the_str)
        matching_seqs.sort(key=lambda tup: tup[0])
        return matching_seqs

    def set_to_listbox(self, data: list) -> None:
        """param data: list of tuples (id:str, description:str)"""
        self._list_box.delete(0, END)
        if len(data) == 0: return
        for item in data:
            self._list_box.insert(END, f"{item[0]} -- {item[1]}")

    def on_cancel_event(self) -> None:
        self.on_closing()

    def on_ok_event(self) -> None:
        """A sequence has been selected. Propagate the choice."""
        user_input = self._entry_field.get()
        if len(user_input) != 0:
            seq_id = user_input.split()[0]
            self.main_window.on_seq_selection(seq_id)
        self.on_closing()

    def on_closing(self):
        self.grab_release()
        self.destroy()


# ------------------------------------------------------------------------

class FloatSpinbox(ctk.CTkFrame):
    """Spinbox for float values."""

    def __init__(self, *args, width: int = 120, height: int = 32, step_size: Union[int, float] = 0.1,
                 minval: float = 0.0, maxval: float = 1.0, **kwargs):
        super().__init__(*args, width=width, height=height, **kwargs)

        self.step_size = step_size
        self.minval = minval
        self.maxval = maxval
        self.highlight_val = maxval

        self.configure(fg_color=("gray78", "gray28"))  # set frame color

        self.grid_columnconfigure((0, 2), weight=0)  # buttons don't expand
        self.grid_columnconfigure(1, weight=1)  # entry expands

        self.subtract_button = ctk.CTkButton(self, text="-", width=height - 6, height=height - 6,
                                             command=self.on_subtract_event)
        self.subtract_button.grid(row=0, column=0, padx=(3, 0), pady=3)

        self.entry = ctk.CTkEntry(self, width=width - (2 * height), height=height - 6, border_width=0)
        self.entry.grid(row=0, column=1, columnspan=1, padx=3, pady=3, sticky="ew")
        self.entry.configure(state="disabled")

        self.add_button = ctk.CTkButton(self, text="+", width=height - 6, height=height - 6, command=self.on_add_event)
        self.add_button.grid(row=0, column=2, padx=(0, 3), pady=3)

        self.entry.insert(0, "0.0")

    def on_add_event(self):
        try:
            currval = round(float(self.entry.get()), 2)
            if currval + self.step_size >= self.maxval:
                value = self.maxval
            else:
                value = round(currval + self.step_size, 2)
            self.set(value)
        except ValueError:
            return

    def on_subtract_event(self):
        try:
            currval = round(float(self.entry.get()), 2)
            if currval - self.step_size <= self.minval:
                value = self.minval
            else:
                value = round(currval - self.step_size, 2)
            self.set(value)
        except ValueError:
            return

    def get(self) -> Union[float, None]:
        try:
            return float(self.entry.get())
        except ValueError:
            return None

    def set(self, value: float):
        self.entry.configure(state="normal")
        self.entry.delete(0, "end")
        self.entry.insert(0, str(float(value)))
        if value == self.highlight_val:
            self.entry.configure(text_color="#0059b3")  # blue
        else:
            self.entry.configure(text_color="black")
        self.entry.configure(state="disabled")

    def set_highlighted(self, highlight_val: float):
        self.highlight_val = highlight_val
        self.set(highlight_val)

    def enable(self, enable: float = True):
        enable_str = "normal" if enable else "disabled"
        self.add_button.configure(state=enable_str)
        self.subtract_button.configure(state=enable_str)
