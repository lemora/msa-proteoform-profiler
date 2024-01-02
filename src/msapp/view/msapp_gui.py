from tkinter import filedialog
import customtkinter as ctki

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from typing import Union, Callable

ctki.set_appearance_mode("System")
ctki.set_default_color_theme("blue")


class App(ctki.CTk):
    def __init__(self, controller):
        super().__init__()
        self.controller = controller

        # ------ configure window

        self.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.bind("<Configure>", self.on_resizing)

        self.title("MSA Proteoform Profiler")
        self.geometry(f"{1400}x{900}")
        self.minsize(width=1200, height=700)

        # configure grid layout (4x4)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure((2, 3), weight=0)
        self.grid_rowconfigure((0, 1, 2, 3), weight=1)

        # ------ left sidebar frame

        self.sidebar_frame = ctki.CTkFrame(self, width=140, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
        self.sidebar_frame.grid_rowconfigure(4, weight=1)
        self.logo_label = ctki.CTkLabel(self.sidebar_frame, text="MSAPP",
                                                 font=ctki.CTkFont(size=20, weight="bold"))
        self.logo_label.grid(row=0, column=0, padx=20, pady=(20, 10))
        self.sidebar_button_1 = ctki.CTkButton(self.sidebar_frame, text='Load FASTA',
                                                        command=self.on_load_file)
        self.sidebar_button_1.grid(row=1, column=0, padx=20, pady=10)
        self.sidebar_button_save = ctki.CTkButton(self.sidebar_frame, text="Save", command=self.on_save)
        self.sidebar_button_save.grid(row=2, column=0, padx=20, pady=10)
        self.sidebar_button_close = ctki.CTkButton(self.sidebar_frame, text="Quit", command=self.on_closing)
        self.sidebar_button_close.grid(row=3, column=0, padx=20, pady=10)
        self.appearance_mode_label = ctki.CTkLabel(self.sidebar_frame, text="Appearance Mode:", anchor="w")
        self.appearance_mode_label.grid(row=5, column=0, padx=20, pady=(10, 0))
        self.appearance_mode_optionemenu = ctki.CTkOptionMenu(self.sidebar_frame,
                                                                       values=["Light", "Dark", "System"],
                                                                       command=self.on_change_appearance_mode_event)
        self.appearance_mode_optionemenu.grid(row=6, column=0, padx=20, pady=(10, 10))

        # ------ CENTRAL VIEW

        # ------ MAT DISPLAY FRAME (upper middle)

        self.tabview_mat = ctki.CTkTabview(self)
        self.tabview_mat.grid(row=0, rowspan=2, column=1, padx=(20, 0), pady=(20, 0), sticky="nsew")
        self.tabview_mat.grid_columnconfigure(0, weight=1)
        self.tabview_mat.add("MSA")
        self.tabview_mat.add("Consensus")
        self.tabview_mat.tab("MSA").grid_columnconfigure(0, weight=1)
        self.tabview_mat.tab("MSA").grid_rowconfigure(0, weight=1)
        self.tabview_mat.tab("Consensus").grid_columnconfigure(0, weight=1)
        self.tabview_mat.tab("Consensus").grid_rowconfigure(0, weight=1)

        # ------ mat view

        self.mat_display_frame = ctki.CTkFrame(self.tabview_mat.tab("MSA"), fg_color="transparent")
        self.mat_display_frame.grid(row=0, column=0, padx=(0, 0), pady=(0, 0), sticky="nsew")
        self.mat_display_frame.grid_columnconfigure(0, weight=1)
        self.mat_display_frame.grid_rowconfigure(0, weight=1)
        self.mat_label = ctki.CTkLabel(self.mat_display_frame, text="MSA Matrix",
                                                font=ctki.CTkFont(size=15))
        self.mat_label.grid(row=0, column=0, padx=20, pady=(0, 0), sticky="n")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fig.set_figheight(4)
        ax.axis("off")
        self.canvas_mat = FigureCanvasTkAgg(fig, self.mat_display_frame)
        self.canvas_mat.draw()
        self.canvas_mat.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

        # ------ consensus view

        self.consensus_frame = ctki.CTkFrame(self.tabview_mat.tab("Consensus"), fg_color="transparent")
        self.consensus_frame.grid(row=0, column=0, padx=(0, 0), pady=(0, 0), sticky="nsew")
        self.consensus_frame.grid_columnconfigure(0, weight=1)
        self.consensus_frame.grid_rowconfigure(0, weight=1)
        self.mat_label = ctki.CTkLabel(self.consensus_frame, text="MSA Consensus",
                                                font=ctki.CTkFont(size=15))
        self.mat_label.grid(row=0, column=0, padx=20, pady=(0, 0), sticky="n")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fig.set_figheight(4)
        ax.axis("off")
        self.canvas_consensus = FigureCanvasTkAgg(fig, self.consensus_frame)
        self.canvas_consensus.draw()
        self.canvas_consensus.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

        # ------ DENDROGRAM DISPLAY FRAME (lower middle)

        self.tabview_dendro = ctki.CTkTabview(self)
        self.tabview_dendro.grid(row=2, rowspan=2, column=1, padx=(20, 0), pady=(20, 20), sticky="nsew")
        self.tabview_dendro.grid_columnconfigure(0, weight=1)
        self.tabview_dendro.add("Dendrogram")
        self.tabview_dendro.add("Dendro-Clusters")
        self.tabview_dendro.tab("Dendrogram").grid_columnconfigure(0, weight=1)
        self.tabview_dendro.tab("Dendrogram").grid_rowconfigure(0, weight=1)
        self.tabview_dendro.tab("Dendro-Clusters").grid_columnconfigure(0, weight=1)
        self.tabview_dendro.tab("Dendro-Clusters").grid_rowconfigure(0, weight=1)

        # ------ dendrogram view

        self.dendrogram_display_frame = ctki.CTkFrame(self.tabview_dendro.tab("Dendrogram"), fg_color="transparent")
        self.dendrogram_display_frame.grid(row=0, column=0, padx=(0, 0), pady=(0, 0), sticky="nsew")
        self.dendrogram_display_frame.grid_columnconfigure(0, weight=1)
        self.dendrogram_display_frame.grid_rowconfigure(0, weight=1)
        self.dendrogram_label = ctki.CTkLabel(self.dendrogram_display_frame, text="Dendrogram",
                                                       font=ctki.CTkFont(size=15))
        self.dendrogram_label.grid(row=0, column=0, padx=20, pady=(0, 0), sticky="n")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fig.set_figheight(4)
        ax.axis("off")
        self.canvas_dendro = FigureCanvasTkAgg(fig, self.dendrogram_display_frame)
        self.canvas_dendro.draw()
        self.canvas_dendro.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

        # ------ dendrogram height cluster view

        self.dendrogram_clustercount_frame = ctki.CTkFrame(self.tabview_dendro.tab("Dendro-Clusters"), fg_color="transparent")
        self.dendrogram_clustercount_frame.grid(row=0, column=0, padx=(0, 0), pady=(0, 0), sticky="nsew")
        self.dendrogram_clustercount_frame.grid_columnconfigure(0, weight=1)
        self.dendrogram_clustercount_frame.grid_rowconfigure(0, weight=1)
        self.dendro_height_label = ctki.CTkLabel(self.dendrogram_clustercount_frame, text="Cluster Count Per Height",
                                                       font=ctki.CTkFont(size=15))
        self.dendro_height_label.grid(row=0, column=0, padx=20, pady=(0, 0), sticky="n")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fig.set_figheight(4)
        ax.axis("off")
        self.canvas_dendro_clustercount = FigureCanvasTkAgg(fig, self.dendrogram_clustercount_frame)
        self.canvas_dendro_clustercount.draw()
        self.canvas_dendro_clustercount.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

        # ------ options/operations frame

        self.tabview = ctki.CTkTabview(self, width=280)
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

        self.filter_msa_selector = ctki.CTkOptionMenu(self.tabview.tab("Global"), width=110,
                                                     dynamic_resizing=False, values=["Standard", "Aggressive"])
        self.filter_msa_selector.grid(row=0, column=0, padx=(10, 5), pady=(15, 0), sticky="nw")

        self.button_filter_msa = ctki.CTkButton(self.tabview.tab("Global"), width=120,
                                                         text="Filter MSA", command=self.on_filter_msa)
        self.button_filter_msa.grid(row=0, column=1, padx=(0, 10), pady=(15, 0))

        # dendro cutoff spinbox
        self.dendro_cutoff_spinbox = FloatSpinbox(self.tabview.tab("Global"), width=110, step_size=0.05,
                                                                                    minval=0.25, maxval=1.0)
        self.dendro_cutoff_spinbox.grid(row=1, column=0, padx=(10, 5), pady=(10, 0), sticky="nw")

        self.button_dendro_hcutoff = ctki.CTkButton(self.tabview.tab("Global"), width=120, text="Use dendro cutoff",
                                                         command=self.on_dendro_height_change)
        self.button_dendro_hcutoff.grid(row=1, column=1, padx=(0, 10), pady=(10, 0))

        # visualization checkboxes
        self.checkbox_var_hide_empty_cols = ctki.BooleanVar()
        self.checkbox_hide_empty_cols = ctki.CTkCheckBox(self.tabview.tab("Global"),
                                                                  text="Hide empty columns",
                                                                  command=self.on_hide_empty_cols_switch,
                                                                  variable=self.checkbox_var_hide_empty_cols,
                                                                  onvalue=True, offvalue=False)
        self.checkbox_hide_empty_cols.grid(row=2, column=0, columnspan=2, pady=(20, 0), padx=(10, 10), sticky="nw")

        self.checkbox_var_reorder_mat_rows = ctki.BooleanVar()
        self.checkbox_reorder_mat_rows = ctki.CTkCheckBox(self.tabview.tab("Global"),
                                                                   text="Sort sequences",
                                                                   command=self.on_reorder_mat_rows_switch,
                                                                   variable=self.checkbox_var_reorder_mat_rows,
                                                                   onvalue=True, offvalue=False)
        self.checkbox_reorder_mat_rows.grid(row=3, column=0, columnspan=2, pady=(10, 0), padx=(10, 10), sticky="nw")


        # Tab 2: SingleSeq

        self.seq_selector = ctki.CTkOptionMenu(self.tabview.tab("SingleSeq"), width=110,
                                                        dynamic_resizing=False, values=["-"])
        self.seq_selector.grid(row=0, column=0, padx=(10, 5), pady=(15, 0), sticky="nw")

        self.button_choose_seq = ctki.CTkButton(self.tabview.tab("SingleSeq"), width=120, text="Select Sequence",
                                                         command=self.on_select_seq)
        self.button_choose_seq.grid(row=0, column=1, padx=(0, 10), pady=(15, 0))

        # selected sequence
        self.selected_seq_label = ctki.CTkLabel(self.tabview.tab("SingleSeq"), text="Details:",
                                                         font=ctki.CTkFont(size=15))
        self.selected_seq_label.grid(row=1, column=0, columnspan=2, padx=(10, 0), pady=(10, 0), sticky="nw")
        self.selected_seq_value = ctki.CTkLabel(self.tabview.tab("SingleSeq"), text="",
                                                         font=ctki.CTkFont(size=15))
        self.selected_seq_value.grid(row=3, column=0, padx=(10, 0), pady=(10, 0), sticky="nw")

        # ------ info frame

        self.scrollable_frame = ctki.CTkScrollableFrame(self, label_text="Info")
        self.scrollable_frame.grid(row=0, rowspan=2, column=3, padx=(20, 20), pady=(37, 0), sticky="nsew")
        self.scrollable_frame.grid_columnconfigure(0, weight=1)

        # filter score
        self.seq_count_label = ctki.CTkLabel(self.scrollable_frame, text="Sequence count:",
                                                         font=ctki.CTkFont(size=15))
        self.seq_count_label.grid(row=0, column=0, padx=(10, 0), pady=(10, 0), sticky="nw")
        self.seq_count_value = ctki.CTkLabel(self.scrollable_frame, text="",
                                                         font=ctki.CTkFont(size=15))
        self.seq_count_value.grid(row=0, column=1, padx=(0, 10), pady=(10, 0), sticky="nw")

        self.cluster_count_label = ctki.CTkLabel(self.scrollable_frame, text="Cluster count:",
                                                      font=ctki.CTkFont(size=15))
        self.cluster_count_label.grid(row=1, column=0, padx=(10, 0), pady=(5, 0), sticky="nw")
        self.cluster_count_value = ctki.CTkLabel(self.scrollable_frame, text="",
                                                      font=ctki.CTkFont(size=15))
        self.cluster_count_value.grid(row=1, column=1, padx=(0, 10), pady=(5, 0), sticky="nw")


        # ------ textbox

        self.logging_frame = ctki.CTkFrame(self)
        self.logging_frame.grid(row=2, rowspan=2, column=2, columnspan=2, padx=(20, 20), pady=(20, 20), sticky="nsew")
        self.logging_frame.grid_columnconfigure((0, 1), weight=1)
        self.logging_frame.grid_rowconfigure(0, weight=1)
        self.textbox = ctki.CTkTextbox(self.logging_frame)
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

        self.seq_selector.set("-")
        self.seq_selector.configure(state="disabled")
        self.button_choose_seq.configure(state="disabled")

        self.textbox.configure(state="disabled")
        self.last_mat_frame_hwratio = 0

    # ------ general

    def on_closing(self):
        self.quit()

    def on_change_appearance_mode_event(self, new_appearance_mode: str):
        ctki.set_appearance_mode(new_appearance_mode)

    def on_resizing(self, event):
        if event.widget == self:
            if getattr(self, "_after_id", None):
                self.after_cancel(self._after_id)
            self._after_id = self.after(500, lambda: self.eval_significant_resize_event())

    def eval_significant_resize_event(self):
        curr_mfratio = self.get_mat_frame_wh_ratio()
        if (abs(self.last_mat_frame_hwratio - curr_mfratio) > 0.2):
            self.controller.on_show_msa_mat(force=True),
            self.controller.on_show_dendrogram(force=True)
            self.last_mat_frame_hwratio = curr_mfratio

    # ------ calling the controller

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

        self.controller.on_show_msa_mat()
        self.controller.on_show_dendrogram()
        self.controller.on_show_consensus()
        self.controller.on_show_dendro_clustercount()

        self.button_filter_msa.configure(state="normal")
        self.filter_msa_selector.configure(state="normal")
        self.button_dendro_hcutoff.configure(state="normal")
        self.dendro_cutoff_spinbox.enable(True)
        self.seq_count_value.configure(text=self.controller.msa.nrows)
        # self.seq_selector.configure(state="normal", values=self.controller.msa.seq_names)

    def on_filter_msa(self):
        msa_filter_type = self.filter_msa_selector.get()
        self.controller.run_filtering_pipeline(msa_filter_type == "Aggressive")
        self.controller.on_show_msa_mat()
        self.controller.on_show_dendrogram()
        self.controller.on_show_consensus()
        self.controller.on_show_dendro_clustercount()
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

    def on_dendro_height_change(self):
        val = self.dendro_cutoff_spinbox.get()
        highlight_val = self.dendro_cutoff_spinbox.highlight_val
        if highlight_val == val: return
        self.dendro_cutoff_spinbox.set_highlighted(val)
        self.add_to_textbox(f"Dendrogram height cutoff set to {val}")
        self.controller.set_dendro_hcutoff(val)
        self.controller.on_show_dendrogram()
        self.controller.on_show_consensus()
        self.controller.on_show_dendro_clustercount()

    def on_select_seq(self):
        print("on select seq clicked")


    # ------ called by controller

    def show_matrix(self, mat_fig: Figure):
        plt.close()
        ax = mat_fig.subplots()
        ax.axis("off")
        mat_fig.set_figheight(4)
        mat_fig.set_tight_layout(True)
        self.canvas_mat = FigureCanvasTkAgg(mat_fig, self.mat_display_frame)
        self.canvas_mat.draw()
        self.canvas_mat.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

    def show_dendrogram(self, dendro_fig: Figure):
        plt.close()
        ax = dendro_fig.subplots()
        ax.axis("off")
        dendro_fig.set_figheight(4)
        dendro_fig.set_tight_layout(True)
        self.canvas_dendro = FigureCanvasTkAgg(dendro_fig, self.dendrogram_display_frame)
        self.canvas_dendro.draw()
        self.canvas_dendro.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="new")

    def show_consensus(self, consensus_fig: Figure):
        plt.close()
        ax = consensus_fig.subplots()
        ax.axis("off")
        consensus_fig.set_figheight(4)
        consensus_fig.set_tight_layout(True)
        self.canvas_consensus = FigureCanvasTkAgg(consensus_fig, self.consensus_frame)
        self.canvas_consensus.draw()
        self.canvas_consensus.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="new")

    def show_dendro_clustercount(self, dc_fig: Figure):
        plt.close()
        ax = dc_fig.subplots()
        ax.axis("off")
        dc_fig.set_figheight(4)
        dc_fig.set_tight_layout(True)
        self.canvas_dendro_clustercount = FigureCanvasTkAgg(dc_fig, self.dendrogram_clustercount_frame)
        self.canvas_dendro_clustercount.draw()
        self.canvas_dendro_clustercount.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="new")

    def add_to_textbox(self, text_to_add: str):
        self.textbox.configure(state="normal")
        self.textbox.insert(0.0, f"{text_to_add}\n")
        self.textbox.configure(state="disabled")

    def set_cluster_count(self, cluster_count: int):
        self.cluster_count_value.configure(text=cluster_count)

    def get_mat_frame_wh_ratio(self):
        """Get the current matrix frame width to height ratio with reduced height (for plot title + x-axis label)."""
        w = self.mat_display_frame._current_width
        h = self.mat_display_frame._current_height - 150
        wh_ratio = w / h
        return wh_ratio


# ------------------------------------------------------------------------

class FloatSpinbox(ctki.CTkFrame):
    def __init__(self, *args,
                 width: int = 120,
                 height: int = 32,
                 step_size: Union[int, float] = 0.1,
                 minval: float = 0.0,
                 maxval: float = 1.0,
                 command: Callable = None,
                 **kwargs):
        super().__init__(*args, width=width, height=height, **kwargs)

        self.step_size = step_size
        self.minval = minval
        self.maxval = maxval
        self.highlight_val = maxval
        self.command = command

        self.configure(fg_color=("gray78", "gray28"))  # set frame color

        self.grid_columnconfigure((0, 2), weight=0)  # buttons don't expand
        self.grid_columnconfigure(1, weight=1)  # entry expands

        self.subtract_button = ctki.CTkButton(self, text="-", width=height-6, height=height-6,
                                                       command=self.subtract_button_callback)
        self.subtract_button.grid(row=0, column=0, padx=(3, 0), pady=3)

        self.entry = ctki.CTkEntry(self, width=width-(2*height), height=height-6, border_width=0)
        self.entry.grid(row=0, column=1, columnspan=1, padx=3, pady=3, sticky="ew")
        self.entry.configure(state="disabled")

        self.add_button = ctki.CTkButton(self, text="+", width=height-6, height=height-6,
                                                  command=self.add_button_callback)
        self.add_button.grid(row=0, column=2, padx=(0, 3), pady=3)

        # default value
        self.entry.insert(0, "0.0")

    def add_button_callback(self):
        if self.command is not None:
            self.command()
        try:
            currval = round(float(self.entry.get()), 2)
            if currval + self.step_size >= self.maxval:
                value = self.maxval
            else:
                value = round(currval + self.step_size, 2)
            self.set(value)
        except ValueError:
            return

    def subtract_button_callback(self):
        if self.command is not None:
            self.command()
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
            self.entry.configure(text_color="#0059b3")
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