from tkinter import filedialog
import customtkinter

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

customtkinter.set_appearance_mode("System")
customtkinter.set_default_color_theme("blue")


class App(customtkinter.CTk):
    def __init__(self, controller):
        super().__init__()
        self.controller = controller

        # ------ configure window

        self.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.bind("<Configure>", self.on_resizing)

        self.title("MSA Proteoform Profiler")
        self.geometry(f"{1200}x{700}")
        self.minsize(width=1200, height=700)

        # configure grid layout (4x4)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure((2, 3), weight=0)
        self.grid_rowconfigure((0, 1, 2, 3), weight=1)

        # ------ left sidebar frame

        self.sidebar_frame = customtkinter.CTkFrame(self, width=140, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
        self.sidebar_frame.grid_rowconfigure(4, weight=1)
        self.logo_label = customtkinter.CTkLabel(self.sidebar_frame, text="MSAPP",
                                                 font=customtkinter.CTkFont(size=20, weight="bold"))
        self.logo_label.grid(row=0, column=0, padx=20, pady=(20, 10))
        self.sidebar_button_1 = customtkinter.CTkButton(self.sidebar_frame, text='Load FASTA',
                                                        command=self.on_load_file)
        self.sidebar_button_1.grid(row=1, column=0, padx=20, pady=10)
        self.sidebar_button_save = customtkinter.CTkButton(self.sidebar_frame, text="Save", command=self.on_save)
        self.sidebar_button_save.grid(row=2, column=0, padx=20, pady=10)
        self.sidebar_button_close = customtkinter.CTkButton(self.sidebar_frame, text="Quit", command=self.on_closing)
        self.sidebar_button_close.grid(row=3, column=0, padx=20, pady=10)
        self.appearance_mode_label = customtkinter.CTkLabel(self.sidebar_frame, text="Appearance Mode:", anchor="w")
        self.appearance_mode_label.grid(row=5, column=0, padx=20, pady=(10, 0))
        self.appearance_mode_optionemenu = customtkinter.CTkOptionMenu(self.sidebar_frame,
                                                                       values=["Light", "Dark", "System"],
                                                                       command=self.on_change_appearance_mode_event)
        self.appearance_mode_optionemenu.grid(row=6, column=0, padx=20, pady=(10, 10))

        # ------ CENTRAL VIEW

        # ------ MAT DISPLAY FRAME (upper middle)

        self.tabview_mat = customtkinter.CTkTabview(self)
        self.tabview_mat.grid(row=0, rowspan=2, column=1, padx=(20, 0), pady=(20, 0), sticky="nsew")
        self.tabview_mat.grid_columnconfigure(0, weight=1)
        self.tabview_mat.add("MSA")
        self.tabview_mat.add("Consensus")
        self.tabview_mat.tab("MSA").grid_columnconfigure(0, weight=1)
        self.tabview_mat.tab("MSA").grid_rowconfigure(0, weight=1)
        self.tabview_mat.tab("Consensus").grid_columnconfigure(0, weight=1)
        self.tabview_mat.tab("Consensus").grid_rowconfigure(0, weight=1)

        # ------ mat view

        self.mat_display_frame = customtkinter.CTkFrame(self.tabview_mat.tab("MSA"), fg_color="transparent")
        self.mat_display_frame.grid(row=0, column=0, padx=(0, 0), pady=(0, 0), sticky="nsew")
        self.mat_display_frame.grid_columnconfigure(0, weight=1)
        self.mat_display_frame.grid_rowconfigure(0, weight=1)
        self.mat_label = customtkinter.CTkLabel(self.mat_display_frame, text="MSA Matrix",
                                                font=customtkinter.CTkFont(size=15))
        self.mat_label.grid(row=0, column=0, padx=20, pady=(0, 0), sticky="n")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fig.set_figheight(4)
        ax.axis("off")
        self.canvas_mat = FigureCanvasTkAgg(fig, master=self.mat_display_frame)
        self.canvas_mat.draw()
        self.canvas_mat.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

        # ------ consensus view

        self.consensus_frame = customtkinter.CTkFrame(self.tabview_mat.tab("Consensus"), fg_color="transparent")
        self.consensus_frame.grid(row=0, column=0, padx=(0, 0), pady=(0, 0), sticky="nsew")
        self.consensus_frame.grid_columnconfigure(0, weight=1)
        self.consensus_frame.grid_rowconfigure(0, weight=1)
        self.mat_label = customtkinter.CTkLabel(self.consensus_frame, text="MSA Consensus",
                                                font=customtkinter.CTkFont(size=15))
        self.mat_label.grid(row=0, column=0, padx=20, pady=(0, 0), sticky="n")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fig.set_figheight(4)
        ax.axis("off")
        self.canvas_consensus = FigureCanvasTkAgg(fig, master=self.consensus_frame)
        self.canvas_consensus.draw()
        self.canvas_consensus.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

        # ------ DENDROGRAM DISPLAY FRAME (lower middle)

        self.tabview_dendro = customtkinter.CTkTabview(self)
        self.tabview_dendro.grid(row=2, rowspan=2, column=1, padx=(20, 0), pady=(20, 20), sticky="nsew")
        self.tabview_dendro.grid_columnconfigure(0, weight=1)
        self.tabview_dendro.add("Dendrogram")
        self.tabview_dendro.add("Dendro-Clusters")
        self.tabview_dendro.tab("Dendrogram").grid_columnconfigure(0, weight=1)
        self.tabview_dendro.tab("Dendrogram").grid_rowconfigure(0, weight=1)
        self.tabview_dendro.tab("Dendro-Clusters").grid_columnconfigure(0, weight=1)
        self.tabview_dendro.tab("Dendro-Clusters").grid_rowconfigure(0, weight=1)

        # ------ dendrogram view

        self.dendrogram_display_frame = customtkinter.CTkFrame(self.tabview_dendro.tab("Dendrogram"), fg_color="transparent")
        self.dendrogram_display_frame.grid(row=0, column=0, padx=(0, 0), pady=(0, 0), sticky="nsew")
        self.dendrogram_display_frame.grid_columnconfigure(0, weight=1)
        self.dendrogram_display_frame.grid_rowconfigure(0, weight=1)
        self.dendrogram_label = customtkinter.CTkLabel(self.dendrogram_display_frame, text="Dendrogram",
                                                       font=customtkinter.CTkFont(size=15))
        self.dendrogram_label.grid(row=0, column=0, padx=20, pady=(0, 0), sticky="n")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fig.set_figheight(4)
        ax.axis("off")
        self.canvas_dendro = FigureCanvasTkAgg(fig, master=self.dendrogram_display_frame)
        self.canvas_dendro.draw()
        self.canvas_dendro.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

        # ------ dendrogram height cluster view

        self.dendrogram_clustercount_frame = customtkinter.CTkFrame(self.tabview_dendro.tab("Dendro-Clusters"), fg_color="transparent")
        self.dendrogram_clustercount_frame.grid(row=0, column=0, padx=(0, 0), pady=(0, 0), sticky="nsew")
        self.dendrogram_clustercount_frame.grid_columnconfigure(0, weight=1)
        self.dendrogram_clustercount_frame.grid_rowconfigure(0, weight=1)
        self.dendro_height_label = customtkinter.CTkLabel(self.dendrogram_clustercount_frame, text="Cluster Count Per Height",
                                                       font=customtkinter.CTkFont(size=15))
        self.dendro_height_label.grid(row=0, column=0, padx=20, pady=(0, 0), sticky="n")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fig.set_figheight(4)
        ax.axis("off")
        self.canvas_dendro_clustercount = FigureCanvasTkAgg(fig, master=self.dendrogram_clustercount_frame)
        self.canvas_dendro_clustercount.draw()
        self.canvas_dendro_clustercount.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

        # ------ operations frame

        self.tabview = customtkinter.CTkTabview(self, width=250)
        self.tabview.grid(row=0, rowspan=2, column=2, padx=(20, 0), pady=(20, 0), sticky="nsew")
        self.tabview.add("Basic Operations")
        self.tabview.add("Tab 2")
        self.tabview.tab("Basic Operations").grid_columnconfigure(0, weight=1)
        self.tabview.tab("Tab 2").grid_columnconfigure(0, weight=1)

        # Tab 1
        self.button_filter_msa = customtkinter.CTkButton(self.tabview.tab("Basic Operations"), text="Filter MSA",
                                                         command=self.on_filter_msa)
        self.button_filter_msa.grid(row=1, column=0, padx=20, pady=(10, 10))

        # Tab 2
        self.optionmenu_2 = customtkinter.CTkOptionMenu(self.tabview.tab("Tab 2"), dynamic_resizing=False,
                                                        values=["Value 1", "Value 2", "Value 3"])
        self.optionmenu_2.grid(row=0, column=0, padx=20, pady=(20, 10))

        # ------ info frame

        self.scrollable_frame = customtkinter.CTkScrollableFrame(self, label_text="Visualization")
        self.scrollable_frame.grid(row=0, rowspan=2, column=3, padx=(20, 20), pady=(20, 0), sticky="nsew")
        self.scrollable_frame.grid_columnconfigure(0, weight=1)
        self.checkbox_var_hide_empty_cols = customtkinter.BooleanVar()
        self.checkbox_hide_empty_cols = customtkinter.CTkCheckBox(master=self.scrollable_frame,
                                                                  text="Hide empty columns",
                                                                  command=self.on_hide_empty_cols_switch,
                                                                  variable=self.checkbox_var_hide_empty_cols,
                                                                  onvalue=True, offvalue=False)
        self.checkbox_hide_empty_cols.grid(row=1, column=0, pady=10, padx=(10, 10), sticky="nw")

        self.checkbox_var_reorder_mat_rows = customtkinter.BooleanVar()
        self.checkbox_reorder_mat_rows = customtkinter.CTkCheckBox(master=self.scrollable_frame,
                                                                   text="Sort sequences",
                                                                   command=self.on_reorder_mat_rows_switch,
                                                                   variable=self.checkbox_var_reorder_mat_rows,
                                                                   onvalue=True, offvalue=False)
        self.checkbox_reorder_mat_rows.grid(row=2, column=0, pady=10, padx=(10, 10), sticky="nw")

        # ------ textbox

        self.logging_frame = customtkinter.CTkFrame(self)
        self.logging_frame.grid(row=2, rowspan=2, column=2, columnspan=2, padx=(20, 20), pady=(20, 20), sticky="nsew")
        self.logging_frame.grid_columnconfigure((0, 1), weight=1)
        self.logging_frame.grid_rowconfigure(0, weight=1)
        self.textbox = customtkinter.CTkTextbox(master=self.logging_frame)
        self.textbox.grid(row=0, column=0, columnspan=2, padx=(10, 10), pady=(10, 10), sticky="nsew")

        # ------ set default values

        self.appearance_mode_optionemenu.set("System")
        self.optionmenu_2.set("Option Menu")

        self.sidebar_button_save.configure(state="disabled")

        self.button_filter_msa.configure(state="disabled")
        self.checkbox_hide_empty_cols.deselect()
        self.controller.hide_empty_cols = False
        self.textbox.configure(state="disabled")

        self.last_mat_frame_hwratio = 0

    # ------ general

    def on_closing(self):
        self.quit()

    def on_change_appearance_mode_event(self, new_appearance_mode: str):
        customtkinter.set_appearance_mode(new_appearance_mode)

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
        if len(filename) != 0:
            success: bool = self.controller.initialize_from_file(filename)
            if success:
                self.controller.on_show_msa_mat()
                self.controller.on_show_dendrogram()
                self.controller.on_show_consensus()
                self.controller.on_show_dendro_clustercount()
                self.button_filter_msa.configure(state="normal")

    def on_filter_msa(self):
        self.controller.run_filtering_pipeline()
        self.controller.on_show_msa_mat()
        self.controller.on_show_dendrogram()
        self.controller.on_show_consensus()
        self.controller.on_show_dendro_clustercount()
        self.button_filter_msa.configure(state="disabled")

    def on_hide_empty_cols_switch(self):
        hide = self.checkbox_var_hide_empty_cols.get()
        self.controller.toggle_hide_empty_cols(hide)
        self.controller.on_show_msa_mat()

    def on_reorder_mat_rows_switch(self):
        reorder = self.checkbox_var_reorder_mat_rows.get()
        self.controller.toggle_reorder_mat_rows(reorder)
        self.controller.on_show_msa_mat()

    # ------ called by controller

    def show_matrix(self, mat_fig: Figure):
        plt.close()
        ax = mat_fig.subplots()
        ax.axis("off")
        mat_fig.set_figheight(4)
        mat_fig.set_tight_layout(True)
        self.canvas_mat = FigureCanvasTkAgg(mat_fig, master=self.mat_display_frame)
        self.canvas_mat.draw()
        self.canvas_mat.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="nsew")

    def show_dendrogram(self, dendro_fig: Figure):
        plt.close()
        ax = dendro_fig.subplots()
        ax.axis("off")
        dendro_fig.set_figheight(4)
        dendro_fig.set_tight_layout(True)
        self.canvas_dendro = FigureCanvasTkAgg(dendro_fig, master=self.dendrogram_display_frame)
        self.canvas_dendro.draw()
        self.canvas_dendro.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="new")

    def show_consensus(self, consensus_fig: Figure):
        plt.close()
        ax = consensus_fig.subplots()
        ax.axis("off")
        consensus_fig.set_figheight(4)
        consensus_fig.set_tight_layout(True)
        self.canvas_consensus = FigureCanvasTkAgg(consensus_fig, master=self.consensus_frame)
        self.canvas_consensus.draw()
        self.canvas_consensus.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="new")

    def show_dendro_clustercount(self, dc_fig: Figure):
        plt.close()
        ax = dc_fig.subplots()
        ax.axis("off")
        dc_fig.set_figheight(4)
        dc_fig.set_tight_layout(True)
        self.canvas_dendro_clustercount = FigureCanvasTkAgg(dc_fig, master=self.dendrogram_clustercount_frame)
        self.canvas_dendro_clustercount.draw()
        self.canvas_dendro_clustercount.get_tk_widget().grid(row=0, column=0, padx=(10, 10), pady=(25, 10), sticky="new")

    def add_to_textbox(self, text_to_add: str):
        self.textbox.configure(state="normal")
        self.textbox.insert(0.0, f"{text_to_add}\n")
        self.textbox.configure(state="disabled")

    def get_mat_frame_wh_ratio(self):
        """Get the current matrix frame width to height ratio with reduced height (for plot title + x-axis label)."""
        w = self.mat_display_frame._current_width
        h = self.mat_display_frame._current_height - 150
        wh_ratio = w / h
        return wh_ratio
