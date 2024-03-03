import matplotlib.pyplot as plt
from datetime import datetime
import os

import msapp.gconst as gc
from msapp.view.visualization import color_clusters, save_figure, show_as_subimages, visualize_dendrogram, \
    visualize_domains


def show_save_results(msa, dcutoff: float, dom_calc_mode: str, save: bool = True, out_dir="out"):
    res_dir = f'msapp-{datetime.now():%Y-%m-%d-%H:%M}'
    dir_path = f"{out_dir}/{res_dir}"
    if save:
        try:
            if not os.path.isdir(dir_path):
                os.makedirs(dir_path)
        except PermissionError:
            err_msg = f"Cannot write to directory '{out_dir}'."
            raise PermissionError(err_msg)
        if gc.VERBOSE: print(f"Writing to output directory: '{dir_path}'")

    fig_dendro = visualize_dendrogram(msa.get_linkage_mat(), dcutoff)
    if save:
        if gc.VERBOSE: print("Saving dendrogram...")
        save_figure(fig_dendro, "dendrogram", dir_path)

    mat = msa.get_mat(hide_empty_cols=True, reorder_rows=True)
    cluster_labels = msa.get_cluster_labels(perc_threshold=dcutoff, dendro_ordered=True)
    img = color_clusters(mat, cluster_labels)
    fig_mat = show_as_subimages(img, "MSA: Filtered, Reordered and Colored by Dendrogram")
    if save:
        if gc.VERBOSE: print("Saving colored MSA...")
        save_figure(fig_mat, "msa_colored_reordered", dir_path)
    plt.close(fig_mat)

    domains = msa.calculate_domains(dcutoff, dom_calc_mode)
    fig_domains = visualize_domains(domains)
    if save:
        if gc.VERBOSE: print("Saving domains...\n")
        save_figure(fig_domains, "detected_domains", dir_path)

    if gc.VERBOSE: print(f"Detected {len(domains)} proteoforms, with {[len(pf) for pf in domains]} domains.")
