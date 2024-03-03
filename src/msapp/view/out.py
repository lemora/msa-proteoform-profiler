import matplotlib.pyplot as plt
from datetime import datetime
import os

from msapp.view.visualization import color_clusters, save_figure, show_as_subimages, visualize_dendrogram, \
    visualize_domains


def show_save_results(msa, dcutoff: float, dom_calc_mode: str, save: bool = True):
    the_dir = f'msapp-{datetime.now():%Y-%m-%d-%H:%M}'
    if save:
        pdir = f"out/{the_dir}"
        print(f"Writing to output directory: '{pdir}'")
        if not os.path.isdir(pdir):
            os.makedirs(pdir)

    fig_dendro = visualize_dendrogram(msa.get_linkage_mat(), dcutoff)
    if save:
        print("Saving dendrogram...")
        save_figure(fig_dendro, f"{the_dir}/dendrogram")

    mat = msa.get_mat(hide_empty_cols=True, reorder_rows=True)
    cluster_labels = msa.get_cluster_labels(perc_threshold=dcutoff, dendro_ordered=True)
    img = color_clusters(mat, cluster_labels)
    fig_mat = show_as_subimages(img, "MSA: Filtered, Reordered and Colored by Dendrogram")
    if save:
        print("Saving colored MSA...")
        save_figure(fig_mat, f"{the_dir}/msa_colored_reordered")
    plt.close(fig_mat)

    domains = msa.calculate_domains(dcutoff, dom_calc_mode)
    fig_domains = visualize_domains(domains)
    if save:
        print("Saving domains...\n")
        save_figure(fig_domains, f"{the_dir}/detected_domains")

    print(f"Detected {len(domains)} proteoforms, with {[len(pf) for pf in domains]} domains.")
