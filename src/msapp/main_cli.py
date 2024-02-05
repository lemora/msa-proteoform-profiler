import argparse
from datetime import datetime
import os
import sys

import msapp.gconst as gc
from msapp.model.msa import MultiSeqAlignment
from msapp.view.visualization import (color_clusters, save_figure, show, show_as_subimages, visualize_dendrogram,
                                      visualize_domains)

def parse_command_line(argv) -> argparse.ArgumentParser:
    """Process command line arguments."""
    p = argparse.ArgumentParser()
    p.add_argument('msafile', metavar='msafilename', type=str,
                   help='name of the msa file containing a multiple sequence alignment in fasta format')
    p.add_argument('-v', '--verbose', action='store_true', default=False, help='enable verbose mode')
    p.add_argument('-d', '--display', action='store_true', default=False,
                   help='display the visualized intermediate results')
    p.add_argument('-s', '--save', action='store_true', default=False, help='save the results at the end')
    p.add_argument('-r', '--reorder', action='store_true', default=False,
                   help='reorder the rows in an alignment matrix')
    p.add_argument('-c', '--dcutoff', type=float, default=0.75,
                   help='Dendrogram cut height to determine cluster count (1: root, single cluster; 0: leaves, '
                        'one cluster per distinct sequence)')
    p.add_argument('-f', '--filter', type=str, default="standard", choices=["mild", "standard", "aggressive"],
                   help="How aggressively to filter the alignment matrix.")
    args = p.parse_args(argv)
    return args


def run() -> None:
    print("Welcome to the msa proteoform profiler.")
    print("------------------------------------------")
    argv = sys.argv[1:]
    p = parse_command_line(argv)
    gc.DISPLAY = p.display
    gc.VERBOSE = p.verbose
    hide_empty_cols = True
    reorder_rows = p.reorder

    if p.dcutoff < 0.0 or p.dcutoff > 1.0:
        print("ERR: the dendrogram cutoff must be a float value between 0.0 and 1.0. Exiting.")
        quit()

    print("Settings:")
    print(f"- Verbose: {gc.VERBOSE}")
    print(f"- Display: {gc.DISPLAY}")
    print(f"- Hide empty columns: {hide_empty_cols}")
    print(f"- Reorder rows: {reorder_rows}")
    print(f"- Dendrogram cutoff: {p.dcutoff}")
    print(f"- Filter mode: {p.filter}")
    print("-------------------------------------------")

    # --- attempt to load MSA

    try:
        msa: MultiSeqAlignment = MultiSeqAlignment(p.msafile)
    except Exception as e:
        print(f"ERR: Failed to load MSA from file '{p.msafile}'. Cause: {str(e)}")
        print("Exiting.")
        quit()

    if gc.VERBOSE: print()

    if gc.DISPLAY:
        show(msa.get_mat(hide_empty_cols, reorder_rows), "MSA after loading")

    # --- run filtering pipeline to remove noise

    msa.run_filtering_pipeline(filter_type=p.filter)
    if gc.VERBOSE: print()

    if gc.DISPLAY:
        show(msa.get_mat(hide_empty_cols, reorder_rows), "After image processing")

    # --- optionally show/save results

    the_dir = f'msapp-{datetime.now():%Y-%m-%d-%H:%M}'
    if p.save:
        pdir = f"out/{the_dir}"
        print(f"Writing to output directory: '{pdir}'")
        if not os.path.isdir(pdir):
            os.makedirs(pdir)

    fig_dendro = visualize_dendrogram(msa.get_linkage_mat(), p.dcutoff)
    if p.save:
        print("Saving dendrogram...")
        save_figure(fig_dendro, f"{the_dir}/dendrogram")

    mat = msa.get_mat(hide_empty_cols=True, reorder_rows=True)
    cluster_labels = msa.get_cluster_labels(perc_threshold=p.dcutoff, dendro_ordered=True)
    img = color_clusters(mat, cluster_labels)
    fig_mat = show_as_subimages(img, "")
    if p.save:
        print("Saving colored MSA...")
        save_figure(fig_mat, f"{the_dir}/msa_colored_reordered")

    domains = msa.retrieve_domains_via_dendrogram(p.dcutoff)
    fig_domains = visualize_domains(domains)
    if p.save:
        print("Saving domains...\n")
        save_figure(fig_domains, f"{the_dir}/detected_domains")

    print(f"Detected {len(domains)} proteoforms, with {[len(pf) for pf in domains]} domains.")
