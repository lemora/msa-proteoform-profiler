import argparse
import sys

import msapp.gconst as gc
from msapp.msa import MultiSeqAlignment
from msapp.imgmagic import cross_convolve


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
    p.add_argument('-f', '--filter', type=int, default=-1,
                   help='filter the msa based on the sequence with this index (default: 0)')
    p.add_argument('-c', '--cut_height_dendogram', type=float, default=0.5,
                   help='Dendogram cut height to determine cluster count (1: root, single cluster; 0: leaves, '
                        'one cluster per distinct sequence)')
    args = p.parse_args(argv)
    return args


def run() -> None:
    print("Welcome to the msa proteoform profiler")
    print("------------------------------------------")
    argv = sys.argv[1:]
    p = parse_command_line(argv)
    gc.DISPLAY = p.display
    gc.VERBOSE = p.verbose
    gc.DENDOGRAM_CUT_HEIGHT = p.cut_height_dendogram
    if gc.VERBOSE:
        print("Finished parsing command line arguments")
        print(f"verbose:{gc.VERBOSE}; display:{gc.DISPLAY}")

    msa: MultiSeqAlignment = MultiSeqAlignment(p.msafile)

    # --- remove sequences that are much too long; > 3 sigma

    sort_by = lambda row: sum(row)
    msa.filter_by_length_statistic()
    msa.remove_empty_cols(show=gc.DISPLAY)

    # --- image processing to remove noise

    msa.img_process(cross_convolve, col_size=17, row_size=7)

    col_size = msa.nrows // 10
    col_size = col_size if col_size % 2 == 1 else col_size + 1
    msa.img_process(cross_convolve, col_size=col_size, row_size=3)

    # cleaning step, optional. Loses horizontal fragmentation
    msa.remove_empty_cols(show=gc.DISPLAY)

    # --- clustering

    if p.filter > 0:
        msa.filter_by_reference(p.filter)

    if p.reorder: msa.sort_by_metric(sorting_metric=sort_by)

    gc.DISPLAY = True

    msa.linkage_cluster()

    if p.save: msa.save_to_file("proteoform-img")
