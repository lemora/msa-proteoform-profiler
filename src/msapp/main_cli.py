import argparse
import sys

import msapp.gconst as gc
from msapp.model.msa import MultiSeqAlignment
from msapp.view.out import show_save_results
from msapp.view.visualization import show


def parse_command_line(argv) -> argparse.ArgumentParser:
    """Process command line arguments."""
    p = argparse.ArgumentParser()
    p.add_argument('msafile', metavar='msafilename', type=str,
                   help='name of the msa file containing a multiple sequence alignment in fasta format')
    p.add_argument('-v', '--verbose', action='store_true', default=False, help='enable verbose mode')
    p.add_argument('-d', '--display', action='store_true', default=False,
                   help='display the visualized intermediate results')
    p.add_argument('-s', '--save', action='store_true', default=False, help='save the results at the end')
    p.add_argument('-c', '--dcutoff', type=float, default=0.75,
                   help='Dendrogram cut height to determine cluster count (1: root, single cluster; 0: leaves, '
                        'one cluster per distinct sequence)')
    p.add_argument('-f', '--filter', type=str, default="standard", choices=["mild", "standard", "aggressive"],
                   help="How aggressively to filter the alignment matrix.")
    p.add_argument('-m', '--domains', type=str, default="quick", choices=["quick", "thorough"],
                   help="How thoroughly to perform the domains calculation.")
    p.add_argument('-o', '--out', type=str, default="out",
                   help="The directory in which to store the results.")
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

    if p.dcutoff < 0.0 or p.dcutoff > 1.0:
        print("ERR: the dendrogram cutoff must be a float value between 0.0 and 1.0. Exiting.")
        quit()

    print("Settings:")
    print(f"- Verbose: {gc.VERBOSE}")
    print(f"- Display: {gc.DISPLAY}")
    print(f"- Dendrogram cutoff: {p.dcutoff}")
    print(f"- Filter mode: {p.filter}")
    print(f"- Domain calculation: {p.domains}")
    print(f"- Save: {p.save}")
    if p.save:
        print(f"- Output directory: {p.out}")
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
        show(msa.get_mat(hide_empty_cols, reorder_rows=False), "Binary MSA After Loading")

    # --- run filtering pipeline to remove noise

    msa.run_filtering_pipeline(filter_type=p.filter)
    if gc.VERBOSE: print()

    if gc.DISPLAY:
        show(msa.get_mat(hide_empty_cols, reorder_rows=False), "MSA After Image Processing")

    # --- optionally show/save results

    try:
        gc.VERBOSE = True
        show_save_results(msa, p.dcutoff, p.domains, p.save, p.out)
    except PermissionError as e:
        err_msg = str(e)
        if len(err_msg) > 0:
            print(f"ERR: {err_msg}")
        return
