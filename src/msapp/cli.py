import argparse
import sys

import msapp.gconst as gc
from msapp.msa import MultiSeqAlignment
from msapp.imgmagic import convolve, dilate_erode, gaussian_blur, median_blur

def parse_command_line(argv) -> argparse.ArgumentParser:
  """Process command line arguments."""
  p = argparse.ArgumentParser()
  p.add_argument('msafile', metavar='msafilename', type=str,
                 help='name of the msa file containing a multiple sequence alignment in fasta format')
  p.add_argument('-v', '--verbose', action='store_true', default=False,
                 help='enable verbose mode')
  p.add_argument('-d', '--display', action='store_true', default=False,
                 help='display the visualized intermediate results')
  p.add_argument('-s', '--save', action='store_true', default=False,
                 help='save the results at the end')
  p.add_argument('-r', '--reorder', action='store_true', default=False,
                 help='reorder the rows in an alignment matrix')
  p.add_argument('-f', '--filter', type=int, default=0,
                 help='filter the msa based on the sequence with this index (default: 0)')
  args = p.parse_args(argv)
  return args


def run():
  print("Welcome to the msa proteoform profiler")
  argv = sys.argv[1:]
  p = parse_command_line(argv)
  gc.DISPLAY = p.display
  gc.VERBOSE = p.verbose
  if gc.VERBOSE:
    print("Finished parsing command line arguments")
    print(f"verbose:{gc.VERBOSE}; display:{gc.DISPLAY}")

  msa: MultiSeqAlignment = MultiSeqAlignment(p.msafile)

  # --- image processing pipeline to detect regions of interest

  # filtering as an early step just for now in order to see results of img processing better
  msa.filter_by_reference(p.filter)
  if p.reorder:
    sort_by = lambda row: sum(row)
    msa.sort_by_metric(sorting_metric=sort_by)

  msa.img_process(convolve, col_size=17, row_size=7)

  col_size = msa.nrows // 10
  col_size = col_size if col_size % 2 == 1 else col_size + 1
  # filter with a large emphasis on columns
  msa.img_process(convolve, col_size=col_size, row_size=3)

  if p.save: msa.save_to_file("proteoform-img")


