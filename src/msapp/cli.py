import argparse
import sys

import msapp.gconst as gc
from msapp.msa import MultiSeqAlignment

def parse_command_line(argv) -> argparse.ArgumentParser:
  """Process command line arguments."""
  p = argparse.ArgumentParser()
  p.add_argument('msafile', metavar='msafilename', type=str,
                 help='name of the msa file containing a multiple sequence alignment in fasta format')
  p.add_argument('-v', '--verbose', action='store_true', default=False,
                 help='enable verbose mode')
  p.add_argument('-d', '--display', action='store_true', default=False,
                 help='displays the visualized results')
  p.add_argument('-s', '--save', action='store_true', default=False,
                 help='saves the results at the end')
  p.add_argument('-r', '--reference', type=int, default=0,
                 help='maps the general msa to the sequence with this index in the msa as a first step')
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

  # image processing pipeline to detect regions of interest

  msa.dilate_erode(ksize=3)
  msa.gaussian_blur(ksize=3)

  msa.filter_by_reference(p.reference)  # 500?

  if p.save: msa.save_to_file("proteoform-img")


