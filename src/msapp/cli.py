import argparse
import sys

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
  p.add_argument('-r', '--reference', type=int, default=0,
                 help='maps the general msa to the sequence with this index in the msa as a first step')
  args = p.parse_args(argv)
  return args


def run():
  print("Welcome to the msa proteoform profiler")
  argv = sys.argv[1:]
  p = parse_command_line(argv)
  if p.verbose: print("Finished parsing command line arguments")

  msa: MultiSeqAlignment = MultiSeqAlignment(p.msafile)
  if p.verbose: msa.print_msa(1)
  if p.display: msa.visualize()

  msa.filter_by_reference(p.reference) # 500?
  if p.display: msa.visualize()

  msa.dilate_erode(show=p.display)
  msa.blur(show=p.display)


