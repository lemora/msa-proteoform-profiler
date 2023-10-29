import argparse
import sys

def parse_command_line(argv) -> argparse.ArgumentParser:
  """Process command line arguments."""
  p = argparse.ArgumentParser()
  p.add_argument('mfile', metavar='msafilename', type=str,
                 help='name of the msa file containing a multiple sequence alignment in fasta format')
  p.add_argument('-v', '--verbose', action='store_true', default=False,
                 help='enable verbose mode')
  p.add_argument('-d', '--display', type=int, default=3,
                 help='displays the visualized results')
  args = p.parse_args(argv)
  return args


def run():
  print("Welcome to the msa proteoform profiler")
  argv = sys.argv[1:]
  p = parse_command_line(argv)
  if p.verbose: print("Finished parsing command line arguments")
