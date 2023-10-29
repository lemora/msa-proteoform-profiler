import copy
import numpy as np
from pysam import FastaFile, FastxFile

from msapp.visualize import show, imgsave
from msapp.imgmagic import blur, dilate_erode, gaussian_blur


class MultiSeqAlignment():
  """Class that contains the state of a multiple sequence alignment that is being processed."""

  def __init__(self, filename):
    self._filename = filename
    self._showsteps = False

    # preprocessing
    msa = FastaFile(filename)
    self.nrows = len(msa.lengths)
    self.ncols = np.max(msa.lengths)
    if np.min(msa.lengths) != self.ncols:
      raise ValueError("ERR: The MSA contains sequences of different lengths!")
    msa.close()

    self._to_binary_matrix()

    self._removed_cols = [] # stores the rows that were removed during filtering
    self._msa_mat_filtered = None
    print(f"Successfully loaded MSA ({self.nrows} sequences of size {self.ncols})")


  def _to_binary_matrix(self):
    """Creates a binary matrix from the alignment, where white (1) is a gap and black (0) an aligned position."""
    msa_mat = np.ones((self.nrows, self.ncols), dtype=np.uint8)
    with FastxFile(self._filename) as fh:
      for i, row in enumerate(fh):
        for j, letter in enumerate(row.sequence):
          if letter.isalpha():
            msa_mat[i, j] = 0 # black, aligned position
      fh.close()
    self._msa_mat = msa_mat


  ############### filter operations

  def filter_by_reference(self, reference_idx, force=False):
    """Filters the alignment matrix by a given row in that MSA."""
    if self._msa_mat_filtered is not None:
      if not force:
        print("WARN: There is already a filtered version. Use force to override")
        return
      else:
        print("INFO: Forcefully overriding stored filtered version")

    self._refseq = reference_idx
    refseq = self._msa_mat[reference_idx]
    ones = np.sum([1 for i in refseq if i == 0])
    removed_cols = [i for i, val in enumerate(refseq) if val == 1]
    # print("Removed cols:", len(removed_cols), removed_cols)

    msa_mat_filtered = np.delete(self._msa_mat, removed_cols, 1)  # delete all empty cols
    self._msa_mat_filtered = msa_mat_filtered
    self._removed_cols = removed_cols

    count_ones = sum([sum(array) for array in msa_mat_filtered])


  ############### image processing

  def blur(self, ksize=9, show=False):
    self._msa_mat_filtered = blur(self._msa_mat_filtered, ksize, show)

  def gaussian_blur(self, ksize=5, show=False):
    self._msa_mat_filtered = gaussian_blur(self._msa_mat_filtered, ksize, show)

  def dilate_erode(self, ksize=5, show=False):
    self._msa_mat_filtered = dilate_erode(self._msa_mat_filtered, ksize, show)


  ############### helper/debug

  def print_msa(self, n):
    if n < 1: return
    i = 0
    with FastxFile(self._filename) as fh:
      for entry in fh:
        print(entry.name)
        print(entry.sequence)
        print(entry.comment)
        print(entry.quality)
        i += 1
        if i >= n: break
      fh.close()


  def visualize(self):
    """Shows a binary image of the alignment."""
    if self._msa_mat_filtered is None:
      show(self._msa_mat, "Original Alignment (Row: Sequence)")
    else:
      show(self._msa_mat_filtered, "Filtered Alignment (Row: Sequence)")

  def save_to_file(self, filename):
    """Saves the final alignment image as well as the identified proteoform positions."""
    imgsave(self._msa_mat_filtered, filename)
    print(f"Wrote filtered alignment image to file out/{filename}.png")
    # TODO: invent file format and save proteoforms + meta information


  ############### getter/setter

  def get_msa_mat(self):
    return self._msa_mat

  def get_filtered_msa(self):
    return self._removed_cols, self._msa_mat_filtered
