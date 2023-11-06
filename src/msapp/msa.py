import copy
import numpy as np
from pysam import FastaFile, FastxFile
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist


import msapp.gconst as gc
from msapp.visualize import show, imgsave, visualize_clusters



class MultiSeqAlignment():
  """Class that contains the state of a multiple sequence alignment that is being processed."""

  def __init__(self, filename):
    self._filename = filename

    # preprocessing
    msa = FastaFile(filename)
    self.nrows = len(msa.lengths)
    self.ncols = np.max(msa.lengths)
    if np.min(msa.lengths) != self.ncols:
      raise ValueError("ERR: The MSA contains sequences of different lengths!")
    msa.close()

    self._to_binary_matrix()

    self._msa_mat_filtered = None
    self._filter_idx = None  # index of the sequence by which the MSA has been filtered
    self.filtered = False

    print(f"Successfully loaded MSA ({self.nrows} sequences of length {self.ncols})")

    if gc.VERBOSE: self.print_msa(1)
    if gc.DISPLAY: self.visualize()
    self._post_op()


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

  def filter_by_reference(self, idx, force=False):
    """Filters the alignment matrix by a given row in that MSA."""
    if self.filtered:
      if not force:
        print("WARN: There is already a filtered version. Use force to override")
        return
      else:
        print("INFO: Forcefully overriding stored filtered version")

    if idx >= self.nrows:
      print(f"ERR: Reference index {idx} is too large. Aborting.")
      return

    mat = copy.deepcopy(self._get_curr_mat())

    print(f"OP: Filtering by reference with index {idx}")
    refseq = self._msa_mat[idx]
    removed_cols = [i for i, val in enumerate(refseq) if val == 1]
    filtered = np.delete(mat, removed_cols, 1)  # delete all empty cols

    self.nrows = filtered.shape[0]
    self.ncols = filtered.shape[1]
    self._msa_mat_filtered = filtered
    self._filter_idx = idx
    self.filtered = True
    if gc.DISPLAY: self.visualize()
    self._post_op()


  ############### sort operations

  def sort_by_metric(self, sorting_metric):
    """Sorts a MSA binary matrix by the given function that works on a binary list."""
    print("OP: sort binary MSA rows")
    mat = copy.deepcopy(self._get_curr_mat())
    sorting_metrics = np.apply_along_axis(sorting_metric, axis=1, arr=mat)
    sorted_indices = np.argsort(sorting_metrics)
    sorted_matrix = mat[sorted_indices]

    self._msa_mat_filtered = sorted_matrix
    if gc.DISPLAY: self.visualize()
    self._post_op()


  ############### cluster operations, constructs implicit phylogenetic tree

  def cluster(self):
    """Clusters a MSA by some method.
    The Result is a number of clusters and their assigned sequences? Or positions, to directly identify domains?
    """
    # TODO: work in progress...

    mat = copy.deepcopy(self._get_curr_mat())

    # Calculate similarity matrix
    distance_matrix = pdist(mat, metric='hamming')
    linkage_mat = linkage(distance_matrix, method='complete')
    dendrogram(linkage_mat)

    visualize_clusters(mat, linkage_mat)
    self._post_op()


  ############### image processing

  def img_process(self, img_fun, *args, **kwargs):
    """Process the alignment by passing the current alignment matrix/image into the given image processing function."""
    # TODO: maybe create enum for all possible operations so as not to expose actual implementations to caller
    mat = img_fun(img=self._get_curr_mat(), *args, **kwargs)
    self._msa_mat_filtered = mat
    self._post_op()


  ############### analysis/metrics

  def calc_average_row_transitions(self):
    """Calculates the average number of value transitions per row as a metric for the noisiness of the MSA."""
    mat = self._get_curr_mat()
    assert mat.shape[0] == self.nrows and mat.shape[1] == self.ncols

    transitions = 0
    for row in range(self.nrows):
      current_state = mat[row, 0]
      for col in range(1, self.ncols):
        if mat[row, col] != current_state:
          transitions += 1
          current_state = mat[row, col]

    avg_row_transitions = float(transitions) / float(self.nrows)
    return avg_row_transitions


  def calc_average_col_transitions(self):
    """Calculates the average number of value transitions per column as a metric for the noisiness of the MSA."""
    mat = self._get_curr_mat()
    assert mat.shape[0] == self.nrows and mat.shape[1] == self.ncols

    transitions = 0
    for col in range(self.ncols):
      current_state = mat[0, col]
      for row in range(1, self.nrows):
        if mat[row, col] != current_state:
          transitions += 1
          current_state = mat[row, col]

    avg_col_transitions = float(transitions) / float(self.ncols)
    return avg_col_transitions


  def calc_cr_metric(self, verbose=None):
    rt = self.calc_average_row_transitions()
    ct = self.calc_average_col_transitions()
    prod = rt * ct
    if verbose if verbose is not None else gc.VERBOSE:
      print(f"Avg transitions. col: {ct:1.2f}, row: {rt:1.2f} -> prod: {prod:1.2f}")
    self.cr_score = prod
    return rt, ct


  def _post_op(self):
    """Operations to perform at the end of an MSA processing step."""
    self.calc_cr_metric(verbose=True)

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
      show(self._msa_mat, "Original Alignment (1 seq/row, white=gap)")
    else:
      aligned_str = "" if self._filter_idx is None else f", ref. row {self._filter_idx}"
      show(self._msa_mat_filtered, f"Filtered Alignment (1 seq/row, white=gap{aligned_str})")


  def save_to_file(self, filename):
    """Saves the final alignment image as well as the identified proteoform positions."""
    imgsave(self._msa_mat_filtered, filename)
    print(f"Wrote filtered alignment image to file out/{filename}.png")
    # TODO: invent file format and save proteoforms + meta information


  ############### getter/setter

  def _get_curr_mat(self):
    return self._msa_mat if self._msa_mat_filtered is None else self._msa_mat_filtered
