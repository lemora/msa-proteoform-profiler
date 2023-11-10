import copy
import numpy as np
from pysam import FastaFile, FastxFile
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

import msapp.gconst as gc
from msapp.visualize import show, imgsave, visualize_clusters, hcluster, show_hist



class MultiSeqAlignment():
  """Class that contains the state of a multiple sequence alignment that is being processed."""

  def __init__(self, filename) -> None:
    """Constructor, loads msa, checks state and initializes binary matrix."""
    self._filename: str = filename
    msa_file = FastaFile(filename)
    self.nrows: int = len(msa_file.lengths)
    self.ncols: int = np.max(msa_file.lengths)
    if np.min(msa_file.lengths) != self.ncols:
      raise ValueError("ERR: The MSA contains sequences of different lengths!")
    msa_file.close()

    self._to_binary_matrix()
    self._filter_idx: int = None  # index of the sequence by which the MSA has been filtered
    self._filtered_by_reference: bool = False  # filtered

    print(f"Successfully loaded MSA ({self.nrows} sequences of length {self.ncols})")
    if gc.VERBOSE: self.print_msa(1)
    if gc.DISPLAY: self.visualize("Original")
    self._post_op()


  def _to_binary_matrix(self) -> None:
    """Creates a binary matrix from the alignment, where white (1) is a gap and black (0) an aligned position."""
    msa_mat = np.ones((self.nrows, self.ncols), dtype=np.uint8)
    with FastxFile(self._filename) as fh:
      for i, row in enumerate(fh):
        for j, letter in enumerate(row.sequence):
          if letter.isalpha():
            msa_mat[i, j] = 0 # black, aligned position
      fh.close()
    self._mat = msa_mat


  ############### filter operations

  def filter_by_length_statistic(self) -> None:
    """Creates a histogram of sequence lengths, then removes all longer than 3 std deviations above the median."""
    seq_lengths = [self.ncols - sum(row) for row in self._mat]
    median = np.median(seq_lengths)
    std = np.std(seq_lengths)
    print(f"Median: {median}, std: {std}")
    if gc.DISPLAY: show_hist(seq_lengths, nbins=100)

    three_std = median + 3 * std
    # print(f"3 std: {three_std}")
    seqs_over_three_std = [idx for idx, l in enumerate(seq_lengths) if l > three_std]
    self.remove_seqs_from_alignment(idx_list=seqs_over_three_std, cols=False)
    # print(f"Seqs over 3 std: {seqs_over_three_std}")
    self.remove_empty_cols()
    if gc.DISPLAY: self.visualize(rf"Removed {len(seqs_over_three_std)} seqs > 3 $\sigma$ long")
    self._post_op()


  def filter_by_reference(self, idx, force=False) -> None:
    """Filters the alignment matrix by a given row in that MSA. Final step, visualization purposes."""
    assert idx < self.nrows
    if self._filtered_by_reference:
      if not force:
        print("WARN: There is already a filtered version. Use force to override")
        return
      else:
        print("INFO: Forcefully overriding stored filtered version")

    mat = copy.deepcopy(self._mat)
    print(f"OP: Filtering by reference with index {idx}")
    refseq = mat[idx]
    cols_to_remove = [i for i, val in enumerate(refseq) if val == 1]
    self.remove_seqs_from_alignment(cols_to_remove, cols=True)
    self._filter_idx = idx
    self._filtered_by_reference = True

    if gc.DISPLAY: self.visualize(f"Filtered by reference {idx}")
    self._post_op()


  def remove_empty_cols(self, show: bool = False):
    """Removes all columns that are empty from the matrix, meaning they only contain the value 1."""
    empty_columns = np.where(np.all(self._mat == 1, axis=0))[0]
    self.remove_seqs_from_alignment(empty_columns, cols=True)
    if show: self.visualize(f"Removed {len(empty_columns)} empty columns")


  def remove_seqs_from_alignment(self, idx_list: list[int], cols: bool = True) -> None:
    """Removes the columns (else rows) which have the given indices.
    param cols: should remove columns, else remove rows."""
    if len(idx_list) == 0: return
    max_idx = self.ncols if cols else self.nrows
    assert min(idx_list) >= 0 and max(idx_list) < max_idx

    mat = copy.deepcopy(self._mat)
    filtered = np.delete(mat, idx_list, axis=(1 if cols else 0))
    self.nrows = filtered.shape[0]
    self.ncols = filtered.shape[1]
    self._mat = filtered


  ############### sort operations

  def sort_by_metric(self, sorting_metric) -> None:
    """Sorts a MSA binary matrix by the given function that works on a binary list."""
    print("OP: sort binary MSA rows")
    mat = copy.deepcopy(self._mat)
    sorting_metrics = np.apply_along_axis(sorting_metric, axis=1, arr=mat)
    sorted_indices = np.argsort(sorting_metrics)
    sorted_matrix = mat[sorted_indices]

    self._mat = sorted_matrix
    if gc.DISPLAY: self.visualize("Rows sorted")
    self._post_op()


  ############### image processing

  def img_process(self, img_fun, *args, **kwargs) -> None:
    """Process the alignment by passing the current alignment matrix/image into the given image processing function."""
    # TODO: maybe create enum for all possible operations so as not to expose actual implementations to caller
    mat = img_fun(img=self._mat, *args, **kwargs)
    self._mat = mat
    self._post_op()


  ############### cluster operations, constructs implicit phylogenetic tree

  def linkage_cluster(self, cmethod: str = "complete") -> None:
    """Clusters a MSA by some method.
    The Result is a number of clusters and their assigned sequences? Or positions, to directly identify domains?
    """
    # Calculate similarity matrix
    distance_matrix = pdist(self._mat, metric='hamming')
    linkage_mat = linkage(distance_matrix, method=cmethod)

    visualize_clusters(self._mat, linkage_mat)
    hcluster(self._mat, linkage_mat)

    self._post_op()


  ############### analysis/metrics

  def calc_average_row_transitions(self) -> float:
    """Calculates the average number of value transitions per row as a metric for the noisiness of the MSA."""
    assert self._mat.shape[0] == self.nrows and self._mat.shape[1] == self.ncols

    transitions = 0
    for row in range(self.nrows):
      current_state = self._mat[row, 0]
      for col in range(1, self.ncols):
        if self._mat[row, col] != current_state:
          transitions += 1
          current_state = self._mat[row, col]

    avg_row_transitions = float(transitions) / float(self.nrows)
    return avg_row_transitions


  def calc_average_col_transitions(self) -> float:
    """Calculates the average number of value transitions per column as a metric for the noisiness of the MSA."""
    assert self._mat.shape[0] == self.nrows and self._mat.shape[1] == self.ncols

    transitions = 0
    for col in range(self.ncols):
      current_state = self._mat[0, col]
      for row in range(1, self.nrows):
        if self._mat[row, col] != current_state:
          transitions += 1
          current_state = self._mat[row, col]

    avg_col_transitions = float(transitions) / float(self.ncols)
    return avg_col_transitions


  def calc_cr_metric(self, verbose=None) -> (float, float):
    """Calculates a column-row noisiness metric, that is the product of the average number of column value transitions
    times the average number of row value transitions."""
    rt = self.calc_average_row_transitions()
    ct = self.calc_average_col_transitions()
    prod = rt * ct
    if (verbose if verbose is not None else gc.VERBOSE):
      print(f"Avg transitions. col: {ct:1.2f}, row: {rt:1.2f} -> prod: {prod:1.2f}")
    self.cr_score = prod
    return rt, ct


  def _post_op(self) -> None:
    """Operations to perform at the end of an MSA processing step."""
    if gc.VERBOSE: self.calc_cr_metric(verbose=True)


  ############### helper/debug

  def print_msa(self, n) -> None:
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


  def visualize(self, title_addon: str = "") -> None:
    """Shows a binary image of the alignment."""
    aligned_str = "" if self._filter_idx is None else f", ref. row {self._filter_idx}"
    addstr = "" if title_addon == "" else f": {title_addon}"
    show(self._mat, f"MSA{addstr} (1 seq/row, white=gap{aligned_str})")


  def save_to_file(self, filename: str) -> None:
    """Saves the final alignment image as well as the identified proteoform positions."""
    imgsave(self._mat, filename)
    print(f"Wrote filtered alignment image to file out/{filename}.png")
    # TODO: save annotated image to file
    # TODO: invent file format and save proteoforms + meta information

