import copy
import numpy as np
from pysam import FastxFile
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import mode

import msapp.gconst as gc
from msapp.visualize import imgsave, visualize_clusters, visualize_cluster_consensuses, show, show_hist


class MultiSeqAlignment():
  """Class that contains the state of a multiple sequence alignment that is being processed."""

  def __init__(self, filename: str = "") -> None:
    """Constructor. If filename given, loads msa, checks state and initializes binary matrix."""
    # TODO: track indices of rows and cols to be able to reconstruct sequences/indices
    self.initialized = False
    self._filename = ""
    self.nrows: int = -1
    self.ncols: int = -1
    self._ridx: np.array = []  # list of row indices
    self._cidx: np.array = []  # list of col indices
    self._filtered_by_reference: bool = False  # filtered
    self._filter_idx: int = -1  # index of the sequence by which the MSA has been filtered
    self.nclusters: int = 3  # likely number of discovered isoforms. maybe later list ranked by likelihood?

    if filename is not None and filename != "":
      self.init_from_file(filename)


  def init_from_mat(self, binary_mat: np.array):
    """Initializes a MultiSeqAlignment object with a given binary matrix. Mainly for testing.
    arg binary_mat: two-dimensional np.array with np.uint8 values"""
    if self.initialized: raise ValueError("The MSA object has already been initialized.")
    assert binary_mat is not None

    self._mat = binary_mat
    self.nrows: int = binary_mat.shape[0]
    self.ncols: int = binary_mat.shape[1]
    self._ridx = np.array([i for i in range(self.nrows)])
    self._cidx = np.array([i for i in range(self.ncols)])
    self.initialized = True
    print(f"Successfully loaded MSA ({self.nrows} sequences of length {self.ncols})")


  def init_from_file(self, filename: str):
    """Initializes a MultiSeqAlignment object from a msa fasta file."""
    if self.initialized: raise ValueError("The MSA object has already been initialized.")
    self._filename: str = filename

    # Create a binary matrix from the multiple sequence alignment file
    msa_mat = []
    with FastxFile(self._filename) as fh:
      for row in fh:
        the_row = []
        for letter in row.sequence:
          the_row.append(0 if letter.isalpha() else 1)
        msa_mat.append(the_row)
      fh.close()

    msa_mat = np.array(msa_mat, dtype=np.uint8)
    self._mat = msa_mat
    self.nrows = msa_mat.shape[0]
    self.ncols = msa_mat.shape[1]
    self._ridx = np.array([i for i in range(self.nrows)])
    self._cidx = np.array([i for i in range(self.ncols)])
    self.initialized = True

    print(f"Successfully loaded MSA ({self.nrows} sequences of length {self.ncols})\n")
    if gc.VERBOSE: self.print_msa(1)
    if gc.DISPLAY: self.visualize("Original")
    self._post_op()


  ############### filter operations

  def filter_by_length_statistic(self) -> None:
    """Creates a histogram of sequence lengths, then removes all longer than 3 std deviations above the median."""
    print("\n-- OP: Filtering by length statistic > 3 sigma.")

    seq_lengths = [self.ncols - sum(row) for row in self._mat]
    median = np.median(seq_lengths)
    std = np.std(seq_lengths)
    if gc.DISPLAY: show_hist(seq_lengths, nbins=100)

    over_three_std = median + 3 * std
    seqs_over_three_std = np.array([idx for idx, l in enumerate(seq_lengths) if l > over_three_std])
    self.remove_seqs_from_alignment(idx_list=seqs_over_three_std, cols=False)

    if gc.DISPLAY: self.visualize(rf"Removed {len(seqs_over_three_std)} seqs > 3 $\sigma$ length")
    self._post_op()


  def filter_by_reference(self, idx: int, force=False) -> None:
    """Filters the alignment matrix by a given row in that MSA. Final step, visualization purposes."""
    if idx >= self.nrows: raise ValueError("The index needs to be smaller than the number of rows.")

    if self._filtered_by_reference:
      if not force:
        print("WARN: There is already a filtered version. Use force to override")
        return
      else:
        print("INFO: Forcefully overriding stored filtered version")

    print(f"\n-- OP: Filter by reference with index {idx}")
    refseq = self._mat[idx]
    cols_to_remove = np.array([i for i, val in enumerate(refseq) if val == 1])
    self.remove_seqs_from_alignment(cols_to_remove, cols=True)
    self._filter_idx = idx
    self._filtered_by_reference = True

    if gc.DISPLAY: self.visualize(f"Filtered by reference {idx}")
    self._post_op()


  def remove_empty_cols(self, show: bool = False):
    """Removes all columns that are empty from the matrix, meaning they only contain the value 1."""
    print("\n-- OP: Removing empty columns.")
    empty_columns = np.where(np.all(self._mat == 1, axis=0))[0]
    self.remove_seqs_from_alignment(empty_columns, cols=True)
    if show: self.visualize(f"Removed {len(empty_columns)} empty columns")
    self._post_op()


  def remove_seqs_from_alignment(self, idx_list: np.ndarray[int], cols: bool = True) -> None:
    """Removes the columns (else rows) which have the given indices.
    param cols: should remove columns, else remove rows."""
    if len(idx_list) == 0: return
    max_idx = self.ncols if cols else self.nrows
    if min(idx_list) < 0 or max(idx_list) >= max_idx:
      raise ValueError("The indices o delete must be between 0 and max row or col count.")

    # update index lists
    if cols:
      self._cidx = np.delete(self._cidx, idx_list)
    else:
      self._ridx = np.delete(self._ridx, idx_list)

    mat = copy.deepcopy(self._mat)
    filtered = np.delete(mat, idx_list, axis=(1 if cols else 0))
    self.nrows = filtered.shape[0]
    self.ncols = filtered.shape[1]
    self._mat = filtered


  ############### sort operations

  def sort_by_metric(self, sorting_metric = lambda row: sum(row)) -> None:
    """Sorts a MSA binary matrix by the given function that works on a binary list."""
    print("\n-- OP: sorting MSA rows.")
    mat = copy.deepcopy(self._mat)
    sorting_metrics = np.apply_along_axis(sorting_metric, axis=1, arr=mat)
    
    # print("Sorting metrics (maybe idx odering?): ", sorting_metrics)
    # TODO: update index lists!

    sorted_indices = np.argsort(sorting_metrics)
    sorted_matrix = mat[sorted_indices]


    self._mat = sorted_matrix
    if gc.DISPLAY: self.visualize("Rows sorted")
    self._post_op()


  ############### image processing

  def img_process(self, img_fun, *args, **kwargs) -> None:
    """Process the alignment by passing the current alignment matrix/image into the given image processing function."""
    # TODO: maybe create enum for all possible operations so as not to expose actual implementations to caller?
    mat = img_fun(img=self._mat, *args, **kwargs)
    self._mat = mat
    self._post_op()


  ############### cluster operations, constructs implicit phylogenetic tree

  def linkage_cluster(self, cmethod: str = "complete") -> None:
    """Clusters a MSA by some method.
    The Result is a number of clusters and their assigned sequences? Or positions, to directly identify domains?
    """
    print("\n-- OP: Linkage clustering.")
    distance_matrix = pdist(self._mat, metric='hamming')
    linkage_mat = linkage(distance_matrix, method=cmethod)
    if gc.DISPLAY: visualize_clusters(self._mat, linkage_mat)

    self.calc_consensus_custers(linkage_mat, self.nclusters)
    self._post_op()


  def calc_consensus_custers(self, linkage_mat, nclusters: int = 3):
    """Calculate clusters and a consensus sequence (average) per cluster."""
    mat = copy.deepcopy(self._mat)
    cluster_labels = fcluster(linkage_mat, nclusters, criterion='maxclust')
    consensus_list = []
    for i in range(1, nclusters + 1):
      cluster_indices = np.where(cluster_labels == i)[0]
      cluster_data = mat[cluster_indices]
      consensus_sequence = mode(cluster_data, axis=0).mode
      consensus_list.append(list(consensus_sequence))

    if gc.DISPLAY: visualize_cluster_consensuses(consensus_list)


  ############### analysis/metrics

  def calc_avg_row_col_transitions(self) -> (float, float):
    """Calculates the average number of value transitions per row as a metric for the noisiness of the MSA."""
    row_prev = self._mat[0]
    col_trans = 0
    row_trans = 0
    for row in range(self.nrows):
      row_curr = self._mat[row]
      if row > 0:
        xor_result = row_prev ^ row_curr
        col_trans += sum(xor_result)
      row_prev = row_curr

      current_state = row_curr[0]
      for col in range(1, self.ncols):
        if row_curr[col] != current_state:
          row_trans += 1
          current_state = self._mat[row, col]

    avg_row_trans = float(row_trans) / float(self.nrows)
    avg_col_trans = float(col_trans) / float(self.ncols)
    return avg_row_trans, avg_col_trans


  def calc_cr_metric(self, verbose=None) -> (float, float):
    """Calculates a column-row noisiness metric, that is the product of the average number of column value transitions
    times the average number of row value transitions."""
    rt, ct = self.calc_avg_row_col_transitions()
    prod = rt * ct
    if (verbose if verbose is not None else gc.VERBOSE):
      print(f"CR transitions. c:{ct:1.2f}, r:{rt:1.2f} -> c*r:{prod:1.2f}")
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


  ############### out

  def visualize(self, title_addon: str = "") -> None:
    """Shows a binary image of the alignment."""
    aligned_str = "" if self._filter_idx == -1 else f", ref. row {self._filter_idx}"
    addstr = "" if title_addon == "" else f": {title_addon}"
    show(self._mat, f"MSA{addstr} (1 seq/row, white=gap{aligned_str})")


  def save_to_file(self, filename: str) -> None:
    """Saves the final alignment image as well as the identified proteoform positions."""
    imgsave(self._mat, filename)
    print(f"Wrote filtered alignment image to file out/{filename}.png")
    # TODO: save annotated image to file
    # TODO: invent file format and save proteoforms + meta information

