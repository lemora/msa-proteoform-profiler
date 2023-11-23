import copy
from matplotlib.figure import Figure
import numpy as np
from pysam import FastxFile
import re
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import mode

import msapp.gconst as gc
from msapp.model.mat_manipulation import remove_empty_cols, sort_by_metric
from msapp.view.visualization import imgsave, visualize_clusters, create_cluster_consensus_visualization, show, show_hist


class MultiSeqAlignment:
    """Class that contains the state of a multiple sequence alignment that is being processed."""

    def __init__(self, filename: str = "") -> None:
        """Constructor. If filename given, loads msa, checks state and initializes binary matrix."""
        self.initialized = False
        self._filename = ""
        self._mat = None
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
        if binary_mat is None: raise ValueError("The MSA matrix cannot be None.")

        self._mat = binary_mat
        self.nrows: int = binary_mat.shape[0]
        self.ncols: int = binary_mat.shape[1]
        self._ridx = np.array([i for i in range(self.nrows)])
        self._cidx = np.array([i for i in range(self.ncols)])
        self.initialized = True
        print(f"Successfully loaded MSA ({self.nrows} sequences of length {self.ncols})")
        return True

    def init_from_file(self, filename: str) -> bool:
        """Initializes a MultiSeqAlignment object from a msa fasta file."""
        if self.initialized: raise ValueError("The MSA object has already been initialized.\n")

        # Create a binary matrix from the multiple sequence alignment file
        msa_mat = []
        min_rlen = -1
        max_rlen = -1
        with FastxFile(filename) as fh:
            for entry in fh:
                seq = entry.sequence
                if not re.match("^[A-Za-z-]*$", seq):
                    err_msg = f"Could not load MSA; a line does not match the expected pattern.\n"
                    print(err_msg)
                    raise ValueError(err_msg)
                the_row = []
                rlen = len(seq)
                if min_rlen == -1 or rlen < min_rlen:
                    min_rlen = rlen
                if max_rlen == -1 or rlen > max_rlen:
                    max_rlen = rlen
                for letter in seq:
                    the_row.append(0 if letter.isalpha() else 1)  # the_row = np.array(the_row, dtype=np.uint8)
                msa_mat.append(the_row)
            fh.close()

        if msa_mat == None or len(msa_mat) == 0:
            print("Could not initialize MSA object from the given fasta file.\n")
            return False

        print(f"MSA mat. min len: {min_rlen}, max len: {max_rlen}")
        if min_rlen != max_rlen:
            err_msg = "The fasta sequences are not equally long, this in not a valid alignment file.\n"
            print(err_msg)
            raise ValueError(err_msg)

        self._filename: str = filename
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
        return True

    # ------ filter operations

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
        # TODO: now redundant. remove?
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
        # TODO: now redundant. remove?
        print("\n-- OP: Removing empty columns.")
        empty_columns = np.where(np.all(self._mat == 1, axis=0))[0]
        self.remove_seqs_from_alignment(empty_columns, cols=True)
        if show: self.visualize(f"Removed {len(empty_columns)} empty columns")
        self._post_op()

    def remove_seqs_from_alignment(self, idx_list: np.ndarray[int], cols: bool = True) -> None:
        """Removes the columns (else rows) which have the given indices.
        param cols: should remove columns, else remove rows."""
        # TODO: now redundant. remove?
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

    # ------ sort operations

    def sort_by_metric(self, sorting_metric=lambda row: sum(row)) -> None:
        """Sorts a MSA binary matrix by the given function that works on a binary list."""
        # TODO: now redundant. remove?
        print("\n-- OP: sorting MSA rows.")
        mat = copy.deepcopy(self._mat)
        sorting_metrics = np.apply_along_axis(sorting_metric, axis=1, arr=mat)

        # TODO: update index lists!
        sorted_indices = np.argsort(sorting_metrics)
        sorted_matrix = mat[sorted_indices]

        self._mat = sorted_matrix
        if gc.DISPLAY: self.visualize("Rows sorted")
        self._post_op()

    # ------ image processing

    def img_process(self, img_fun, *args, **kwargs) -> None:
        """Process the alignment by passing the current alignment matrix/image into the given image processing
        function."""
        mat = img_fun(img=self._mat, *args, **kwargs)
        self._mat = mat
        self._post_op()

    # ------ cluster operations, constructs implicit phylogenetic tree

    def get_linkage_mat(self, cmethod: str = "complete") -> np.array:
        distance_matrix = pdist(self._mat, metric='hamming')
        linkage_mat = linkage(distance_matrix, method=cmethod)
        return linkage_mat

    def linkage_cluster(self, dendogram_cut_height: float = 0.5) -> None:
        """Clusters a MSA by some method.
        The Result is a number of clusters and their assigned sequences? Or positions, to directly identify domains?
        """
        print("\n-- OP: Linkage clustering.")
        linkage_mat = self.get_linkage_mat()
        if gc.DISPLAY: visualize_clusters(self._mat, linkage_mat)

        self.calc_consensus_clusters(linkage_mat, perc_threshold=dendogram_cut_height)
        self._post_op()

    def calc_consensus_clusters(self, linkage_mat: np.array, perc_threshold: float):
        """Calculate clusters based on relative similarity in dendrogram and a consensus sequence (average) per cluster.
        arg perc_threshold: relative dendrogram height to cut at (0 to 1), where 1 is the root (one cluster),
        0 the leaves (one cluster for every distinct row)"""
        if perc_threshold < 0 or perc_threshold > 1: raise ValueError(
            "The dendrogram cut height must be between 0 and 1.")

        mat = copy.deepcopy(self._mat)
        mat = remove_empty_cols(mat)
        max_val = max(linkage_mat[:, 2])
        dist_threshold = perc_threshold * max_val
        if gc.VERBOSE: print(
            f"Normalized dist threshold ({perc_threshold:1.2f} * {max_val:1.2f}): {dist_threshold:1.2f}")

        cluster_labels = fcluster(linkage_mat, t=dist_threshold, criterion='distance')
        nclusters = len(set(cluster_labels))
        if gc.VERBOSE: print(f"Clustering, n discovered: {nclusters}")

        consensus_list = []
        for i in range(1, nclusters + 1):
            cluster_indices = np.where(cluster_labels == i)[0]
            cluster_data = mat[cluster_indices]
            consensus_sequence = mode(cluster_data, axis=0).mode
            consensus_list.append(list(consensus_sequence))

        if gc.DISPLAY: create_cluster_consensus_visualization(consensus_list)
        return consensus_list

    # -------- analysis/metrics

    def calc_cr_metric(self, verbose=None) -> (float, float):
        """Calculates a column-row noisiness metric, that is the product of the average number of column value
        transitions times the average number of row value transitions."""
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

        rt = float(row_trans) / float(self.nrows)
        ct = float(col_trans) / float(self.ncols)
        prod = rt * ct

        if (verbose if verbose is not None else gc.VERBOSE):
            print(f"CR transitions. c:{ct:1.2f}, r:{rt:1.2f} -> c*r:{prod:1.2f}")
        self.cr_score = prod
        return rt, ct

    def _post_op(self) -> None:
        """Operations to perform at the end of an MSA processing step."""
        if gc.VERBOSE: self.calc_cr_metric(verbose=True)

    # -------- helper/debug

    def print_msa(self, n) -> None:
        if self._filename == "": return
        if n < 1: return
        i = 0
        with FastxFile(self._filename) as fh:
            for entry in fh:
                print(f"name: {entry.name}")
                print(f"seq: {entry.sequence}")
                print(f"comment: {entry.comment}")
                print(f"quality: {entry.quality}\n")
                i += 1
                if i >= n: break
            fh.close()

    # -------- output

    def visualize(self, title_addon: str = "") -> None:
        """Shows a binary image of the alignment."""
        if not self.initialized or self._mat is None:
            print("Matrix not initialized; cannot visualize.")
            return
        aligned_str = "" if self._filter_idx == -1 else f", ref. row {self._filter_idx}"
        addstr = "" if title_addon == "" else f": {title_addon}"
        show(self._mat, f"MSA{addstr} (1 seq/row, white=gap{aligned_str})")

    def get_mat_visualization(self, hide_empty_cols: bool = False, reorder_rows: bool = False,
                              max_row_width: int = -1) -> Figure:
        if not self.initialized or self._mat is None:
            print("Matrix not initialized; cannot visualize.")
            return
        mat = self._mat
        if hide_empty_cols:
            mat = copy.deepcopy(self._mat)
            mat = remove_empty_cols(mat)
        if reorder_rows:
            mat = sort_by_metric(mat)

        return show(mat, "", max_row_width)

    def get_mat(self, hide_empty_cols: bool = False, reorder_rows: bool = False) -> np.array:
        if not self.initialized or self._mat is None:
            print("Matrix not initialized; cannot visualize.")
            return
        mat = copy.deepcopy(self._mat)
        if hide_empty_cols:
            mat = remove_empty_cols(mat)
        if reorder_rows:
            print("Reordering rows..")
            mat = sort_by_metric(mat)
        return mat

    def save_to_file(self, filename: str) -> None:
        """Saves the final alignment image as well as the identified proteoform positions."""
        imgsave(self._mat, filename)
        print(
            f"Wrote filtered alignment image to file out/{filename}.png")  # TODO: save annotated image to file  #
        # TODO: invent file format and save proteoforms + meta information
