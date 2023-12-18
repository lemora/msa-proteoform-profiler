import copy
import numpy as np

from Bio.SeqIO.FastaIO import SimpleFastaParser
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import mode

import msapp.gconst as gc
from msapp.model.mat_manipulation import remove_empty_cols, remove_seqs_from_alignment, sort_by_metric
from msapp.view.visualization import create_cluster_consensus_visualization, imgsave, visualize_clusters, show, show_hist


class MultiSeqAlignment:
    """Class that contains the state of a multiple sequence alignment that is being processed."""

    def __init__(self, filename: str = "") -> None:
        """Constructor. If filename given, loads msa, checks state and initializes binary matrix."""
        self.initialized = False
        self._filename = ""
        self._mat = None
        self.seq_names = None
        self.nrows: int = -1
        self.ncols: int = -1
        self._ridx: np.array = []  # list of row indices
        self._cidx: np.array = []  # list of col indices
        self._csort_indices = None
        self.linkage_mat = LinkageMat()

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
        msa_mat = np.array([])
        seq_names = []
        entry_length = -1
        with open(filename) as fh:
            for entry in SimpleFastaParser(handle=fh):
                seq_names.append(entry[0])
                seq = entry[1]
                rlen = len(seq)
                if entry_length == -1:
                    entry_length = rlen
                elif rlen != entry_length:
                    err_msg = "The fasta sequences are not equally long, this in not a valid alignment file.\n"
                    print(err_msg)
                    raise ValueError(err_msg)

                the_row = np.fromiter((0 if letter.isalpha() else 1 for letter in seq), dtype=np.uint8)
                if len(the_row) > 0:
                    msa_mat = np.vstack((msa_mat, the_row)) if msa_mat.size else the_row

        if not msa_mat.size or len(msa_mat) == 0:
            print("Could not initialize MSA object from the given fasta file.\n")
            return False

        self._filename: str = filename
        self.seq_names = seq_names
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
        self._mat = remove_seqs_from_alignment(self._mat, idx_list=seqs_over_three_std, cols=False)
        self.nrows = self._mat.shape[0]
        self.ncols = self._mat.shape[1]
        self.linkage_mat.mat_changed()
        self._csort_indices = None

        if gc.DISPLAY: self.visualize(rf"Removed {len(seqs_over_three_std)} seqs > 3 $\sigma$ length")
        self._post_op()

    def remove_empty_cols(self, show: bool = False):
        """Removes all columns that are empty from the matrix, meaning they only contain the value 1."""
        self._mat = remove_empty_cols(self._mat)
        self.nrows = self._mat.shape[0]
        self.ncols = self._mat.shape[1]
        if show: self.visualize(f"Removed empty columns")
        self._post_op()


    # ------ sort operations

    def sort_by_metric(self, sorting_metric=lambda row: sum(row)) -> None:
        """Sorts the rows of an MSA binary matrix by the given binary list-comparing sorting function."""
        print("\n-- OP: sorting MSA rows.")
        # # TODO: update index lists!
        self._mat = sort_by_metric(self._mat, sorting_metric)
        if gc.DISPLAY: self.visualize("Rows sorted")
        self._post_op()

    # ------ image processing

    def img_process(self, img_fun, *args, **kwargs) -> None:
        """Process the alignment by passing the current alignment matrix into the given image processing function."""
        mat = img_fun(img=self._mat, *args, **kwargs)
        self._mat = mat
        self.linkage_mat.mat_changed()
        self._post_op()

    # ------ cluster operations, constructs implicit phylogenetic tree

    def get_linkage_mat(self, cmethod: str = "complete") -> np.array:
        """Returns a linkage matrix created by means of the given distance metric."""
        return self.linkage_mat.get(self._mat, cmethod)

    def calc_consensus_clusters(self, perc_threshold: float = 0.75):
        """Calculate clusters based on relative similarity in dendrogram and a consensus sequence (average) per cluster.
        arg perc_threshold: relative dendrogram height to cut at (0 to 1), where 1 is the root (one cluster),
        0 the leaves (one cluster for every distinct row)"""
        if perc_threshold < 0 or perc_threshold > 1: raise ValueError(
            "The dendrogram cut height must be between 0 and 1.")

        linkage_mat = self.get_linkage_mat()
        if gc.DISPLAY: visualize_clusters(self._mat, linkage_mat)

        mat = self._mat
        max_val = max(linkage_mat[:, 2])
        dist_threshold = perc_threshold * max_val
        if gc.VERBOSE: print(
            f"Normalized dist threshold ({perc_threshold:1.2f} * {max_val:1.2f}): {dist_threshold:1.2f}")

        cluster_labels = fcluster(linkage_mat, t=dist_threshold, criterion='distance')
        nclusters = len(set(cluster_labels))
        if gc.VERBOSE: print(f"Clustering, n discovered: {nclusters}")

        consensus_list = np.array([])
        for i in range(1, nclusters + 1):
            cluster_indices = np.where(cluster_labels == i)[0]
            cluster_data = mat[cluster_indices]
            cseq = mode(cluster_data, axis=0).mode
            consensus_list = np.vstack((consensus_list, cseq)) if consensus_list.size else cseq

        if gc.DISPLAY: create_cluster_consensus_visualization(consensus_list)
        if nclusters == 1:
            consensus_list = [consensus_list]
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
        """Fetches n entries from the MSA file and prints certain values."""
        if self._filename == "": return
        if n < 1: return
        i = 0
        with open(self._filename) as fh:
            for entry in SimpleFastaParser(handle=fh):
                print("info:", entry[0])
                print("seq:", entry[1])
                i += 1
                if i >= n: break

    # -------- output

    def visualize(self, title_addon: str = "") -> None:
        """Shows a binary image of the alignment."""
        if not self.initialized or self._mat is None:
            print("Matrix not initialized; cannot visualize.")
            return
        addstr = "" if title_addon == "" else f": {title_addon}"
        show(self._mat, f"MSA{addstr} (1 seq/row, white=gap)")


    def get_mat(self, hide_empty_cols: bool = False, reorder_rows: bool = False) -> np.array:
        """Creates a copy of the matrix, hides empty columns and reorders rows if specified, then returns it."""
        if not self.initialized or self._mat is None:
            print("Matrix not initialized; cannot visualize.")
            return
        mat = copy.deepcopy(self._mat)
        if hide_empty_cols:
            mat = remove_empty_cols(mat)
        if reorder_rows:
            if self._csort_indices is None:
                cluster_labels = fcluster(self.get_linkage_mat(), t=0, criterion='distance')
                self._csort_indices = np.argsort(cluster_labels)
            mat = mat[self._csort_indices]
        return mat

    def save_to_file(self, filename: str) -> None:
        """Saves the final alignment image as well as the identified proteoform information."""
        imgsave(self._mat, filename)
        print(
            f"Wrote filtered alignment image to file out/{filename}.png")  # TODO: save annotated image to file  #
        # TODO: invent file format and save proteoforms + meta information


class LinkageMat:
    """Class that caches a linkage matrix of given type to avoid unnecessary recalculations."""

    def __init__(self) -> None:
        """Constructor."""
        self.dist_mat = None
        self.dmat_changed = True
        self.link_mat = None
        self.link_cmethod = ""

    def mat_changed(self):
        self.dmat_changed = True

    def get(self, mat: np.ndarray, cmethod: str):
        """Returns the distance matrix for the given mat."""
        lmat_changed = False
        if self.dmat_changed:
            if gc.VERBOSE: print("Calculating distance matrix")
            self.dst_mat = pdist(mat, metric='hamming')
            self.dmat_changed =  False
            lmat_changed = True

        if lmat_changed or self.link_cmethod != cmethod:
            if gc.VERBOSE: print("Calculating linkage matrix")
            self.link_mat = linkage(self.dst_mat, method=cmethod)
            self.link_cmethod =  cmethod

        return self.link_mat

