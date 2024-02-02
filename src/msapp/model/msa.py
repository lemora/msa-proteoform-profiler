import copy
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import mode

import msapp.gconst as gc
from msapp.model.mat_manipulation import clear_seqs_in_alignment, cross_convolve, gaussian_blur, remove_empty_cols
from msapp.view.visualization import imgsave


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
        self.linkage_mat = LinkageMat()
        self.seq_indexer = SequenceIndexer()
        self.cr_score = None

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
                    err_msg = "The fasta sequences are not equally long, this in not a valid alignment file."
                    raise ValueError(err_msg)

                the_row = np.fromiter((0 if letter.isalpha() else 1 for letter in seq), dtype=np.uint8)
                if len(the_row) > 0:
                    msa_mat = np.vstack((msa_mat, the_row)) if msa_mat.size else the_row

        if not msa_mat.size or len(msa_mat) == 0:
            err_msg = "Could not initialize MSA object from the given fasta file."
            raise ValueError(err_msg)

        self._filename: str = filename
        self.seq_names = seq_names
        self.seq_indexer.init(seq_names)
        msa_mat = np.array(msa_mat, dtype=np.uint8)
        self._mat = msa_mat
        self.nrows = msa_mat.shape[0]
        self.ncols = msa_mat.shape[1]
        self._ridx = np.array([i for i in range(self.nrows)])
        self._cidx = np.array([i for i in range(self.ncols)])
        self.initialized = True

        if gc.VERBOSE:
            print(f"Successfully loaded MSA from file '{filename}' (sequences: {self.nrows}, length: {self.ncols})\n")
            print("First entry:")
            self.print_msa(1)
        return True

    def get_seq_indexer(self):
        return self.seq_indexer

    # ------ filter operations

    def filter_by_length_statistic(self) -> None:
        """Creates a histogram of sequence lengths, then removes all longer than 3 std deviations above the median."""
        if gc.VERBOSE: print("-- OP: Filtering by length statistic > 3 sigma.")

        seq_lengths = [self.ncols - sum(row) for row in self._mat]
        median = np.median(seq_lengths)
        std = np.std(seq_lengths)
        # if gc.DISPLAY: show_hist(seq_lengths, nbins=100)

        over_three_std = median + 3 * std
        seqs_over_three_std = np.array([idx for idx, l in enumerate(seq_lengths) if l > over_three_std])
        if gc.VERBOSE and len(seqs_over_three_std) > 0:
            print(f"Removed {len(seqs_over_three_std)} sequences > 3 sigma.")

        self._mat = clear_seqs_in_alignment(self._mat, idx_list=seqs_over_three_std)
        self.seq_indexer.indices_dendro_changed()
        self.linkage_mat.mat_changed()

    def remove_isolated_connections(self):
        if gc.VERBOSE: print("-- OP: Removing isolated connections.")

        black_pixels_count = np.sum(self._mat == 0, axis=0)
        columns_to_change = np.where(black_pixels_count <= 1)[0]
        self._mat[:, columns_to_change] = 1

    # ------ image processing

    def img_process(self, img_fun, *args, **kwargs) -> None:
        """Process the alignment by passing the current alignment matrix into the given image processing function."""
        mat = img_fun(img=self._mat, *args, **kwargs)
        self._mat = mat
        self.seq_indexer.indices_dendro_changed()
        self.linkage_mat.mat_changed()

    def run_filtering_pipeline(self, filter_type: str = "standard"):
        filter_type = filter_type.lower()
        filter_type = "standard" if filter_type not in ["mild", "standard", "aggressive"] else filter_type
        if gc.VERBOSE: print(f"Running {filter_type} filtering pipeline.")

        self.filter_by_length_statistic()
        self.remove_isolated_connections()
        self.img_process(gaussian_blur, ksize=5)

        if filter_type.lower() == "standard" or filter_type.lower() == "aggressive":
            self.img_process(cross_convolve, col_size=5)
            self.img_process(cross_convolve, row_size=7, col_size=17)

        if filter_type.lower() == "aggressive":
            vsize = self.nrows // 10
            vsize = vsize if vsize % 2 == 1 else vsize + 1
            self.img_process(cross_convolve, col_size=vsize, row_size=3)

    # ------ cluster operations, constructs implicit phylogenetic tree

    def get_linkage_mat(self, cmethod: str = "complete") -> np.array:
        """Returns a linkage matrix created by means of the given distance metric."""
        return self.linkage_mat.get(self._mat, cmethod)

    def get_cluster_labels(self, perc_threshold: float = 0.75, dendro_ordered = False):
        linkage_mat = self.get_linkage_mat()
        max_val = max(linkage_mat[:, 2])
        dist_threshold = perc_threshold * max_val
        cluster_labels = fcluster(linkage_mat, t=dist_threshold, criterion='distance')
        if dendro_ordered:
            indices_dendro = self.get_seq_indexer().get_indices_dendro()
            cluster_labels = [cluster_labels[i] for i in indices_dendro]
        return cluster_labels

    def retrieve_domains_via_dendrogram(self, perc_threshold: float = 0.75):
        """Calculates and returns domains that are recognizable in the MSA.
                Returns a list of lists, where each list corresponds to a cluster and
                each list contains tuples, one per domain: (start, end) between 0.0 and 10.0."""
        cluster_labels = self.get_cluster_labels(perc_threshold)
        nclusters = len(set(cluster_labels))

        consensus_list = np.array([])
        cluster_sizes = []
        for i in range(nclusters, 0, -1):
            cluster_indices = np.where(cluster_labels == i)[0]
            cluster_sizes.append(len(cluster_indices))
            #cluster_data = self._mat[cluster_indices]
            mat = self.get_mat(hide_empty_cols=True)
            cluster_data = mat[cluster_indices]
            cseq = mode(cluster_data, axis=0).mode
            consensus_list = np.vstack((consensus_list, cseq)) if consensus_list.size else cseq

        if nclusters == 1:
            consensus_list = [consensus_list]
        domains = self.calculate_black_regions(consensus_list)
        # print(f"calculated domains: {domains}")
        return domains

    def calculate_black_regions(self, consensus_list):
        black_regions_list = []
        min_dist = self.ncols / 700
        min_width = 6

        for consensus_seq in consensus_list:
            black_regions = []
            current_start = None

            for i, value in enumerate(consensus_seq):
                if value == 0:
                    if current_start is None:
                        current_start = i
                elif current_start is not None:
                    # white found
                    current_end = i - 1
                    merged = False
                    if len(black_regions) > 0:
                        prev = black_regions[-1]
                        if (current_start - prev[1] < min_dist \
                            or (prev[1] - prev[0] < min_width and current_start - prev[1] < min_dist*2)):
                            black_regions[-1] = (prev[0], current_end)
                            merged = True
                    if not merged:
                        # if current_end - current_start >= min_width:
                        black_regions.append((current_start, current_end))
                    current_start = None
            if current_start is not None:
                current_end = len(consensus_seq) - 1
                black_regions.append((current_start, current_end))
            black_regions_list.append(black_regions)

        return black_regions_list

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

        if verbose if verbose is not None else gc.VERBOSE:
            print(f"CR transitions. c:{ct:1.2f}, r:{rt:1.2f} -> c*r:{prod:1.2f}")
        self.cr_score = prod
        return rt, ct

    # -------- getter

    def get_mat(self, hide_empty_cols: bool = False, reorder_rows: bool = False) -> np.array:
        """Creates a copy of the matrix, hides empty columns and reorders rows if specified, then returns it."""
        if not self.initialized or self._mat is None:
            print("Matrix not initialized; cannot visualize.")
            return
        mat = copy.deepcopy(self._mat)
        if hide_empty_cols:
            mat = remove_empty_cols(mat)
        if reorder_rows:
            self._refresh_dendro_indices()
            indices_dendro = self.seq_indexer.get_indices_dendro()
            mat = mat[indices_dendro]
        return mat

    def _refresh_dendro_indices(self):
        if self.seq_indexer.get_indices_dendro() is None:
            cluster_labels = fcluster(self.get_linkage_mat(), t=0, criterion='distance')
            dendro_indices = np.argsort(cluster_labels)
            self.seq_indexer.set_indices_dendro(dendro_indices)

    # -------- output

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

    def save_to_file(self, filename: str) -> None:
        """Saves the final alignment image as well as the identified proteoform information."""
        mat = self.get_mat(hide_empty_cols=True, reorder_rows=True)
        imgsave(mat, filename)
        print(f"Wrote filtered alignment image to file out/{filename}.png")
        # TODO: save annotated image to file
        # TODO: invent file format and save proteoforms + meta information


# ------------------------------------------------------------------------

class LinkageMat:
    """Class that caches a linkage matrix of given type to avoid unnecessary recalculations."""

    def __init__(self) -> None:
        self.dist_mat = None
        self.link_mat = None
        self.link_cmethod = ""

    def mat_changed(self):
        self.dist_mat = None
        self.link_mat = None

    def _update_if_needed(self, mat: np.ndarray, cmethod: str = "complete") -> None:
        if self.dist_mat is None:
            self.dist_mat = pdist(mat, metric='hamming')
            self.link_mat = None

        if self.link_mat is None or self.link_cmethod != cmethod:
            self.link_mat = linkage(self.dist_mat, method=cmethod, metric="hamming")
            self.link_cmethod = cmethod

    def get(self, mat: np.ndarray, cmethod: str):
        """Returns the distance matrix for the given mat."""
        self._update_if_needed(mat, cmethod)
        return self.link_mat


# ------------------------------------------------------------------------

class SequenceIndexer:
    """Stores a list of sequence names and their order. Provides access methods via index, id and name."""

    def __init__(self):
        self.num_entries = 0
        self.seqid_dict = None  # key: seq_id, value: tuple (text:str, idx_mat:int, idx_dendro:int)
        self.seqid_list = None
        self.indices_dendro = None

    def init(self, id_name_list: np.array) -> None:
        """Initialization.
        param id_name_list: list of strings, each consisting of at least two words (seq id and description)."""
        self.num_entries = len(id_name_list)
        seqid_descr_list = [(s.split()[0], ' '.join(s.split()[1:])) for s in id_name_list]
        self.seqid_list = [entry[0] for entry in seqid_descr_list]
        self.seqid_dict = {sid: (text, idx, -1) for idx, (sid, text) in enumerate(seqid_descr_list)}

    def get_seq_infos_containing_string(self, search_string) -> list:
        """Returns a list of (seq_id:str,description:str) tuples for sequences where one of the two fields
        contains the given search_string (case agnostic)."""
        search_string = search_string.lower()
        filtered_tuples = [(sid, val[0]) for sid, val in self.seqid_dict.items() if
                           search_string in sid.lower() or search_string in val[0].lower()]
        return filtered_tuples

    def get_infos_for_seq_id(self, seq_id) -> tuple:
        """For the given seq_id, returns the following tuple: (description:str, idx_mat:int, idx_dendro:int)"""
        if not seq_id in self.seqid_dict:
            return None
        info = self.seqid_dict[seq_id]
        return info
        # return (seq_id, info[0])

    # --- seq id to index

    def get_matidx_from_seqid(self, seq_id: int) -> int:
        """Returns the row index in the msa matrix for the given sequence ID."""
        if self.seqid_dict is None:
            return None
        return self.seqid_dict[seq_id][1]

    def get_dendroidx_from_seqid(self, seq_id: int) -> int:
        """Returns the row index in the dendro-sorted matrix for the given sequence ID."""
        if self.seqid_dict is None:
            return None
        return self.seqid_dict[seq_id][2]

    def get_ith_seqid_mat(self, i: int) -> str:
        """Returns the ith sequence ID in from the standard matrix sorting indices"""
        if i is None or i == -1:
            return self.seqid_list[0]
        return self.seqid_list[i]

    def get_ith_seqid_dendro(self, i: int) -> str:
        """Returns the ith sequence ID in from the by-dendrogram matrix sorting indices"""
        if i >= self.num_entries:
            return None
        seq_id = [seq_id for seq_id, value in self.seqid_dict.items() if value[2] == i][0]
        return seq_id

    # --- index to seq id

    def get_seqid_from_matidx(self, idx: int):
        return self.seqid_list[idx]

    def indices_dendro_changed(self):
        self.indices_dendro = None

    # --- getter

    def get_indices_dendro(self):
        return self.indices_dendro

    # --- setter

    def set_indices_dendro(self, dendro_indices):
        self.indices_dendro = dendro_indices
        for dendro_idx, mat_idx in enumerate(dendro_indices):
            seq_id = self.get_seqid_from_matidx(mat_idx)
            seq_tup = self.seqid_dict[seq_id]
            self.seqid_dict[seq_id] = (seq_tup[0], seq_tup[1], dendro_idx)
