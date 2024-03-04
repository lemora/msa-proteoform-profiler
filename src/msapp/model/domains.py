import numpy as np
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.neighbors import KNeighborsClassifier

from msapp.model.mat_manipulation import remove_empty_cols


def calculate_domains(mat, cluster_labels: list, mode: str = "Quick", hide_empty_cols=True) -> list:
    """
    From a matrix and a list of labels for each row in that matrix, calculate a consensus for each label value from
    all associated sequences.
    :param mat: a binary matrix (list of lists with values [0,1])
    :param cluster_labels: a list of labels, one per row of mat
    :param mode: how thoroughly to calculate a consensus, one of ["quick", "thorough"]
    :return: list of lists of tuples, with each tuple  (start, end) denoting a black region (domain) in that consensus
    """
    mode = mode.lower()
    nclusters = len(set(cluster_labels))

    consensus_lists = np.array([])
    cluster_sizes = []
    for i in range(nclusters, 0, -1):
        cluster_indices = np.where(cluster_labels == i)[0]
        cluster_sizes.append(len(cluster_indices))
        cluster_data = mat[cluster_indices]

        if mode == "quick":
            cseq = create_average_sequence_median(cluster_data)
        else:
            cseq = create_average_sequence_kmeans(cluster_data)

        cseq = filter_average_sequence(cseq, min_width=6, min_gap_len=3)
        consensus_lists = np.vstack((consensus_lists, cseq)) if consensus_lists.size else cseq

    if nclusters == 1:
        consensus_lists = [consensus_lists]

    if hide_empty_cols:
        consensus_lists = remove_empty_cols(consensus_lists)
    return mat_to_tuples(consensus_lists)


def filter_average_sequence(seq, min_width: int = 3, min_gap_len: int = 2):
    """
    Gets a list with binary values. Fuses black sections that are less than max_gap apart and removes those black
    regions that are smaller than min_width.
    """
    filtered_seq = seq.copy()
    black_stretch_length = 0
    white_stretch_length = 0
    slen = len(seq)
    for i in range(slen):
        if seq[i] == 0:
            black_stretch_length += 1
            if i == slen - 1 and black_stretch_length < min_width:
                filtered_seq[i - black_stretch_length:i+1] = 1
            elif white_stretch_length < min_gap_len and i > white_stretch_length:
                filtered_seq[i - white_stretch_length:i] = 0
                black_stretch_length = get_black_region_size_before_idx(filtered_seq, i)
            white_stretch_length = 0
        else:
            white_stretch_length += 1
            if black_stretch_length < min_width:
                filtered_seq[i - black_stretch_length:i] = 1
            black_stretch_length = 0
    return np.array(filtered_seq)

def get_black_region_size_before_idx(seq, index):
    count = 0
    for i in range(index, -1, -1):
        if seq[i] == 0:
            count += 1
        else:
            break
    return count

# ----------------- to tuples

def mat_to_tuples(consensus_lists: list):
    """For a binary matrix, turns it into a list of lists, one per row, containing black regions as tuples:
    (start, end)."""
    tuples_lists = []
    for row in consensus_lists:
        tuples_lists.append(blist_to_region_tuples(row))
    return tuples_lists


def blist_to_region_tuples(binary_list: list):
    """Turns a list with binary values into a list of tuples denoting black regions (start, end)."""
    black_regions = []
    start = None

    for i, value in enumerate(binary_list):
        if value == 0 and start is None:
            start = i
        elif value == 1 and start is not None:
            black_regions.append((start, i - 1))
            start = None

    if start is not None:
        black_regions.append((start, len(binary_list) - 1))

    return black_regions


# ----------------- consensus methods

def create_average_sequence_median(binary_sequences: list):
    """Gets a list of binary-valued lists, creates one consensus list by taking the median."""
    binary_sequences = np.array(binary_sequences)
    average_sequence = np.median(binary_sequences, axis=0)
    return np.round(average_sequence).astype(int)


def create_average_sequence_kmeans(binary_sequences: list):
    """Gets a list of binary-valued lists, creates one consensus list by averaging the rounded centroids calculated
    by kmeans."""
    if len(binary_sequences) < 3:
        return create_average_sequence_median(binary_sequences)

    binary_sequences = np.array(binary_sequences)
    kmeans = KMeans(n_clusters=4, random_state=0)
    kmeans.fit(binary_sequences)
    centroids = kmeans.cluster_centers_
    average_sequence = np.round(centroids).astype(int)
    consensus_sequence = np.mean(average_sequence, axis=0)
    consensus_sequence = np.round(consensus_sequence).astype(int)
    return consensus_sequence


def create_average_sequence_agg(binary_sequences: list, n_clusters: int = 2):
    """Gets a list of binary-valued lists, creates one consensus list by averaging the cluster averages calculated
    by agglomerative clustering."""
    if len(binary_sequences) < n_clusters:
        return create_average_sequence_median(binary_sequences)

    binary_sequences = np.array(binary_sequences)
    clustering = AgglomerativeClustering(n_clusters, metric="hamming", linkage="complete")
    cluster_labels = clustering.fit_predict(binary_sequences)

    cluster_counts = np.zeros((n_clusters, binary_sequences.shape[1]), dtype=float)
    cluster_weights = np.zeros((n_clusters, binary_sequences.shape[1]), dtype=float)
    for i in range(binary_sequences.shape[0]):
        cluster_counts[cluster_labels[i]] += binary_sequences[i]
        cluster_weights[cluster_labels[i]] += 1

    cluster_averages = cluster_counts / (cluster_weights + 1e-10)
    average_sequence = np.round(np.median(cluster_averages, axis=0)).astype(int)

    return average_sequence


def create_average_sequence_knn(binary_sequences: list, n_neighbors: int = 3):
    if len(binary_sequences) < n_neighbors:
        return create_average_sequence_median(binary_sequences)

    sequences_array = np.array(binary_sequences)
    average_sequence = np.zeros_like(sequences_array[0])

    for i in range(sequences_array.shape[1]):
        X = np.arange(len(binary_sequences))[:, np.newaxis]
        y = sequences_array[:, i]
        knn = KNeighborsClassifier(n_neighbors=n_neighbors, weights='distance')
        knn.fit(X, y)
        average_value = knn.predict([[len(binary_sequences) // 2]])
        average_sequence[i] = average_value

    return average_sequence
