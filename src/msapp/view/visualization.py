import copy
import cv2
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
import numpy as np
from scipy.cluster.hierarchy import dendrogram
import seaborn as sns

import msapp.gconst as gc


# ----------------- mat manipulation as plotting preparation

def filter_by_length_statistic(mat) -> None:
    """Creates a histogram of sequence lengths, then removes all longer than 3 std deviations above the median."""
    print("\n-- OP: Filtering by length statistic > 3 sigma.")
    ncols = mat.shape[1]
    seq_lengths = [ncols - sum(row) for row in mat]
    median = np.median(seq_lengths)
    std = np.std(seq_lengths)
    # if gc.DISPLAY: show_hist(seq_lengths, nbins=100)

    over_three_std = median + 3 * std
    seqs_over_three_std = np.array([idx for idx, l in enumerate(seq_lengths) if l > over_three_std])
    remove_seqs_from_alignment(mat, idx_list=seqs_over_three_std, cols=False)


def filter_by_reference(mat, idx: int) -> None:
    """Filters the alignment matrix by a given row in that MSA. Final step, visualization purposes."""
    if idx >= mat.shape[0]: raise ValueError("The index needs to be smaller than the number of rows.")

    print(f"\n-- OP: Filter by reference with index {idx}")
    refseq = mat[idx]
    cols_to_remove = np.array([i for i, val in enumerate(refseq) if val == 1])
    filtered = remove_seqs_from_alignment(mat, cols_to_remove, cols=True)
    return filtered


def remove_empty_cols(mat):
    """Removes all columns that are empty from the matrix, meaning they only contain the value 1."""
    print("\n-- OP: Removing empty columns.")
    empty_columns = np.where(np.all(mat == 1, axis=0))[0]
    filtered = remove_seqs_from_alignment(mat, empty_columns, cols=True)
    return filtered


def remove_seqs_from_alignment(mat, idx_list: np.ndarray[int], cols: bool = True) -> np.array:
    """Removes the columns (else rows) which have the given indices.
    param cols: should remove columns, else remove rows."""
    if len(idx_list) == 0: return mat
    max_idx = mat.shape[1] if cols else mat.shape[0]
    if min(idx_list) < 0 or max(idx_list) >= max_idx:
        raise ValueError("The indices o delete must be between 0 and max row or col count.")

    # mat = copy.deepcopy(mat)
    filtered = np.delete(mat, idx_list, axis=(1 if cols else 0))
    return filtered


def create_subimages_mat(mat, max_row_width: int = -1):
    binary_image = np.array(mat, dtype=np.uint8) * 255
    # Split the image into equal columns
    height, width = binary_image.shape
    if max_row_width != -1:
        print(f"width: {width}, max row width: {max_row_width}")
        splits = max(round(width / max_row_width), 2)
    else:
        splits = 3 if mat.shape[1] > 8000 else 2

    subimage_width = width // splits
    concat_img = cv2.cvtColor(binary_image, cv2.COLOR_GRAY2BGR)

    separator = np.zeros((12, subimage_width, 3), dtype=np.uint8)
    separator[:, :] = (150, 150, 0)  # colourful border

    subimages = []
    for i in range(splits):
        start_col = i * subimage_width
        end_col = (i + 1) * subimage_width
        subimage = binary_image[:, start_col:end_col]
        subimage = cv2.cvtColor(subimage, cv2.COLOR_GRAY2BGR)
        subimages.append(subimage)
        if i < splits - 1:
            subimages.append(separator)

    return np.vstack(subimages)


# ----------------- general plotting


def show_pre_post(pre, post, title: str) -> None:
    if pre.shape[1] > 3000:
        if gc.VERBOSE: print("INFO: Image is too large to show pre/post. Just showing post version.")
        show_as_subimages(post, title)
        return

    plt.subplot(121)
    plt.imshow(pre, cmap="gray"), plt.title(f'Before [{pre.shape[0]}x{pre.shape[1]}]')
    plt.xlabel("Position in aligned sequence")
    plt.ylabel("Sequence number")

    fig = plt.subplot(122)
    plt.imshow(post, cmap="gray"), plt.title(f'After [{post.shape[0]}x{post.shape[1]}]')
    plt.xlabel("Position in aligned sequence")
    plt.ylabel("Sequence number")

    plt.suptitle(f"{title}")
    plt.show()
    if gc.DISPLAY: plt.show()
    return fig


def show(msa_mat, title: str, max_row_width: int = -1) -> Figure:
    if msa_mat.shape[1] > 3000 or max_row_width != -1:
        return show_as_subimages(msa_mat, title, max_row_width)
    else:
        return show_as_one(msa_mat, title)


def show_as_one(mat, title: str) -> Figure:
    """Shows the alignment as a binary image."""
    img = np.array(mat, dtype=np.uint8) * 255

    # fig = plt.subplot()
    figure = plt.figure(figsize=(8, 8))
    figure.subplots()
    plt.imshow(img, cmap="gray"), plt.title(f"{title} [{mat.shape[0]}x{mat.shape[1]}]")
    plt.xlabel("Position in aligned sequence")
    plt.ylabel("Sequence number")
    if gc.DISPLAY: plt.show()
    return img


def show_as_subimages(mat, title: str, max_row_width: int = -1) -> Figure:
    """Shows the alignmet as a binary image split over several rows."""
    concat_img = create_subimages_mat(mat, max_row_width)
    # scale down image: otherwise too large to properly display. mainly a cv2 problem
    # concat_img = cv2.resize(concat_img, (concat_img.shape[1] // 2, concat_img.shape[0] // 2))

    figure = plt.figure(figsize=(8, 8))
    figure.subplots()
    # fig = plt.subplot()
    plt.imshow(concat_img, cmap="gray"), plt.title(f"{title} [{mat.shape[0]}x{mat.shape[1]}]")
    plt.xticks([]), plt.yticks([])
    plt.xlabel("Position in aligned sequence")
    plt.ylabel("Sequence number")

    if gc.DISPLAY: plt.show()
    return figure


# ----------------- clustering

def visualize_cluster_consensuses(consensus_sequences: list[list]):
    """Shows the consensus of a number of sequence clusters based on similarity."""
    sorted_cseqs = consensus_sequences
    num_blocks = len(sorted_cseqs)
    rows_per_block = 30
    if num_blocks > 100: rows_per_block = 1

    subimages = []
    for cseq in sorted_cseqs:
        for j in range(rows_per_block):
            subimages.append(cseq)

    concat_img = np.vstack(subimages)

    plt.subplot(), plt.imshow(concat_img, cmap="gray"), plt.title(f"Consensuses of {num_blocks} sequence clusters")
    plt.xlabel("Position in aligned sequence")
    plt.ylabel("Sequence number")
    plt.show()


def visualize_clusters(mat, linkage_mat) -> None:
    """Plot the original matrix with highlighted clusters in the form of a dendrogram."""
    dendrogram(linkage_mat)

    sns.set(style="white")
    sns.clustermap(mat, row_linkage=linkage_mat, col_cluster=False, method='complete')

    plt.show()


def create_clustermap(mat, linkage_mat) -> Figure:
    sns.set(style="white")
    g = sns.clustermap(mat, row_linkage=linkage_mat, col_cluster=False, method='complete', figsize=(4, 4))
    return g.fig


def create_dendogram(linkage_mat) -> Figure:
    fig, ax = plt.subplots()
    dendrogram(linkage_mat, ax=ax)
    plt.xticks([])
    ytick_len = len(ax.get_yticklabels())
    ldist = 1.1/ytick_len
    ytick_labels = [f"{ldist*i:1.1f}" for i in range(ytick_len)]
    ax.set_yticklabels(ytick_labels)
    plt.xlabel("Aligned sequences reordered by similarity")
    plt.ylabel("Branching by similarity")
    return fig

# ----------------- statistics

def show_hist(counts: list, nbins: int = 50) -> None:
    plt.subplot(), plt.title(f"Histogram of sequence lengths in the alignment")
    plt.grid()
    plt.hist(counts, bins=nbins, color="black")
    plt.xlabel("Sequence length")
    plt.ylabel("Count")
    plt.show()


# ----------------- saving

def imgsave(img, filename="proteoform-img") -> None:
    fig, ax = plt.subplots(nrows=1, ncols=1)
    plt.imshow(img, cmap="gray")
    plt.xticks([]), plt.yticks([])
    fig.savefig(f"out/{filename}.png", bbox_inches='tight')
    plt.close(fig)
