import cv2
import math
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster
import seaborn as sns

import msapp.gconst as gc


# ----------------- general plotting

def show_pre_post(pre, post, title: str) -> None:
    """Shows two matrix images next to each other if they are not too wide. Intended for before/after some change was
    applied to the same matrix."""
    if pre.shape[1] > 3000:
        if gc.VERBOSE: print("INFO: Image is too large to show pre/post. Just showing post version.")
        si_mat = create_colsplit_subimages_mat(post, -1)
        show_as_subimages(si_mat, title)
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


def show(msa_mat, title: str, splits: int = -1) -> Figure:
    """Fetches a matrix manipulation with a hard-coded split if the width exceeds a certain size."""
    if msa_mat.shape[1] > 3000 or splits != -1:
        si_mat = create_colsplit_subimages_mat(msa_mat, splits)
        return show_as_subimages(si_mat, title)
    return show_as_one(msa_mat, title)


def show_as_one(mat, title: str) -> Figure:
    """Shows the alignment as a binary image."""
    img = np.array(mat, dtype=np.uint8) * 255
    figure = plt.figure(figsize=(8, 8))
    figure.subplots()
    plt.imshow(img, cmap="gray"), plt.title(f"{title} [{mat.shape[0]}x{mat.shape[1]}]")
    plt.xlabel("Position in aligned sequence")
    plt.ylabel("Sequence number")
    if gc.DISPLAY: plt.show()
    return figure


def show_as_subimages(mat, title: str) -> Figure:
    """Shows the alignmet as a binary image split over several rows. Expects a subimaged mat to be passed in."""
    # scale down image: otherwise too large to properly display. mainly a cv2 problem
    # concat_img = cv2.resize(concat_img, (concat_img.shape[1] // 2, concat_img.shape[0] // 2))

    figure = plt.figure(figsize=(8, 8))
    figure.subplots()
    plt.imshow(mat, cmap="gray") #, plt.title(f"{title} [{mat.shape[0]}x{mat.shape[1]}]")
    plt.xticks([]), plt.yticks([])
    plt.xlabel("Position in aligned sequence")
    plt.ylabel("Sequence number")

    if gc.DISPLAY: plt.show()
    return figure


def create_resized_mat_visualization(mat: np.array, target_ratio: float):
    """Based on a given width to height ratio, this method calculates the optimal number of blocks the matrix needs to
    be split either by row or by column, so that if the blocks are concatenated in the other axis direction,
    they best approximate the given target aspect ratio.
    Another method is called that then does the matrix splitting and concatenation."""
    original_height, original_width = mat.shape

    # Calculate the potential splits for both row and column
    splits_by_row = math.ceil(math.sqrt(original_height * target_ratio / original_width))
    splits_by_column = math.ceil(math.sqrt(original_width / (original_height * target_ratio)))

    new_ratio_row = (original_width * splits_by_row) / (original_height / splits_by_row)
    new_ratio_column = (original_width / splits_by_column) / (original_height * splits_by_column)

    if abs(new_ratio_row - target_ratio) < abs(new_ratio_column - target_ratio):
        split_count = splits_by_row
        split_by_row =  True
    else:
        split_count = splits_by_column
        split_by_row =  False

    # if splitting by row is more optimal, transpose matrix before and after calling the splitting by col function
    if split_count == 1: return mat
    if split_by_row:
        mat = cv2.transpose(mat)
    reordered_mat = create_colsplit_subimages_mat(mat, split_count)
    if split_by_row:
        return cv2.transpose(reordered_mat)
    return reordered_mat


def create_colsplit_subimages_mat(mat, splits: int = -1):
    """Splits a matrix into several ('splits' parameter) row blocks and appends them as column blocks.
    The blocks are divided by a colourful border."""
    binary_image = np.array(mat, dtype=np.uint8) * 255
    # Split the image into equal columns
    height, width = binary_image.shape
    if splits == -1:
        splits = 3 if mat.shape[1] > 8000 else 2

    subimage_width = width // splits
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


# ----------------- clustering

def create_cluster_consensus_visualization(consensus_sequences: list[list]):
    """Shows the consensus of a number of sequence clusters based on similarity."""
    sorted_cseqs = consensus_sequences
    num_blocks = len(sorted_cseqs)
    row_len = len(sorted_cseqs[0])
    rows_per_block = int(row_len/ (3.5 * num_blocks))

    subimages = []
    for cseq in sorted_cseqs:
        for j in range(rows_per_block):
            subimages.append(cseq)

    concat_img = np.vstack(subimages)

    fig, ax = plt.subplots()
    plt.subplot(), plt.imshow(concat_img, cmap="gray")#, plt.title(f"Consensus of {num_blocks} sequence clusters")
    plt.xticks([]), plt.yticks([])
    plt.xlabel("Position in aligned sequence")
    plt.ylabel("Clusters")
    if gc.DISPLAY: plt.show()
    return fig


def visualize_clusters(mat, linkage_mat) -> None:
    """Plot the original matrix with highlighted clusters in the form of a dendrogram."""
    dendrogram(linkage_mat)
    sns.set(style="white")
    sns.clustermap(mat, row_linkage=linkage_mat, col_cluster=False, method='complete')
    plt.show()


def visualize_clustermap(mat, linkage_mat) -> Figure:
    """Creates a clustermap + dendrogram visualization of a matrix and a corresponding linkage matrix."""
    sns.set(style="white")
    g = sns.clustermap(mat, row_linkage=linkage_mat, col_cluster=False, method='complete', figsize=(4, 4))
    if gc.DISPLAY: plt.show()
    return g.fig


def visualize_dendrogram(linkage_mat, dheight: float = 0.75) -> Figure:
    """Creates a dendrogram visualization of a matrix and a corresponding linkage matrix. Colors the trees differently
    that start below a certain significant height (0<=dheight<=1) + adds a horizontal line for clarity."""
    if dheight < 0.0 or dheight > 1.0:
        raise ValueError("The dendrogram height at which to color differently needs to be between 0 and 1.")
    fig, ax = plt.subplots()
    dheight = dheight * max(linkage_mat[:,2])
    dendrogram(linkage_mat, ax=ax, color_threshold=dheight)
    plt.axhline(y=dheight, linestyle='--', color='gray', label='Desired Clusters')

    y = [0]
    ylabels = ['0.0']
    ax.set_yticks(y, labels=ylabels)
    plt.xticks([])
    plt.xlabel("Aligned sequences reordered by similarity")
    plt.ylabel("Branching by sequence similarity")
    if gc.DISPLAY: plt.show()
    return fig


def create_dendrogram_height_cluster_count_plot(linkage_mat, dheight: float = 0.75):
    """Creates a line plot based on a given linkage matrix and corresponding dendrogram, which shows the number of
    clusters per normalized height ([0: leaves, 1: root])."""
    max_distance = linkage_mat[:, 2].max()
    distances = np.linspace(0, max_distance, num=200)
    num_clusters_at_distance = [len(np.unique(fcluster(linkage_mat, t=d, criterion='distance'))) for d in distances]
    normalized_distances = distances / max_distance

    fig, ax = plt.subplots()
    plt.plot(normalized_distances, num_clusters_at_distance, color='black')
    plt.axvline(x=dheight, linestyle='--', color='gray', label='Desired Clusters')

    xticks = np.linspace(0, 10, num=11)
    xlabels = ["leaves" if i == 0 else "root" if i == 10 else f"{(i/10):1.1f}" for i in xticks]
    ax.set_xticks(xticks/10, labels=xlabels)
    plt.xlabel('Normalized distance in dendogram')
    plt.ylabel('Number of clusters')
    # plt.yscale('log')
    if gc.DISPLAY: plt.show()
    return fig


def get_empty_plot():
    """Creates an empty plot. Mainly for resetting GUI canvases."""
    fig, ax = plt.subplots()
    ax.axis("off")
    plt.xticks([]), plt.yticks([])
    return fig


# ----------------- statistics

def show_hist(counts: list, nbins: int = 50) -> None:
    fig, ax = plt.subplots()
    plt.title(f"Histogram of sequence lengths in the alignment")
    plt.grid()
    plt.hist(counts, bins=nbins, color="black")
    plt.xlabel("Sequence length")
    plt.ylabel("Count")
    if gc.DISPLAY: plt.show()
    return fig


# ----------------- saving

def imgsave(img, filename="proteoform-img") -> None:
    fig, ax = plt.subplots(nrows=1, ncols=1)
    plt.imshow(img, cmap="gray")
    plt.xticks([]), plt.yticks([])
    plt.set_tight_layout(True)
    fig.savefig(f"out/{filename}.png", bbox_inches='tight')
    plt.close(fig)
