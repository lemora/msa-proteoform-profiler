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
    if pre.shape[1] > 3000:
        if gc.VERBOSE: print("INFO: Image is too large to show pre/post. Just showing post version.")
        si_mat = create_subimages_mat(post, -1)
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
    if msa_mat.shape[1] > 3000 or splits != -1:
        si_mat = create_subimages_mat(msa_mat, splits)
        return show_as_subimages(si_mat, title)
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
    return figure


def show_as_subimages(mat, title: str, splits: int = -1) -> Figure:
    """Shows the alignmet as a binary image split over several rows. Expects a subimaged mat to be passed in."""
    # concat_img = create_subimages_mat(mat, splits)
    concat_img = mat
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


def create_subimages_mat(mat, splits: int = -1):
    """Splits a matrix into several row blocks and appends them as columns with a separator.
    For visualoization."""
    # TODO: instead use hight to width ratio for mat formatting!
    # TODO: create row or col subimages, depending on dominating axis!
    binary_image = np.array(mat, dtype=np.uint8) * 255
    # Split the image into equal columns
    height, width = binary_image.shape
    if splits == -1:
        splits = 3 if mat.shape[1] > 8000 else 2

    subimage_width = width // splits
    # concat_img = cv2.cvtColor(binary_image, cv2.COLOR_GRAY2BGR)

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


def create_resized_mat_visualization(mat: np.array, target_ratio: float):
    original_height, original_width = mat.shape
    original_ratio = original_width / original_height

    # Calculate the potential splits for both row and column
    splits_by_row = math.ceil(math.sqrt(original_height * target_ratio / original_width))
    splits_by_column = math.ceil(math.sqrt(original_width / (original_height * target_ratio)))

    # Calculate the resulting aspect ratios for both cases
    new_ratio_row = (original_width * splits_by_row) / (original_height / splits_by_row)
    new_ratio_column = (original_width / splits_by_column) / (original_height * splits_by_column)

    # Choose the splits that brings us closer to the target ratio
    if abs(new_ratio_row - target_ratio) < abs(new_ratio_column - target_ratio):
        split_count = splits_by_row
        split_by_row =  True
    else:
        split_count = splits_by_column
        split_by_row =  False

    # if splitting by row is more optimal, transpose matrix pre and post calling the splitting by col function
    if split_count == 1: return mat
    if split_by_row:
        mat = cv2.transpose(mat)
    subimaged = create_subimages_mat(mat, split_count)
    if split_by_row:
        return cv2.transpose(subimaged)
    return subimaged


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
    # desired_num_clusters = 0.75
    # plt.axhline(y=desired_num_clusters, linestyle='--', color='red', label='Desired Clusters')

    yticks_curr = ax.get_yticklabels()
    ytick_len = len(yticks_curr)
    ldist = 1/(ytick_len-2)
    # ytick_labels = [f"{ldist*i:1.1f}" for i in range(ytick_len)]
    # ytick_labels = [f"{ldist*i:1.1f}" for i in range(ytick_len-1)]
    # ax.set_yticklabels(ytick_labels)
    y = [0]
    ylabels = ['0.0']
    ax.set_yticks(y, labels=ylabels)
    plt.xticks([])
    # plt.yticks([])
    plt.xlabel("Aligned sequences reordered by similarity")
    plt.ylabel("Branching by similarity")
    return fig


def create_dendogram_height_cluster_count_plot(linkage_mat):
    max_val = max(linkage_mat[:, 2])
    xvals = np.arange(0.0, 1.1, 0.05)
    d = max_val/len(xvals)
    yvals = [len(set(fcluster(linkage_mat, t=(d*i), criterion='distance'))) for i in xvals]

    fig, ax = plt.subplots()
    ax.plot(xvals, yvals, linewidth=2.0)
    plt.xlabel("Dendogram tree height, from leaves to root")
    plt.ylabel("Number of clusters")
    return fig

    # dendrogram_data = dendrogram(linkage_mat, no_plot=True)
    # heights = dendrogram_data['dcoord'][0]
    # clusters_at_height = [len(np.unique(fcluster(linkage_mat, t=h, criterion='distance'))) for h in heights]

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
