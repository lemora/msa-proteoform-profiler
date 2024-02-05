import cv2
import math
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster
import seaborn as sns

import msapp.gconst as gc


# ----------------- general plotting


def get_empty_plot() -> Figure:
    """Returns an empty plot. Mainly for resetting GUI canvases."""
    fig, ax = plt.subplots()
    ax.axis("off")
    plt.xticks([]), plt.yticks([])
    return fig


def show_pre_post(pre, post, title: str) -> Figure:
    """Shows two matrix images next to each other if they are not too wide. Intended for before/after some change was
    applied to the same matrix."""
    if pre.shape[1] > 3000:
        if gc.VERBOSE: print("INFO: Image is too large to show pre/post. Just showing post version.")
        si_mat = create_colsplit_subimages_mat(post, -1)
        return show_as_subimages(si_mat, title)

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
    """Shows the alignment as a binary image split over several rows. Expects a subimaged mat to be passed in."""
    # scale down image: otherwise too large to properly display. mainly a cv2 problem
    # concat_img = cv2.resize(concat_img, (concat_img.shape[1] // 2, concat_img.shape[0] // 2))
    height = len(mat)
    width = len(mat[0])

    # figure = plt.figure(figsize=(8, 8))
    figure = plt.figure(figsize=(width/100, height/100))
    figure.subplots()
    plt.imshow(mat, cmap="gray")
    plt.xticks([]), plt.yticks([])
    plt.xlabel(title if title != "" else "Position in aligned sequence")
    plt.ylabel("Sequence number")

    if gc.DISPLAY: plt.show()
    return figure


def create_resized_mat_visualization(image: np.array, target_ratio: float, color_sep: bool = False) -> np.array:
    """Based on a given width to height ratio, this method calculates the optimal number of blocks the matrix needs to
    be split either by row or by column, so that if the blocks are concatenated in the other axis direction,
    they best approximate the given target aspect ratio.
    It also highlights a selected row.
    Another method is called that then does the matrix splitting and concatenation."""
    original_height = len(image)
    original_width = len(image[0])

    # Calculate the potential splits for both row and column
    splits_by_row = math.ceil(math.sqrt(original_height * target_ratio / original_width))
    splits_by_column = math.ceil(math.sqrt(original_width / (original_height * target_ratio)))

    new_ratio_row = (original_width * splits_by_row) / (original_height / splits_by_row)
    new_ratio_column = (original_width / splits_by_column) / (original_height * splits_by_column)

    diff = abs(new_ratio_row - target_ratio) - abs(new_ratio_column - target_ratio)
    if abs(diff) <= 0.2:
        split_count = 1
        should_split_by_row = False
    elif diff < 0:
        split_count = splits_by_row
        should_split_by_row = True
    else:
        split_count = splits_by_column
        should_split_by_row = False

    # if splitting by row is more optimal, transpose matrix before and after calling the splitting by col function
    # image = highlight_row(image, highlight_row_idx)
    if should_split_by_row:
        image = cv2.transpose(image)
    reordered_mat = create_colsplit_subimages_mat(image, split_count, color_sep)
    if should_split_by_row:
        return cv2.transpose(reordered_mat)
    return reordered_mat


def create_colsplit_subimages_mat(imgmat, splits: int = -1, color_sep: bool = False) -> np.array:
    """Splits an image (3 or 1 value per entry) into several ('splits' parameter) row blocks and appends them as
    column blocks. The blocks are divided by a gray border.
    An RGB image with 3 values per entry is returned."""
    if isinstance(imgmat[0][0], np.uint8):
        # create image with three channels
        imgmat = np.array(imgmat, dtype=np.uint8) * 255
        imgmat = cv2.cvtColor(imgmat, cv2.COLOR_GRAY2BGR)

    width = len(imgmat[0])
    # Split the image into equal columns
    if splits == 1:
        return imgmat
    elif splits == -1:
        splits = 3 if width > 8000 else 2

    subimage_width = width // splits
    separator = np.zeros((12, subimage_width, 3), dtype=np.uint8)
    sep_color = (150, 150, 0) if color_sep else (128, 128, 128)
    separator[:, :] = sep_color

    subimages = []
    for i in range(splits):
        start_col = i * subimage_width
        end_col = (i + 1) * subimage_width
        subimage = imgmat[:, start_col:end_col]
        subimages.append(subimage)
        if i < splits - 1:
            subimages.append(separator)

    return np.vstack(subimages)


def highlight_row(imgmat: np.array, row_idx: int) -> np.array:
    if isinstance(imgmat[0][0], np.uint8):
        # to image with three channels
        imgmat = np.array(imgmat, dtype=np.uint8) * 255
        imgmat = cv2.cvtColor(imgmat, cv2.COLOR_GRAY2BGR)
    height = len(imgmat)

    # highlight selected row
    white_highlight = (255, 0, 0)  # red
    black_highlight = (0, 89, 178)  # blue
    if 0 <= row_idx < height:
        # colored_row = np.array([white_highlight if (i == 255).all() else black_highlight for i in the_row])
        colored_row = np.where((np.all(imgmat[row_idx] == [255, 255, 255], axis=1)[:, None]), white_highlight, black_highlight)
        imgmat[row_idx] = colored_row

    return imgmat


def color_clusters(mat: np.array, cluster_labels: np.ndarray) -> np.array:
    height, width = mat.shape

    # to image with three channels
    image = np.array(mat, dtype=np.uint8) * 255
    image = cv2.cvtColor(image, cv2.COLOR_GRAY2BGR)

    if height != len(cluster_labels):
        return image

    nclusters = len(set(cluster_labels))
    palette = sns.color_palette(None, nclusters + 1)[1:]
    palette_255 = np.array([(int(color[0] * 255), int(color[1] * 255), int(color[2] * 255)) for color in palette])
    white_highlight = (255, 255, 255)
    # colour clusters
    for i in range(height):
        lbl = cluster_labels[i] - 1
        cluster_col = palette_255[lbl]
        colored_row = np.where((np.all(image[i] == [0, 0, 0], axis=1)[:, None]), cluster_col, white_highlight)
        image[i] = colored_row

    return image


# ----------------- clusterings

def create_cluster_consensus_visualization(seqs: list[list], cluster_sizes: list = None) -> Figure:
    """Shows the consensus of a number of sequence clusters based on similarity."""
    num_blocks = len(seqs)
    row_len = len(seqs[0])
    goal_height = int(row_len / 3)
    min_blocksize = int(goal_height / 20)
    if cluster_sizes is None:
        rows_per_block = int(row_len / (3 * num_blocks))
        cluster_sizes = [rows_per_block for _ in range(num_blocks)]
    else:
        curr_height = sum(cluster_sizes)
        block_factor = int(goal_height / curr_height)
        block_divisor = int(curr_height / goal_height)
        if block_factor > 1:
            cluster_sizes = [max(i * block_factor, min_blocksize) for i in cluster_sizes]
        elif block_divisor > 1:
            cluster_sizes = [max(int(i / block_divisor), min_blocksize) for i in cluster_sizes]

    height = sum(cluster_sizes)

    image = np.zeros([height, row_len, 3], dtype=np.uint8) * 255

    white = (255, 255, 255)
    black = (0, 0, 0)
    # sepcolor = (0, 89, 178)
    sepcolor = (0, 51, 102)
    pos = 0
    for i, cseq in enumerate(seqs):
        col = black if i % 2 == 0 else sepcolor
        colored_row = np.array([col if i == 0 else white for i in cseq])
        for j in range(cluster_sizes[i]):
            image[pos] = colored_row
            pos += 1

    fig, ax = plt.subplots()
    plt.subplot()
    plt.imshow(image, cmap="gray")
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
    maxval = max(linkage_mat[:, 2])
    dheight = dheight * maxval
    for i in np.arange(0.25, 1.25, 0.25):
        plt.axhline(y=i*maxval, color='#f4f4f4')
    dendrogram(linkage_mat, ax=ax, color_threshold=dheight, above_threshold_color='#ccc9ca')
    plt.axhline(y=dheight, linestyle='--', color='gray', label='Desired Clusters')

    y = [0, 0.5 * maxval, 1 * maxval]
    ylabels = ['0.0', '0.5', '1.0']
    ax.set_yticks(y, labels=ylabels)
    plt.xticks([])
    plt.xlabel("Aligned sequences reordered by similarity")
    plt.ylabel("Branching by sequence similarity")
    if gc.DISPLAY: plt.show()
    return fig


def create_dendrogram_height_cluster_count_plot(linkage_mat, dheight: float = 0.75) -> Figure:
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
    xlabels = ["leaves" if i == 0 else "root" if i == 10 else f"{(i / 10):1.1f}" for i in xticks]
    ax.set_xticks(xticks / 10, labels=xlabels)
    plt.xlabel('Normalized distance in dendogram')
    plt.ylabel('Number of clusters')
    # plt.yscale('log')
    if gc.DISPLAY: plt.show()
    return fig


# ----------------- domain visualization

def visualize_domains(domains_lists: np.array) -> Figure:
    """Gets a list of lists of tuples with start and end values between 0 and 10 (domains)."""
    msalen = max((t[1] for sublist in domains_lists for t in sublist), default=None)
    num_lines = len(domains_lists)
    image_width = 10
    image_height = num_lines
    fig, ax = plt.subplots(figsize=(image_width, image_height))
    ax.set_xlim(0, image_width)
    ax.set_ylim(0, image_height)

    normalization_factor = 10.0 / float(msalen)

    # Draw horizontal lines with filled vertical rectangles
    line_height = 0.8
    line_color = "#000000"
    for i, cluster in enumerate(domains_lists):
        line_y = i + 0.5
        ax.hlines(line_y, 0, image_width, color=line_color, linewidth=2)
        for domain in cluster:
            start = domain[0] * normalization_factor
            end = domain[1] * normalization_factor
            ax.fill_betweenx(y=[line_y - line_height / 2, line_y + line_height / 2],
                             x1=start, x2=end, color=line_color)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.xlabel("Domains")
    plt.ylabel("Proteoforms")
    if gc.DISPLAY: plt.show()
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
    # plt.set_tight_layout(True)
    fig.savefig(f"out/{filename}.png", bbox_inches='tight')
    plt.close(fig)

def save_figure(fig, filename, pad=True) -> None:
    if pad:
        fig.savefig(f'out/{filename}.png', bbox_inches='tight', dpi='figure')
    else:
        fig.savefig(f'out/{filename}.png', bbox_inches='tight', pad_inches=0, dpi='figure')
