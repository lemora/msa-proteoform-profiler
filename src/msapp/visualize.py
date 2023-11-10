import cv2
from matplotlib import pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster
from scipy.stats import mode
import seaborn as sns

import msapp.gconst as gc


# ----------------- general plotting


def show_pre_post(pre, post, title: str) -> None:
  if pre.shape[1] > 3000:
    if gc.VERBOSE: print("INFO: Image is too large to show pre/post. Just showing post version.")
    show_as_subimages(post, title)
    return

  plt.subplot(121), plt.imshow(pre, cmap="gray"), plt.title(f'Before [{pre.shape[0]}x{pre.shape[1]}]')
  plt.xlabel("Position in aligned sequence")
  plt.ylabel("Sequence number")

  plt.subplot(122), plt.imshow(post, cmap="gray"), plt.title(f'After [{post.shape[0]}x{post.shape[1]}]')
  plt.xlabel("Position in aligned sequence")
  plt.ylabel("Sequence number")

  mng = plt.get_current_fig_manager()
  mng.resize(*mng.window.maxsize())
  plt.suptitle(f"{title}")
  plt.show()


def show(msa_mat, title: str) -> None:
  if msa_mat.shape[1] > 3000:
    show_as_subimages(msa_mat, title)
  else:
    show_as_one(msa_mat, title)


def show_as_one(mat, title: str) -> None:
  """Shows the alignment as a binary image."""
  img = np.array(mat, dtype=np.uint8) * 255

  plt.subplot(), plt.imshow(img, cmap="gray"), plt.title(f"{title} [{mat.shape[0]}x{mat.shape[1]}]")
  plt.xlabel("Position in aligned sequence")
  plt.ylabel("Sequence number")
  mng = plt.get_current_fig_manager()
  mng.resize(*mng.window.maxsize())
  plt.show()


def show_as_subimages(mat, title: str) -> None:
  """Shows the alignmet as a binary image split over several rows."""
  binary_image = np.array(mat, dtype=np.uint8) * 255

  # Split the image into equal columns
  splits = 3
  height, width = binary_image.shape
  subimage_width = width // splits
  concat_img = cv2.cvtColor(binary_image, cv2.COLOR_GRAY2BGR)

  separator = np.zeros((12, subimage_width, 3), dtype=np.uint8)
  separator[:, :] = (150, 150, 0) # colourful border

  subimages = []
  for i in range(splits):
    start_col = i * subimage_width
    end_col = (i + 1) * subimage_width
    subimage = binary_image[:, start_col:end_col]
    subimage = cv2.cvtColor(subimage, cv2.COLOR_GRAY2BGR)
    subimages.append(subimage)
    if i < splits-1:
      subimages.append(separator)

  concat_img = np.vstack(subimages)
  # scale down image: otherwise too large to properly display. mainly a cv2 problem
  # concat_img = cv2.resize(concat_img, (concat_img.shape[1] // 2, concat_img.shape[0] // 2))

  plt.subplot(), plt.imshow(concat_img, cmap="gray"), plt.title(f"{title} [{mat.shape[0]}x{mat.shape[1]}]")
  # plt.xticks([])
  plt.xlabel("Position in aligned sequence")
  plt.ylabel("Sequence number")
  
  mng = plt.get_current_fig_manager()
  mng.resize(*mng.window.maxsize())
  plt.show()


# ----------------- clustering


def hcluster(mat, linkage_mat):
  """Perform hierarchical clustering"""

  # Define the number of clusters you want to identify
  n_clusters = 4

  # Get the cluster labels for each row
  cluster_labels = fcluster(linkage_mat, n_clusters, criterion='maxclust')

  # Identify consensus sequences within each cluster
  consensus_list = []
  for i in range(1, n_clusters+1):
    cluster_indices = np.where(cluster_labels == i)[0]
    cluster_data = mat[cluster_indices]
    # print(f"cluster {i} data:", cluster_data)
    consensus_sequence = mode(cluster_data, axis=0).mode
    consensus_list.append(list(consensus_sequence))

  # Print consensus sequences
  # for i, consensus in enumerate(consensus_list):
  #   print(f"Cluster {i + 1} Consensus Sequence: {consensus}")

  visualize_cluster_consensuses(consensus_list)


def visualize_cluster_consensuses(consensus_sequences: list[list]):
  """Shows the consensus of a number of sequence clusters based on similarity."""
  sorted_cseqs = sorted(consensus_sequences, key=lambda x: x.count(1))

  num_blocks = len(sorted_cseqs)
  rows_per_block = 30

  subimages = []
  for cseq in sorted_cseqs:
    for j in range(rows_per_block):
      subimages.append(cseq)

  concat_img = np.vstack(subimages)

  plt.subplot(), plt.imshow(concat_img, cmap="gray"), plt.title(f"Consensuses of {num_blocks} sequence clusters")
  plt.xlabel("Position in aligned sequence")
  plt.ylabel("Sequence number")

  mng = plt.get_current_fig_manager()
  mng.resize(*mng.window.maxsize())
  plt.show()


def visualize_clusters(mat, linkage_mat) -> None:
    """Plot the original matrix with highlighted clusters in the form of a dendrogram."""
    dendrogram(linkage_mat)

    sns.set(style="white")
    sns.clustermap(mat, row_linkage=linkage_mat, col_cluster=False, method='complete')

    plt.show()


# ----------------- statistics

def show_hist(counts: list, nbins: int = 50) -> None:
  plt.subplot(), plt.title(f"Histogram of sequence lengths in the alignment")
  plt.grid()
  plt.hist(counts, bins = nbins, color="black")
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
