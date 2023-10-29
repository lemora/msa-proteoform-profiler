import cv2
from matplotlib import pyplot as plt
import numpy as np


def show_pre_post(pre, post, title:str):
  plt.subplot(121), plt.imshow(pre, cmap="gray"), plt.title(f'Before ({pre.shape[0]} rows, {pre.shape[1]} cols)')
  plt.xticks([]), plt.yticks([])
  plt.subplot(122), plt.imshow(post, cmap="gray"), plt.title(f'After ({post.shape[0]} rows, {post.shape[1]} cols)')
  plt.xticks([]), plt.yticks([])
  mng = plt.get_current_fig_manager()
  mng.resize(*mng.window.maxsize())
  mng.window.title(f"{title} ({pre.shape[0]} rows, {pre.shape[1]} cols)")
  plt.show()


def show(msa_mat, title: str):
  if msa_mat.shape[1] > 3000:
    show_as_subimages(msa_mat, title)
  else:
    show_as_one(msa_mat, title)


def show_as_one(mat, title: str):
  """Shows the alignmet as a binary image."""
  img = np.array(mat, dtype=np.uint8) * 255

  plt.subplot(), plt.imshow(img, cmap="gray"), plt.title(f"{title} ({mat.shape[0]} rows, {mat.shape[1]} cols)")
  plt.xticks([]), plt.yticks([])
  mng = plt.get_current_fig_manager()
  mng.resize(*mng.window.maxsize())
  plt.show()


def show_as_subimages(mat, title: str):
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

  plt.subplot(), plt.imshow(concat_img, cmap="gray"), plt.title(f"{title} ({mat.shape[0]} rows, {mat.shape[1]} cols)")
  plt.xticks([]), plt.yticks([])
  mng = plt.get_current_fig_manager()
  mng.resize(*mng.window.maxsize())
  plt.show()