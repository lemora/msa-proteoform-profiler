import cv2
import numpy as np

import msapp.gconst as gc
from msapp.visualize import show_pre_post

def blur(bimg, ksize=9):
  blur = cv2.blur(bimg, (ksize, ksize))
  if gc.DISPLAY: show_pre_post(bimg, blur, f"Blurring, kernel size={ksize}")
  return to_binary_matrix(blur)


def gaussian_blur(bimg, ksize=5):
  blur = cv2.GaussianBlur(bimg, (ksize, ksize), cv2.BORDER_DEFAULT)
  if gc.DISPLAY: show_pre_post(bimg, blur, f"Gaussian Blur, kernel size={ksize}")
  return to_binary_matrix(blur)


def dilate_erode(bimg, ksize=5):
  kernel = np.ones((ksize, ksize), np.uint8)
  morphed = cv2.morphologyEx(bimg, cv2.MORPH_OPEN, kernel)
  if gc.DISPLAY: show_pre_post(bimg, morphed, f"Morphing (Dilate, then Erode), kernel size={ksize}")
  return to_binary_matrix(morphed)


############### pre/post formatting

def to_binary_matrix(img):
  """Makes sure that after image processing the contained values are binary again: 0 or 1."""
  maxval = np.max(img)
  m = img / maxval
  m = np.round(m)
  threshold = 0.5
  m = (m >= threshold).astype(int)
  mat = np.array(img, dtype=np.uint8)
  return mat