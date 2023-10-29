import cv2 as cv
from matplotlib import pyplot as plt
import numpy as np

from msapp.visualize import show_pre_post

def blur(bimg, show=False):
  blur = cv.blur(bimg, (9, 9))
  if show: show_pre_post(bimg, blur, "Blurring")
  return to_binary_matrix(blur)


def dilate_erode(bimg, show=False):
  kernel = np.ones((5, 5), np.uint8)
  morphed = cv.morphologyEx(bimg, cv.MORPH_OPEN, kernel)
  if show: show_pre_post(bimg, morphed, "Morphing (Dilate, then Erode)")
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