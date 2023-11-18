import copy

import cv2
import numpy as np

import msapp.gconst as gc
from msapp.view.visualization import show_pre_post


def blur(img: np.array, ksize=9, show=None):
    """Blur image."""
    print(f"\n-- OP: Blur ({ksize}x{ksize} kernel)")
    show = show if show is not None else gc.DISPLAY
    processed = cv2.blur(img, (ksize, ksize))
    processed = to_binary_matrix(processed)
    if show: show_pre_post(img, processed, f"Blurring, {ksize}x{ksize} kernel")
    return processed


def gaussian_blur(img: np.array, ksize=5, show=None):
    """Gaussian blur image."""
    print(f"\n-- OP: Gaussian blur ({ksize}x{ksize} kernel)")
    show = show if show is not None else gc.DISPLAY
    processed = cv2.GaussianBlur(img, (ksize, ksize), cv2.BORDER_DEFAULT)
    processed = to_binary_matrix(processed)
    if show: show_pre_post(img, processed, f"Gaussian Blur, {ksize}x{ksize} kernel")
    return processed


def median_blur(img: np.array, ksize=5, show=None):
    """Median blur image."""
    print(f"\n-- OP: Median blur ({ksize}x{ksize} kernel)")
    show = show if show is not None else gc.DISPLAY
    processed = cv2.medianBlur(img, ksize)
    processed = to_binary_matrix(processed)
    if show: show_pre_post(img, processed, f"Median Blur, {ksize}x{ksize} kernel")
    return processed


def dilate_erode(img: np.array, ksize=5, show=None):
    """Dilate, then erode image."""
    print(f"\n-- OP: Dilate/erode ({ksize}x{ksize} kernel)")
    show = show if show is not None else gc.DISPLAY
    kernel = np.ones((ksize, ksize), np.uint8)
    processed = cv2.morphologyEx(img, cv2.MORPH_OPEN, kernel)
    processed = to_binary_matrix(processed)
    if show: show_pre_post(img, processed, f"Morphing (Dilate, then Erode), {ksize}x{ksize} kernel")
    return processed


def cross_convolve(img: np.array, row_size: int = 1, col_size: int = 3, show=None):
    """Convolves with a cross-shaped kernel with row_size ones on the center row and col_size 1s on the center
    column."""
    ksize = max(row_size, col_size)
    print(f"\n-- OP: Convolving with {ksize}x{ksize} cross-kernel (r:{row_size}, c:{col_size})")

    kernel = create_cross_kernel(row_size=row_size, col_size=col_size)
    processed = cv2.filter2D(src=img, ddepth=-1, kernel=kernel)
    processed = to_binary_matrix(processed)

    if show if show is not None else gc.DISPLAY:
        show_pre_post(img, processed,
                      f"Convolving with {ksize}x{ksize} cross-kernel (1s on row:{row_size}, on col:{col_size})")
    return processed


def create_cross_kernel(row_size: int = 1, col_size: int = 3):
    """Creates a kernel of size max(row_size, col_size), that is zero everywhere except for the center column and row.
    param row_size: hom many centered 1s to have on the center row
    param col_size: hom many centered 1s to have on the center column
    """
    assert row_size > 0 and col_size > 0
    if row_size % 2 == 0 or col_size % 2 == 0:
        raise ValueError("Kernel size must be an odd number.")

    if gc.VERBOSE: print(f"Cross kernel with row size: {row_size}, col size: {col_size}")
    ksize = max(col_size, row_size)
    kernel = np.zeros((ksize, ksize), dtype=int)

    middle = ksize // 2

    # middle col
    if col_size > 1:
        start_row = middle - (col_size // 2)
        end_row = start_row + col_size
        kernel[start_row:end_row, middle] = 1

    # middle row
    if row_size > 1:
        start_col = middle - (row_size // 2)
        end_col = start_col + row_size
        kernel[middle, start_col:end_col] = 1

    return kernel


############### pre/post formatting

def to_binary_matrix(img: np.array):
    """Makes sure that after image processing the contained values are binary again: 0 or 1."""
    mat = copy.deepcopy(img)
    minval = np.min(mat)
    # print(f"Min matrix val: {minval}")
    if minval < 0:
        print(f"Min matrix entry: {minval} < 0!")
        mat = mat - minval

    maxval = np.max(mat)
    # print(f"Max matrix val: {maxval}")
    if maxval <= 0:
        print(f"Max matrix entry {maxval} <= 0!")
    elif maxval != 1:
        mat = mat / maxval
    mat = np.round(mat)
    mat = np.array(mat, dtype=np.uint8)

    if gc.VERBOSE:
        print(f"After mat binarization, min: {np.min(mat)}, max: {np.max(mat)}")

    # print(mat)
    # print(f"Sum of values in mat: {sum(mat)}")
    return mat
