import copy
import cv2
import numpy as np

import msapp.gconst as gc


# ----------------- mat manipulation and plotting preparation

def filter_by_reference(mat, idx: int) -> np.array:
    """Filters the alignment matrix by a given row in that MSA for visualization purposes."""
    if idx >= mat.shape[0]: raise ValueError("The index needs to be smaller than the number of rows.")
    if gc.VERBOSE: print(f"-- OP: Filter by reference with index {idx}")

    refseq = mat[idx]
    cols_to_remove = np.array([i for i, val in enumerate(refseq) if val == 1])
    filtered = remove_seqs_from_alignment(mat, cols_to_remove, cols=True)
    return filtered


def remove_empty_cols(mat) -> np.array:
    """Removes all columns that are empty from the matrix, meaning they only contain the value 1."""
    empty_columns = np.where(np.all(mat == 1, axis=0))[0]
    filtered = remove_seqs_from_alignment(mat, empty_columns, cols=True)
    return filtered


def remove_seqs_from_alignment(mat, idx_list: np.ndarray[int], cols: bool = True) -> np.array:
    """Removes the columns (else rows) which have the given indices.
    param cols: should remove columns, else remove rows."""
    if len(idx_list) == 0: return mat
    max_idx = mat.shape[1] if cols else mat.shape[0]
    if min(idx_list) < 0 or max(idx_list) >= max_idx:
        raise ValueError("The indices to delete must be between 0 and max row or col count.")

    filtered = np.delete(mat, idx_list, axis=(1 if cols else 0))
    return filtered


def clear_seqs_in_alignment(mat, idx_list: np.ndarray[int]) -> np.array:
    """Removes the columns (else rows) which have the given indices.
    param cols: should remove columns, else remove rows."""
    if len(idx_list) == 0: return mat
    if min(idx_list) < 0 or max(idx_list) >= mat.shape[0]:
        raise ValueError("The indices to delete must be between 0 and max row or col count.")
    mat[idx_list] = 1
    return mat


def sort_by_metric(mat, sorting_metric=lambda row: sum(row)) -> np.array:
    """Sorts a MSA binary matrix by the given function that works on a binary list."""
    sorting_metrics = np.apply_along_axis(sorting_metric, axis=1, arr=mat)
    sorted_indices = np.argsort(sorting_metrics)
    sorted_matrix = mat[sorted_indices]
    return sorted_matrix


# ----------------- classic image processing, mainly convolution

def blur(img: np.array, ksize=9) -> np.array:
    """Blur image."""
    if gc.VERBOSE: print(f"-- OP: Blur ({ksize}x{ksize} kernel)")
    processed = cv2.blur(img, (ksize, ksize))
    processed = to_binary_matrix(processed)
    return processed


def gaussian_blur(img: np.array, ksize=5) -> np.array:
    """Gaussian blur image."""
    if gc.VERBOSE: print(f"-- OP: Gaussian blur ({ksize}x{ksize} kernel)")
    processed = cv2.GaussianBlur(img, (ksize, ksize), cv2.BORDER_DEFAULT)
    processed = to_binary_matrix(processed)
    return processed


def median_blur(img: np.array, ksize=5) -> np.array:
    """Median blur image."""
    if gc.VERBOSE: print(f"-- OP: Median blur ({ksize}x{ksize} kernel)")
    processed = cv2.medianBlur(img, ksize)
    processed = to_binary_matrix(processed)
    return processed


def dilate_erode(img: np.array, ksize=5) -> np.array:
    """Dilate, then erode image."""
    if gc.VERBOSE: print(f"-- OP: Dilate/erode ({ksize}x{ksize} kernel)")
    kernel = np.ones((ksize, ksize), np.uint8)
    processed = cv2.morphologyEx(img, cv2.MORPH_OPEN, kernel)
    processed = to_binary_matrix(processed)
    return processed


def cross_convolve(img: np.array, row_size: int = 1, col_size: int = 3) -> np.array:
    """Convolves with a cross-shaped kernel with row_size ones on the center row and col_size 1s on the center
    column."""
    ksize = max(row_size, col_size)
    if gc.VERBOSE: print(f"-- OP: Convolving with {ksize}x{ksize} cross-kernel (r:{row_size}, c:{col_size})")

    kernel = create_cross_kernel(row_size=row_size, col_size=col_size)
    processed = cv2.filter2D(src=img, ddepth=-1, kernel=kernel)
    processed = to_binary_matrix(processed)

    return processed


def create_cross_kernel(row_size: int = 1, col_size: int = 3) -> np.array:
    """Creates a kernel of size max(row_size, col_size), that is zero everywhere except for the center column and row.
    param row_size: hom many centered 1s to have on the center row
    param col_size: hom many centered 1s to have on the center column
    """
    assert row_size > 0 and col_size > 0
    if row_size % 2 == 0 or col_size % 2 == 0:
        raise ValueError("Kernel size must be an odd number.")

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


# ----------------- pre/post formatting

def to_binary_matrix(img: np.array) -> np.array:
    """Makes sure that after image processing the contained values are binary again: 0 or 1."""
    mat = copy.deepcopy(img)
    minval = np.min(mat)
    if minval < 0:
        mat = mat - minval

    maxval = np.max(mat)
    if maxval <= 0:
        raise ValueError(f"Max matrix entry {maxval} <= 0!")
    if maxval != 1:
        mat = mat / maxval
    mat = np.round(mat)
    mat = np.array(mat, dtype=np.uint8)

    return mat
