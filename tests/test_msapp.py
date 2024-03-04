import copy
import numpy as np
from numpy.testing import assert_array_equal, assert_approx_equal
import pytest

import msapp.gconst as gc
from msapp.model.msa import MultiSeqAlignment
from msapp.model.mat_manipulation import filter_by_reference, remove_empty_cols, sort_by_metric
from msapp.model.domains import filter_average_sequence

# ------------------------------
gc.DISPLAY = False
gc.VERBOSE = False


############### msa metrics: row column changes

def test_should_report_correct_col_row_changes() -> None:
    """Should return the correct number of average column and average row changes."""
    msa = MultiSeqAlignment()
    mat = np.array([[0, 1, 1],
                    [1, 1, 0],
                    [0, 1, 1]], dtype=np.uint8)
    msa.init_from_mat(mat)

    rt, ct = msa.calc_cr_metric()
    ct_should_be = 4/3  # 4/3
    rt_should_be = 1  # 3/3
    assert_approx_equal(ct, ct_should_be)
    assert_approx_equal(rt, rt_should_be)



############### mat: sort by metric

def test_msa_should_sort_rows_correctly() -> None:
    """Should reorder rows incrementally by the number of 1s they contain."""
    mat = np.array([[0, 1, 1],
                    [1, 1, 1],
                    [0, 1, 0]], dtype=np.uint8)
    result = sort_by_metric(mat)
    should_be = np.array([[0, 1, 0],
                          [0, 1, 1],
                          [1, 1, 1]], dtype=np.uint8)
    assert_array_equal(result, should_be)


############### mat: remove empty cols

def test_should_remove_empty_cols() -> None:
    """Should remove columns that only contain 1s."""
    mat = np.array([[0, 1, 1],
                    [0, 1, 0],
                    [0, 1, 1]], dtype=np.uint8)
    result = remove_empty_cols(mat)
    should_be = np.array([[0, 1],
                          [0, 0],
                          [0, 1]], dtype=np.uint8)
    assert_array_equal(result, should_be)


def test_should_not_remove_nonempty_cols() -> None:
    """Should not remove any columns if there are none that only contain 1s."""
    mat = np.array([[0, 0, 1],
                    [0, 1, 0],
                    [0, 1, 1]], dtype=np.uint8)
    should_be = copy.deepcopy(mat)
    result = remove_empty_cols(mat)
    assert_array_equal(result, should_be)


############### mat: filter by reference

def test_filter_by_full_row_should_be_noop() -> None:
    """Should not change the matrix if the row that is filtered by is 0 everywhere."""
    mat = np.array([[0, 1, 1],
                    [1, 1, 1],
                    [0, 0, 0]], dtype=np.uint8)
    should_be = copy.deepcopy(mat)
    result = filter_by_reference(mat, 2)
    assert_array_equal(result, should_be)


def test_should_filter_by_correct_idx() -> None:
    """Should filter by the row values with the given valid index."""
    mat = np.array([[0, 1, 0],
                    [1, 1, 1],
                    [0, 1, 1]], dtype=np.uint8)
    result = filter_by_reference(mat, 0)
    should_be = np.array([[0, 0],
                          [1, 1],
                          [0, 1]], dtype=np.uint8)
    assert_array_equal(result, should_be)


def test_should_throw_exception_when_trying_to_filter_by_bad_idx() -> None:
    """Should throw a ValueError exception if trying to filter by a bad index."""
    mat = np.array([[0, 1, 1],
                    [1, 1, 1],
                    [0, 1, 0]], dtype=np.uint8)

    with pytest.raises(ValueError, match=r'The index needs to be smaller than the number of rows.'):
        filter_by_reference(mat, 3)

############### domains

# connecting regions

def test_should_connect_close_black_regions() -> None:
    before = np.array([0, 0, 0, 1, 0, 0, 1, 1, 0, 0])
    min_width = 1
    min_gap_len = 3
    result = filter_average_sequence(before.copy(), min_width, min_gap_len)
    should_be = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    assert_array_equal(result, should_be)

def test_should_not_connect_beginning() -> None:
    before = np.array([1, 1, 0, 0, 0, 0, 0])
    min_width = 1
    min_gap_len = 3
    result = filter_average_sequence(before.copy(), min_width, min_gap_len)
    should_be = np.array([1, 1, 0, 0, 0, 0, 0])
    assert_array_equal(result, should_be)

def test_should_not_connect_far_apart_black_regions() -> None:
    before = np.array([1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0])
    min_width = 1
    min_gap_len = 3
    result = filter_average_sequence(before.copy(), min_width, min_gap_len)
    should_be = np.array([1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0])
    assert_array_equal(result, should_be)

# removing regions

def test_should_remove_small_black_regions() -> None:
    before = np.array([1, 1, 0, 1, 0, 0, 1, 1, 1])
    min_width = 3
    min_gap_len = 1
    result = filter_average_sequence(before.copy(), min_width, min_gap_len)
    should_be = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1])
    assert_array_equal(result, should_be)

def test_should_remove_small_black_end_regions() -> None:
    before = np.array([1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0])
    min_width = 2
    min_gap_len = 1
    result = filter_average_sequence(before.copy(), min_width, min_gap_len)
    should_be = np.array([1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1])
    assert_array_equal(result, should_be)

def test_should_keep_large_black_regions() -> None:
    before = np.array([0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1])
    min_width = 3
    min_gap_len = 1
    result = filter_average_sequence(before.copy(), min_width, min_gap_len)
    should_be = np.array([0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1])
    assert_array_equal(result, should_be)

# connecting and removing regions

def test_fuse_keep_remove_black_regions() -> None:
    before = np.array([0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1])
    min_width = 3
    min_gap_len = 3
    result = filter_average_sequence(before.copy(), min_width, min_gap_len)
    should_be = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1])
    assert_array_equal(result, should_be)