import copy
import numpy as np
from numpy.testing import assert_array_equal, assert_approx_equal
import pytest

import msapp.gconst as gc
from msapp.model.msa import MultiSeqAlignment

# ------------------------------
gc.DISPLAY = False
gc.VERBOSE = False


############### msa: remove empty cols

def test_msa_should_remove_empty_cols():
    """Should remove columns that only contain 1s."""
    msa = MultiSeqAlignment()
    mat = np.array([[0, 1, 1],
                    [0, 1, 0],
                    [0, 1, 1]], dtype=np.uint8)
    msa.init_from_mat(mat)
    msa.remove_empty_cols()

    should_be = np.array([[0, 1],
                          [0, 0],
                          [0, 1]], dtype=np.uint8)
    assert_array_equal(msa._mat, should_be)


def test_msa_should_not_remove_nonempty_cols():
    """Should not remove any columns if there are none that only contain 1s."""
    msa = MultiSeqAlignment()
    mat = np.array([[0, 0, 1],
                    [0, 1, 0],
                    [0, 1, 1]], dtype=np.uint8)
    should_be = copy.deepcopy(mat)
    msa.init_from_mat(mat)

    msa.remove_empty_cols()
    assert_array_equal(msa._mat, should_be)


############### msa: sort by metric

def test_msa_should_sort_rows_correctly():
    """Should reorder rows incrementally by the number of 1s they contain."""
    msa = MultiSeqAlignment()
    mat = np.array([[0, 1, 1],
                    [1, 1, 1],
                    [0, 1, 0]], dtype=np.uint8)
    msa.init_from_mat(mat)

    msa.sort_by_metric()
    should_be = np.array([[0, 1, 0],
                          [0, 1, 1],
                          [1, 1, 1]], dtype=np.uint8)
    assert_array_equal(msa._mat, should_be)


############### msa: filter by reference

def test_msa_filter_by_full_row_should_be_noop():
    """Should not change the matrix if the row that is filtered by is 0 eveywhere."""
    msa = MultiSeqAlignment()
    mat = np.array([[0, 1, 1],
                    [1, 1, 1],
                    [0, 0, 0]], dtype=np.uint8)
    should_be = copy.deepcopy(mat)
    msa.init_from_mat(mat)

    msa.filter_by_reference(2)
    assert_array_equal(msa._mat, should_be)


def test_msa_should_filter_by_correct_idx():
    """Should filter by the row values with the given valid index."""
    msa = MultiSeqAlignment()
    mat = np.array([[0, 1, 0],
                    [1, 1, 1],
                    [0, 1, 1]], dtype=np.uint8)
    msa.init_from_mat(mat)

    msa.filter_by_reference(0)
    should_be = np.array([[0, 0],
                          [1, 1],
                          [0, 1]], dtype=np.uint8)
    assert_array_equal(msa._mat, should_be)


def test_msa_should_throw_exception_when_trying_to_filter_by_bad_idx():
    """Should throw a ValueError exception if trying to filter by a bad index."""
    msa = MultiSeqAlignment()
    mat = np.array([[0, 1, 1],
                    [1, 1, 1],
                    [0, 1, 0]], dtype=np.uint8)
    msa.init_from_mat(mat)

    with pytest.raises(ValueError, match=r'The index needs to be smaller than the number of rows.'):
        msa.filter_by_reference(3)


############### msa: row column changes

def test_should_report_correct_col_row_changes():
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
