import numpy as np
import pytest

from csd import embedding as emb


def test_sum_accumulate_uses_complete_blocks_and_ignores_nan():
    out = emb.sum_accumulate([1, 2, 3, np.nan, 5, 6, 7], 3)
    np.testing.assert_allclose(out, [6.0, 11.0])  # last value dropped, NaN skipped


def test_sum_accumulate_all_nan_block_is_nan():
    out = emb.sum_accumulate([np.nan, np.nan, 1.0, 1.0], 2)
    assert np.isnan(out[0]) and out[1] == 2.0


def test_delay_vectors_column_order_matches_matlab():
    v = emb.delay_vectors(np.arange(10.0), 3, 2)
    assert v.shape == (6, 3)
    np.testing.assert_array_equal(v[0], [4, 2, 0])  # most recent value first
    np.testing.assert_array_equal(v[-1], [9, 7, 5])


def test_series_from_vectors_inverts_delay_vectors():
    x = np.random.default_rng(0).random(30)
    v = emb.delay_vectors(x, 4, 3)
    np.testing.assert_allclose(emb.series_from_vectors(v, 4, 3), x)


def test_delay_vectors_rejects_short_series():
    with pytest.raises(ValueError):
        emb.delay_vectors(np.arange(3.0), 3, 2)


def test_paired_rows_equal_embeddings_reduce_to_matlab_truncation():
    gcm_rows, sta_rows = emb.paired_rows(100, 100, 2, 2, lag=4)
    assert (gcm_rows.start, gcm_rows.stop) == (0, 96)
    assert (sta_rows.start, sta_rows.stop) == (4, 100)


def test_paired_rows_pair_by_time_when_embeddings_differ():
    n = 50
    sta = emb.delay_vectors(np.arange(float(n)), 3, 2)  # first coordinate = time of row + 4
    gcm = emb.delay_vectors(np.arange(float(n)), 2, 1)  # first coordinate = time of row + 1
    g, s = emb.paired_rows(len(sta), len(gcm), 4, 1, lag=0)
    np.testing.assert_array_equal(sta[list(s)][:, 0], gcm[list(g)][:, 0])
