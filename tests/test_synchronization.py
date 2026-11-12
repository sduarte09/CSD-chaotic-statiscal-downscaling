import numpy as np
import pytest

from csd.calibration import calibrate
from csd.synchronization import rmse, synchronize


def _run(sta, mod, n_cal=500, dims=(2, 2), taus=(1, 1)):
    c = calibrate(sta[:n_cal], mod[:n_cal], dims[0], taus[0], dims[1], taus[1], step=5)
    return c, synchronize(sta[:n_cal], mod[:n_cal], mod[n_cal:], c, dims[0], taus[0], dims[1], taus[1])


def test_identical_station_and_model_are_reproduced(chaotic_pair):
    sta, _ = chaotic_pair
    _, result = _run(sta, sta.copy())
    assert len(result.series) == len(sta) - 500
    # chaos amplifies small errors, so only the first steps can match
    np.testing.assert_allclose(result.series[:8], sta[500:508], rtol=0.02, atol=0.5)


def test_output_has_one_value_per_application_step_and_is_finite(chaotic_pair):
    sta, mod = chaotic_pair
    for dims, taus in [((2, 2), (1, 1)), ((3, 2), (2, 1))]:
        c, result = _run(sta, mod, dims=dims, taus=taus)
        assert len(result.series) == len(mod) - 500
        assert np.isfinite(result.series).all()
        assert result.series.min() >= 0  # roots below zero are discarded


def test_downscaled_series_keeps_the_station_distribution(chaotic_pair):
    sta, mod = chaotic_pair  # the model has a different mean and spread
    _, result = _run(sta, mod)
    assert result.series.mean() == pytest.approx(sta.mean(), rel=0.1)
    assert result.series.std() == pytest.approx(sta.std(), rel=0.1)


def test_synchronization_is_deterministic(chaotic_pair):
    sta, mod = chaotic_pair
    a = _run(sta, mod)[1].series
    b = _run(sta, mod)[1].series
    np.testing.assert_array_equal(a, b)


def test_calibration_of_a_different_length_is_rejected(chaotic_pair):
    sta, mod = chaotic_pair
    c = calibrate(sta[:400], mod[:400], 2, 1, 2, 1, step=5)
    with pytest.raises(ValueError):
        synchronize(sta[:500], mod[:500], mod[500:], c, 2, 1, 2, 1)


def test_rmse_ignores_nan_pairs_and_uses_common_length():
    assert rmse([1.0, 2.0, np.nan, 9.0], [1.0, 4.0, 3.0]) == pytest.approx(np.sqrt(2.0))
