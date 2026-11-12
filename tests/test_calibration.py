import numpy as np
import pytest

from csd import calibration as cal


def test_amplitude_grid_starts_at_a_tenth_and_ends_near_the_larger_ratio():
    sta, gcm = np.array([2.0, 8.0, 4.0]), np.array([1.0, 2.0, 2.0])
    grid = cal.amplitude_grid(sta, gcm)
    assert grid[0] == pytest.approx(0.1)
    assert grid[-1] <= 4.0 + 1e-9 and grid[-1] > 3.5  # max ratio 4, mean ratio 2.14
    assert np.all(np.diff(grid) > 0)


def test_calibration_recovers_lag_and_amplitude(chaotic_pair):
    sta, _ = chaotic_pair
    model = np.roll(sta, -3) * 0.5  # model leads the station by 3 steps and is half as large
    result = cal.calibrate(sta[:500], model[:500], 2, 1, 2, 1, step=5)
    assert result.lag == 3
    assert result.amplitude == pytest.approx(2.0, rel=0.1)
    assert result.sg.shape == (result.amplitudes.size, result.lags.size)
    assert result.lags[-1] == 30 // 5 - 1


def test_calibration_ss_has_one_value_per_paired_vector(chaotic_pair):
    sta, mod = chaotic_pair
    result = cal.calibrate(sta[:300], mod[:300], 2, 1, 2, 1, step=10)
    assert len(result.ss) == len(result.u)
    assert result.u[0] == pytest.approx(result.ss[0])
