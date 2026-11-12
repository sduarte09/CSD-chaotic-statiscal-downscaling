import numpy as np
import pandas as pd
import pytest

from csd.pipeline import downscale_station


def _series(values, start="2000-01-01"):
    return pd.Series(values, index=pd.date_range(start, periods=len(values), freq="D"))


def test_downscale_station_end_to_end_on_synthetic_chaos():
    from conftest import logistic

    n = 2400
    station = _series(logistic(0.3, n) * 10)
    model = _series(logistic(0.37, n) * 6 + 1)
    res = downscale_station(station, model, step=5, dims=(2, 2), delays=(1, 1),
                            cal_end="2004-12-31", val_end="2005-12-31", app_end="2006-06-30")
    assert res.calibration.lag >= 0
    assert len(res.validation) == len(res.observed_validation) == 365 // 5
    assert res.validation_rmse >= 0
    assert res.downscaled.index[0] == pd.Timestamp("2004-12-30")  # day after the last complete block
    assert np.isfinite(res.downscaled.values).all()
    assert res.downscaled.index.is_monotonic_increasing
    gaps = res.downscaled.index.to_series().diff().dropna().dt.days
    assert (gaps == 5).all()


def test_short_calibration_is_rejected():
    s = _series(np.random.default_rng(0).random(60))
    with pytest.raises(ValueError):
        downscale_station(s, s, step=5, dims=(3, 3), delays=(2, 2), cal_end="2000-01-20")
