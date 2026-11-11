import numpy as np
import pandas as pd
import pytest

from csd import accumulation as acc


def test_mean_accumulate_uses_complete_blocks_only():
    x = np.arange(1.0, 12.0)  # 11 values
    out = acc.mean_accumulate(x, 5)
    assert list(out) == [3.0, 8.0]  # means of 1..5 and 6..10; the 11 is dropped


def test_mean_accumulate_ignores_nan_inside_a_block_and_keeps_all_nan_blocks_as_nan():
    x = np.array([1.0, np.nan, 3.0, np.nan, np.nan, np.nan])
    out = acc.mean_accumulate(x, 3)
    assert out[0] == pytest.approx(2.0)
    assert np.isnan(out[1])


def test_mean_accumulate_short_series_and_bad_step():
    assert acc.mean_accumulate([1.0, 2.0], 5).size == 0
    with pytest.raises(ValueError):
        acc.mean_accumulate([1.0, 2.0], 0)


def test_disaggregate_preserves_block_totals():
    x = np.array([2.0, 0.5, 4.0])
    out = acc.disaggregate(x, 5, rng=0)
    assert out.shape == (15,)
    np.testing.assert_allclose(out.reshape(3, 5).sum(axis=1), x * 5)
    np.testing.assert_allclose(acc.mean_accumulate(out, 5), x)


def test_disaggregate_is_reproducible_with_a_seed():
    x = np.array([1.0, 2.0, 3.0])
    np.testing.assert_array_equal(acc.disaggregate(x, 7, rng=42), acc.disaggregate(x, 7, rng=42))
    assert not np.array_equal(acc.disaggregate(x, 7, rng=1), acc.disaggregate(x, 7, rng=2))


def _daily(start, end):
    dates = pd.date_range(start, end, freq="D")
    return dates, np.ones(len(dates))


def test_monthly_series_sums_days_and_drops_partial_years():
    dates, values = _daily("1999-06-15", "2002-03-10")  # complete years: 2000 and 2001
    monthly = acc.monthly_series(values, dates)
    assert len(monthly) == 24
    assert monthly.loc[(2000, 2)] == 29  # leap year
    assert monthly.loc[(2001, 2)] == 28
    assert monthly.loc[(2001, 12)] == 31
    assert (1999, 12) not in monthly.index


def test_monthly_climatology_averages_over_years():
    dates, values = _daily("2000-01-01", "2002-12-31")
    clim = acc.monthly_climatology(values, dates)
    assert list(clim.index) == list(range(1, 13))
    assert clim.loc[2] == pytest.approx((29 + 28 + 28) / 3)
    assert clim.loc[1] == 31


def test_monthly_mean_for_non_precipitation_variables():
    dates, _ = _daily("2000-01-01", "2000-12-31")
    temperature = np.full(len(dates), 15.0)
    monthly = acc.monthly_series(temperature, dates, how="mean")
    assert (monthly == 15.0).all()


def test_annual_series_totals():
    dates, values = _daily("2000-01-01", "2002-12-31")
    annual = acc.annual_series(values, dates)
    assert list(annual.index) == [2000, 2001, 2002]
    assert list(annual.values) == [366, 365, 365]
