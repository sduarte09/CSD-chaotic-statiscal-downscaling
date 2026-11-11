"""Accumulation of daily series over several days, and aggregation to months and years."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd


def mean_accumulate(x, step):
    """Mean of consecutive blocks of ``step`` values (NaN ignored within a block).

    Only complete blocks are returned, so the result has ``len(x) // step`` values.
    """
    x = np.asarray(x, dtype="float64")
    if step < 1:
        raise ValueError("step must be at least 1")
    n_blocks = x.size // step
    if n_blocks == 0:
        return np.array([], dtype="float64")
    blocks = x[: n_blocks * step].reshape(n_blocks, step)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)  # blocks that are all NaN
        return np.nanmean(blocks, axis=1)


def disaggregate(x, step, rng=None):
    """Spread each block mean over ``step`` days with random weights.

    Each block keeps its total (``x[k] * step``), so ``mean_accumulate`` of the
    result recovers ``x``. Pass an integer seed or a ``numpy.random.Generator``
    as ``rng`` to make the result reproducible.
    """
    x = np.asarray(x, dtype="float64")
    rng = np.random.default_rng(rng)
    weights = rng.random((x.size, step))
    return (x[:, None] * step * weights / weights.sum(axis=1, keepdims=True)).ravel()


def _complete_years(values, dates):
    series = pd.Series(np.asarray(values, dtype="float64"), index=pd.DatetimeIndex(dates))
    first, last = series.index[0], series.index[-1]
    first_year = first.year if (first.month, first.day) == (1, 1) else first.year + 1
    last_year = last.year if (last.month, last.day) == (12, 31) else last.year - 1
    return series[(series.index.year >= first_year) & (series.index.year <= last_year)]


def monthly_series(values, dates, how="sum"):
    """Monthly values for the complete calendar years in a daily series.

    ``how`` is ``"sum"`` for precipitation (mm per month) or ``"mean"`` for
    variables such as temperature. Returns a Series indexed by (year, month).
    """
    series = _complete_years(values, dates)
    return series.groupby([series.index.year, series.index.month]).agg(how)


def monthly_climatology(values, dates, how="sum"):
    """Mean over the years of each calendar month (12 values, January first)."""
    monthly = monthly_series(values, dates, how)
    return monthly.groupby(level=1).mean()


def annual_series(values, dates, how="sum"):
    """Annual values for the complete calendar years in a daily series."""
    series = _complete_years(values, dates)
    return series.groupby(series.index.year).agg(how)
