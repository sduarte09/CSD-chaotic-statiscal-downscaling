"""Delay vectors and block sums used by the calibration, validation and synchronization."""
from __future__ import annotations

import numpy as np


def sum_accumulate(x, step):
    """Sum of consecutive complete blocks of ``step`` values (NaN ignored within a block).

    A block that is entirely NaN gives NaN. Returns ``len(x) // step`` values.
    """
    x = np.asarray(x, dtype="float64")
    if step < 1:
        raise ValueError("step must be at least 1")
    n_blocks = x.size // step
    blocks = x[: n_blocks * step].reshape(n_blocks, step)
    out = np.nansum(blocks, axis=1)
    out[np.isnan(blocks).all(axis=1)] = np.nan
    return out


def delay_vectors(series, dimension, delay):
    """Delay vectors in the column order used by the MATLAB code.

    Row ``r`` is ``[x[r + (D-1)*tau], ..., x[r + tau], x[r]]``: the first column is the
    most recent value and the last column the oldest. The result has
    ``len(series) - (D-1)*tau`` rows.
    """
    series = np.asarray(series, dtype="float64")
    n_rows = series.size - (dimension - 1) * delay
    if n_rows <= 0:
        raise ValueError("Series too short for this dimension and delay")
    return np.column_stack(
        [series[delay * (dimension - 1 - i) : delay * (dimension - 1 - i) + n_rows] for i in range(dimension)]
    )


def series_from_vectors(vectors, dimension, delay):
    """Invert :func:`delay_vectors`: last column plus the tail of the first column."""
    vectors = np.asarray(vectors, dtype="float64")
    return np.concatenate([vectors[:, -1], vectors[len(vectors) - delay * (dimension - 1) :, 0]])


def paired_rows(n_station, n_model, offset_station, offset_model, lag):
    """Row ranges that pair a model vector at time ``t`` with a station vector at ``t + lag``.

    The first coordinate of station row ``r`` is the series value at time
    ``r + offset_station`` (``offset = (D - 1) * tau``) and likewise for the model.
    Returns ``(model_rows, station_rows)`` as ``range`` objects of equal length. Pairing by
    time keeps the two attractors aligned even when their embeddings differ; when the
    embeddings are equal it reduces to the MATLAB ``STACALVEC(1+LAG:end)`` /
    ``GCMCALVEC(1:end-LAG)`` truncation.
    """
    shift = offset_model + lag - offset_station  # station row = model row + shift
    lo = max(0, -shift)
    hi = min(n_model, n_station - shift)
    if hi <= lo:
        raise ValueError("No overlapping vectors for this lag, dimension and delay")
    return range(lo, hi), range(lo + shift, hi + shift)
