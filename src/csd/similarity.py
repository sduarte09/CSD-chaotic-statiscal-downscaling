"""Calibration of the synchronization model: amplitude, lag and mutual false neighbours."""
from __future__ import annotations

import numpy as np


def amplitude_ratio(station_cal, model_cal):
    """Ratio of the calibration-period means, used to rescale the model series."""
    return float(np.mean(station_cal) / np.mean(model_cal))


def similarity_function(station, model, lag, amplitude=1.0):
    """Generalized similarity function Sg between station and lagged model series.

    Smaller is more similar. Both series are cut and centred as in the original code,
    so that station values from position ``lag`` onwards are compared with the model values
    from the start (lag 0 compares the series as they are). ``amplitude`` rescales the centred model series; as in the MATLAB
    expression it enters the normalisation once (not squared). With ``amplitude=1`` the
    result is the plain similarity function.
    """
    station = np.asarray(station, dtype="float64")
    model = np.asarray(model, dtype="float64")
    n = station.size
    shifted = station[lag:]
    base = station[: n - lag]
    other = model[: n - lag]
    if shifted.size == 0:
        return np.nan
    diff = (shifted - shifted.mean()) - amplitude * (other - other.mean())
    scale = np.sqrt(np.mean((base - base.mean()) ** 2) * amplitude * np.mean((other - other.mean()) ** 2))
    return float(np.sqrt(np.mean(diff**2) / scale))


def best_lag(station, model, step, max_days=180):
    """Lag (in accumulation steps) that minimises the similarity function.

    Lags from 0 to ``max_days // step - 1`` are tried; the first minimum is returned.
    Returns ``(lag, similarity_values)``.
    """
    n_lags = max_days // step
    values = np.array([similarity_function(station, model, lag) for lag in range(n_lags)])
    return int(np.nanargmin(values)), values


def _nearest_other(vectors, t):
    """Index of and distance to the nearest vector at non-zero distance from vector t."""
    dist = np.linalg.norm(vectors - vectors[t], axis=1)
    dist[dist == 0] = np.inf  # excludes the vector itself and exact duplicates
    j = int(np.argmin(dist))
    return (j, float(dist[j])) if np.isfinite(dist[j]) else (None, np.nan)


def mutual_false_neighbors(station_vectors, model_vectors):
    """Mutual false nearest neighbour ratios ``ss`` and their running mean ``u``.

    For each time t, take the nearest neighbour of the model vector and of the station
    vector (ignoring identical vectors) and compare the distances measured in the other
    space: ``ss = (d_sta(t, nn_model) * d_model(t, nn_sta)) / (d_model(t, nn_model) *
    d_sta(t, nn_sta))``. ``u[t]`` is the mean of ``ss[: t + 1]``.

    This follows the MATLAB implementation. The Python 2 script measured one distance
    in the wrong space and indexed neighbours after deleting the vector itself; see the
    README.
    """
    station_vectors = np.asarray(station_vectors, dtype="float64")
    model_vectors = np.asarray(model_vectors, dtype="float64")
    n = len(station_vectors)
    if len(model_vectors) != n:
        raise ValueError("Station and model vectors must have the same length")
    ss = np.full(n, np.nan)
    for t in range(n):
        j_model, d_model_model = _nearest_other(model_vectors, t)
        j_sta, d_sta_sta = _nearest_other(station_vectors, t)
        if j_model is None or j_sta is None:
            continue
        d_sta_model = np.linalg.norm(station_vectors[t] - station_vectors[j_model])
        d_model_sta = np.linalg.norm(model_vectors[t] - model_vectors[j_sta])
        ss[t] = (d_sta_model * d_model_sta) / (d_model_model * d_sta_sta)
    u = np.cumsum(ss) / np.arange(1, n + 1)
    return ss, u
