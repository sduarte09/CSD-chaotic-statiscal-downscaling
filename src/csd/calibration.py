"""Calibration: lag and amplitude from the generalized similarity function, then ss and u."""
from __future__ import annotations

from typing import NamedTuple

import numpy as np

from .embedding import delay_vectors, paired_rows
from .similarity import amplitude_ratio, mutual_false_neighbors, similarity_function


class Calibration(NamedTuple):
    lag: int  # in accumulation steps
    amplitude: float
    ss: np.ndarray  # mutual false neighbour ratio over the calibration vectors
    u: np.ndarray  # running mean of ss
    sg: np.ndarray  # similarity surface, shape (n_amplitudes, n_lags)
    amplitudes: np.ndarray
    lags: np.ndarray


def amplitude_grid(station_cum, model_cum, n_points=100, start=0.1):
    """Amplitudes tried by the calibration: ``start`` up to the larger of the max-ratio and
    mean-ratio of station to model, in about ``n_points`` steps (MATLAB ``a:d:b`` range)."""
    ratios = [np.max(station_cum) / np.max(model_cum), amplitude_ratio(station_cum, model_cum)]
    upper, lower = max(ratios), min(ratios)
    if upper <= start:
        return np.array([upper])
    d = (upper - lower) / (n_points - 1)
    if d == 0:
        return np.array([upper])
    n = int(np.floor((upper - start) / d + 1e-9)) + 1
    return start + d * np.arange(n)


def calibrate(station_cum, model_cum, dim_station, delay_station, dim_model, delay_model, step,
              max_days=30, n_amplitudes=100):
    """Calibrate lag, amplitude, ``ss`` and ``u`` on accumulated calibration series.

    ``station_cum`` and ``model_cum`` are the series after accumulation with ``step``
    days. The lag is searched from 0 to ``max_days // step - 1`` accumulation steps.
    The similarity function is evaluated on the accumulated series (the MATLAB code does
    this through linear indexing, which amounts to the first coordinate of the delay
    vectors). Vectors are paired by time, see :func:`csd.embedding.paired_rows`.
    """
    station_cum = np.asarray(station_cum, dtype="float64")
    model_cum = np.asarray(model_cum, dtype="float64")
    n = min(len(station_cum), len(model_cum))

    lags = np.arange(max(max_days // step, 1))
    amps = amplitude_grid(station_cum[:n], model_cum[:n], n_amplitudes)
    sg = np.full((amps.size, lags.size), np.nan)
    for j, lag in enumerate(lags):
        for i, amp in enumerate(amps):
            sg[i, j] = similarity_function(station_cum[:n], model_cum[:n], int(lag), float(amp))
    row, col = np.unravel_index(np.nanargmin(sg), sg.shape)
    lag, amp = int(lags[col]), float(amps[row])

    sta_vec = delay_vectors(station_cum, dim_station, delay_station)
    gcm_vec = delay_vectors(model_cum, dim_model, delay_model)
    gcm_rows, sta_rows = paired_rows(len(sta_vec), len(gcm_vec), (dim_station - 1) * delay_station,
                                     (dim_model - 1) * delay_model, lag)
    ss, u = mutual_false_neighbors(sta_vec[sta_rows.start : sta_rows.stop],
                                   gcm_vec[gcm_rows.start : gcm_rows.stop] * amp)
    return Calibration(lag, amp, ss, u, sg, amps, lags)
