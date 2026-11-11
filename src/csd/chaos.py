"""Chaos analysis of rainfall series: delay, embedding dimension and Lyapunov exponents.

The phase-space tools come from the ``snlda`` package; this module adds the choices
specific to CSD (delay definition, false-neighbour threshold, Lyapunov settings).
"""
from __future__ import annotations

import math
from typing import NamedTuple, Optional

import numpy as np
from snlda import dynamics

from .accumulation import mean_accumulate


class ChaosResult(NamedTuple):
    lag: int
    dimension: Optional[int]  # None if no dimension passes the false-neighbour threshold
    lyapunov: Optional[np.ndarray]
    is_chaotic: bool


class StepSelection(NamedTuple):
    step: int
    station: ChaosResult
    model: ChaosResult


def _import_nolds():
    """Import ``nolds`` lazily, with a clear message if it cannot be imported."""
    try:
        import nolds
    except (ImportError, TypeError) as error:
        raise ImportError(
            "nolds could not be imported in this environment "
            f"({type(error).__name__}: {error}). Use Python 3.12 or newer."
        ) from error
    return nolds


def csd_delay(x):
    """Time delay used by CSD: the smaller of the autocorrelation and mutual-information
    delays, plus one (as in the original Python implementation)."""
    return min(dynamics.autocorrelation_delay(x), dynamics.mutual_information_delay(x)) + 1


def classify_lyapunov(lyapunov):
    """Type of motion from a Lyapunov spectrum (MATLAB implementation rules)."""
    lyapunov = np.asarray(lyapunov, dtype="float64")
    total, largest = lyapunov.sum(), lyapunov.max()
    if total < 0 and largest > 0:
        return "Deterministic_Chaos"
    if total > 0 and largest > 0:
        return "Noise_(Random)"
    if total < 0 and largest == 0:
        return "Stable_Fixed_Point"
    return "Stable_Limit_Cycle"


def chaos_analysis(x, fnn_threshold=0.15):
    """Delay, embedding dimension and Lyapunov spectrum of a series.

    The embedding dimension is the smallest dimension whose false-nearest-neighbour
    fraction is at or below ``fnn_threshold``. The series is chaotic if the sum of the
    Lyapunov exponents is negative and the largest is positive.
    """
    nolds = _import_nolds()
    x = np.asarray(x, dtype="float64")
    lag = csd_delay(x)
    max_dim = 2 * int(math.floor(2 * math.log10(x.size))) + 1
    _, fnn = dynamics.false_nearest_neighbors(x, lag, min_dims=1, max_dims=max_dim)
    passing = np.flatnonzero(fnn <= fnn_threshold)
    if passing.size == 0:
        return ChaosResult(int(lag), None, None, False)
    dim = int(passing[0]) + 1
    lyapunov = np.asarray(nolds.lyap_e(x, emb_dim=2 * dim - 1, matrix_dim=dim, tau=lag))
    chaotic = bool(lyapunov.sum() < 0 and lyapunov.max() > 0)
    return ChaosResult(int(lag), dim, lyapunov, chaotic)


def find_common_chaotic_step(station, model, steps, model_cache=None, **kwargs):
    """First accumulation step at which both series are deterministic chaos.

    ``station`` and ``model`` are daily series. Steps are tried in the order given.
    The model is only analysed when the station is chaotic at that step. Pass a dict as
    ``model_cache`` to reuse model results between stations that share a model cell.
    Returns a :class:`StepSelection`, or ``None`` if no step qualifies.
    """
    cache = {} if model_cache is None else model_cache
    for step in steps:
        station_result = chaos_analysis(mean_accumulate(station, step), **kwargs)
        if not station_result.is_chaotic:
            continue
        if step not in cache:
            cache[step] = chaos_analysis(mean_accumulate(model, step), **kwargs)
        if cache[step].is_chaotic:
            return StepSelection(step, station_result, cache[step])
    return None
