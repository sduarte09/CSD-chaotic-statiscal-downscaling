"""Synchronization of the model with the station attractor (validation and application).

For every new model vector ``Xn`` the method looks for the station vector ``Yn`` whose
first coordinate is the unknown of a quadratic equation, obtained by asking that the
mutual false neighbour ratio equals a value ``FS`` taken around its calibrated mean.
The same routine is used for validation (against observations) and for the application
period (future scenario). It follows ``Synchronization.m`` / ``Validation.m``.
"""
from __future__ import annotations

from typing import NamedTuple, Optional

import numpy as np

from .calibration import Calibration
from .embedding import delay_vectors, paired_rows, series_from_vectors


class SyncResult(NamedTuple):
    series: np.ndarray  # downscaled accumulated series for the application period
    ss: np.ndarray
    u: np.ndarray


def _fs_grid(ss, n_points):
    half = ss[len(ss) // 2 - 1 :] if len(ss) >= 2 else ss
    sso = np.nanmean(half)
    eps = np.nanstd(half, ddof=1) if np.sum(~np.isnan(half)) > 1 else 0.0
    fs_min = sso - eps if sso - eps > 0 else 0.0
    return np.linspace(fs_min, sso + eps, n_points + 1)


def _nanmin(a):
    return np.nan if np.isnan(a).all() else np.nanmin(a)


def _solve_step(x_new, gcm, sta, delay_sta, fs_values):
    """One synchronization step: returns Yn plus the neighbours used by the ss update."""
    n, dim = sta.shape
    dist_to_all = np.linalg.norm(gcm[: n + 1] - x_new, axis=1)  # last row is Xn itself
    dist_to_all[dist_to_all == 0] = np.inf
    j = int(np.argmin(dist_to_all))
    xnnd, ynnd = gcm[j], sta[j]
    d_xnnd = np.linalg.norm(x_new - xnnd)

    with np.errstate(divide="ignore", invalid="ignore"):
        disx = d_xnnd / np.linalg.norm(gcm[:n] - x_new, axis=1)
    base = np.zeros(dim)
    base[1:] = sta[n - delay_sta, : dim - 1]
    rest = base[1:]
    c2d = -2.0 * rest @ ynnd[1:] + np.sum(ynnd**2)  # same for every candidate row
    sta_sq = np.sum(sta**2, axis=1)
    sta_rest_dot = sta[:, 1:] @ rest
    rest_sq = np.sum(rest**2)

    best_dist, best_c, best_eps = 1e5, None, None
    with np.errstate(divide="ignore", invalid="ignore"):
        for fs in fs_values:
            factor = (fs * disx) ** 2
            a2 = 1.0 - factor
            b2 = 2.0 * factor * sta[:, 0] - 2.0 * ynnd[0]
            c2 = c2d + 2.0 * factor * sta_rest_dot - factor * sta_sq + (1.0 - factor) * rest_sq
            disc = b2**2 - 4.0 * a2 * c2
            root = np.sqrt(np.where(disc >= 0, disc, np.nan))
            cands = []
            for sign in (1.0, -1.0):
                c1 = (-b2 + sign * root) / (2.0 * a2)
                c1 = np.where(c1 < 0, np.nan, c1)
                diff = sta[:, 1:] - rest
                d = np.sqrt((sta[:, 0] - c1) ** 2 + np.sum(diff**2, axis=1))
                cands.append((_nanmin(d), c1, d))
            (m1, c_a, d_a), (m2, c_b, d_b) = cands
            m = np.nanmin([m1, m2]) if not (np.isnan(m1) and np.isnan(m2)) else np.nan
            chosen_c, chosen_d = (c_a, d_a) if m1 < m2 else (c_b, d_b)
            if m < best_dist:
                best_dist, best_c, best_eps = m, chosen_c, chosen_d
    if best_c is None:  # no admissible root for any FS: keep the station coordinate
        best_c = sta[:, 0].copy()
        diff = sta[:, 1:] - rest
        best_eps = np.sqrt(np.sum(diff**2, axis=1))
    order = np.argsort(best_eps, kind="stable")
    nonzero = [k for k in order if best_eps[k] != 0]
    k = int(nonzero[0] if nonzero else order[0])
    yn = np.concatenate([[best_c[k]], rest])
    return yn, xnnd, ynnd, k


def synchronize(station_cum, model_cal_cum, model_app_cum, calibration: Calibration,
                dim_station, delay_station, dim_model, delay_model, n_fs=25,
                scale_application=True):
    """Downscale the accumulated model series of the application period.

    Parameters
    ----------
    station_cum : accumulated station series of the calibration period
    model_cal_cum : accumulated model series over the same period
    model_app_cum : accumulated model series to downscale (continues ``model_cal_cum``)
    calibration : result of :func:`csd.calibration.calibrate`
    n_fs : number of intervals of the ``FS`` grid (25 in ``Synchronization.m``, 30 in
        ``Validation.m``)
    scale_application : multiply the application vectors by the calibrated amplitude.
        ``Synchronization.m`` does not (only the calibration vectors are scaled), while the
        Python 2 script scales the whole model series. Scaling is the consistent choice.

    Returns a :class:`SyncResult` whose ``series`` has one value per accumulation step.
    """
    lag, amp = calibration.lag, calibration.amplitude
    station_cum = np.asarray(station_cum, dtype="float64")
    n_station, n_model_cal = len(station_cum), len(model_cal_cum)
    off_s, off_g = (dim_station - 1) * delay_station, (dim_model - 1) * delay_model
    sta_all = delay_vectors(station_cum, dim_station, delay_station)
    # vectors over calibration + application so that none are lost at the junction
    gcm_all = delay_vectors(np.concatenate([model_cal_cum, model_app_cum]), dim_model, delay_model)

    gcm_rows, sta_rows = paired_rows(len(sta_all), n_model_cal - off_g, off_s, off_g, lag)
    n0 = len(sta_rows)
    if n0 != len(calibration.ss):
        raise ValueError("Calibration was made with a different series length or embedding")
    sta_lo = sta_rows.start
    n_app = len(model_app_cum)
    # station time n_station + v is paired with the model vector at time n_station + v - lag
    first = n_station - lag - off_g
    if first < 0:
        raise ValueError("Lag larger than the calibration record")
    gcm_app = gcm_all[first : first + n_app]
    if len(gcm_app) < n_app:
        raise ValueError("Model series too short for the application period")
    if scale_application:
        gcm_app = gcm_app * amp

    sta_buf = np.empty((n0 + n_app, dim_station))
    sta_buf[:n0] = sta_all[sta_rows.start : sta_rows.stop]
    gcm_buf = np.empty((n0 + n_app, dim_model))
    gcm_buf[:n0] = gcm_all[gcm_rows.start : gcm_rows.stop] * amp
    ss = np.concatenate([calibration.ss, np.full(n_app, np.nan)])
    u = np.concatenate([calibration.u, np.full(n_app, np.nan)])
    t = n0

    for v in range(n_app):
        n = n0 + v
        x_new = gcm_app[v]
        gcm_buf[n] = x_new
        fs_values = _fs_grid(ss[: t + v], n_fs)
        yn, xnnd, ynnd, k = _solve_step(x_new, gcm_buf[: n + 1], sta_buf[:n], delay_station, fs_values)
        xnnr, ynnr = gcm_buf[k], sta_buf[k]
        den = np.linalg.norm(yn - ynnr) * np.linalg.norm(x_new - xnnd)
        with np.errstate(divide="ignore", invalid="ignore"):
            ss[t + v] = (np.linalg.norm(yn - ynnd) * np.linalg.norm(x_new - xnnr)) / den
        u[t + v] = np.nansum(ss[: t + v + 1]) / (t + v + 1)
        sta_buf[n] = yn

    full = series_from_vectors(sta_buf, dim_station, delay_station)
    # full[i] is the station series at time sta_lo + i; the first new value is at n_station
    series = full[n_station - sta_lo :]
    return SyncResult(series, ss, u)


def rmse(simulated, observed):
    """Root mean squared error over the common length (NaN pairs ignored)."""
    n = min(len(simulated), len(observed))
    d = np.asarray(simulated[:n], dtype=float) - np.asarray(observed[:n], dtype=float)
    return float(np.sqrt(np.nanmean(d**2)))
