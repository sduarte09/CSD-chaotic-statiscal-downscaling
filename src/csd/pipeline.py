"""End-to-end downscaling of one station: calibration, validation and application."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

from .calibration import Calibration, calibrate
from .embedding import sum_accumulate
from .synchronization import rmse, synchronize


@dataclass
class DownscalingResult:
    step: int
    calibration: Calibration
    validation_rmse: Optional[float]
    validation: pd.Series  # downscaled accumulated series, validation period
    observed_validation: pd.Series  # station accumulated over the same blocks
    downscaled: pd.Series  # application period, indexed by the first day of each block


def _block_starts(index, n_blocks, step):
    return index[: n_blocks * step : step]


def downscale_station(station, model, step, dims, delays, cal_end, val_end=None, app_end=None,
                      n_fs_validation=30, n_fs_application=25):
    """Run the CSD method for one station.

    Parameters
    ----------
    station, model : daily ``pandas.Series`` (precipitation, same units) with DatetimeIndex
    step : accumulation step in days
    dims, delays : ``((dim_station, dim_model), (delay_station, delay_model))`` as given by
        :func:`csd.chaos.chaos_analysis`
    cal_end : last day of the calibration period. Application starts the day after the
        last complete accumulation block that ends on or before this date.
    val_end : last day of the validation period (defaults to the end of the station record)
    app_end : last day of the series to downscale (defaults to the end of the model series)

    Notes
    -----
    Station and model are cut to their common start. Accumulation uses block sums on
    complete blocks only. Validation runs the same synchronization over the period for
    which observations exist and compares it with the *observed* accumulated series;
    ``Validation.m`` compared with the calibration series, which looks like an error.
    """
    (dim_s, dim_m), (tau_s, tau_m) = dims, delays
    start = max(station.index[0], model.index[0])
    station, model = station[start:], model[start:]
    cal_end = pd.Timestamp(cal_end)
    n_cal_days = int((station.index <= cal_end).sum())
    n_cal_blocks = n_cal_days // step
    n_cal_days = n_cal_blocks * step
    if n_cal_blocks <= (max(dim_s, dim_m) - 1) * max(tau_s, tau_m) + 2:
        raise ValueError("Calibration period too short for this step, dimension and delay")

    sta_cal = sum_accumulate(station.values[:n_cal_days], step)
    gcm_cal = sum_accumulate(model.reindex(station.index[:n_cal_days]).values, step)
    cal = calibrate(sta_cal, gcm_cal, dim_s, tau_s, dim_m, tau_m, step)

    # validation: observed period right after calibration
    val_end = pd.Timestamp(val_end) if val_end is not None else station.index[-1]
    obs = station.iloc[n_cal_days:][:val_end]
    n_val_blocks = len(obs) // step
    validation = pd.Series(dtype=float)
    observed = pd.Series(dtype=float)
    err = None
    if n_val_blocks > 0:
        n_val_days = n_val_blocks * step
        dates = obs.index[:n_val_days]
        model_val = sum_accumulate(model.reindex(dates).values, step)
        sim = synchronize(sta_cal, gcm_cal, model_val, cal, dim_s, tau_s, dim_m, tau_m, n_fs_validation)
        starts = dates[::step]
        observed = pd.Series(sum_accumulate(obs.values[:n_val_days], step), index=starts)
        validation = pd.Series(sim.series[:n_val_blocks], index=starts[: len(sim.series)])
        err = rmse(validation.values, observed.values)

    # application: model from the end of calibration up to app_end
    app_end = pd.Timestamp(app_end) if app_end is not None else model.index[-1]
    app_dates = model.index[(model.index > station.index[n_cal_days - 1]) & (model.index <= app_end)]
    n_app_blocks = len(app_dates) // step
    app_dates = app_dates[: n_app_blocks * step]
    model_app = sum_accumulate(model.reindex(app_dates).values, step)
    app = synchronize(sta_cal, gcm_cal, model_app, cal, dim_s, tau_s, dim_m, tau_m, n_fs_application)
    starts = app_dates[::step]
    downscaled = pd.Series(app.series[: len(starts)], index=starts[: len(app.series)], name="precipitation")
    return DownscalingResult(step, cal, err, validation, observed, downscaled)
