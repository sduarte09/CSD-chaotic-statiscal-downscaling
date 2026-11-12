"""NetCDF input and output for CSD.

Reading uses ``netCDF4`` when installed and falls back to ``scipy.io.netcdf_file``
(classic NetCDF3 files such as the CMIP5 daily files). Writing always produces a
classic NetCDF file with scipy, which every NetCDF tool can read.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

SECONDS_PER_DAY = 86400.0  # kg m-2 s-1 -> mm day-1


class _Reader:
    """Minimal common interface over netCDF4.Dataset and scipy's netcdf_file."""

    def __init__(self, path):
        try:
            import netCDF4

            self._ds = netCDF4.Dataset(path)
            self._kind = "netcdf4"
        except ImportError:
            from scipy.io import netcdf_file

            self._ds = netcdf_file(path, mmap=False)
            self._kind = "scipy"

    def array(self, name, index=None):
        var = self._ds.variables[name]
        data = var[:] if index is None else var[index]
        return np.ma.filled(np.ma.asarray(data, dtype="float64"), np.nan) if self._kind == "netcdf4" else np.array(data, dtype="float64")

    def attr(self, name, attr, default=None):
        var = self._ds.variables[name]
        value = getattr(var, attr, default)
        return value.decode() if isinstance(value, bytes) else value

    def fill_value(self, name):
        for key in ("_FillValue", "missing_value"):
            value = getattr(self._ds.variables[name], key, None)
            if value is not None:
                return float(np.ravel(value)[0])
        return None

    def close(self):
        self._ds.close()


def _mask_fill(data, fill):
    """NaN where ``data`` equals the fill value (compared with a tolerance, since the fill
    attribute is often stored in single precision)."""
    if fill is not None:
        data[np.abs(data - fill) <= 1e-6 * abs(fill)] = np.nan
    return data


def decode_time(values, units, calendar="standard", midday_offset=True):
    """Daily NetCDF time values to a DatetimeIndex of calendar days.

    Only calendars with a pandas equivalent are supported (standard, gregorian,
    proleptic_gregorian). Daily CMIP files stamp days at 12:00, so 0.5 day is removed
    as in the MATLAB code (``time - 0.5``).
    """
    if calendar not in ("standard", "gregorian", "proleptic_gregorian"):
        raise ValueError(f"Unsupported calendar {calendar!r}; convert the file to a gregorian calendar first")
    unit, _, origin = units.partition(" since ")
    if unit.strip() != "days":
        raise ValueError(f"Only 'days since ...' time units are supported, got {units!r}")
    values = np.asarray(values, dtype="float64") - (0.5 if midday_offset else 0.0)
    return pd.DatetimeIndex(pd.Timestamp(origin.strip()) + pd.to_timedelta(np.floor(values + 1e-6), unit="D"))


def find_cell(bounds, value):
    """Index of the first interval ``[lo, hi]`` of ``bounds`` (n, 2) that contains ``value``."""
    bounds = np.asarray(bounds, dtype="float64")
    inside = np.flatnonzero((bounds[:, 0] <= value) & (value <= bounds[:, 1]))
    if inside.size == 0:
        raise ValueError(f"{value} lies outside the grid bounds")
    return int(inside[0])


def read_model_cell(path, latitude, longitude, variable="pr", scale=SECONDS_PER_DAY):
    """Daily series of the model grid cell that contains a point.

    Negative longitudes are converted to 0-360. ``scale`` converts the units; the default
    turns precipitation flux (kg m-2 s-1) into mm/day. Returns a ``pandas.Series`` named
    after the variable, plus ``(lat_index, lon_index)`` and the cell bounds in ``.attrs``.
    """
    longitude = longitude + 360.0 if longitude < 0 else longitude
    reader = _Reader(path)
    try:
        i = find_cell(reader.array("lat_bnds"), latitude)
        j = find_cell(reader.array("lon_bnds"), longitude)
        data = reader.array(variable, (slice(None), i, j))
        data = _mask_fill(data, reader.fill_value(variable))
        index = decode_time(reader.array("time"), reader.attr("time", "units"),
                            reader.attr("time", "calendar", "standard"))
        lat_b, lon_b = reader.array("lat_bnds")[i], reader.array("lon_bnds")[j]
    finally:
        reader.close()
    series = pd.Series(data * scale, index=index, name=variable)
    series.attrs.update(cell=(i, j), lat_bounds=tuple(lat_b), lon_bounds=tuple(lon_b))
    return series


def concat_model(historic, future):
    """Join a historic and a scenario series; the scenario wins where dates overlap."""
    joined = pd.concat([historic, future])
    return joined[~joined.index.duplicated(keep="last")].sort_index()


def write_netcdf(path, series, name="pr", units="mm", long_name="CSD downscaled precipitation",
                 history="Created with csd"):
    """Write a series with a DatetimeIndex to a classic NetCDF file (``days since 1850-01-01``)."""
    from scipy.io import netcdf_file

    index = pd.DatetimeIndex(series.index)
    days = (index - pd.Timestamp("1850-01-01")) / pd.Timedelta(days=1)
    with netcdf_file(path, "w") as f:
        f.createDimension("time", len(series))
        t = f.createVariable("time", "d", ("time",))
        t[:] = np.asarray(days, dtype="float64")
        t.units, t.calendar, t.standard_name = "days since 1850-01-01 00:00:00", "proleptic_gregorian", "time"
        v = f.createVariable(name, "d", ("time",))
        v[:] = np.nan_to_num(np.asarray(series, dtype="float64"), nan=1e20)
        v.units, v.long_name, v._FillValue = units, long_name, 1e20
        f.history = history


def read_netcdf_series(path, name="pr"):
    """Read back a series written by :func:`write_netcdf`."""
    reader = _Reader(path)
    try:
        values = reader.array(name)
        values = _mask_fill(values, reader.fill_value(name))
        index = decode_time(reader.array("time"), reader.attr("time", "units"),
                            reader.attr("time", "calendar", "standard"), midday_offset=False)
    finally:
        reader.close()
    return pd.Series(values, index=index, name=name)


def _read_table(path, **kwargs):
    path = str(path)
    if path.lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(path, **kwargs)  # needs openpyxl
    return pd.read_csv(path, **kwargs)


def read_station_table(path, date_column=0, dayfirst=True):
    """Daily station table: one date column (``dd/mm/yyyy`` by default) and one column per station.

    Returns a DataFrame indexed by date, sorted, with the station names as columns. CSV and
    Excel files are accepted.
    """
    table = _read_table(path)
    dates = pd.to_datetime(table.iloc[:, date_column], dayfirst=dayfirst)
    table = table.drop(columns=table.columns[date_column])
    table.index = pd.DatetimeIndex(dates, name="date")
    return table.apply(pd.to_numeric, errors="coerce").sort_index()


def read_catalog(path):
    """Station catalog: station name, latitude, longitude (in that order, or found by name).

    Returns a DataFrame indexed by station name with ``latitude`` and ``longitude`` columns.
    """
    table = _read_table(path)
    names = table.columns[0]
    lower = {c: str(c).lower() for c in table.columns}
    lat = next((c for c, l in lower.items() if l.startswith("lat")), table.columns[1])
    lon = next((c for c, l in lower.items() if l.startswith("lon")), table.columns[2])
    out = pd.DataFrame({"latitude": pd.to_numeric(table[lat]), "longitude": pd.to_numeric(table[lon])})
    out.index = table[names].astype(str).str.strip()
    return out
