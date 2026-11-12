import numpy as np
import pandas as pd
import pytest
from scipy.io import netcdf_file

from csd import io


def _write_fake_model(path, n_days=10, fill=1e20):
    lat_b = np.array([[-90, 0], [0, 90.0]])
    lon_b = np.array([[-0.9375, 0.9375], [0.9375, 2.8125], [2.8125, 4.6875]])
    with netcdf_file(path, "w") as f:
        f.createDimension("time", n_days); f.createDimension("lat", 2)
        f.createDimension("lon", 3); f.createDimension("bnds", 2)
        t = f.createVariable("time", "d", ("time",)); t[:] = np.arange(n_days) + 56978.5
        t.units, t.calendar = "days since 1850-1-1 00:00:00", "proleptic_gregorian"
        for name, b in (("lat_bnds", lat_b), ("lon_bnds", lon_b)):
            v = f.createVariable(name, "d", (name[:3], "bnds")); v[:] = b
        pr = f.createVariable("pr", "f", ("time", "lat", "lon"))
        data = np.full((n_days, 2, 3), 1e-5, dtype="float32")
        data[:, 1, 1] = 2e-5
        data[3, 1, 1] = fill
        pr[:] = data
        pr._FillValue = np.float32(fill)


def test_find_cell_uses_bounds():
    bounds = np.array([[-90, 0], [0, 90.0]])
    assert io.find_cell(bounds, 10) == 1
    with pytest.raises(ValueError):
        io.find_cell(bounds, 100)


def test_decode_time_removes_the_midday_offset():
    idx = io.decode_time([56978.5, 56979.5], "days since 1850-1-1 00:00:00", "proleptic_gregorian")
    assert list(idx) == [pd.Timestamp("2006-01-01"), pd.Timestamp("2006-01-02")]


def test_decode_time_rejects_unsupported_calendars():
    with pytest.raises(ValueError):
        io.decode_time([0.5], "days since 1850-1-1", "360_day")


def test_read_model_cell_converts_units_wraps_longitude_and_masks_fill(tmp_path):
    path = tmp_path / "pr.nc"
    _write_fake_model(path)
    s = io.read_model_cell(str(path), latitude=10.0, longitude=1.5)
    assert s.attrs["cell"] == (1, 1)
    assert s.iloc[0] == pytest.approx(2e-5 * 86400, rel=1e-4)  # kg m-2 s-1 -> mm/day
    assert np.isnan(s.iloc[3])
    assert s.index[0] == pd.Timestamp("2006-01-01")
    west = io.read_model_cell(str(path), latitude=10.0, longitude=-358.5)  # 1.5 degrees east
    assert west.attrs["cell"] == (1, 1)


def test_concat_model_prefers_the_scenario_where_dates_overlap():
    a = pd.Series([1.0, 2.0], index=pd.to_datetime(["2005-12-30", "2005-12-31"]))
    b = pd.Series([20.0, 3.0], index=pd.to_datetime(["2005-12-31", "2006-01-01"]))
    out = io.concat_model(a, b)
    assert list(out) == [1.0, 20.0, 3.0]


def test_netcdf_roundtrip_keeps_dates_values_and_gaps(tmp_path):
    s = pd.Series([1.5, np.nan, 3.0], index=pd.date_range("2030-01-01", periods=3, freq="5D"), name="pr")
    io.write_netcdf(str(tmp_path / "o.nc"), s)
    back = io.read_netcdf_series(str(tmp_path / "o.nc"))
    assert (back.index == s.index).all()
    np.testing.assert_allclose(back.values, s.values)


def test_read_station_table_and_catalog(tmp_path):
    (tmp_path / "st.csv").write_text("Fecha,A,B\n01/02/2000,1.5,\n02/02/2000,0,3\n", encoding="utf-8")
    t = io.read_station_table(tmp_path / "st.csv")
    assert t.index[0] == pd.Timestamp("2000-02-01") and np.isnan(t.loc["2000-02-01", "B"])
    assert list(t.columns) == ["A", "B"] and t.loc["2000-02-02", "B"] == 3

    (tmp_path / "cat.csv").write_text("Name,Latitude,Longitude,Elev\nA ,4.6,-74.1,2600\nB,5,-75,10\n", encoding="utf-8")
    c = io.read_catalog(tmp_path / "cat.csv")
    assert c.loc["A", "latitude"] == 4.6 and c.loc["B", "longitude"] == -75


def test_run_notebook_is_valid_json():
    import json
    from pathlib import Path

    nb = json.loads((Path(__file__).parents[1] / "notebooks" / "run_csd.ipynb").read_text(encoding="utf-8"))
    assert nb["nbformat"] == 4 and any(c["cell_type"] == "code" for c in nb["cells"])
