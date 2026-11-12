![tests](https://github.com/sduarte09/chaotic-statiscal-downscaling/actions/workflows/tests.yml/badge.svg)

# CSD: Chaotic Statistical Downscaling

CSD is a statistical downscaling method that links a climate model and a local rain gauge through three steps: phase-space reconstruction for different accumulation time steps, identification of deterministic chaos, and a predictive model based on the general synchronization of two dynamical systems. It was tested in the Bogotá River Basin (Colombia) with 47 rain gauges and the MPI-ESM-MR climate model.

> Duarte, F., Corzo, G., Santos, G., Hernández, O. (2018). Chaotic Statistical Downscaling (CSD): Application and Comparison in the Bogotá River Basin. In *HIC 2018. 13th International Conference on Hydroinformatics*, EPiC Series in Engineering, vol. 3, pp. 626-634. [https://doi.org/10.29007/wkcx](https://doi.org/10.29007/wkcx)

## Status

This repository holds the original MATLAB / Python 2 code (`legacy/`) and a Python 3 package (`src/csd`) that covers the whole chain from NetCDF input to NetCDF output.

| Part | Status |
|---|---|
| Accumulation, disaggregation, monthly/annual aggregates (`csd.accumulation`) | tested |
| Delay, embedding dimension, Lyapunov, choice of accumulation step (`csd.chaos`) | tested |
| Amplitude and lag calibration, mutual false nearest neighbours (`csd.similarity`) | tested |
| Calibration: lag, amplitude, mutual false neighbours (`csd.calibration`) | tested on synthetic data |
| Synchronization / prediction (`csd.synchronization`) | tested on synthetic data |
| NetCDF model cell reader and writer, station table and catalogue readers (`csd.io`) | tested; reader also checked on the MPI-ESM-MR file |
| Per-station pipeline: calibration, validation, application (`csd.pipeline`) | tested on synthetic data |
| Run notebook (`notebooks/run_csd.ipynb`) | runs end to end in demo mode |

**What has not been verified.** There is no reference output from the MATLAB code in this repository, so the Python port is checked against synthetic chaotic series (identical series are recovered; independent series keep the station's mean and spread), not against the paper. The paper's RMSE of 73.37 has not been reproduced yet; that needs the Bogotá station data run through both implementations. `csd.chaos` needs `snlda` and `nolds` and is tested in CI, not in every environment.

See `docs/LEGACY_NOTES.md` for what was found in the original code and which choices the package makes.

## Installation

CSD needs **Python 3.12 or newer**.

```bash
git clone https://github.com/sduarte09/chaotic-statiscal-downscaling.git
cd chaotic-statiscal-downscaling
pip install -e ".[test]"
```

This also installs [SNLDA](https://github.com/sduarte09/SNLDA), which provides the phase-space tools.

## Quick start

```python
import numpy as np
from csd import accumulation, chaos, similarity

# daily series -> 10-day means
accumulated = accumulation.mean_accumulate(daily_precipitation, 10)

# first step at which station and climate model are both deterministic chaos
selection = chaos.find_common_chaotic_step(station_daily, model_daily, steps=[5, 7, 10, 15, 30])

# lag that best aligns station and model, and the mutual false neighbour ratios
lag, sg = similarity.best_lag(station_accumulated, model_accumulated, step=selection.step)
ss, u = similarity.mutual_false_neighbors(station_vectors, model_vectors)
```

## Run a station

Open `notebooks/run_csd.ipynb`, set the paths to the model NetCDF file(s) and to the station table and catalogue, and run all cells. Without a station file it runs a clearly labelled demo with a synthetic station. The same steps in code:

```python
from csd import io, pipeline

model = io.read_model_cell("pr_day_MPI-ESM-MR_rcp45_r1i1p1_20060101-20091231.nc", latitude=4.6, longitude=-74.1)
station = io.read_station_table("data/stations_daily.csv")["station_name"].dropna()

result = pipeline.downscale_station(station, model, step=10, dims=(3, 3), delays=(2, 2), cal_end="2005-12-31")
print(result.calibration.lag, result.calibration.amplitude, result.validation_rmse)
io.write_netcdf("downscaled.nc", result.downscaled, name="pr", units="mm per 10 days")
```

Model precipitation is converted from kg m-2 s-1 to mm/day. Classic NetCDF files are read with SciPy; install the `netcdf` extra (`pip install -e ".[netcdf]"`) for NetCDF4/HDF5 files.

## Repository layout

```
src/csd/      Python 3 package
notebooks/    run_csd.ipynb
tests/        pytest tests on small synthetic examples
legacy/       original MATLAB and Python 2 implementations, unchanged
data/         Bogotá station data and catalogue
docs/         notes on the legacy code
```

## Tests

```bash
pytest
```

## License

MIT, see `LICENSE`.
