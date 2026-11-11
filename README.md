![tests](https://github.com/sduarte09/chaotic-statiscal-downscaling/actions/workflows/tests.yml/badge.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24133815.svg)](https://doi.org/10.5281/zenodo.24133815)

# CSD: Chaotic Statistical Downscaling

CSD is a statistical downscaling method that links a climate model and a local rain gauge through three steps: phase-space reconstruction for different accumulation time steps, identification of deterministic chaos, and a predictive model based on the general synchronization of two dynamical systems. It was tested in the Bogotá River Basin (Colombia) with 47 rain gauges and the MPI-ESM-MR climate model.

> Duarte, F., Corzo, G., Santos, G., Hernández, O. (2018). Chaotic Statistical Downscaling (CSD): Application and Comparison in the Bogotá River Basin. In *HIC 2018. 13th International Conference on Hydroinformatics*, EPiC Series in Engineering, vol. 3, pp. 626-634. [https://doi.org/10.29007/wkcx](https://doi.org/10.29007/wkcx)

## Status

This repository holds the original code and a new Python 3 package that is being built step by step.

| Part | Status |
|---|---|
| Accumulation, disaggregation, monthly/annual aggregates (`csd.accumulation`) | tested |
| Delay, embedding dimension, Lyapunov, choice of accumulation step (`csd.chaos`) | tested |
| Amplitude and lag calibration, mutual false nearest neighbours (`csd.similarity`) | tested |
| Synchronization (prediction) and the full per-station pipeline | not yet ported; use `legacy/` |

The unported parts need the MPI-ESM-MR NetCDF files, which are not in this repository. See `docs/LEGACY_NOTES.md` for what was found in the original code and which choices the package makes.

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

## Repository layout

```
src/csd/      Python 3 package
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
