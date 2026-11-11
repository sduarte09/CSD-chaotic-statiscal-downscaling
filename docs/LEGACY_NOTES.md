# Notes on the legacy code

Found while reading the original MATLAB and Python 2 implementations. Nothing in `legacy/`
was changed. The package in `src/csd` makes the choices described below; check them
against your own results before relying on them.

## Differences between the two legacy versions

| | MATLAB | Python 2 |
|---|---|---|
| Accumulation of daily data | sum (`nansum`) | mean (`prom_data`) |
| Time delay | `min(autocorrelation, mutual information)` | the same, plus 1 |
| False-neighbour threshold | 0.1 | 0.15 |
| Lyapunov exponents | TISEAN `lyap_spec` | `nolds.lyap_e` |
| Lag search range | `floor(30 / step)` | `floor(180 / step)` |
| Amplitude correction | grid search over a range of ratios | ratio of calibration means |

The package follows the Python version for the chaos analysis (it uses `nolds`) and the
MATLAB version for the mutual false nearest neighbours (see below).

## Issues in `CSD.py` (Python 2)

- Python 2 syntax: print statement, `open(..., 'rb')` for CSV, `keys()[-1]`, integer
  division in a slice, a generator passed to `np.vstack`. It cannot run on Python 3.
- It runs at import (no functions to call from tests) and has hard-coded paths.
- `deprom_data` draws weights with an unseeded `random`, so results are not reproducible.
  `csd.accumulation.disaggregate` takes a seed.
- Mutual false nearest neighbours: `dist_sta_sta` is taken from `distanceA1`, which is the
  distance in the model space, where the MATLAB version uses the station space. After
  `np.delete(distanceA, t)` the neighbour index is also shifted by one for neighbours after
  position `t`. `csd.similarity.mutual_false_neighbors` follows the MATLAB version
  (nearest neighbour at non-zero distance, correct indices).
- The synchronization loop uses `t`, left over from the previous loop.
- `next` statements do nothing (`continue` was probably intended) and `his.close` /
  `rcp.close` are missing parentheses.

## Issues in the MATLAB code

- `Dimension.m`: the `else` branch uses `mmax`, which is not defined (the variable is `MMAX`).
- `CSD.m`: `for z=2:2` processes only one station (looks like a debugging leftover);
  `station_longitude(lon_change)=station_longitude+360` assigns a vector of the wrong size.
- External dependencies not in the repository (see `legacy/README.md`). TISEAN is, to my
  knowledge, GPL-licensed; check before redistributing its executables.

## Reproduction target

The paper reports RMSE = 73.37 for the Bogotá case. Running the legacy code on the data
in `data/` and checking that this value is reproduced is the natural regression test before
porting the synchronization step.
