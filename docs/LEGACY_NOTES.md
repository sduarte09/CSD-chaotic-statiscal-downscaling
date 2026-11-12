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

## Choices in the Python 3 port of Calibration, Validation and Synchronization

- **Delay vectors** follow the MATLAB order: first column is the most recent value, row `r` is `[x[r+(D-1)τ], ..., x[r]]`.
- **Time pairing.** Station and model vectors are paired by the time of their first coordinate plus the calibrated lag (`csd.embedding.paired_rows`), not by row number. Row-number pairing drifts when the two embeddings differ in dimension or delay.
- **Similarity function.** The lag shifts the station series against the model series; the amplitude enters the normalisation once, as in `Calibration.m`. The earlier Python 2 version was off by one in the lag alignment.
- **Application vectors** are built from the concatenated calibration + application model series, so no vector straddles a gap at the junction. The calibration amplitude is applied to the model vectors.
- **Synchronization** solves the quadratic in the first coordinate over the FS grid (mean and standard deviation, ddof=1, of the second half of `ss`), picks the candidate with minimum distance, updates `ss` and `u`, and appends the new value. If no root is admissible, the nearest-neighbour value is used.
- **Validation** compares the downscaled series with the *observed* accumulated series over the validation period. `Validation.m` compares against the calibration series, which looks like an error.
- **Accumulation** uses complete blocks only (all-NaN blocks stay NaN). The application period starts the day after the last complete block of the calibration period.
- **NetCDF.** `pr` is multiplied by 86400 (mm/day); model times are at 12:00 and are shifted to the calendar day; the 1e20 fill value is masked with a tolerance because it is stored as float32.
- **Chaos is not reproducible pointwise.** For independent chaotic series only the distribution (mean, spread) is expected to be preserved; pointwise recovery happens only when station and model series are identical.

## Reproduction target

The paper reports RMSE = 73.37 for the Bogotá case. Running the legacy code and the Python port on the same Bogotá
data and comparing the results is the missing regression test; it has not been done.
