# Legacy implementations

The original code behind the HIC 2018 paper, kept unchanged for reproducibility.

- `matlab/` - MATLAB implementation (`CSD.m` is the entry point). It needs files that are
  **not** in this repository: the auxiliary functions (`Fill_Data.m`, `cao_deneme_mod.m`),
  the TISEAN executables used through `system(...)` calls (`mutual`, `d2`, `false_nearest`,
  `lyap_spec`, `delay`, `stp`, `poincare`), `nctoolbox` and the hydrographic zoning
  shapefile. Paths use Windows separators.
- `CSD_python2.py` - Python 2 version. It runs as a script (no functions to import), has
  hard-coded Windows paths and needs the MPI-ESM-MR NetCDF files, which are not included.

See `docs/LEGACY_NOTES.md` for what was found when reading them.
