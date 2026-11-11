"""CSD: Chaotic Statistical Downscaling.

Tested building blocks of the CSD method (Duarte et al., 2018):

* :mod:`csd.accumulation` - mean accumulation of daily series, random disaggregation,
  monthly and annual aggregates
* :mod:`csd.chaos`        - time delay, embedding dimension, Lyapunov exponents and the
  search for an accumulation step at which predictor and predictand are both chaotic
* :mod:`csd.similarity`   - amplitude and lag calibration (generalized similarity
  function) and the mutual false nearest neighbour parameters

The original MATLAB and Python 2 implementations are kept unchanged in ``legacy/``.
"""
from . import accumulation, chaos, similarity

__version__ = "0.1.0"

__all__ = ["accumulation", "chaos", "similarity", "__version__"]
