import numpy as np
import pytest


def logistic(x0, n, r=3.9):
    x = np.empty(n)
    x[0] = x0
    for i in range(n - 1):
        x[i + 1] = r * x[i] * (1 - x[i])
    return x


@pytest.fixture
def chaotic_pair():
    """Two independent trajectories of the same chaotic map, with different scale."""
    return logistic(0.3, 900) * 100, logistic(0.37, 900) * 60 + 5
