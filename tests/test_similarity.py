import numpy as np
import pytest

from csd import similarity as sim


def test_amplitude_ratio():
    assert sim.amplitude_ratio([2.0, 4.0], [1.0, 2.0]) == pytest.approx(2.0)


def test_similarity_function_matches_the_original_expression():
    rng = np.random.default_rng(0)
    station = rng.random(50)
    model = rng.random(50)
    lag = 3
    n = len(station)
    # Expression as written in the original script, with its extra (n - lag) factors
    a = station[1 + lag : n + 1] - np.mean(station[1 + lag : n + 1])
    b = station[: n - lag - 1] - np.mean(station[: n - lag - 1])
    c = model[: n - lag - 1] - np.mean(model[: n - lag - 1])
    original = (
        (np.mean((a - c) ** 2) / (n - lag))
        / (np.mean(b**2) * np.mean(c**2) / (n - lag) ** 2) ** 0.5
    ) ** 0.5
    assert sim.similarity_function(station, model, lag) == pytest.approx(original)


def test_similarity_is_zero_for_identical_series_at_lag_zero():
    x = np.random.default_rng(1).random(40)
    # at lag 0 the compared windows are x[1:] and x[:-1], so use a constant-step ramp
    ramp = np.arange(40.0)
    assert sim.similarity_function(ramp, ramp, 0) == pytest.approx(0.0, abs=1e-12)
    assert sim.similarity_function(x, x, 0) > 0  # random series shifted by one step differ


def test_best_lag_recovers_a_known_shift():
    t = np.arange(400)
    station = np.sin(2 * np.pi * t / 40)
    model = np.roll(station, -6)  # model[i] = station[i + 6]
    # the comparison pairs station[j + 1 + lag] with model[j], so the match is at lag = 6 - 1
    lag, values = sim.best_lag(station, model, step=5, max_days=180)
    assert lag == 5
    assert np.nanmin(values) == values[lag]


def test_best_lag_number_of_candidate_lags():
    station = np.random.default_rng(2).random(100)
    _, values = sim.best_lag(station, station, step=10, max_days=180)
    assert len(values) == 18


def test_mutual_false_neighbors_hand_computed_case():
    model = np.array([[0.0], [1.0], [10.0]])
    station = np.array([[0.0], [5.0], [1.0]])
    ss, u = sim.mutual_false_neighbors(station, model)
    # t=0: model neighbour j=1 (distance 1), station distance to it |0-5| = 5;
    #      station neighbour j=2 (distance 1), model distance to it |0-10| = 10
    assert ss[0] == pytest.approx((5 * 10) / (1 * 1))
    assert u[0] == pytest.approx(ss[0])
    assert u[-1] == pytest.approx(np.mean(ss))


def test_identical_spaces_give_ratio_one():
    x = np.random.default_rng(3).random((30, 3))
    ss, u = sim.mutual_false_neighbors(x, x.copy())
    np.testing.assert_allclose(ss, 1.0)
    np.testing.assert_allclose(u, 1.0)


def test_mutual_false_neighbors_ignores_duplicate_vectors():
    model = np.array([[0.0], [0.0], [2.0], [5.0]])  # vectors 0 and 1 are identical
    station = np.array([[0.0], [3.0], [7.0], [1.0]])
    ss, _ = sim.mutual_false_neighbors(station, model)
    # t=0: nearest different model vector is index 2 (distance 2), not the duplicate at index 1;
    #      nearest station vector is index 3 (distance 1)
    #      ss = (|sta0 - sta2| * |mod0 - mod3|) / (2 * 1) = (7 * 5) / 2
    assert ss[0] == pytest.approx(17.5)


def test_mutual_false_neighbors_checks_lengths():
    with pytest.raises(ValueError):
        sim.mutual_false_neighbors(np.ones((3, 2)), np.ones((4, 2)))
