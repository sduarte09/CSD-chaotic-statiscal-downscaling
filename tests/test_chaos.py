import sys
from unittest import mock

import numpy as np
import pytest

from csd import chaos


def _sine(n=1500, period=50):
    return np.sin(2 * np.pi * np.arange(n) / period)


def _fake_nolds(spectrum):
    fake = mock.MagicMock()
    fake.lyap_e.return_value = np.array(spectrum)
    return fake


def test_classify_lyapunov_covers_the_four_cases():
    assert chaos.classify_lyapunov([0.5, -2.0]) == "Deterministic_Chaos"
    assert chaos.classify_lyapunov([0.5, 0.2]) == "Noise_(Random)"
    assert chaos.classify_lyapunov([0.0, -1.0]) == "Stable_Fixed_Point"
    assert chaos.classify_lyapunov([-0.1, -1.0]) == "Stable_Limit_Cycle"


def test_csd_delay_is_the_snlda_delay_plus_one():
    from snlda import dynamics

    x = _sine()
    expected = min(dynamics.autocorrelation_delay(x), dynamics.mutual_information_delay(x)) + 1
    assert chaos.csd_delay(x) == expected


def test_chaos_analysis_uses_the_first_passing_dimension_and_the_lyapunov_settings():
    fake = _fake_nolds([0.3, -0.9])
    x = _sine()
    with mock.patch.dict(sys.modules, {"nolds": fake}):
        res = chaos.chaos_analysis(x)
    assert res.is_chaotic
    assert res.dimension >= 1
    kwargs = fake.lyap_e.call_args.kwargs
    assert kwargs["matrix_dim"] == res.dimension
    assert kwargs["emb_dim"] == 2 * res.dimension - 1
    assert kwargs["tau"] == res.lag


def test_chaos_analysis_not_chaotic_when_the_spectrum_does_not_qualify():
    fake = _fake_nolds([0.3, 0.1])  # positive sum: noise
    with mock.patch.dict(sys.modules, {"nolds": fake}):
        res = chaos.chaos_analysis(_sine())
    assert not res.is_chaotic


def test_chaos_analysis_without_a_passing_dimension():
    fake = _fake_nolds([0.3, -0.9])
    noise = np.random.default_rng(0).random(400)
    with mock.patch.dict(sys.modules, {"nolds": fake}):
        res = chaos.chaos_analysis(noise, fnn_threshold=-1.0)
    assert res.dimension is None and res.lyapunov is None and not res.is_chaotic
    fake.lyap_e.assert_not_called()


def test_clear_error_when_nolds_cannot_be_imported():
    with mock.patch.dict(sys.modules, {"nolds": None}):
        with pytest.raises(ImportError):
            chaos.chaos_analysis(_sine(600, 40))


def _result(chaotic):
    return chaos.ChaosResult(lag=2, dimension=3, lyapunov=np.array([0.1, -0.5]), is_chaotic=chaotic)


def test_find_common_chaotic_step_returns_the_first_step_where_both_are_chaotic():
    calls = []

    def fake_analysis(series, **kwargs):
        calls.append(len(series))
        # station series are the long ones (identical lengths per step), model has 60 extra days
        return _result(chaotic=(len(series) in station_ok or len(series) in model_ok))

    station = np.ones(300)
    model = np.ones(360)
    station_ok = {300 // 5, 300 // 10}  # station is chaotic at steps 5 and 10
    model_ok = {360 // 10}  # model only at step 10
    with mock.patch.object(chaos, "chaos_analysis", side_effect=fake_analysis):
        selection = chaos.find_common_chaotic_step(station, model, [3, 5, 10, 15])
    assert selection.step == 10
    assert selection.station.is_chaotic and selection.model.is_chaotic


def test_find_common_chaotic_step_returns_none_and_reuses_the_cache():
    cache = {}
    with mock.patch.object(chaos, "chaos_analysis", return_value=_result(False)):
        assert chaos.find_common_chaotic_step(np.ones(100), np.ones(100), [5, 10], cache) is None
    assert cache == {}  # the model is never analysed when the station is not chaotic
    with mock.patch.object(chaos, "chaos_analysis", return_value=_result(True)) as analysis:
        first = chaos.find_common_chaotic_step(np.ones(100), np.ones(100), [5, 10], cache)
        assert 5 in cache and first.step == 5
        analysis.reset_mock()
        chaos.find_common_chaotic_step(np.ones(100), np.ones(100), [5, 10], cache)
        assert analysis.call_count == 1  # only the station; the model result came from the cache
