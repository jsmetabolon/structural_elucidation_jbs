"""Parity and cache tests for the parallel verifier scoring path.

The verifier optimization must be output-identical: worker count and cache
state may change wall-clock only, never scores, sources, or row order.
"""

from __future__ import annotations

import pandas as pd
import pytest

import structural_elucidation.nce.verifier_v2 as verifier_v2
from structural_elucidation.nce.score_cache import VerifierScoreCache
from structural_elucidation.nce.verifier_v2 import annotate_predicted_spectrum_scores_v2


@pytest.fixture()
def stub_scorer(monkeypatch):
    """Deterministic stand-in for the NCE matcher; counts invocations."""
    calls: list[tuple] = []

    def fake_score(smiles, ion_bins, loss_bins, *, lib_id=209, fallback_on_error=True):
        calls.append((smiles, ion_bins, loss_bins, lib_id))
        return round((len(str(smiles)) * 0.01 + int(lib_id) * 0.001) % 1.0, 6)

    monkeypatch.setattr(verifier_v2, "predicted_spectrum_score_v2", fake_score)
    return calls


def _candidates() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "control_id": ["X-1", "X-1", "X-2", "X-3", "X-1"],
            "predicted_smiles": ["CCO", "CCO", "CCN", "CC(=O)O", "CCC"],
            "lib_id": [209, 209, 209, 209, 209],
        }
    )


def _observed() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "control_id": ["X-1", "X-2"],
            "observed_ion_bins_json": ["[100, 200]", "[110]"],
            "observed_loss_bins_json": ["[50]", "[60, 70]"],
            "spectrum_source": ["srcA", "srcB"],
            "lib_id": [209, 209],
        }
    )


def test_serial_and_parallel_outputs_identical(stub_scorer):
    serial = annotate_predicted_spectrum_scores_v2(_candidates(), _observed(), workers=1)
    parallel = annotate_predicted_spectrum_scores_v2(_candidates(), _observed(), workers=3)
    pd.testing.assert_frame_equal(serial, parallel)


def test_duplicate_keys_scored_once(stub_scorer):
    annotate_predicted_spectrum_scores_v2(_candidates(), _observed(), workers=1)
    # X-1/CCO appears twice but shares one key; X-3 has no observed spectrum.
    assert len(stub_scorer) == 3


def test_missing_observed_rows_marked(stub_scorer):
    out = annotate_predicted_spectrum_scores_v2(_candidates(), _observed(), workers=1)
    missing = out[out["control_id"] == "X-3"]
    assert (missing["predicted_spectrum_score_v2"] == 0.0).all()
    assert (missing["predicted_spectrum_score_v2_source"] == "missing_observed_ms2").all()


def test_score_cache_roundtrip(tmp_path):
    path = tmp_path / "scores.csv"
    cache = VerifierScoreCache(path)
    key = ("CCO", "[100]", "[50]", 209)
    assert cache.get(key) is None
    cache.put(key, 0.4321)
    cache.save()
    reloaded = VerifierScoreCache(path)
    assert reloaded.get(key) == pytest.approx(0.4321)
    assert len(reloaded) == 1


def test_warm_cache_skips_scoring_and_matches_cold(stub_scorer, tmp_path):
    path = tmp_path / "scores.csv"
    cold_cache = VerifierScoreCache(path)
    cold = annotate_predicted_spectrum_scores_v2(
        _candidates(), _observed(), workers=1, score_cache=cold_cache
    )
    cold_cache.save()
    n_cold_calls = len(stub_scorer)
    assert n_cold_calls == 3

    warm_cache = VerifierScoreCache(path)
    warm = annotate_predicted_spectrum_scores_v2(
        _candidates(), _observed(), workers=1, score_cache=warm_cache
    )
    assert len(stub_scorer) == n_cold_calls  # zero new scorer calls
    pd.testing.assert_frame_equal(cold, warm)


def test_corrupt_cache_file_is_ignored(tmp_path):
    path = tmp_path / "scores.csv"
    path.write_text("not,a,valid\ncache")
    cache = VerifierScoreCache(path)
    assert len(cache) == 0
    cache.put(("C", "[1]", "[2]", 209), 0.5)
    cache.save()
    assert VerifierScoreCache(path).get(("C", "[1]", "[2]", 209)) == pytest.approx(0.5)
