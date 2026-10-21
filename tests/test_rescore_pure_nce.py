"""Pure matcher scores distinguish a valid zero from missing or failed evidence."""

from types import SimpleNamespace

import pandas as pd
import pytest

from structural_elucidation.nce import runtime
from tools import rescore_pure_nce as N


@pytest.fixture
def matcher(monkeypatch):
    scorer = SimpleNamespace(model_dict={"ion_100": object()})
    scorer.score_a_smiles_for_spectrum = lambda *_: 0.25
    monkeypatch.setattr(N, "_valid_smiles", lambda smiles: bool(smiles))
    monkeypatch.setattr(runtime, "get_predicted_spectrum_matcher", lambda *a, **kw: scorer)
    return scorer


@pytest.mark.parametrize("value,status", [(0.25, "scored"), (0.0, "zero")])
def test_successful_scores_keep_their_value(matcher, value, status):
    matcher.score_a_smiles_for_spectrum = lambda *_: value
    assert N.score_key(("CCO", "[100]", "[]", 400)) == (value, status)


def test_matcher_exception_is_not_a_successful_zero(matcher):
    def fail(*_):
        raise RuntimeError("model failure")
    matcher.score_a_smiles_for_spectrum = fail
    assert N.score_key(("CCO", "[100]", "[]", 400)) == (0.0, "error")


def test_matcher_load_exception_is_unavailable(monkeypatch):
    def fail(*_a, **_kw):
        raise FileNotFoundError("model absent")
    monkeypatch.setattr(N, "_valid_smiles", lambda smiles: True)
    monkeypatch.setattr(runtime, "get_predicted_spectrum_matcher", fail)
    assert N.score_key(("CCO", "[100]", "[]", 400)) == (0.0, "error")


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -0.1, 1.1])
def test_invalid_matcher_values_are_unavailable(matcher, value):
    matcher.score_a_smiles_for_spectrum = lambda *_: value
    assert N.score_key(("CCO", "[100]", "[]", 400)) == (0.0, "error")


@pytest.mark.parametrize("ions,losses", [("[]", "[]"), ("[0]", "[]"), ("invalid", "[]")])
def test_missing_bins_do_not_load_the_matcher(monkeypatch, ions, losses):
    def unexpected(*_a, **_kw):
        raise AssertionError("no model is needed without observations")
    monkeypatch.setattr(runtime, "get_predicted_spectrum_matcher", unexpected)
    assert N.score_key(("CCO", ions, losses, 400)) == (0.0, "no_observed_bins")


def test_observed_lookup_keeps_library_precedence_and_missing_losses(tmp_path):
    for lib, ion, bins in [(400, "positive", "[200]"), (209, "negative", "[100]")]:
        pd.DataFrame([{"control_id": "X-1", "observed_ion_bins_json": bins}]).to_csv(
            tmp_path / f"lib{lib}_{ion}_observed_ms2_lookup.csv", index=False
        )
    row = N.load_observed(tmp_path, {"X-1"}).iloc[0]
    assert (row.lib_id, row.ion_json, row.loss_json) == (209, "[100]", "[]")
    assert N.load_observed(tmp_path, {"X-2"}).empty


def test_empty_observed_lookup_retains_merge_columns(tmp_path):
    assert list(N.load_observed(tmp_path, {"X-1"}).columns) == [
        "control_id", "ion_json", "loss_json", "lib_id"
    ]
