"""Asset routing uses selected model data and retains canonical checkout code."""

import os
import pickle
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

from structural_elucidation.nce import runtime as R
from structural_elucidation.nce import verifier_v2


@pytest.fixture(autouse=True)
def isolated_runtime(monkeypatch):
    monkeypatch.delenv("NCE_REFERENCE_DIR", raising=False)
    monkeypatch.setenv("REFERENCE_DIR", "pretest")
    monkeypatch.setattr(R, "install_local_nce_shims", lambda: None)
    monkeypatch.setattr(sys, "path", list(sys.path))
    R._get_predicted_spectrum_matcher.cache_clear()
    R._get_cached_mlims_service.cache_clear()
    R._get_cached_mlims_service_from_pickle.cache_clear()
    yield
    R._get_predicted_spectrum_matcher.cache_clear()
    R._get_cached_mlims_service.cache_clear()
    R._get_cached_mlims_service_from_pickle.cache_clear()


def test_old_default_and_explicit_precedence(monkeypatch, tmp_path):
    assert R.resolve_reference_dir() == R.DEFAULT_REFERENCE_DIR
    monkeypatch.setenv("NCE_REFERENCE_DIR", str(tmp_path / "environment"))
    assert R.resolve_reference_dir() == tmp_path / "environment"
    assert R.resolve_reference_dir(tmp_path / "explicit") == tmp_path / "explicit"
    monkeypatch.setenv("NCE_REFERENCE_DIR", "   ")
    assert R.resolve_reference_dir() == R.DEFAULT_REFERENCE_DIR


def test_prepare_runtime_selects_data_without_redirecting_code(monkeypatch, tmp_path):
    monkeypatch.setenv("NCE_REFERENCE_DIR", str(tmp_path / "models"))
    R.prepare_nce_runtime()
    assert os.environ["REFERENCE_DIR"] == str(tmp_path / "models")
    assert str(R.DEFAULT_NCE_ROOT) in sys.path
    assert str(tmp_path / "vendor/nce") not in sys.path


def test_matcher_cache_is_keyed_by_resolved_model_directory(monkeypatch, tmp_path):
    calls = []
    module = types.ModuleType("elucidation.predictedSpectrumMatcher")
    class FakeMatcher:
        @staticmethod
        def get_instance(lib_id):
            value = (lib_id, os.environ["REFERENCE_DIR"])
            calls.append(value)
            return value
    module.PredictedSpectrumMatcher = FakeMatcher
    monkeypatch.setitem(sys.modules, module.__name__, module)
    monkeypatch.setenv("NCE_REFERENCE_DIR", str(tmp_path / "first"))
    assert R.get_predicted_spectrum_matcher(400) == (400, str(tmp_path / "first"))
    assert R.get_predicted_spectrum_matcher(400) == (400, str(tmp_path / "first"))
    monkeypatch.setenv("NCE_REFERENCE_DIR", str(tmp_path / "second"))
    assert R.get_predicted_spectrum_matcher(400) == (400, str(tmp_path / "second"))
    assert R.get_predicted_spectrum_matcher(400, reference_dir=tmp_path / "explicit") == (400, str(tmp_path / "explicit"))
    assert len(calls) == 3


def test_actual_verifier_call_passes_selected_data_to_matcher(monkeypatch, tmp_path):
    calls = []
    matcher = types.SimpleNamespace(model_dict={"bin_100": object()},
                                    score_a_smiles_for_spectrum=lambda *args: 0.42)
    def get_matcher(lib_id, **kwargs):
        calls.append((lib_id, kwargs))
        return matcher
    monkeypatch.setattr(R, "get_predicted_spectrum_matcher", get_matcher)
    monkeypatch.setenv("NCE_REFERENCE_DIR", str(tmp_path / "external_models"))
    assert verifier_v2.predicted_spectrum_score_v2("CCO", "[100]", "[]", lib_id=400) == 0.42
    assert calls == [(400, {"nce_root": str(R.DEFAULT_NCE_ROOT),
                           "reference_dir": str(tmp_path / "external_models")})]
    verifier_v2.predicted_spectrum_score_v2("CCO", "[100]", "[]", reference_dir=str(tmp_path / "explicit"))
    assert calls[-1][1]["reference_dir"] == str(tmp_path / "explicit")


def test_mlims_service_uses_selected_directory(monkeypatch, tmp_path):
    module = types.ModuleType("elucidation.library_services.cachedMlimsLibraryService")
    class FakeService:
        @staticmethod
        def get_instance():
            return os.environ["REFERENCE_DIR"]
    module.CachedMlimsLibraryService = FakeService
    monkeypatch.setitem(sys.modules, module.__name__, module)
    monkeypatch.setenv("NCE_REFERENCE_DIR", str(tmp_path / "environment"))
    assert R.get_cached_mlims_service() == str(tmp_path / "environment")
    assert R.get_cached_mlims_service(reference_dir=tmp_path / "explicit") == str(tmp_path / "explicit")


def test_pickle_service_resolves_directory_before_caching(monkeypatch, tmp_path):
    for name in ("first", "second"):
        folder = tmp_path / name
        folder.mkdir()
        with (folder / "mlims_tables.pkl").open("wb") as handle:
            pickle.dump({"chro_comp_entry": pd.DataFrame([{"comp_id": 1, "name": name}])}, handle)
        monkeypatch.setenv("NCE_REFERENCE_DIR", str(folder))
        assert R.get_cached_mlims_service_from_pickle().get_name_dict() == {1: name}
