"""Unit tests for the retention-index feasibility audit."""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


RT = load_tool("audit_rt_feasibility")

from structural_elucidation.utils.domain import ri_from_feature_key  # noqa: E402


def test_ri_is_recovered_from_the_feature_key():
    key = "feature_X_209_2021-10-14-1545_mz100.0404_ri1248_X-23917"
    assert ri_from_feature_key(key) == 1248.0


def test_a_zero_retention_index_means_not_measured():
    # The upstream tables write 0 for "no index", which is not an elution at zero.
    assert math.isnan(ri_from_feature_key("feature_X_400_mz100.0_ri0_X-1"))
    assert math.isnan(ri_from_feature_key("feature_with_no_index_X-1"))
    assert math.isnan(ri_from_feature_key(None))


def test_linear_fit_recovers_a_known_line():
    xs = [1.0, 2.0, 3.0, 4.0]
    ys = [3.0, 5.0, 7.0, 9.0]
    slope, intercept = RT.fit_linear(xs, ys)
    assert slope == pytest.approx(2.0)
    assert intercept == pytest.approx(1.0)
    assert RT.fit_linear([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) is None


def test_pearson_matches_the_obvious_cases():
    assert RT.pearson([1.0, 2.0, 3.0], [2.0, 4.0, 6.0]) == pytest.approx(1.0)
    assert RT.pearson([1.0, 2.0, 3.0], [6.0, 4.0, 2.0]) == pytest.approx(-1.0)
    assert RT.pearson([1.0, 2.0], [1.0, 2.0]) is None


def test_auc_is_half_when_the_residual_carries_no_signal():
    assert RT.auc([1.0, 1.0, 1.0, 1.0], [True, False, True, False]) == 0.5
    assert RT.auc([2.0, 1.0], [True, False]) == 1.0


def _library(**over) -> dict:
    base = dict(
        n=30, ion_mode="negative", sufficient_controls=True,
        recalibrated_loo_rmse=500.0, fit={"slope": 3.0, "intercept": 100.0},
    )
    base.update(over)
    return base


def test_a_library_qualifies_only_on_all_three_conditions():
    calibration = {
        "209": _library(),
        "305": _library(fit={"slope": -3.9, "intercept": 3200.0}),
        "402": _library(n=2, sufficient_controls=False, fit=None),
        "999": _library(),
    }
    decoy = {
        "by_library": {
            "209": {"auc": 0.68},
            "305": {"auc": 0.72},  # good AUC cannot rescue a backwards fit
            "402": {"auc": None},
            "999": {"auc": 0.55},  # plausible fit, but no discrimination
        }
    }
    result = RT.qualifying_libraries(calibration, decoy)
    assert result["209"]["qualifies"] is True
    assert result["305"]["qualifies"] is False
    assert any("negative" in reason for reason in result["305"]["excluded_because"])
    assert result["402"]["qualifies"] is False
    assert result["999"]["qualifies"] is False
    assert any("below the" in reason for reason in result["999"]["excluded_because"])


def test_library_phrase_reads_as_english():
    assert RT._library_phrase([]) == "no library"
    assert RT._library_phrase(["209"]) == "library 209"
    assert RT._library_phrase(["209", "400"]) == "libraries 209 and 400"


def test_decoy_scoring_prefers_the_structure_that_fits_the_observed_index():
    records = [
        {
            "control_id": "X-1", "lib_id": 209, "ion_mode": "negative",
            "observed_ri": 1000.0, "exact_mass": 200.0, "recalibrated_ri": 1000.0,
        },
        {
            "control_id": "X-2", "lib_id": 209, "ion_mode": "negative",
            "observed_ri": 3000.0, "exact_mass": 210.0, "recalibrated_ri": 3000.0,
        },
    ]
    result = RT.decoy_discrimination(records, mass_tolerance=50.0)
    # Each true structure predicts its own index exactly and the other's badly.
    assert result["auc"] == 1.0
    assert result["n_controls"] == 2
    assert result["controls_beating_median_decoy"] == 2


def test_decoys_are_drawn_only_from_the_same_library():
    records = [
        {
            "control_id": "X-1", "lib_id": 209, "ion_mode": "negative",
            "observed_ri": 1000.0, "exact_mass": 200.0, "recalibrated_ri": 1000.0,
        },
        {
            "control_id": "X-2", "lib_id": 400, "ion_mode": "positive",
            "observed_ri": 3000.0, "exact_mass": 205.0, "recalibrated_ri": 3000.0,
        },
    ]
    # Different columns are different methods, so neither compound is a decoy for
    # the other and there is nothing to score.
    assert RT.decoy_discrimination(records, mass_tolerance=50.0)["n_controls"] == 0


def test_committed_audit_still_reflects_the_committed_inputs():
    pytest.importorskip("rdkit", reason="the retention proxy needs rdkit")
    committed = json.loads(
        (ROOT / "docs/paper/data/retention/audit_rt_feasibility.json").read_text()
    )
    records = RT.build_records(
        RT.load_truth_structures(ROOT / "data/known_X_all_available_structures_simple.csv"),
        RT.load_observed_ri(ROOT / "data/feature_rows"),
    )
    # calibration_by_library writes the leave-one-out fit onto the records, so the
    # decoy check has to follow it, as it does in the tool.
    assert RT.calibration_by_library(records) == committed["calibration"]
    assert (
        RT.decoy_discrimination(records, use_recalibrated=True) == committed["decoy_recalibrated"]
    )
    assert (
        RT.decoy_discrimination(records, use_recalibrated=False, keep_per_control=False)
        == committed["decoy_raw"]
    )


def test_the_raw_proxy_is_not_usable_and_recalibration_is_what_rescues_it():
    """Guard the reason the axis loads a fitted calibration rather than constants.

    The published proxy constants score exactly at chance against mass-matched
    decoys. If that ever changes the axis could read the proxy directly.
    """

    committed = json.loads(
        (ROOT / "docs/paper/data/retention/audit_rt_feasibility.json").read_text()
    )
    assert committed["decoy_raw"]["auc"] < 0.55
    assert committed["decoy_recalibrated"]["auc"] >= RT.DECOY_AUC_GATE


def test_only_audited_libraries_appear_in_the_calibration_the_model_loads():
    calibration = json.loads((ROOT / "docs/paper/data/retention/ri_calibration.json").read_text())
    audit = json.loads(
        (ROOT / "docs/paper/data/retention/audit_rt_feasibility.json").read_text()
    )
    qualified = {lib for lib, block in audit["libraries"].items() if block["qualifies"]}
    assert set(calibration["libraries"]) == qualified
    assert all(block["slope"] > 0 for block in calibration["libraries"].values())
