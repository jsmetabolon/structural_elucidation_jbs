"""Unit tests for the threshold-measurement harness."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


T = load_tool("audit_thresholds")


# --- per-source calibration ---------------------------------------------------


def _control(raw, correct, source):
    return {"raw_full": str(raw), "top1_correct_0p70": "True" if correct else "False", "source": source}


def test_per_source_calibration_needs_a_source_column():
    result = T.per_source_calibration([_control(0.5, True, "")])
    assert result["status"] == "not_evaluated"


def test_per_source_calibration_measures_each_source_separately():
    rows = (
        [_control(0.2 + 0.05 * i, i > 4, "exact") for i in range(10)]
        + [_control(0.3 + 0.04 * i, i > 4, "predicted") for i in range(10)]
    )
    result = T.per_source_calibration(rows)
    assert set(result["strata"]) == {"exact", "predicted"}
    for block in result["strata"].values():
        assert block["n"] == 10
        assert block["auc_raw"] is not None
        assert block["auc_per_source_loo"] is not None


def test_a_perfectly_ordering_raw_score_reads_as_such():
    rows = [_control(0.1 * i, i >= 5, "exact") for i in range(10)]
    result = T.per_source_calibration(rows)
    assert result["strata"]["exact"]["auc_raw"] == 1.0


def test_the_recommendation_names_the_raw_score_when_calibration_does_not_help():
    """The committed cohort's answer, and the one the tool has to be able to give."""
    committed = json.loads(
        (ROOT / "docs" / "paper" / "data" / "thresholds" / "audit_thresholds_summary.json").read_text()
    )
    block = committed["source_stratification"]["per_source_calibration"]
    assert block["per_source_beats_pooled"] is True
    assert block["per_source_beats_raw"] is False
    assert "raw evidence score" in block["recommendation"]


# --- retune traps -------------------------------------------------------------


def test_the_sweep_grid_always_contains_the_tier_cutoffs():
    """A retuned cutoff off the hand-written grid used to vanish from the artifact."""
    for value in (T.C.TIER_LOW, T.C.TIER_MODERATE, T.C.TIER_HIGH):
        assert round(value, 4) in T.SCORE_GRID


def test_named_operating_points_are_keyed_by_tier_not_by_value():
    """Keys embedding the number would assert the old cutoff after a retune."""
    for name in T.NAMED_SCORE_POINTS:
        if name.startswith("tier_"):
            assert not any(character.isdigit() for character in name), name


def test_unpicking_the_grid_union_fails_at_import_rather_than_silently():
    original = dict(T.NAMED_SCORE_POINTS)
    T.NAMED_SCORE_POINTS["tier_high"] = 0.6123
    try:
        with pytest.raises(AssertionError, match="off the sweep grid"):
            T._assert_named_points_are_on_the_grid()
    finally:
        T.NAMED_SCORE_POINTS.clear()
        T.NAMED_SCORE_POINTS.update(original)


def test_a_deliberately_narrow_grid_may_still_omit_a_named_point():
    """Unit callers pass narrow grids on purpose; only the default grid is guarded."""
    assert T._named_points([{"threshold": 0.9}], (0.9,)) == {}


def test_every_named_point_survives_a_sweep_on_the_default_grid():
    scores = [0.05 * i for i in range(20)]
    labels = [i % 2 == 0 for i in range(20)]
    result = T.score_cutoff_sweep(scores, labels)
    assert set(result["named_operating_points"]) == set(T.NAMED_SCORE_POINTS)


def test_degeneracy_check_needs_levels_to_say_anything():
    assert T.positive_class_degeneracy([_control(0.5, True, "exact")], None)["status"] == "not_evaluated"


def test_a_cohort_whose_positives_are_all_level_1_is_flagged():
    rows = [_control(0.9, True, "exact"), _control(0.2, False, "exact")]
    rows[0]["control_id"], rows[1]["control_id"] = "X-1", "X-2"
    result = T.positive_class_degeneracy(rows, {"X-1": "1", "X-2": "3"})
    assert result["positive_class_is_level_1"] is True
    assert result["usable_for_non_level_1_reports"] is False
    assert result["correct_outside_level_1"] == 0


def test_a_cohort_with_correct_calls_outside_level_1_is_usable():
    rows = [_control(0.9, True, "exact"), _control(0.8, True, "exact"), _control(0.2, False, "exact")]
    for row, cid in zip(rows, ("X-1", "X-2", "X-3")):
        row["control_id"] = cid
    result = T.positive_class_degeneracy(rows, {"X-1": "1", "X-2": "3", "X-3": "3"})
    assert result["positive_class_is_level_1"] is False
    assert result["usable_for_non_level_1_reports"] is True
    assert result["correct_outside_level_1"] == 1


def test_the_committed_cohort_cannot_validate_a_tier_for_unconfirmed_reports():
    """The finding that stops a High tier being set from this data.

    If this ever flips, the truth cohort has grown correct calls outside Level 1
    and the High question genuinely reopens - so the assertion is the trigger,
    not an invariant to preserve.
    """
    committed = json.loads(
        (ROOT / "docs" / "paper" / "data" / "thresholds" / "audit_thresholds_summary.json").read_text()
    )
    block = committed["source_stratification"]["positive_class_degeneracy"]
    assert block["positive_class_is_level_1"] is True
    assert block["correct_outside_level_1"] == 0
    assert block["usable_for_non_level_1_reports"] is False


def test_committed_per_source_calibration_keeps_raw_above_both_calibrations():
    """Guards the A3 conclusion: bands belong on the raw score, within a source."""
    committed = json.loads(
        (ROOT / "docs" / "paper" / "data" / "thresholds" / "audit_thresholds_summary.json").read_text()
    )
    for source, block in committed["source_stratification"]["per_source_calibration"]["strata"].items():
        assert block["auc_raw"] > block["auc_pooled_loo"], source
        assert block["auc_raw"] > block["auc_per_source_loo"], source
        assert block["auc_per_source_loo"] > block["auc_pooled_loo"], source


def test_auc_is_one_for_a_perfectly_ordering_score():
    assert T.auc([0.9, 0.8, 0.2, 0.1], [True, True, False, False]) == 1.0


def test_auc_is_half_for_a_score_that_carries_no_signal():
    assert T.auc([0.5, 0.5, 0.5, 0.5], [True, False, True, False]) == 0.5
    assert T.auc([0.9, 0.1, 0.8, 0.2], [True, True, False, False]) == 0.5


def test_auc_is_undefined_without_both_classes():
    assert T.auc([0.9, 0.8], [True, True]) is None


def test_score_sweep_reports_coverage_and_lift():
    scores = [0.9, 0.7, 0.3, 0.1]
    labels = [True, True, False, False]
    result = T.score_cutoff_sweep(scores, labels, grid=(0.0, 0.5))
    assert result["base_rate"] == 0.5
    everything, confident = result["grid"]
    assert everything["coverage"] == 1.0
    assert everything["lift_over_base_rate"] == 0.0
    assert confident["n"] == 2
    assert confident["accuracy"] == 1.0
    assert confident["lift_over_base_rate"] == 0.5


def test_score_sweep_tolerates_a_threshold_that_keeps_nothing():
    result = T.score_cutoff_sweep([0.1, 0.2], [True, False], grid=(0.9,))
    point = result["grid"][0]
    assert point["n"] == 0
    assert point["accuracy"] is None
    assert point["lift_over_base_rate"] is None


def _depth_rows(triples: list[tuple[int, int, int]]) -> list[dict[str, str]]:
    return [
        {
            "top1_correct_0p70": str(one),
            "top3_correct_0p70": str(three),
            "top10_correct_0p70": str(ten),
        }
        for one, three, ten in triples
    ]


def test_candidate_depth_counts_incorrect_structures_displayed():
    # Two of four controls have the correct structure at rank 1; all four have it
    # somewhere in the top 3.
    rows = _depth_rows([(1, 1, 1), (1, 1, 1), (0, 1, 1), (0, 1, 1)])
    points = T.candidate_depth_sweep(rows)["points"]
    top1, top3, top10 = points
    assert top1["hit_rate"] == 0.5
    assert top1["incorrect_shown_per_compound"] == 0.5
    assert top3["hit_rate"] == 1.0
    assert top3["incorrect_shown_per_compound"] == 2.0
    # Half a hit gained for 1.5 more incorrect structures on the page.
    assert top3["incorrect_shown_per_additional_hit"] == 3.0
    # A deeper list that gains nothing has no exchange rate to report.
    assert top10["incorrect_shown_per_additional_hit"] is None


def _axis_record(**over) -> dict:
    base = dict(
        structure=0.9,
        mass=0.95,
        ms2=0.5,
        observed_mz_available=True,
        source="predicted",
        correct=True,
    )
    base.update(over)
    return base


def test_msi_replay_matches_the_delivered_rule():
    strong = _axis_record()
    assert T._msi_level(strong, T.C.DEFAULT_MSI_THRESHOLDS) == 2
    # No observed MS2 spectrum drops the same evidence to a hypothesis.
    assert T._msi_level(_axis_record(ms2=None), T.C.DEFAULT_MSI_THRESHOLDS) == 3
    # A candidate whose mass disagrees cannot claim a structure at all.
    assert T._msi_level(_axis_record(mass=0.0), T.C.DEFAULT_MSI_THRESHOLDS) == 4


def test_msi_sweep_tightening_the_bar_names_fewer_compounds():
    records = [
        _axis_record(structure=0.95, correct=True),
        _axis_record(structure=0.75, correct=False),
        _axis_record(structure=0.65, correct=False),
    ]
    sweep = T.msi_rule_sweep(records)["sweep"]
    by_rule = {
        (point["struct_moderate"], point["ms2_present"], point["mass_strong"]): point
        for point in sweep
        if point["ms2_probable_exact"] == 0.50
    }
    loose = by_rule[(0.60, 0.30, 0.90)]
    tight = by_rule[(0.85, 0.30, 0.90)]
    assert loose["n"] == 3
    assert loose["accuracy"] == round(1 / 3, 4)
    assert tight["n"] == 1
    assert tight["accuracy"] == 1.0
    assert tight["incorrect_names"] == 0


def test_msi_sweep_never_offers_a_catalog_hit_an_easier_bar():
    sweep = T.msi_rule_sweep([_axis_record()])["sweep"]
    assert sweep
    assert all(point["ms2_probable_exact"] >= point["ms2_present"] for point in sweep)


def test_source_stratification_finds_the_pooled_curve_is_confounded():
    """Guard the finding behind the exact-source MSI bar.

    A catalog hit scores high and is right less often than a de novo prediction,
    so pooling the two inverts the signal. If that ever stops being true the
    exact-source bar should be revisited.
    """

    import json

    committed = json.loads(
        (ROOT / "docs/paper/data/thresholds/audit_thresholds_summary.json").read_text()
    )["source_stratification"]
    assert committed["confounded"] is True
    exact = committed["strata"]["exact"]["raw_full"]
    predicted = committed["strata"]["predicted"]["raw_full"]
    # Scores higher, correct less often.
    assert exact["median_score"] > predicted["median_score"]
    assert exact["base_rate"] < predicted["base_rate"]
    # Yet within each source the raw score does order correctness.
    assert min(exact["auc"], predicted["auc"]) > committed["pooled_auc"]


def test_axis_records_are_skipped_when_the_controls_file_predates_them():
    assert T.axis_records_from_controls([{"control_id": "X-1", "top1_correct_0p70": "1"}]) == []


def test_committed_threshold_audit_matches_a_fresh_computation():
    """The committed audit must still follow from the committed per-control data."""

    import csv
    import json

    controls = ROOT / "docs/paper/data/structure/audit_structure_controls.csv"
    with controls.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    committed = json.loads(
        (ROOT / "docs/paper/data/thresholds/audit_thresholds_summary.json").read_text()
    )

    scored = [row for row in rows if T._float(row.get("loo_probability_full")) is not None]
    scores = [float(row["loo_probability_full"]) for row in scored]
    labels = [T._truthy(row.get("top1_correct_0p70", "")) for row in scored]

    levels = T.load_msi_levels(ROOT / "docs/paper/data/exec_summary/exec_summary_by_compound.csv")

    assert T.score_cutoff_sweep(scores, labels) == committed["score_cutoff"]
    assert T.candidate_depth_sweep(rows) == committed["candidate_depth"]
    assert T.source_stratification(rows, levels) == committed["source_stratification"]


def test_committed_audit_shows_the_score_does_not_order_correctness():
    """Guard the finding the cutoff decisions rest on.

    If the evidence score ever becomes discriminative, this test fails and the
    naming and candidate-count rules should be revisited: they are rule-based
    precisely because a score-based cutoff cannot be justified today.
    """

    import json

    committed = json.loads(
        (ROOT / "docs/paper/data/thresholds/audit_thresholds_summary.json").read_text()
    )
    interval = committed["score_cutoff"]["auc_ci95"]
    assert interval[0] < 0.5 < interval[1]
