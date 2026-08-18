"""Unit tests for the MSI sub-level confidence-band audit."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


T = load_tool("audit_sublevel_bands")
C = T.C


def row(score, *, level="3", correct=False, truth=True, sublevel=""):
    return {
        "x_compound_id": f"X-{abs(hash((score, level, correct))) % 100000:05d}",
        "msi_level": level,
        "metabolon_sublevel": sublevel,
        "rank1_candidate_confidence": str(score),
        "rank1_evidence_confidence": str(score),
        "truth_basis": "library" if truth else "",
        "reference_answer": "CCO" if truth else "",
        "matched_at_rank_1": "True" if correct else "False",
    }


# --- reachability -------------------------------------------------------------


def test_platt_map_bounds_the_evidence_score_below_one():
    reach = T.platt_reachability()
    assert 0.0 < reach["evidence_score_min"] < reach["evidence_score_max"] < 1.0


def test_tier_high_is_not_reachable_by_the_evidence_model():
    """A cutoff above the Platt ceiling can only be entered by a policy override."""
    reach = T.platt_reachability()
    assert reach["cutoffs"]["TIER_HIGH"]["reachable"] is False
    assert reach["cutoffs"]["TIER_HIGH"]["raw_evidence_required"] > 1.0
    assert reach["unreachable_cutoffs"] == ["TIER_HIGH"]


def test_the_lower_two_tier_cutoffs_are_reachable():
    reach = T.platt_reachability()
    assert reach["cutoffs"]["TIER_MODERATE"]["reachable"] is True
    assert reach["cutoffs"]["TIER_LOW"]["reachable"] is True


# --- banding ------------------------------------------------------------------


def test_bands_are_assigned_from_the_top_down():
    cuts = (0.65, 0.48, 0.38)
    assert T._band_of(0.90, cuts) == "a"
    assert T._band_of(0.65, cuts) == "a"
    assert T._band_of(0.50, cuts) == "b"
    assert T._band_of(0.40, cuts) == "c"
    assert T._band_of(0.10, cuts) == "below"
    assert T._band_of(None, cuts) == "missing"


def test_occupancy_names_the_bands_nobody_can_enter():
    rows = [row(0.40), row(0.42), row(0.44)]
    result = T.band_occupancy(rows, "rank1_candidate_confidence", (0.65, 0.48, 0.38))
    assert result["counts"]["c"] == 3
    assert result["empty_bands"] == ["a", "b"]


# --- ordering -----------------------------------------------------------------


def test_a_correctly_ordered_score_reads_as_monotone():
    rows = (
        [row(0.90, correct=True) for _ in range(9)] + [row(0.90)]
        + [row(0.55, correct=True) for _ in range(5)] + [row(0.55) for _ in range(5)]
        + [row(0.40, correct=True)] + [row(0.40) for _ in range(9)]
    )
    result = T._precision(rows, "rank1_candidate_confidence", (0.65, 0.48, 0.38))
    assert result["monotone"] is True
    assert result["inversions"] == []
    assert result["bands"]["a"]["p_correct"] == 0.9


def test_an_inverted_score_is_reported_as_not_ordered():
    """The failure this audit exists to catch: a lower band that is more accurate."""
    rows = (
        [row(0.90, correct=True)] + [row(0.90) for _ in range(9)]
        + [row(0.55, correct=True) for _ in range(9)] + [row(0.55)]
    )
    result = T._precision(rows, "rank1_candidate_confidence", (0.65, 0.48, 0.38))
    assert result["monotone"] is False
    assert "a<b" in result["inversions"]


def test_a_flat_score_is_not_ordered_even_without_an_inversion():
    rows = [row(0.90, correct=i % 2 == 0) for i in range(10)] + [row(0.55, correct=i % 2 == 0) for i in range(10)]
    result = T._precision(rows, "rank1_candidate_confidence", (0.65, 0.48, 0.38))
    assert result["monotone"] is False
    assert result["inversions"] == []


def test_auc_reports_reversed_discrimination_below_half():
    assert T._auc([0.1, 0.2], [0.8, 0.9]) == 0.0
    assert T._auc([0.8, 0.9], [0.1, 0.2]) == 1.0
    assert T._auc([0.5, 0.5], [0.5, 0.5]) == 0.5
    assert T._auc([], [0.5]) is None


# --- curated sub-levels -------------------------------------------------------


def test_fine_sublevel_spread_ignores_the_bare_levels_and_curation_states():
    rows = (
        [row(0.40, sublevel="3A") for _ in range(30)]
        + [row(0.41, sublevel="3B") for _ in range(30)]
        + [row(0.97, sublevel="Released") for _ in range(30)]
    )
    result = T.curated_sublevel_separation(rows, "rank1_candidate_confidence")
    assert result["fine_sublevels"] == ["3A", "3B"]
    assert result["fine_median_spread"] < 0.05
    assert result["median_spread"] > 0.5


def test_a_compound_carrying_two_sublevels_counts_in_both():
    rows = [row(0.40, sublevel="3A; 3B") for _ in range(30)]
    result = T.curated_sublevel_separation(rows, "rank1_candidate_confidence")
    assert result["per_sublevel"]["3A"]["n"] == 30
    assert result["per_sublevel"]["3B"]["n"] == 30


# --- recovering raw_evidence --------------------------------------------------


def test_raw_evidence_round_trips_through_the_platt_map():
    """Recovery has to invert the forward map the confidence model applies."""
    for raw in (0.15, 0.4, 0.6609, 0.95):
        overall = C._sigmoid(C.PLATT_A * raw + C.PLATT_B)
        rows = [row(round(overall, 6))]
        T.raw_evidence_rows(rows)
        assert abs(float(rows[0][T.RAW_EVIDENCE_COLUMN]) - raw) < 1e-4


def test_a_policy_pinned_headline_is_left_unset_rather_than_clamped():
    """Inverting a pinned value gives raw > 1, which no axis average can produce."""
    rows = [row(0.97, level="1"), row(1.0, level="1"), row(0.4322)]
    recovery = T.raw_evidence_rows(rows)
    assert rows[0][T.RAW_EVIDENCE_COLUMN] == ""
    assert rows[1][T.RAW_EVIDENCE_COLUMN] == ""
    assert rows[2][T.RAW_EVIDENCE_COLUMN] != ""
    assert recovery["recovered"] == 1
    assert recovery["unrecoverable_policy_pinned"] == 2


def test_reanchoring_cannot_change_a_rank_statistic():
    """The measurement behind the re-anchor claim: a rescale is not a fix.

    Built so the two scores disagree in scale but not in order, which is what
    inverting a strictly monotone map guarantees.
    """
    rows = []
    for i in range(40):
        raw = 0.2 + 0.02 * i
        overall = C._sigmoid(C.PLATT_A * raw + C.PLATT_B)
        rows.append(row(round(overall, 6), correct=i < 10))
    T.raw_evidence_rows(rows)
    platt = T.discrimination(rows, "rank1_candidate_confidence")
    raw_side = T.discrimination(rows, T.RAW_EVIDENCE_COLUMN)
    assert platt["auc"] == raw_side["auc"]
    assert platt["median_correct"] != raw_side["median_correct"]


# --- the committed artifact ---------------------------------------------------


def test_committed_audit_shows_the_reanchor_fixing_only_reachability():
    payload = json.loads(
        (ROOT / "docs" / "paper" / "data" / "identification" / "audit_sublevel_bands.json").read_text()
    )
    eq = payload["reanchor_equivalence"]
    assert eq["auc_unchanged"] is True, "re-anchoring changed a rank statistic; re-read the findings"
    assert eq["fixes_reachability"] is True
    assert eq["fixes_ordering"] is False


# --- the committed artifact ---------------------------------------------------


def test_committed_audit_records_the_two_blocking_findings():
    """Guards the conclusions docs/msi_sublevel_band_findings.md rests on."""
    payload = json.loads(
        (ROOT / "docs" / "paper" / "data" / "identification" / "audit_sublevel_bands.json").read_text()
    )
    assert payload["platt_reachability"]["unreachable_cutoffs"] == ["TIER_HIGH"]
    for key in ("candidate_confidence", "evidence_confidence", "raw_evidence"):
        schemes = payload["scores"][key]["band_schemes"]
        assert not any(s["precision_model_scored"]["monotone"] for s in schemes.values()), (
            f"{key} now orders its bands; re-read the findings before relying on them"
        )
