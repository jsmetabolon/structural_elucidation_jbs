"""Tests for the feces pathway re-weighting.

Two things here can quietly produce a plausible but wrong answer: comparing
cohorts filtered on different rules, and letting one class report as two because
the source spells it two ways. Both are covered, along with the arithmetic.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


F = load_tool("build_feces_pathway_update")


def metabolite(pathway: str, bucket: str = "Often (50-90%)", rate: float = 0.6):
    return {
        "chemical_id": "1",
        "metabolite": "m",
        "super_pathway": F.canonical_pathway(pathway),
        "sub_pathway": F.canonical_pathway(pathway + " sub"),
        "gdp_bucket": bucket,
        "gdp_rate": rate,
        # build() copies gdp_rate here; compose() is given records that already have it.
        "detection_rate": rate,
    }


def assignment(pathway: str, tier: str = "high"):
    return {
        "super_pathway": pathway,
        "sub_pathway": pathway + " sub",
        "confidence_tier": tier,
        "assignment_status": "strong_similarity_vote",
    }


# --- vocabulary --------------------------------------------------------------


def test_the_two_spellings_of_secondary_metabolism_collapse_into_one():
    """Left apart, the sheet reports one class as two and both shares are wrong."""

    assert F.canonical_pathway("Secondary metabolism") == F.canonical_pathway(
        "Secondary Metabolism"
    )


def test_a_blank_pathway_is_labelled_rather_than_dropped():
    assert F.canonical_pathway(None) == "(unassigned)"
    assert F.canonical_pathway("  ") == "(unassigned)"


# --- composition arithmetic --------------------------------------------------


def test_shares_are_recomputed_against_the_grown_denominator():
    baseline = [metabolite("Lipid"), metabolite("Lipid"), metabolite("Amino Acid")]
    added = [{"super_pathway": "Amino Acid", "sub_pathway": "s", "detection_rate": 0.5}]
    rows = {row["chemical_pathway"]: row for row in F.compose(baseline, added, "super", False)}
    assert rows["Lipid"]["baseline_pct"] == 66.67
    assert rows["Lipid"]["updated_pct"] == 50.0
    assert rows["Amino Acid"]["updated_pct"] == 50.0
    assert rows["Amino Acid"]["delta_pct_points"] == 16.67


def test_a_class_the_pipeline_never_predicts_is_flagged():
    """Its share falls purely because the denominator grew, which is not a finding."""

    baseline = [metabolite("Hormone metabolism"), metabolite("Lipid")]
    added = [{"super_pathway": "Lipid", "sub_pathway": "s", "detection_rate": 1.0}]
    rows = {row["chemical_pathway"]: row for row in F.compose(baseline, added, "super", False)}
    assert rows["Hormone metabolism"]["outside_predicted_vocabulary"] == "TRUE"
    assert rows["Hormone metabolism"]["delta_pct_points"] < 0
    assert rows["Lipid"]["outside_predicted_vocabulary"] == "FALSE"


def test_weighting_counts_a_rarely_seen_compound_for_less():
    baseline = [metabolite("Lipid", rate=1.0)]
    added = [{"super_pathway": "Amino Acid", "sub_pathway": "s", "detection_rate": 0.02}]
    unweighted = {r["chemical_pathway"]: r for r in F.compose(baseline, added, "super", False)}
    weighted = {r["chemical_pathway"]: r for r in F.compose(baseline, added, "super", True)}
    assert unweighted["Amino Acid"]["updated_pct"] == 50.0
    assert weighted["Amino Acid"]["updated_pct"] < 2.0


def test_a_record_with_no_detection_rate_is_reported_not_silently_dropped(capsys):
    """Weighting turns a missing rate into a zero, which would shrink a class quietly."""

    baseline = [metabolite("Lipid", rate=1.0)]
    added = [{"super_pathway": "Amino Acid", "sub_pathway": "s", "detection_rate": None}]
    F.compose(baseline, added, "super", True)
    assert "no detection rate" in capsys.readouterr().err


# --- cohort filtering --------------------------------------------------------


def test_never_detected_metabolites_are_not_in_the_baseline():
    panel = [metabolite("Lipid", bucket="Never"), metabolite("Lipid")]
    result = F.build(panel, {"X-1": 0.9}, {"X-1": assignment("Lipid")})
    assert result["n_identified_detected"] == 1


def test_the_common_only_view_applies_the_same_bar_to_both_sides():
    """Filtering the named metabolites but not the X-compounds inflates the shift.

    The rarely detected X-compound must drop out alongside the rarely detected
    metabolite, or the comparison is between unlike cohorts.
    """

    panel = [metabolite("Lipid"), metabolite("Lipid", bucket="Rarely (<10%)", rate=0.05)]
    detection = {"X-common": 0.9, "X-rare": 0.05}
    assignments = {"X-common": assignment("Amino Acid"), "X-rare": assignment("Amino Acid")}

    everything = F.build(panel, detection, assignments)
    common = F.build(panel, detection, assignments, common_only=True)

    assert everything["n_identified_detected"] == 2
    assert everything["n_x_compounds_added"] == 2
    assert common["n_identified_detected"] == 1
    assert common["n_x_compounds_added"] == 1


def test_stratifying_by_tier_keeps_only_that_tier():
    panel = [metabolite("Lipid")]
    detection = {"X-1": 0.9, "X-2": 0.9}
    assignments = {"X-1": assignment("Lipid", "high"), "X-2": assignment("Lipid", "low")}
    assert F.build(panel, detection, assignments, tier="high")["n_x_compounds_added"] == 1


def test_an_x_compound_not_detected_in_feces_is_not_added():
    """The bundle covers more compounds than the feces detection history does."""

    panel = [metabolite("Lipid")]
    assignments = {"X-1": assignment("Lipid"), "X-elsewhere": assignment("Lipid")}
    result = F.build(panel, {"X-1": 0.9}, assignments)
    assert result["n_x_compounds_added"] == 1
