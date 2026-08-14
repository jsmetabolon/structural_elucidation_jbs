"""Tests for the feces pathway re-weighting.

Two failure modes matter here and neither is loud. A vocabulary the pipeline
cannot emit will lose share purely because the denominator grew, which reads as
a finding if it is not labelled. And restricting the named metabolites to the
frequently detected ones while adding every X-compound compares unlike cohorts,
which inflates the shift in the direction the exercise is looking for.
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


def panel_row(pathway, bucket="Very Often (>90%)", rate=0.9, sub=None):
    return {
        "chemical_id": "1",
        "metabolite": "m",
        "super_pathway": F.canonical_pathway(pathway),
        "sub_pathway": F.canonical_pathway(sub or pathway),
        "gdp_bucket": bucket,
        "gdp_rate": rate,
    }


def assignment(pathway, tier="high", sub=None):
    return {
        "super_pathway": pathway,
        "sub_pathway": sub or pathway,
        "confidence_tier": tier,
        "assignment_status": "strong_similarity_vote",
    }


# --- vocabulary --------------------------------------------------------------


def test_the_two_spellings_of_secondary_metabolism_are_one_class():
    """The sheet carries both. Left apart, one class reports as two."""

    assert F.canonical_pathway("Secondary metabolism") == F.canonical_pathway(
        "Secondary Metabolism"
    )


def test_a_missing_pathway_is_named_rather_than_blank():
    assert F.canonical_pathway(None) == "(unassigned)"
    assert F.canonical_pathway("  ") == "(unassigned)"


# --- composition -------------------------------------------------------------


def test_shares_are_recomputed_against_the_grown_denominator():
    baseline = [panel_row("Lipid"), panel_row("Lipid"), panel_row("Amino Acid")]
    added = [{"super_pathway": "Amino Acid", "sub_pathway": "Amino Acid", "detection_rate": 0.5}]
    rows = {row["chemical_pathway"]: row for row in F.compose(baseline, added, "super", False)}
    assert rows["Lipid"]["baseline_pct"] == 66.67
    assert rows["Lipid"]["updated_pct"] == 50.0
    assert rows["Amino Acid"]["updated_pct"] == 50.0
    assert rows["Amino Acid"]["x_compounds_added"] == 1


def test_a_class_the_pipeline_never_predicts_is_flagged():
    """Its share can only fall here, so the movement is an artefact of the denominator."""

    baseline = [panel_row("Hormone metabolism"), panel_row("Lipid")]
    added = [{"super_pathway": "Lipid", "sub_pathway": "Lipid", "detection_rate": 1.0}]
    rows = {row["chemical_pathway"]: row for row in F.compose(baseline, added, "super", False)}
    assert rows["Hormone metabolism"]["outside_predicted_vocabulary"] == "TRUE"
    assert rows["Hormone metabolism"]["delta_pct_points"] < 0
    assert rows["Lipid"]["outside_predicted_vocabulary"] == "FALSE"


def test_weighting_counts_a_rarely_seen_compound_for_less():
    """Driven through build() so the detection rates are wired the way the tool wires them."""

    panel = [panel_row("Lipid", rate=1.0)]
    detection = {"X-1": 0.02}
    assignments = {"X-1": assignment("Peptide")}
    unweighted = {
        r["chemical_pathway"]: r for r in F.build(panel, detection, assignments)["super_pathway"]
    }
    weighted = {
        r["chemical_pathway"]: r
        for r in F.build(panel, detection, assignments, weight_by_detection=True)["super_pathway"]
    }
    assert unweighted["Peptide"]["updated_pct"] == 50.0
    assert weighted["Peptide"]["updated_pct"] < 5.0


def test_a_record_with_no_detection_rate_fails_loudly_when_weighting():
    """Silently weighting it zero would delete that side of the comparison."""

    import pytest

    with pytest.raises(KeyError):
        F.compose([{"super_pathway": "Lipid", "sub_pathway": "Lipid"}], [], "super", True)


def test_build_does_not_mutate_the_panel_it_is_given():
    """The panel is reused across every variant the tool reports."""

    panel = [panel_row("Lipid")]
    before = dict(panel[0])
    F.build(panel, {"X-1": 0.5}, {"X-1": assignment("Peptide")})
    assert panel[0] == before


def test_subpathway_level_aggregates_separately():
    baseline = [panel_row("Lipid", sub="Sterol"), panel_row("Lipid", sub="Bile Acid")]
    rows = F.compose(baseline, [], "sub", False)
    assert {row["chemical_pathway"] for row in rows} == {"Sterol", "Bile Acid"}


# --- detection thresholds ----------------------------------------------------


def test_only_detected_metabolites_form_the_baseline():
    panel = [panel_row("Lipid"), panel_row("Peptide", bucket=F.NEVER)]
    result = F.build(panel, {}, {})
    assert result["n_identified_detected"] == 1


def test_the_common_only_view_applies_the_same_bar_to_both_sides():
    """A rarely seen X-compound must not join a frequently seen baseline."""

    panel = [panel_row("Lipid", bucket="Very Often (>90%)")]
    detection = {"X-1": 0.90, "X-2": 0.02}
    assignments = {"X-1": assignment("Peptide"), "X-2": assignment("Peptide")}

    everything = F.build(panel, detection, assignments)
    common = F.build(panel, detection, assignments, common_only=True)
    assert everything["n_x_compounds_added"] == 2
    assert common["n_x_compounds_added"] == 1


def test_an_x_compound_absent_from_the_feces_history_is_not_added():
    """The bundle covers more compounds than feces does."""

    panel = [panel_row("Lipid")]
    result = F.build(panel, {"X-1": 0.5}, {"X-1": assignment("Peptide"), "X-9": assignment("Lipid")})
    assert result["n_x_compounds_added"] == 1


def test_stratifying_by_tier_selects_only_that_tier():
    panel = [panel_row("Lipid")]
    detection = {"X-1": 0.5, "X-2": 0.5}
    assignments = {"X-1": assignment("Peptide", "high"), "X-2": assignment("Peptide", "low")}
    assert F.build(panel, detection, assignments, tier="high")["n_x_compounds_added"] == 1
