"""Tests for the isomer-discrimination audit.

The audit exists to distinguish three things a reader of the top-10 export would
otherwise conflate: a truth candidate the pipeline scored and ranked below a
wrong one, a truth candidate it could not score at all, and a truth candidate
that is evidentially indistinguishable from the one above it. Getting those
apart is the whole point, so the tests concentrate there.

The spectral section needs rdkit and is exercised only for its skip contract;
everything else runs on the standard library.
"""

from __future__ import annotations

import importlib.util
import json
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


A = load_tool("audit_isomer_discrimination")


def candidate(rank, *, smiles="C", formula="C7H8O3S", matched=False, basis="emitted", **kwargs):
    row = {
        "x_compound_id": kwargs.get("x_compound_id", "X-1"),
        "cohort": "truth_control",
        "rank": str(rank),
        "candidate_smiles": smiles,
        "candidate_formula": formula,
        "candidate_confidence": kwargs.get("candidate_confidence", "0.5"),
        "evidence_confidence": kwargs.get("evidence_confidence", "0.99"),
        "evidence_confidence_basis": basis,
        "mass_match_score": kwargs.get("mass_match_score", "0.99"),
        "matched_truth": "TRUE" if matched else "FALSE",
        "msi_level": kwargs.get("msi_level", "4"),
    }
    return row


# --- why the truth is not at rank 1 ------------------------------------------


def test_a_truth_at_rank_1_is_not_counted_as_a_failure():
    controls = {"X-1": [candidate(1, matched=True), candidate(2)]}
    result = A.recovery_breakdown(controls)
    assert result["counts"]["truth_at_rank_1"] == 1
    assert result["n_truth_below_rank_1"] == 0


def test_a_missing_truth_is_absent_rather_than_below_rank_1():
    """Not recovered and recovered-but-misranked are different failures."""

    controls = {"X-1": [candidate(1), candidate(2)]}
    result = A.recovery_breakdown(controls)
    assert result["counts"]["truth_absent_from_top10"] == 1
    assert result["n_truth_below_rank_1"] == 0


def test_a_shared_formula_is_reported_as_a_degenerate_axis():
    """Same formula means the mass axis cannot have made the decision."""

    controls = {
        "X-1": [candidate(1, formula="C7H8O3S"), candidate(2, formula="C7H8O3S", matched=True)]
    }
    result = A.recovery_breakdown(controls)
    assert result["counts"]["evidence_degenerate_same_formula"] == 1
    assert result["n_same_formula_as_rank1"] == 1


def test_a_different_formula_is_not_called_degenerate():
    controls = {
        "X-1": [candidate(1, formula="C7H8O3S"), candidate(2, formula="C6H8N2O4", matched=True)]
    }
    result = A.recovery_breakdown(controls)
    assert result["counts"]["different_formula"] == 1
    assert result["n_same_formula_as_rank1"] == 0


def test_an_unscored_truth_row_is_separated_from_a_scored_one():
    """X-23080: the truth row carried no evidence while rank 1 carried a full score."""

    controls = {
        "X-1": [
            candidate(1, basis="emitted"),
            candidate(2, basis="unavailable", matched=True),
        ]
    }
    result = A.recovery_breakdown(controls)
    assert result["counts"]["truth_row_unscored"] == 1
    assert "evidence_degenerate_same_formula" not in result["counts"]


def test_an_unscored_truth_row_still_counts_toward_the_same_formula_total():
    """The buckets are exclusive; the same-formula tally is not."""

    controls = {
        "X-1": [
            candidate(1, formula="C7H8O3S", basis="emitted"),
            candidate(2, formula="C7H8O3S", basis="unavailable", matched=True),
        ]
    }
    result = A.recovery_breakdown(controls)
    assert result["counts"]["truth_row_unscored"] == 1
    assert result["n_same_formula_as_rank1"] == 1


def test_an_unscored_rank_1_is_not_blamed_on_the_truth_row():
    """When rank 1 is the unscored one, the truth row is not the missing measurement."""

    controls = {
        "X-1": [
            candidate(1, basis="unavailable"),
            candidate(2, basis="emitted", matched=True),
        ]
    }
    result = A.recovery_breakdown(controls)
    assert "truth_row_unscored" not in result["counts"]


def test_the_confidence_gap_is_signed_from_rank_1_down():
    controls = {
        "X-1": [
            candidate(1, candidate_confidence="0.3905"),
            candidate(2, candidate_confidence="0.3799", matched=True),
        ]
    }
    detail = A.recovery_breakdown(controls)["detail"][0]
    assert detail["candidate_confidence_gap"] == 0.0106


# --- ties ---------------------------------------------------------------------


def test_identical_displayed_evidence_is_reported_as_a_tie():
    controls = {
        "X-1": [
            candidate(1, candidate_confidence="0.585", evidence_confidence="0.999405"),
            candidate(
                2, candidate_confidence="0.585", evidence_confidence="0.999405", matched=True
            ),
            candidate(3, candidate_confidence="0.585", evidence_confidence="0.999405"),
        ]
    }
    result = A.tie_breakdown(controls)
    assert result["n_rank1_tied"] == 1
    assert result["max_tie_size"] == 3
    assert result["n_truth_inside_tie_below_rank_1"] == 1


def test_a_distinguishable_rank_1_is_not_a_tie():
    controls = {
        "X-1": [
            candidate(1, candidate_confidence="0.60"),
            candidate(2, candidate_confidence="0.40", matched=True),
        ]
    }
    assert A.tie_breakdown(controls)["n_rank1_tied"] == 0


def test_a_tie_below_the_displayed_precision_still_counts_as_a_tie():
    """Two candidates a reader cannot tell apart are tied, whatever the raw floats."""

    controls = {
        "X-1": [
            candidate(1, candidate_confidence="0.58501"),
            candidate(2, candidate_confidence="0.58502", matched=True),
        ]
    }
    assert A.tie_breakdown(controls)["n_rank1_tied"] == 1


# --- evidence basis -----------------------------------------------------------


def test_the_basis_tally_separates_all_rows_from_rank_1_rows():
    controls = {
        "X-1": [candidate(1, basis="unavailable"), candidate(2, basis="emitted")],
        "X-2": [candidate(1, basis="emitted"), candidate(2, basis="backfilled")],
    }
    result = A.evidence_basis_counts(controls)
    assert result["n_rows"] == 4
    assert result["all_rows"] == {"backfilled": 1, "emitted": 2, "unavailable": 1}
    assert result["rank_1_rows"] == {"emitted": 1, "unavailable": 1}


# --- the spectral section is optional ----------------------------------------


def test_the_spectral_section_records_why_it_was_skipped():
    """Without rdkit the audit reports a reason rather than a silent zero."""

    controls = {"X-1": [candidate(1, matched=True)]}
    result = A.spectral_axis_audit(controls, features={})
    assert result["status"] in {"ok", "skipped"}
    if result["status"] == "skipped":
        assert result["reason"]
    else:
        # No feature rows were supplied, so there is nothing to assess.
        assert result["n_assessed"] == 0


def test_a_summary_serialises_to_json():
    """The audit writes JSON, so every value it produces has to be serialisable."""

    controls = {
        "X-1": [
            candidate(1, formula="C7H8O3S"),
            candidate(2, formula="C7H8O3S", basis="unavailable", matched=True),
        ]
    }
    payload = {
        "recovery": A.recovery_breakdown(controls),
        "ties": A.tie_breakdown(controls),
        "basis": A.evidence_basis_counts(controls),
        "spectral": A.spectral_axis_audit(controls, features={}),
    }
    assert json.loads(json.dumps(payload))["recovery"]["n_truth_below_rank_1"] == 1
