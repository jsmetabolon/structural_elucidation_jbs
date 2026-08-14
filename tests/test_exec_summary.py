"""Tests for the executive-summary export.

The column a reader will act on is ``matched_truth``, so the tests concentrate
on the three ways it can be wrong: claiming a match the evidence does not
support, missing one it does, and reporting a false where the question cannot be
asked at all. The formula comparison gets its own set because the references and
the pipeline do not agree on element ordering, and a naive string compare turns
an exact match into a mismatch.
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


E = load_tool("build_exec_summary")


# --- formula canonicalisation ------------------------------------------------


def test_element_order_does_not_change_a_formula():
    """The references write C6H8SO7 where the pipeline writes C6H8O7S."""

    assert E.normalise_formula("C6H8SO7") == E.normalise_formula("C6H8O7S")


def test_an_implicit_one_is_the_same_as_an_explicit_one():
    assert E.normalise_formula("CH4") == E.normalise_formula("C1H4")


def test_two_element_symbols_are_not_split():
    """Cl is chlorine, not carbon and lithium."""

    assert E.normalise_formula("CCl4") != E.normalise_formula("C2L1I1")
    assert E.normalise_formula("CCl4") == E.normalise_formula("Cl4C")


def test_a_different_formula_still_differs():
    assert E.normalise_formula("C6H8O7S") != E.normalise_formula("C6H8O8S")


def test_something_unparseable_does_not_collapse_into_a_match():
    """A fallback that normalised everything to empty would match everything."""

    assert E.normalise_formula("not a formula") != E.normalise_formula("also not one")


# --- cohort assignment -------------------------------------------------------


def candidate(rank: int, tanimoto: float | None = None, formula: str = "C6H12O6"):
    return {
        "rank": rank,
        "smiles": f"C{rank}",
        "formula": formula,
        "candidate_confidence": 0.5,
        "rank_score_relative": 0.1,
        "evidence_confidence": 0.9,
        "mass_match_score": 0.99,
        "abs_mass_error_da": 0.0001,
        "truth_tanimoto": tanimoto,
    }


def report(*candidates, compound="X-1"):
    return {
        "x_compound_id": compound,
        "display_name": "a feature",
        "msi_level": 3,
        "msi_label": "label",
        "candidates": list(candidates),
    }


def test_a_compound_with_a_reference_structure_is_a_truth_control():
    context = E.classify(report(candidate(1, 0.9)), {}, {"X-1": "glucose"})
    assert context["cohort"] == E.COHORT_TRUTH
    assert context["truth_basis"] == E.BASIS_STRUCTURE
    assert context["reference_answer"] == "glucose"


def test_a_2c_compound_with_a_formula_falls_back_to_the_formula():
    levels = {"X-1": {"sublevels": {"2C"}, "formulae": {"C6H12O6"}, "names": {"glucose-ish"}}}
    context = E.classify(report(candidate(1)), levels, {})
    assert context["cohort"] == E.COHORT_2C
    assert context["truth_basis"] == E.BASIS_FORMULA


def test_a_2c_compound_without_a_formula_has_no_basis_at_all():
    """The tentative name alone cannot be checked: no candidate carries a name."""

    levels = {"X-1": {"sublevels": {"2C"}, "formulae": set(), "names": {"N-formylalanine"}}}
    context = E.classify(report(candidate(1)), levels, {})
    assert context["cohort"] == E.COHORT_2C
    assert context["truth_basis"] == E.BASIS_NONE


def test_structural_truth_wins_over_a_2c_annotation():
    """A structure is the stronger reference, so it is the one used."""

    levels = {"X-1": {"sublevels": {"2C"}, "formulae": {"C6H12O6"}, "names": {"x"}}}
    context = E.classify(report(candidate(1, 0.9)), levels, {"X-1": "glucose"})
    assert context["truth_basis"] == E.BASIS_STRUCTURE


def test_a_compound_with_neither_reference_is_not_assessable():
    context = E.classify(report(candidate(1)), {}, {})
    assert context["cohort"] == E.COHORT_OTHER
    assert context["truth_basis"] == E.BASIS_NONE


# --- the match flag ----------------------------------------------------------


STRUCTURE = {"truth_basis": E.BASIS_STRUCTURE, "reference_formula": ""}
FORMULA = {"truth_basis": E.BASIS_FORMULA, "reference_formula": "C6H8SO7"}
NOTHING = {"truth_basis": E.BASIS_NONE, "reference_formula": ""}


def test_the_structural_bar_is_the_published_threshold():
    assert E.match_flag(candidate(1, E.TRUTH_THRESHOLD), STRUCTURE) is True
    assert E.match_flag(candidate(1, E.TRUTH_THRESHOLD - 0.01), STRUCTURE) is False


def test_a_formula_match_ignores_element_ordering():
    assert E.match_flag(candidate(1, formula="C6H8O7S"), FORMULA) is True


def test_a_formula_mismatch_is_false_not_blank():
    assert E.match_flag(candidate(1, formula="C6H8O8S"), FORMULA) is False


def test_an_unassessable_candidate_is_none_rather_than_false():
    """None becomes an empty cell, so nobody reads it as a failed match."""

    assert E.match_flag(candidate(1, 0.99), NOTHING) is None


# --- rollup ------------------------------------------------------------------


def test_the_first_matching_rank_is_recorded_and_marked_once():
    rows, rollup = E.build_rows(
        report(candidate(1, 0.1), candidate(2, 0.95), candidate(3, 0.99)),
        {**STRUCTURE, "cohort": E.COHORT_TRUTH, "reference_answer": ""},
        {},
    )
    assert rollup["first_matching_rank"] == 2
    assert rollup["matched_at_rank_1"] == "FALSE"
    assert rollup["matched_within_top_3"] == "TRUE"
    assert [row["is_best_match_in_top10"] for row in rows] == ["FALSE", "TRUE", "FALSE"]


def test_a_match_deeper_than_three_counts_only_at_ten():
    candidates = [candidate(rank, 0.1) for rank in range(1, 6)] + [candidate(6, 0.95)]
    _, rollup = E.build_rows(
        report(*candidates), {**STRUCTURE, "cohort": E.COHORT_TRUTH, "reference_answer": ""}, {}
    )
    assert rollup["matched_within_top_3"] == "FALSE"
    assert rollup["matched_within_top_10"] == "TRUE"


def test_an_unassessable_compound_reports_blank_not_false():
    """Blank keeps the 3,642 compounds with no reference out of any match rate."""

    _, rollup = E.build_rows(
        report(candidate(1)), {**NOTHING, "cohort": E.COHORT_OTHER, "reference_answer": ""}, {}
    )
    assert rollup["matched_at_rank_1"] == ""
    assert rollup["matched_within_top_10"] == ""


def test_unassessable_compounds_are_excluded_from_the_rates():
    rollups = [
        {"cohort": E.COHORT_TRUTH, "matched_at_rank_1": "TRUE", "truth_basis": E.BASIS_STRUCTURE},
        {"cohort": E.COHORT_TRUTH, "matched_at_rank_1": "FALSE", "truth_basis": E.BASIS_STRUCTURE},
        {"cohort": E.COHORT_OTHER, "matched_at_rank_1": "", "truth_basis": E.BASIS_NONE},
    ]
    for row in rollups:
        row.setdefault("matched_within_top_3", row["matched_at_rank_1"])
        row.setdefault("matched_within_top_10", row["matched_at_rank_1"])
    summary = E.summarise(rollups)
    assert summary["truth_control"]["rank_1"]["n_assessable"] == 2
    assert summary["truth_control"]["rank_1"]["rate"] == 0.5
    assert summary["other"]["n"] == 1
