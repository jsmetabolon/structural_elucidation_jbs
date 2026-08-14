"""Tests for the executive-summary builder.

The column that carries the risk is ``matched_truth``. Getting it wrong in one
direction invents a recovery that did not happen; getting it wrong in the other
buries one that did. The cases below pin the three bases apart, and in
particular pin the difference between "did not match" and "cannot be asked".
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
    """The annotations and the pipeline disagree on ordering, and mean the same thing.

    A reference recorded as C6H8SO7 against a candidate written C6H8O7S is an
    exact match, and comparing the raw strings would report it as a miss.
    """

    assert E.normalise_formula("C6H8SO7") == E.normalise_formula("C6H8O7S")


def test_an_implicit_one_is_the_same_as_an_explicit_one():
    assert E.normalise_formula("CH4") == E.normalise_formula("C1H4")


def test_a_genuinely_different_formula_still_differs():
    assert E.normalise_formula("C6H8O7S") != E.normalise_formula("C6H8O8S")


def test_repeated_elements_are_summed():
    assert E.normalise_formula("CH3CH3") == E.normalise_formula("C2H6")


def test_an_unparseable_formula_does_not_collapse_into_a_match():
    """Junk must not normalise to something another piece of junk also reaches."""

    assert E.normalise_formula("not a formula") != E.normalise_formula("also not one")


# --- cohorts -----------------------------------------------------------------


def candidate(rank, tanimoto=None, formula="C6H12O6"):
    return {
        "rank": rank,
        "smiles": "C",
        "formula": formula,
        "candidate_confidence": 0.5,
        "rank_score_relative": 0.1,
        "evidence_confidence": 0.9,
        "mass_match_score": 0.99,
        "abs_mass_error_da": 0.001,
        "truth_tanimoto": tanimoto,
    }


def report(*candidates, compound="X-1"):
    return {
        "x_compound_id": compound,
        "display_name": "a feature",
        "msi_level": 3,
        "msi_label": "Structure",
        "candidates": list(candidates),
    }


def test_a_compound_with_a_reference_structure_is_a_truth_control():
    context = E.classify(report(candidate(1, 0.9)), {}, {"X-1": "glucose"})
    assert context["cohort"] == E.COHORT_TRUTH
    assert context["truth_basis"] == E.BASIS_STRUCTURE
    assert context["reference_answer"] == "glucose"


def test_a_2c_compound_with_a_formula_falls_back_to_formula():
    levels = {"X-1": {"sublevels": {"2C"}, "names": {"N-formylalanine"}, "formulae": {"C4H7NO3"}}}
    context = E.classify(report(candidate(1)), levels, {})
    assert context["cohort"] == E.COHORT_2C
    assert context["truth_basis"] == E.BASIS_FORMULA


def test_a_2c_compound_without_a_formula_has_no_basis_at_all():
    """Most of them. A tentative name alone cannot be matched against a SMILES."""

    levels = {"X-1": {"sublevels": {"2C"}, "names": {"N-formylalanine"}, "formulae": set()}}
    context = E.classify(report(candidate(1)), levels, {})
    assert context["cohort"] == E.COHORT_2C
    assert context["truth_basis"] == E.BASIS_NONE


def test_structural_truth_wins_over_a_2c_annotation():
    """Only one compound needs this, but the precedence must be deliberate."""

    levels = {"X-1": {"sublevels": {"2C"}, "names": {"x"}, "formulae": {"C4H7NO3"}}}
    context = E.classify(report(candidate(1, 0.95)), levels, {"X-1": "known"})
    assert context["truth_basis"] == E.BASIS_STRUCTURE


def test_an_unannotated_compound_is_not_assessable():
    context = E.classify(report(candidate(1)), {}, {})
    assert context["cohort"] == E.COHORT_OTHER
    assert context["truth_basis"] == E.BASIS_NONE


# --- the match flag ----------------------------------------------------------


def test_a_structural_match_is_scored_at_the_published_threshold():
    context = {"truth_basis": E.BASIS_STRUCTURE, "reference_formula": ""}
    assert E.match_flag(candidate(1, 0.70), context) is True
    assert E.match_flag(candidate(1, 0.699), context) is False


def test_a_formula_match_ignores_element_order():
    context = {"truth_basis": E.BASIS_FORMULA, "reference_formula": "C6H8SO7"}
    assert E.match_flag(candidate(1, formula="C6H8O7S"), context) is True


def test_no_basis_yields_none_rather_than_false():
    """The distinction the CSV depends on: unassessable is not a failed match."""

    context = {"truth_basis": E.BASIS_NONE, "reference_formula": ""}
    assert E.match_flag(candidate(1), context) is None


# --- rows and rollup ---------------------------------------------------------


def test_the_first_matching_rank_is_recorded_and_marked_once():
    context = {
        "cohort": E.COHORT_TRUTH,
        "truth_basis": E.BASIS_STRUCTURE,
        "reference_answer": "glucose",
        "reference_formula": "",
    }
    rows, rollup = build(context, candidate(1, 0.1), candidate(2, 0.95), candidate(3, 0.99))
    assert rollup["first_matching_rank"] == 2
    assert rollup["matched_at_rank_1"] == "FALSE"
    assert rollup["matched_within_top_3"] == "TRUE"
    assert [row["is_best_match_in_top10"] for row in rows] == ["FALSE", "TRUE", "FALSE"]


def test_the_best_tanimoto_is_reported_even_when_nothing_matched():
    context = {
        "cohort": E.COHORT_TRUTH,
        "truth_basis": E.BASIS_STRUCTURE,
        "reference_answer": "",
        "reference_formula": "",
    }
    _, rollup = build(context, candidate(1, 0.4), candidate(2, 0.65))
    assert rollup["matched_within_top_10"] == "FALSE"
    assert rollup["best_truth_tanimoto"] == "0.65"


def test_an_unassessable_compound_leaves_the_match_columns_empty():
    context = {
        "cohort": E.COHORT_OTHER,
        "truth_basis": E.BASIS_NONE,
        "reference_answer": "",
        "reference_formula": "",
    }
    rows, rollup = build(context, candidate(1), candidate(2))
    assert {row["matched_truth"] for row in rows} == {""}
    assert rollup["matched_at_rank_1"] == ""
    assert rollup["matched_within_top_10"] == ""


def build(context, *candidates):
    return E.build_rows(report(*candidates), context, {"levels": {"1"}, "sublevels": {"KNOWN"}})
