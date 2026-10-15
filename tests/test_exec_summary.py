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


def test_the_truncated_sheet_keeps_the_reviewed_columns():
    """The slim top-10 matches the hand-edited executive_summary_top10 sheet."""

    assert E.TRUNCATED_COLUMNS == [
        "x_compound_id",
        "metabolon_level",
        "metabolon_sublevel",
        "known_x_tentative_information",
        "known_x_cleaned_external_info",
        "msi_level",
        "msi_label",
        "feature_name",
        "report_display_name",
        "reported_level",
        "reported_identification",
        "identity_candidate_rank",
        "rank",
        "candidate_smiles",
        "candidate_formula",
        "evidence_confidence",
        "mass_match_score",
        "abs_mass_error_ppm",
        "is_best_match_in_top10",
        "matching_ranks",
    ]
    assert set(E.TRUNCATED_COLUMNS) <= set(E.LONG_COLUMNS)


def test_the_displayed_metabolon_level_is_the_best_across_libraries():
    """Level 3 in one library and Level 5 in another is reported as 3."""

    entry = {
        "level_rows": {("3", "3A"), ("5", "5"), ("", "")},
    }
    assert E.best_metabolon(entry) == ("3", "3A")


def test_a_definite_level_beats_the_same_number_with_a_question_mark():
    assert E.best_metabolon({"level_rows": {("4?", "4?"), ("4", "4")}}) == ("4", "4")


def test_sublevels_at_the_winning_level_are_all_kept():
    assert E.best_metabolon({"level_rows": {("3", "3C"), ("3", "3A"), ("5", "5")}}) == ("3", "3A; 3C")


def test_a_label_only_wins_when_there_is_no_numbered_level():
    assert E.best_metabolon({"level_rows": {("Exclude", "Exclude"), ("Artifact", "Artifact")}})[0] == "Artifact"
    assert E.best_metabolon({"level_rows": set()}) == ("", "")


# --- "Isomer of" labelling for Naming v2 collisions ---------------------------


def v2_fields(display_name: str, *, rank1_reason: str = E.COLLISION_REASON, **extra):
    """A Naming v2 confidence_inputs record whose rank-1 collided with the library."""

    matched = [{"chemical_id": "100000443", "name": "phenylglucopyranoside", "connectivity": "NEZJDVYDSZTRFS"}]
    fields = {
        "naming_policy_version": "v2",
        "display_name": display_name,
        "name_source": "ranked_candidate",
        "exception_basis": "",
        "raw_rank1_candidate": {
            "model_rank": 1,
            "identity_exclusion_reason": rank1_reason,
            "matched_library": matched if rank1_reason == E.COLLISION_REASON else [],
        },
        "excluded_library_matches": [
            {"model_rank": 1, "reason": rank1_reason, "matches": matched if rank1_reason == E.COLLISION_REASON else []},
            {"model_rank": 4, "reason": E.COLLISION_REASON,
             "matches": [{"chemical_id": "9", "name": "arbutin", "connectivity": "AAAA"}]},
        ],
    }
    fields.update(extra)
    return fields


def test_a_withheld_formula_becomes_isomer_of_the_colliding_library_metabolite():
    """X-02269: rank 1 is phenylglucopyranoside, which is already in the library."""

    fields = v2_fields("C12H16O6; structure unresolved")
    assert E.isomer_of_feature_name(fields["display_name"], fields) == "Isomer of phenylglucopyranoside"


def test_an_exact_mass_placeholder_is_also_labelled_as_an_isomer():
    fields = v2_fields("Exact mass 264.1242; unidentified", name_source="formula_or_mass")
    assert E.isomer_of_feature_name(fields["display_name"], fields) == "Isomer of phenylglucopyranoside"


def test_tentative_structure_is_also_labelled_as_an_isomer():
    fields = v2_fields("Tentative Structure")
    assert E.isomer_of_feature_name(fields["display_name"], fields) == "Isomer of phenylglucopyranoside"


def test_a_curated_known_x_hypothesis_is_not_overridden():
    """A ``Possible *`` label from known_X is the report's settled name."""

    fields = v2_fields("Possible 5,6-dihydroxyindole sulfate", name_source="known_x_tentative")
    assert E.isomer_of_feature_name(fields["display_name"], fields) == "Possible 5,6-dihydroxyindole sulfate"


def test_an_authentic_standard_identity_is_not_overridden():
    fields = v2_fields("N-(4-aminobenzoyl)-glutamic acid", name_source="authentic_standard",
                       exception_basis="standard_confirmed")
    assert E.isomer_of_feature_name(fields["display_name"], fields) == "N-(4-aminobenzoyl)-glutamic acid"


def test_a_formula_placeholder_without_a_rank1_collision_stays_a_formula():
    """Deeper-ranked collisions do not make the compound an isomer of anything."""

    fields = v2_fields("C12H16O6; structure unresolved", rank1_reason="")
    assert E.isomer_of_feature_name(fields["display_name"], fields) == "C12H16O6; structure unresolved"
    fields = v2_fields("C12H16O6; structure unresolved", rank1_reason="mass_gate_failed")
    assert E.isomer_of_feature_name(fields["display_name"], fields) == "C12H16O6; structure unresolved"


def test_v1_records_are_left_to_the_sheet_feature_name_rule():
    fields = v2_fields("C12H16O6; structure unresolved", naming_policy_version="v1")
    assert E.isomer_of_feature_name(fields["display_name"], fields) == "C12H16O6; structure unresolved"


def test_stereoisomers_sharing_a_connectivity_are_all_named():
    fields = v2_fields("C26H43NO5; structure unresolved")
    fields["raw_rank1_candidate"]["matched_library"] = [
        {"chemical_id": "1", "name": "glycohyodeoxycholic acid", "connectivity": "K"},
        {"chemical_id": "2", "name": "glycoursodeoxycholic acid", "connectivity": "K"},
        {"chemical_id": "3", "name": "Glycohyodeoxycholic acid", "connectivity": "K"},
    ]
    assert E.isomer_of_feature_name(fields["display_name"], fields) == (
        "Isomer of glycohyodeoxycholic acid / glycoursodeoxycholic acid"
    )


def test_a_long_stereoisomer_list_is_capped_in_the_label():
    """X-10358 collides with six pentonate stereoisomers; the label stays readable."""

    fields = v2_fields("C5H10O6; structure unresolved")
    fields["raw_rank1_candidate"]["matched_library"] = [
        {"chemical_id": str(i), "name": name, "connectivity": "K"}
        for i, name in enumerate(["ribonate", "xylonate", "arabonate", "lyxonate", "arabonate/xylonate"])
    ]
    assert E.isomer_of_feature_name(fields["display_name"], fields) == (
        "Isomer of ribonate / xylonate / arabonate (+2 more)"
    )


def test_excluded_library_reference_lists_every_collision_once():
    fields = v2_fields("C12H16O6; structure unresolved")
    assert E.excluded_library_reference(fields) == "phenylglucopyranoside; arbutin"
    fields["excluded_library_matches"].append(
        {"model_rank": 7, "reason": "mass_gate_failed", "matches": [{"name": "not a collision"}]}
    )
    assert E.excluded_library_reference(fields) == "phenylglucopyranoside; arbutin"


# --- production library panel -------------------------------------------------


PANEL = {
    "ids": {"100000443": "phenylglucopyranoside", "9": "arbutin", "77": "N-(4-aminobenzoyl)-glutamic acid"},
    "names": {"phenylglucopyranoside": "100000443", "arbutin": "9", "n-(4-aminobenzoyl)-glutamic acid": "77"},
}


def test_an_isomer_label_is_marked_by_the_rank1_collision_chemical_id():
    fields = v2_fields("C12H16O6; structure unresolved")
    flag, matches = E.prod_library_marks("Isomer of phenylglucopyranoside", fields, PANEL)
    assert flag == "TRUE"
    assert matches == "phenylglucopyranoside (100000443); arbutin (9)"


def test_an_isomer_label_whose_reference_is_off_panel_is_false():
    fields = v2_fields("C12H16O6; structure unresolved")
    obscure = [{"chemical_id": "555", "name": "obscurine", "connectivity": "Q"}]
    fields["raw_rank1_candidate"]["matched_library"] = obscure
    fields["excluded_library_matches"][0]["matches"] = obscure
    flag, matches = E.prod_library_marks("Isomer of obscurine", fields, PANEL)
    assert flag == "FALSE"
    assert matches == "arbutin (9)"  # the deeper excluded match is still listed


def test_a_bare_or_hedged_name_is_matched_by_name():
    fields = v2_fields("N-(4-aminobenzoyl)-glutamic acid", name_source="authentic_standard",
                       exception_basis="standard_confirmed", rank1_reason="")
    assert E.prod_library_marks("N-(4-aminobenzoyl)-glutamic acid", fields, PANEL)[0] == "TRUE"
    assert E.prod_library_marks("Possible N-(4-aminobenzoyl)-glutamic acid", fields, PANEL)[0] == "TRUE"
    assert E.prod_library_marks("Possible 5,6-dihydroxyindole sulfate", fields, PANEL)[0] == "FALSE"


def test_a_level1_identity_is_matched_by_its_resolved_library_chemical_id():
    """asparaginylasparagine is 'asn-asn' on the panel; the ID bridges the spelling."""

    fields = v2_fields("asparaginylasparagine", name_source="authentic_standard",
                       exception_basis="standard_confirmed", rank1_reason="")
    panel = {"ids": {"100019920": "asn-asn"}, "names": {"asn-asn": "100019920"}}
    flag, matches = E.prod_library_marks("asparaginylasparagine", fields, panel, db_chemical_id="100019920.0")
    assert (flag, matches) == ("TRUE", "asn-asn (100019920)")
    assert E.prod_library_marks("asparaginylasparagine", fields, panel)[0] == "FALSE"


def test_placeholders_leave_the_panel_flag_blank():
    fields = v2_fields("C12H16O6; structure unresolved", rank1_reason="")
    for label in ("C12H16O6; structure unresolved", "Exact mass 153.0197; unidentified", "Tentative Structure"):
        assert E.prod_library_marks(label, fields, PANEL)[0] == ""


def test_no_panel_means_blank_columns():
    fields = v2_fields("C12H16O6; structure unresolved")
    assert E.prod_library_marks("Isomer of phenylglucopyranoside", fields, {"ids": {}, "names": {}}) == ("", "")


def test_capped_isomer_labels_still_resolve_their_listed_names():
    assert E._label_names("Isomer of ribonate / xylonate / arabonate (+2 more)") == [
        "ribonate", "xylonate", "arabonate"
    ]
    assert E._label_names("Possible isobutyl phthalate") == ["isobutyl phthalate"]
    assert E._label_names("C5H10O6; structure unresolved") == []


def test_the_panel_loader_reads_ids_and_names(tmp_path):
    path = tmp_path / "panel.csv"
    path.write_text("CHEMICAL_ID,CHEMICAL_NAME,SUPER_PATHWAY,PATHWAY\n100000443.0,phenylglucopyranoside,Xenobiotics,Food\n")
    panel = E.load_prod_library(path)
    assert panel["ids"] == {"100000443": "phenylglucopyranoside"}
    assert panel["names"] == {"phenylglucopyranoside": "100000443"}
    assert E.load_prod_library(tmp_path / "missing.xlsx") == {"ids": {}, "names": {}}


def test_known_x_annotation_text_is_carried_onto_every_row():
    entry = {
        "levels": {"3"}, "sublevels": {"3A"},
        "names": {"possible isomer of metabolonic lactone sulfate"},
        "cleaned": {"metabolonic lactone sulfate"},
    }
    rows, rollup = E.build_rows(
        report(candidate(1)), {**NOTHING, "cohort": E.COHORT_OTHER, "reference_answer": ""}, entry
    )
    assert rows[0]["known_x_tentative_information"] == "possible isomer of metabolonic lactone sulfate"
    assert rollup["known_x_cleaned_external_info"] == "metabolonic lactone sulfate"


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


def test_ppm_is_da_error_over_exact_mass():
    ppm = E.abs_mass_error_ppm(0.014066, 254.13789)
    assert abs(ppm - 55.35) < 0.01


def test_ppm_is_blank_without_a_mass():
    assert E.abs_mass_error_ppm(0.014066, None) is None
    assert E.abs_mass_error_ppm(0.014066, 0) is None
    assert E.abs_mass_error_ppm(None, 254.13789) is None


def test_mass_reference_prefers_predicted_exact_mass_over_parent_mz():
    assert E._mass_reference({"predicted_exact_mass": "254.13789", "parent_mz": "253.14468"}) == 254.13789
    assert E._mass_reference({"parent_mz": "253.14468"}) == 253.14468
    assert E._mass_reference({}, fallback=253.14468) == 253.14468


# --- cohort assignment -------------------------------------------------------


def candidate(
    rank: int,
    tanimoto: float | None = None,
    formula: str = "C6H12O6",
    evidence: float | None = 0.9,
):
    return {
        "rank": rank,
        "smiles": f"C{rank}",
        "formula": formula,
        "candidate_confidence": 0.5,
        "rank_score_relative": 0.1,
        "evidence_confidence": evidence,
        "mass_match_score": 0.99,
        "abs_mass_error_da": 0.0001,
        "abs_mass_error_ppm": 0.4,
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


# --- evidence_confidence display backfill ------------------------------------


def test_emitted_evidence_is_kept_as_is():
    value, basis = E._display_evidence({"evidence_confidence": "0.5", "formula_prior_score": "1.0"})
    assert (value, basis) == (0.5, E.EVIDENCE_EMITTED)


def test_zero_evidence_is_backfilled_from_the_formula_prior():
    """The truth controls read 0.0 while carrying a formula prior of 1.0."""

    value, basis = E._display_evidence({"evidence_confidence": "0.0", "formula_prior_score": "1.0"})
    assert value == 1.0
    assert basis == E.EVIDENCE_BACKFILLED


def test_no_computable_evidence_is_blank_not_zero():
    value, basis = E._display_evidence({"evidence_confidence": "0.0"})
    assert value is None
    assert basis == E.EVIDENCE_UNAVAILABLE


def test_formula_and_mass_hypothesis_recompute_a_prior():
    """The protected_peptide_permutation rank-1s (e.g. X-21410) carry only a
    formula and a mass hypothesis; the canonical prior must still recover them."""

    value, basis = E._display_evidence(
        {
            "evidence_confidence": "0.0",
            "formula_prior_score": "0.0",
            "predicted_formula": "C19H30O6S",
            "mass_hypothesis_neutral_mass": "386.176916466621",
        }
    )
    assert basis == E.EVIDENCE_BACKFILLED
    assert value > 0.99


def test_ms1_support_adds_the_small_bonus():
    value, _ = E._display_evidence(
        {"evidence_confidence": "", "ms2_match_score": "0.42", "ms1_buddy_supported": "True"}
    )
    assert abs(value - 0.47) < 1e-9


def test_display_evidence_does_not_drift_from_the_canonical_definition():
    """The pure-stdlib display path must match chemistry/evidence.py exactly."""

    import pandas as pd  # noqa: PLC0415
    from structural_elucidation.chemistry.evidence import backfill_evidence_confidence  # noqa: PLC0415

    rows = [
        {"evidence_confidence": "0.0", "formula_prior_score": "1.0"},
        {"evidence_confidence": "0.0", "formula_prior_score": "0.968185"},
        {"evidence_confidence": "0.5", "formula_prior_score": "1.0"},
        {"evidence_confidence": "0.0", "ms2_match_score": "0.42", "ms1_buddy_supported": "True"},
        {"evidence_confidence": "0.0", "predicted_formula": "C19H30O6S", "mass_hypothesis_neutral_mass": "386.1769"},
        {"evidence_confidence": "0.0"},
    ]
    canonical, _ = backfill_evidence_confidence(pd.DataFrame(rows))
    for index, row in enumerate(rows):
        value, _ = E._display_evidence(row)
        assert abs((value or 0.0) - float(canonical.iloc[index])) < 1e-6


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


def test_published_confidence_follows_rank_and_keeps_ties():
    """A later candidate cannot publish a higher confidence than an earlier one."""

    rows, rollup = E.build_rows(
        report(
            candidate(1, evidence=0.60),
            candidate(2, evidence=0.60),
            candidate(3, evidence=0.90),
            candidate(4, evidence=0.40),
        ),
        {**NOTHING, "cohort": E.COHORT_OTHER, "reference_answer": ""},
        {},
    )
    assert [row["evidence_confidence"] for row in rows] == ["0.6", "0.6", "0.6", "0.4"]
    assert rollup["rank1_evidence_confidence"] == "0.6"


def test_the_first_matching_rank_is_recorded_and_marked_once():
    rows, rollup = E.build_rows(
        report(candidate(1, 0.1), candidate(2, 0.95), candidate(3, 0.99)),
        {**STRUCTURE, "cohort": E.COHORT_TRUTH, "reference_answer": ""},
        {},
    )
    assert rollup["first_matching_rank"] == 2
    assert rollup["matching_ranks"] == "2"
    assert rollup["matched_at_rank_1"] == "FALSE"
    assert rollup["matched_within_top_3"] == "TRUE"
    assert [row["is_best_match_in_top10"] for row in rows] == ["FALSE", "TRUE", "FALSE"]
    assert [row["matching_ranks"] for row in rows] == ["2", "2", "2"]
    assert [row["matched_truth"] for row in rows] == ["FALSE", "TRUE", "TRUE"]
    assert rows[0]["abs_mass_error_ppm"] == "0.4"
    assert rollup["rank1_abs_mass_error_ppm"] == "0.4"


def test_a_match_deeper_than_three_counts_only_at_ten():
    candidates = [candidate(rank, 0.1) for rank in range(1, 6)] + [candidate(6, 0.95)]
    _, rollup = E.build_rows(
        report(*candidates), {**STRUCTURE, "cohort": E.COHORT_TRUTH, "reference_answer": ""}, {}
    )
    assert rollup["matched_within_top_3"] == "FALSE"
    assert rollup["matched_within_top_10"] == "TRUE"


def test_an_unassessable_compound_reports_blank_not_false():
    """Blank keeps the 3,642 compounds with no reference out of any match rate."""

    rows, rollup = E.build_rows(
        report(candidate(1)), {**NOTHING, "cohort": E.COHORT_OTHER, "reference_answer": ""}, {}
    )
    assert rollup["matched_at_rank_1"] == ""
    assert rollup["matched_within_top_10"] == ""
    assert rollup["matching_ranks"] == ""
    assert rows[0]["matching_ranks"] == ""


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


# --- white-paper reported identification -------------------------------------


def named_candidate(rank, smiles, formula="C6H12O6", evidence=0.99, mass=0.99):
    row = candidate(rank, formula=formula)
    row["smiles"] = smiles
    row["evidence_confidence"] = evidence
    row["mass_match_score"] = mass
    return row


def test_no_supported_formula_is_level_5():
    assert E.reported_identification("X-1", [named_candidate(1, "CCO", formula="")]) == (
        "Level 5: exact mass of interest"
    )


def test_a_split_formula_is_level_5():
    rows = [
        named_candidate(1, "CCO", formula="C2H6O"),
        named_candidate(2, "CCO", formula="C6H12O6"),
    ]
    assert E.reported_identification("X-1", rows) == "Level 5: exact mass of interest"


def test_an_authentic_standard_keeps_its_name_without_a_formula():
    call = E.reported_identification(
        "X-1",
        [named_candidate(1, "CCO", formula="")],
        {"metabolon_level": "1", "report_display_name": "lactoyl-isoleucine"},
    )
    assert call == "Level 1: lactoyl-isoleucine"


def test_an_authentic_standard_is_level_1():
    call = E.reported_identification(
        "X-1",
        [named_candidate(1, "CCO"), named_candidate(2, "c1ccccc1")],
        {"metabolon_level": "1", "report_display_name": "ethanol"},
    )
    assert call == "Level 1: ethanol"


def test_a_library_collision_is_an_isomer_at_level_2():
    call = E.reported_identification(
        "X-02269",
        [named_candidate(1, "CCO"), named_candidate(2, "c1ccccc1")],
        {"feature_name": "Isomer of phenylglucopyranoside"},
    )
    assert call == "Level 2: isomer of phenylglucopyranoside"


def test_one_connectivity_with_the_reports_own_name_is_level_2_star():
    call = E.reported_identification(
        "X-9",
        [named_candidate(1, "CCO"), named_candidate(2, "CCO")],
        {"report_display_name": "ethanol"},
    )
    assert call == "Level 2: [X-9], Probable ethanol"


def test_one_connectivity_named_only_by_the_annotation_is_level_2():
    call = E.reported_identification(
        "X-9",
        [named_candidate(1, "CCO"), named_candidate(2, "CCO")],
        {
            "report_display_name": "C2H6O; structure unresolved",
            "known_x_tentative_information": "ethanol",
        },
    )
    assert call == "Level 2: [X-9], Probable ethanol"


def test_a_single_structure_that_stands_apart_is_level_2():
    call = E.reported_identification(
        "X-9",
        [
            named_candidate(1, "CCO", evidence=0.90),
            named_candidate(2, "c1ccccc1", evidence=0.20),
        ],
        {"known_x_tentative_information": "ethanol"},
    )
    assert call == "Level 2: [X-9], Probable ethanol"


def test_a_curator_note_is_not_used_as_a_probable_name():
    note = E.reported_identification(
        "X-9",
        [named_candidate(1, "CCO"), named_candidate(2, "CCO")],
        {
            "known_x_tentative_information": "somehow overlapping with X-15956",
            "known_x_cleaned_external_info": "Unknown",
        },
    )
    assert note == "Level 4: C6H12O6"
    glucoside = E.reported_identification(
        "X-9",
        [named_candidate(1, "CCO"), named_candidate(2, "CCO")],
        {
            "known_x_tentative_information": "a glucoside, C10H20O7",
            "known_x_cleaned_external_info": "Glucoside",
        },
    )
    assert glucoside == "Level 3: Glucoside"


def test_a_glucuronide_of_a_formula_is_level_2_star_even_when_candidates_disagree():
    call = E.reported_identification(
        "X-9",
        [
            named_candidate(1, "CCO", evidence=0.99),
            named_candidate(2, "c1ccccc1", evidence=0.99),
        ],
        {"known_x_cleaned_external_info": "glucuronide of C20H26O4"},
    )
    assert call == "Level 2*: glucuronide of C20H26O4*"


def test_library_level_4_uncertain_is_msi_4():
    assert E.msi_from_library_level("4?") == ("4", "Molecular formula")
    assert E.msi_from_library_level("1") == ("1", "Confirmed (authentic standard)")
    assert E.msi_from_library_level("Exclude") == ("", "")


def test_a_class_phrase_is_level_3_when_the_candidates_are_one_structure():
    call = E.reported_identification(
        "X-9",
        [named_candidate(1, "CCO"), named_candidate(2, "CCO")],
        {"known_x_cleaned_external_info": "sulfate of C9H9NO2"},
    )
    assert call == "Level 2*: sulfate of C9H9NO2*"


def test_a_shared_sulfate_is_a_level_3_class():
    pytest = __import__("pytest")
    pytest.importorskip("rdkit")
    call = E.reported_identification(
        "X-9",
        [
            named_candidate(1, "CCCOS(=O)(=O)O", formula="C3H8O4S"),
            named_candidate(2, "CC(C)OS(=O)(=O)O", formula="C3H8O4S"),
        ],
        {},
    )
    assert call == "Level 3: sulfate"


def test_several_unrelated_structures_with_one_formula_are_level_4():
    call = E.reported_identification(
        "X-9",
        [
            named_candidate(1, "CCO", evidence=0.99),
            named_candidate(2, "c1ccccc1", evidence=0.99),
        ],
        {},
    )
    assert call == "Level 4: C6H12O6"


def test_a_name_is_level_2_even_when_the_candidates_disagree():
    """A name is the level. The candidates do not have to pick one structure."""

    call = E.reported_identification(
        "X-9",
        [
            named_candidate(1, "CCO", evidence=0.998),
            named_candidate(2, "c1ccccc1", evidence=0.998),
        ],
        {"known_x_tentative_information": "ethanol"},
    )
    assert call == "Level 2: [X-9], Probable ethanol"


def test_a_named_annotation_is_level_2_without_a_candidate_family():
    call = E.reported_identification(
        "X-9",
        [
            named_candidate(1, "Oc1ccccc1O", evidence=0.99),
            named_candidate(2, "Oc1ccc(O)cc1", evidence=0.99),
        ],
        {"known_x_tentative_information": "dihydroxybenzene"},
    )
    assert call == "Level 2: [X-9], Probable dihydroxybenzene"


def test_build_rows_copies_the_publication_call_onto_every_candidate():
    rows, rollup = E.build_rows(
        report(
            named_candidate(1, "CCO"),
            named_candidate(2, "c1ccccc1"),
            compound="X-02269",
        ),
        {**NOTHING, "cohort": E.COHORT_OTHER, "reference_answer": ""},
        {},
    )
    assert rows[0]["reported_level"] == "4"
    assert rows[0]["reported_identification"] == "Level 4: C6H12O6"
    assert rows[1]["reported_identification"] == rows[0]["reported_identification"]
    assert rollup["reported_identification"] == rows[0]["reported_identification"]
