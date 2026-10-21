"""Compound-level mass-gated pathway vote over the short ranked candidate list."""

from __future__ import annotations

import pandas as pd
import pytest

from structural_elucidation.pipeline.enrichment.pathway import assign_compound_pathways


def _row(**over) -> dict:
    base = dict(
        control_id="X-1",
        rank=1,
        predicted_pathway="Amino Acid",
        predicted_subpathway="Lysine Metabolism",
        pathway_assignment_status="exact_reference_match",
        pathway_vote_share=1.0,
        pathway_best_reference_similarity=1.0,
        formula_prior_score=1.0,
        mass_match_score=1.0,
        abs_mass_error_da=0.0003,
        rerank_score_calibrated=0.9,
        rerank_score=0.9,
    )
    base.update(over)
    return base


def _assign(*rows: dict, **kwargs) -> pd.DataFrame:
    options = dict(policy="exact_first", collapse_amino_acid_peptide=False, apply_name_boost=False)
    options.update(kwargs)
    return assign_compound_pathways(pd.DataFrame(list(rows)), **options)


def test_mass_gate_drops_wrong_mass_exact_hit():
    out = _assign(
        _row(rank=1, predicted_pathway="Amino Acid", predicted_subpathway="Lysine Metabolism"),
        _row(
            rank=3,
            predicted_pathway="Xenobiotics",
            predicted_subpathway="Drug - Other",
            pathway_assignment_status="exact_reference_match",
            formula_prior_score=0.23,
            mass_match_score=0.0,
            abs_mass_error_da=1.056,
            rerank_score_calibrated=0.5,
        ),
    )
    assert out["compound_predicted_pathway"].nunique() == 1
    assert out["compound_predicted_pathway"].iloc[0] == "Amino Acid"
    assert out["compound_predicted_subpathway"].iloc[0] == "Lysine Metabolism"
    assert out["compound_pathway_source"].iloc[0] == "exact_mass_consistent"
    assert int(out["compound_pathway_n_voters"].iloc[0]) == 1


def test_exact_first_ignores_neighbors():
    out = _assign(
        _row(
            rank=1,
            predicted_pathway="Amino Acid",
            predicted_subpathway="Lysine Metabolism",
            rerank_score_calibrated=0.40,
        ),
        _row(
            rank=2,
            predicted_pathway="Lipid",
            predicted_subpathway="Fatty Acid Metabolism",
            pathway_assignment_status="exact_reference_match",
            rerank_score_calibrated=0.80,
        ),
        _row(
            rank=3,
            predicted_pathway="Xenobiotics",
            predicted_subpathway="Chemical",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.99,
            pathway_best_reference_similarity=0.95,
            rerank_score_calibrated=0.99,
        ),
    )
    assert out["compound_predicted_pathway"].iloc[0] == "Lipid"
    assert out["compound_predicted_subpathway"].iloc[0] == "Fatty Acid Metabolism"
    assert out["compound_pathway_source"].iloc[0] == "exact_mass_consistent"
    assert int(out["compound_pathway_n_voters"].iloc[0]) == 2
    assert int(out["compound_pathway_agrees_rank1"].iloc[0]) == 0


def test_weighted_neighbors_can_override_rank1():
    out = _assign(
        _row(
            rank=1,
            predicted_pathway="Amino Acid",
            predicted_subpathway="Lysine Metabolism",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.20,
            pathway_best_reference_similarity=0.40,
            rerank_score_calibrated=1.0,
        ),
        _row(
            rank=2,
            predicted_pathway="Lipid",
            predicted_subpathway="Fatty Acid Metabolism",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.90,
            pathway_best_reference_similarity=0.80,
            rerank_score_calibrated=0.90,
        ),
        _row(
            rank=3,
            predicted_pathway="Lipid",
            predicted_subpathway="Fatty Acid Metabolism",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.80,
            pathway_best_reference_similarity=0.70,
            rerank_score_calibrated=0.80,
        ),
    )
    assert out["compound_predicted_pathway"].iloc[0] == "Lipid"
    assert out["compound_pathway_source"].iloc[0] == "weighted_candidate_consensus"
    assert int(out["compound_pathway_agrees_rank1"].iloc[0]) == 0


def test_subpathway_votes_only_under_winning_super():
    out = _assign(
        _row(
            rank=1,
            predicted_pathway="Lipid",
            predicted_subpathway="Fatty Acid Metabolism",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.70,
            rerank_score_calibrated=1.0,
        ),
        _row(
            rank=2,
            predicted_pathway="Lipid",
            predicted_subpathway="Sterol",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.20,
            rerank_score_calibrated=0.50,
        ),
        _row(
            rank=3,
            predicted_pathway="Amino Acid",
            predicted_subpathway="Lysine Metabolism",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.50,
            rerank_score_calibrated=0.30,
        ),
    )
    assert out["compound_predicted_pathway"].iloc[0] == "Lipid"
    assert out["compound_predicted_subpathway"].iloc[0] == "Fatty Acid Metabolism"


def test_empty_eligibility_leaves_structure_pathway():
    out = _assign(
        _row(
            predicted_pathway="Amino Acid",
            predicted_subpathway="Lysine Metabolism",
            formula_prior_score=0.1,
            mass_match_score=0.1,
            abs_mass_error_da=1.0,
        ),
        _row(
            rank=2,
            predicted_pathway="",
            predicted_subpathway="",
            pathway_assignment_status="no_reference_match",
            formula_prior_score=1.0,
            mass_match_score=1.0,
            abs_mass_error_da=0.0,
        ),
    )
    assert out["predicted_pathway"].iloc[0] == "Amino Acid"
    assert out["compound_predicted_pathway"].iloc[0] == ""
    assert out["compound_pathway_source"].iloc[0] == "none"
    assert out["compound_pathway_status"].iloc[0] == "no_eligible_candidates"
    assert int(out["compound_pathway_n_voters"].iloc[0]) == 0
    assert pd.isna(out["compound_pathway_agrees_rank1"].iloc[0])


def test_broadcasts_compound_columns_across_ranks():
    out = _assign(
        _row(rank=1),
        _row(
            rank=2,
            predicted_pathway="Lipid",
            predicted_subpathway="Fatty Acid Metabolism",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.4,
            formula_prior_score=0.2,
            mass_match_score=0.2,
            abs_mass_error_da=0.9,
        ),
    )
    assert len(out) == 2
    assert out["compound_predicted_pathway"].tolist() == ["Amino Acid", "Amino Acid"]
    assert out["compound_predicted_subpathway"].tolist() == ["Lysine Metabolism", "Lysine Metabolism"]
    assert out["compound_pathway_n_voters"].tolist() == [1, 1]


def test_all_eligible_lets_neighbors_outvote_exact():
    out = _assign(
        _row(
            rank=1,
            predicted_pathway="Lipid",
            predicted_subpathway="Fatty Acid Metabolism",
            pathway_assignment_status="exact_reference_match",
            rerank_score_calibrated=0.20,
        ),
        _row(
            rank=2,
            predicted_pathway="Amino Acid",
            predicted_subpathway="Lysine Metabolism",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.95,
            rerank_score_calibrated=1.0,
        ),
        _row(
            rank=3,
            predicted_pathway="Amino Acid",
            predicted_subpathway="Lysine Metabolism",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.90,
            rerank_score_calibrated=0.90,
        ),
        policy="all_eligible",
    )
    assert out["compound_predicted_pathway"].iloc[0] == "Amino Acid"
    assert int(out["compound_pathway_n_voters"].iloc[0]) == 3


def test_exact_if_neighbors_agree_keeps_exact():
    out = _assign(
        _row(rank=1, predicted_pathway="Lipid", predicted_subpathway="Sterol", rerank_score_calibrated=0.5),
        _row(
            rank=2,
            predicted_pathway="Lipid",
            predicted_subpathway="Fatty Acid Metabolism",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.8,
            rerank_score_calibrated=0.4,
        ),
        policy="exact_if_neighbors_agree",
    )
    assert out["compound_predicted_pathway"].iloc[0] == "Lipid"
    assert int(out["compound_pathway_n_voters"].iloc[0]) == 1
    assert out["compound_pathway_source"].iloc[0] == "exact_mass_consistent"


def test_amino_acid_peptide_collapse_tallies_together():
    out = _assign(
        _row(
            rank=1,
            predicted_pathway="Amino Acid",
            predicted_subpathway="Lysine Metabolism",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.70,
            rerank_score_calibrated=1.0,
        ),
        _row(
            rank=2,
            predicted_pathway="Peptide",
            predicted_subpathway="Dipeptide",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.70,
            rerank_score_calibrated=0.90,
        ),
        _row(
            rank=3,
            predicted_pathway="Lipid",
            predicted_subpathway="Fatty Acid Metabolism",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.90,
            rerank_score_calibrated=1.0,
        ),
        policy="all_eligible",
        collapse_amino_acid_peptide=True,
    )
    assert out["compound_predicted_pathway"].iloc[0] in {"Amino Acid", "Peptide"}
    assert out["compound_predicted_pathway_raw"].iloc[0] == "Lipid"


def test_name_boost_can_flip_subpathway():
    priors = {("Lipid", "cholate"): {"subpathway": "Primary Bile Acid Metabolism", "precision": 1.0, "count": 10}}
    out = _assign(
        _row(
            rank=1,
            predicted_pathway="Lipid",
            predicted_subpathway="Fatty Acid Metabolism",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.51,
            rerank_score_calibrated=1.0,
            candidate_name="cholate sulfate",
        ),
        _row(
            rank=2,
            predicted_pathway="Lipid",
            predicted_subpathway="Primary Bile Acid Metabolism",
            pathway_assignment_status="strong_similarity_vote",
            pathway_vote_share=0.49,
            rerank_score_calibrated=1.0,
            candidate_name="cholate sulfate",
        ),
        policy="all_eligible",
        apply_name_boost=True,
        name_priors=priors,
        name_prior_holdout_frac=0.0,
    )
    assert out["compound_predicted_subpathway"].iloc[0] == "Primary Bile Acid Metabolism"
    assert int(out["compound_pathway_name_boosted"].iloc[0]) == 1


def test_resolve_pathway_display_prefers_compound_over_rank1():
    pytest.importorskip("rdkit")
    from structural_elucidation.reports.compound import resolve_pathway_display

    info = resolve_pathway_display(
        {},
        {
            "predicted_pathway": "Amino Acid",
            "predicted_subpathway": "Lysine Metabolism",
            "compound_predicted_pathway": "Lipid",
            "compound_predicted_subpathway": "Fatty Acid Metabolism",
            "compound_pathway_source": "weighted_candidate_consensus",
            "compound_pathway_status": "weighted_candidate_consensus",
            "compound_pathway_vote_share": 0.72,
            "compound_pathway_n_voters": 3,
            "pathway_confidence_tier": "medium",
        },
    )
    assert info["super"] == "Lipid"
    assert info["sub"] == "Fatty Acid Metabolism"
    assert info["status"] == "weighted_candidate_consensus"
    assert "Compound-level" in info["note"]


def test_resolve_pathway_display_keeps_catalog_exact():
    pytest.importorskip("rdkit")
    from structural_elucidation.reports.compound import resolve_pathway_display

    info = resolve_pathway_display(
        {"metabolon_super": "Peptide", "metabolon_sub": "Dipeptide"},
        {
            "predicted_pathway": "Amino Acid",
            "compound_predicted_pathway": "Lipid",
        },
    )
    assert info["super"] == "Peptide"
    assert info["source"] == "exact"
