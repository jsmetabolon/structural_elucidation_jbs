"""Unit tests for publication audit statistics and cohort semantics."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


P = load_tool("audit_pathway_predictions")
S = load_tool("audit_structure_predictions")


def test_wilson_interval_contains_observed_fraction():
    result = P.wilson(44, 92)
    assert result["accuracy"] == 0.4783
    assert result["ci95"][0] < result["accuracy"] < result["ci95"][1]


def test_assignment_source_separates_lookup_from_neighbor_consensus():
    assert P.assignment_source("exact") == "exact_reference_lookup"
    assert P.assignment_source("predicted") == "neighbor_consensus"
    assert P.assignment_source("") == "neighbor_consensus"


def test_structure_rank_k_is_cumulative():
    records = [
        {"similarities": [0.8, 0.1, 0.1]},
        {"similarities": [0.1, 0.8, 0.1]},
        {"similarities": [0.1, 0.1, 0.8]},
    ]
    result = S.rank_k_table(records)["0.70"]
    assert result["top1"]["correct"] == 1
    assert result["top3"]["correct"] == 3
    assert result["top10"]["correct"] == 3


def test_compound_level_audit_reports_delta_vs_rank1():
    pred = pd.DataFrame(
        [
            {
                "control_id": "X-1",
                "lib": "209",
                "predicted_smiles": "CC",
                "predicted_pathway": "Amino Acid",
                "predicted_subpathway": "Lysine Metabolism",
                "compound_predicted_pathway": "Lipid",
                "compound_predicted_subpathway": "Fatty Acid Metabolism",
                "compound_pathway_source": "weighted_candidate_consensus",
                "compound_pathway_status": "weighted_candidate_consensus",
                "source": "neighbor_consensus",
                "tier": "high",
                "status": "strong_similarity_vote",
                "vote_share": 0.5,
                "best_similarity": 0.8,
            }
        ]
    )
    gt = pd.DataFrame([{"control_id": "X-1", "true_super": "Lipid", "true_sub": "Fatty Acid Metabolism"}])
    rank1, _ = P.build_summary(pred, gt)
    compound = P.compound_level_summary(pred, gt)
    assert rank1["accuracy"]["super_accuracy"] == 0.0
    assert compound["accuracy"]["super_accuracy"] == 1.0
    assert compound["accuracy"]["super_accuracy_aa_peptide"] == 1.0


def test_known_x_hygiene_treats_amino_acid_peptide_and_drops_ask_dts():
    assert P.supers_match("Amino Acid", "Peptide", harmonize_aa_peptide=True)
    assert not P.supers_match("Amino Acid", "Peptide", harmonize_aa_peptide=False)
    assert P.is_unscoreable_sub("Ask DTS")
    assert P.is_unscoreable_sub("Unknown Steroid")
    assert not P.is_unscoreable_sub("Lysine Metabolism")

    pred = pd.DataFrame(
        [
            {
                "control_id": "X-aa",
                "lib": "209",
                "predicted_smiles": "CC",
                "predicted_pathway": "Peptide",
                "predicted_subpathway": "Dipeptide",
                "source": "neighbor_consensus",
                "tier": "high",
                "status": "strong_similarity_vote",
                "vote_share": 0.5,
                "best_similarity": 0.8,
            },
            {
                "control_id": "X-ask",
                "lib": "209",
                "predicted_smiles": "CC",
                "predicted_pathway": "Lipid",
                "predicted_subpathway": "Sterol",
                "source": "neighbor_consensus",
                "tier": "high",
                "status": "strong_similarity_vote",
                "vote_share": 0.5,
                "best_similarity": 0.8,
            },
        ]
    )
    gt = pd.DataFrame(
        [
            {"control_id": "X-aa", "true_super": "Amino Acid", "true_sub": "Lysine Metabolism"},
            {"control_id": "X-ask", "true_super": "Lipid", "true_sub": "Ask DTS"},
        ]
    )
    summary, _ = P.build_summary(pred, gt)
    assert summary["accuracy"]["super_accuracy"] == 0.5
    assert summary["accuracy"]["super_accuracy_aa_peptide"] == 1.0
    assert summary["accuracy"]["n_scoreable_sub"] == 1
    assert summary["accuracy"]["sub_accuracy_scoreable"] == 0.0
    assert summary["accuracy"]["super_ci95"][0] <= 0.5 <= summary["accuracy"]["super_ci95"][1]
    assert P.share_of_ceiling(0.59, 0.625) == 0.944
    assert P.share_of_ceiling(0.59, 0) is None


def test_pathway_axis_ablation_renormalizes_available_axes():
    axes = [
        S.C.Axis("structure", "Structure", 0.8, 0.4, True),
        S.C.Axis("mass", "Mass", 1.0, 0.3, True),
        S.C.Axis("ms2", "MS2", None, 0.2, False),
        S.C.Axis("pathway", "Pathway", 0.2, 0.1, True),
    ]
    expected = (0.8 * 0.4 + 1.0 * 0.3) / (0.4 + 0.3)
    assert abs(S.raw_without_pathway(axes) - expected) < 1e-12

