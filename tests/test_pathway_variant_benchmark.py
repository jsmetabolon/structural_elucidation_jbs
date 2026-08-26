"""Unit tests for the pathway variant benchmark gate and harmonized audit scoring."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


G = load_tool("pathway_variant_benchmark")
A = load_tool("audit_pathway_predictions")


def test_selfcheck_reproduces_committed_audit():
    assert G.run_selfcheck() == 0


def test_harmonize_merges_amino_acid_and_peptide_only():
    for module in (G, A):
        assert module.harmonize("Amino Acid") == module.harmonize("Peptide")
        assert module.harmonize("amino acid") == "amino acid|peptide"
        assert module.harmonize("Lipid") == "lipid"
        assert module.harmonize("Xenobiotics") != module.harmonize("Lipid")


def test_mcnemar_exact_matches_binomial_reference():
    assert G.mcnemar_exact(0, 0) == 1.0
    assert G.mcnemar_exact(5, 5) == 1.0
    # b=1, c=10: two-sided exact = 2 * P(X <= 1 | n=11, p=0.5) = 24/2048.
    assert G.mcnemar_exact(1, 10) == pytest.approx(24 / 2048)
    assert G.mcnemar_exact(10, 1) == G.mcnemar_exact(1, 10)


def _cohort_row(cid: str, true_super: str, baseline: str, tier: str = "high", sim: str = "1.0"):
    return {
        "control_id": cid,
        "predicted_pathway": baseline,
        "true_super": true_super,
        "predicted_subpathway": "",
        "true_sub": "",
        "tier": tier,
        "status": "exact_reference_match" if sim == "1.0" else "strong_similarity_vote",
        "best_similarity": sim,
    }


def test_evaluate_variant_scores_coverage_and_head_to_head():
    cohort = [
        _cohort_row("X-1", "Lipid", "Lipid"),
        _cohort_row("X-2", "Peptide", "Amino Acid"),
        _cohort_row("X-3", "Xenobiotics", "Lipid"),
    ]
    predictions = {
        "X-1": {"predicted_super": "Lipid", "predicted_sub": ""},
        "X-2": {"predicted_super": "Peptide", "predicted_sub": ""},
        # X-3 abstains.
    }
    result = G.evaluate_variant(predictions, cohort, label="test")
    assert result["n_called"] == 2
    assert result["coverage"] == pytest.approx(2 / 3, abs=1e-4)
    assert result["super_accuracy_on_called"]["raw"]["correct"] == 2
    # Baseline got X-2 wrong raw (Amino Acid vs Peptide) but right harmonized.
    assert result["head_to_head_raw"]["variant_only_correct"] == 1
    assert result["head_to_head_harmonized"]["variant_only_correct"] == 0


def test_abstention_policy_filters():
    row_low = _cohort_row("X-1", "Lipid", "Lipid", tier="low", sim="0.4")
    row_high = _cohort_row("X-2", "Lipid", "Lipid", tier="high", sim="1.0")
    assert G._policy_keep("all_assignments", row_low)
    assert not G._policy_keep("drop_low_tier", row_low)
    assert G._policy_keep("drop_low_tier", row_high)
    assert not G._policy_keep("high_tier_only", row_low)
    assert not G._policy_keep("similarity_ge_0.5", row_low)
    assert G._policy_keep("similarity_ge_0.5", row_high)
    assert not G._policy_keep("exact_status_only", row_low)


def test_audit_build_summary_reports_harmonized_metrics():
    import pandas as pd

    pred = pd.DataFrame(
        [
            {
                "control_id": "X-1",
                "predicted_smiles": "",
                "predicted_pathway": "Amino Acid",
                "predicted_subpathway": "",
                "source": "neighbor_consensus",
                "tier": "high",
                "status": "strong_similarity_vote",
                "vote_share": 0.9,
                "best_similarity": 0.8,
            }
        ]
    )
    gt = pd.DataFrame([{"control_id": "X-1", "true_super": "Peptide", "true_sub": ""}])
    summary, ev = A.build_summary(pred, gt)
    assert summary["accuracy"]["super_accuracy"] == 0.0
    assert summary["accuracy"]["super_accuracy_harmonized"] == 1.0
    assert summary["accuracy"]["harmonized_overall_super"]["n"] == 1
    assert bool(ev["correct_super_harmonized"].iloc[0]) is True
