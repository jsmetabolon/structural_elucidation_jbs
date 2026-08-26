"""Unit tests for pathway evidence fusion (T4)."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


F = load_tool("pathway_evidence_fusion")


def test_confident_peptide_detector_outranks_everything():
    fused = F.fuse_one(
        "Lipid",
        "Sterol",
        "high",
        1.0,
        peptide_label="gly-gly dipeptide",
        peptide_confidence=0.9,
        lipid_class="PC",
        mass_super="Xenobiotics",
        mass_share=1.0,
        mass_n_candidates=10,
    )
    assert fused["predicted_super"] == "Peptide"
    assert fused["fusion_source"] == "peptide_detector"
    assert fused["predicted_sub"] == ""  # baseline sub belonged to Lipid, dropped


def test_low_confidence_peptide_detector_does_not_fire():
    fused = F.fuse_one(
        "Amino Acid",
        "",
        "high",
        0.9,
        peptide_label="weak fit",
        peptide_confidence=0.2,
    )
    assert fused["fusion_source"] == "structure_vote"
    assert fused["predicted_super"] == "Amino Acid"


def test_lipid_detector_fires_after_peptide_and_keeps_concordant_sub():
    fused = F.fuse_one(
        "Lipid",
        "Phosphatidylcholine",
        "low",
        0.3,
        lipid_class="PC",
        mass_super="Xenobiotics",
        mass_share=1.0,
        mass_n_candidates=10,
    )
    assert fused["fusion_source"] == "lipid_detector"
    assert fused["predicted_super"] == "Lipid"
    assert fused["predicted_sub"] == "Phosphatidylcholine"  # concordant, kept


def test_mass_prior_backfills_weak_baseline_only():
    weak = F.fuse_one(
        "Amino Acid", "", "low", 0.35, mass_super="Xenobiotics", mass_sub="Benzoate Metabolism",
        mass_share=0.8, mass_n_candidates=5,
    )
    assert weak["fusion_source"] == "mass_prior"
    assert weak["predicted_super"] == "Xenobiotics"
    assert weak["predicted_sub"] == "Benzoate Metabolism"

    strong = F.fuse_one(
        "Amino Acid", "", "high", 0.95, mass_super="Xenobiotics", mass_share=0.8, mass_n_candidates=5
    )
    assert strong["fusion_source"] == "structure_vote"
    assert strong["predicted_super"] == "Amino Acid"


def test_mass_prior_needs_share_and_candidate_thresholds():
    low_share = F.fuse_one(
        "Amino Acid", "", "low", 0.35, mass_super="Xenobiotics", mass_share=0.4, mass_n_candidates=5
    )
    assert low_share["fusion_source"] == "structure_vote"
    few_candidates = F.fuse_one(
        "Amino Acid", "", "low", 0.35, mass_super="Xenobiotics", mass_share=0.9, mass_n_candidates=1
    )
    assert few_candidates["fusion_source"] == "structure_vote"


def test_mass_prior_vote_window_and_tie_break():
    masses = [100.0000, 100.0005, 100.0008, 200.0]
    labels = [("Lipid", "A"), ("Amino Acid", "B"), ("Amino Acid", "B"), ("Xenobiotics", "C")]
    vote = F.mass_prior_vote(100.0004, masses, labels, ppm=10.0)
    assert vote["mass_n_candidates"] == 3  # 200.0 far outside the window
    assert vote["mass_super"] == "Amino Acid"
    assert vote["mass_sub"] == "B"
    assert vote["mass_share"] == pytest.approx(2 / 3, abs=1e-4)

    empty = F.mass_prior_vote(150.0, masses, labels, ppm=10.0)
    assert empty["mass_n_candidates"] == 0
    assert empty["mass_super"] == ""
