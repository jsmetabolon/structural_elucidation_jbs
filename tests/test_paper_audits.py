"""Unit tests for publication audit statistics and cohort semantics."""
from __future__ import annotations

import importlib.util
from pathlib import Path

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


def test_pathway_axis_ablation_renormalizes_available_axes():
    axes = [
        S.C.Axis("structure", "Structure", 0.8, 0.4, True),
        S.C.Axis("mass", "Mass", 1.0, 0.3, True),
        S.C.Axis("ms2", "MS2", None, 0.2, False),
        S.C.Axis("pathway", "Pathway", 0.2, 0.1, True),
    ]
    expected = (0.8 * 0.4 + 1.0 * 0.3) / (0.4 + 0.3)
    assert abs(S.raw_without_pathway(axes) - expected) < 1e-12

