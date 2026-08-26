"""Unit tests for the mummichog X-compound pathway benchmark tooling."""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent

METABOLON_SUPERS = {
    "Amino Acid",
    "Carbohydrate",
    "Cofactors and Vitamins",
    "Energy",
    "Lipid",
    "Nucleotide",
    "Peptide",
    "Xenobiotics",
    "Partially Characterized Molecules",
}


def load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


B = load_tool("mummichog_pathway_benchmark")


def test_crosswalk_loads_unique_rows_with_metabolon_vocabulary():
    mapping = B.load_crosswalk()
    assert len(mapping) >= 100
    unknown_supers = set(mapping.values()) - METABOLON_SUPERS
    assert not unknown_supers, f"crosswalk maps to non-Metabolon supers: {unknown_supers}"


def test_crosswalk_covers_every_model_pathway():
    pytest.importorskip("mummichog")
    from mummichog.JSON_metabolicModels import metabolicModels

    model = metabolicModels["human_model_mfn"]
    mapping = B.load_crosswalk()
    model_names = {p["name"] for p in model["metabolic_pathways"]}
    membership_names = {name for names in model["cpd2pathways"].values() for name in names}
    assert not model_names - set(mapping), "model pathways missing from crosswalk"
    assert not membership_names - set(mapping), "cpd2pathways names missing from crosswalk"


def test_vote_super_splits_weight_across_supers_and_reports_share():
    cpd2pathways = {
        "cpd_a": ["Tryptophan metabolism"],
        "cpd_b": ["Tryptophan metabolism", "Bile acid biosynthesis"],
        "cpd_c": ["Bile acid biosynthesis"],
    }
    crosswalk = {"Tryptophan metabolism": "Amino Acid", "Bile acid biosynthesis": "Lipid"}
    vote = B.vote_super(["cpd_a", "cpd_b", "cpd_c"], cpd2pathways, crosswalk)
    # cpd_a: Amino Acid 1.0; cpd_b: 0.5/0.5; cpd_c: Lipid 1.0 -> tie 1.5 vs 1.5.
    assert vote["tie"] is True
    assert vote["predicted_super"] == "Amino Acid"  # deterministic alphabetical tie-break
    assert vote["vote_share"] == 0.5
    assert vote["n_voting"] == 3


def test_vote_super_ignores_candidates_without_mapped_pathways():
    vote = B.vote_super(["orphan"], {"orphan": []}, {})
    assert vote["predicted_super"] == ""
    assert vote["n_voting"] == 0
    assert math.isnan(vote["vote_share"])


def test_match_mz_respects_ppm_window_and_bin_edges():
    # One theoretical ion at 100.0000; 10 ppm window is +/- 0.001.
    index = [[] for _ in range(2002)]
    index[100].append(("cpd_x", "M-H[-]", 100.0000))
    assert B.match_mz(100.0005, index, ppm=10.0) == [("cpd_x", "M-H[-]", 100.0000)]
    assert B.match_mz(100.0020, index, ppm=10.0) == []
    # Feature whose int() floor is the neighboring bin still finds the ion.
    index[1000].append(("cpd_y", "M-H[-]", 1000.0005))
    assert B.match_mz(999.9999, index, ppm=10.0) == [("cpd_y", "M-H[-]", 1000.0005)]


def test_mcnemar_exact_matches_binomial_reference():
    assert B.mcnemar_exact(0, 0) == 1.0
    assert B.mcnemar_exact(3, 3) == 1.0
    # b=1, c=10: two-sided exact = 2 * P(X <= 1 | n=11, p=0.5) = 24/2048.
    assert B.mcnemar_exact(1, 10) == pytest.approx(24 / 2048)
    assert B.mcnemar_exact(10, 1) == B.mcnemar_exact(1, 10)


def test_evaluable_truth_drops_sentinel_labels():
    import pandas as pd

    truth = pd.DataFrame(
        [
            {"control_id": "X-1", "true_super": "Lipid", "true_sub": ""},
            {"control_id": "X-2", "true_super": "Unknown", "true_sub": ""},
            {"control_id": "X-3", "true_super": "DO NOT INCLUDE", "true_sub": ""},
            {"control_id": "X-4", "true_super": "", "true_sub": ""},
        ]
    )
    kept = B.evaluable_truth(truth)
    assert list(kept["control_id"]) == ["X-1"]
