"""Unit tests for the selective-prediction and taxonomy-harmonization analyses."""
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


SEL = load_tool("analyze_selective_prediction")
HARM = load_tool("analyze_taxonomy_harmonization")


def structure_rows(order: list[tuple[float, str]]) -> list[dict[str, str]]:
    return [
        {"loo_probability_full": str(score), "top1_correct_0p70": correct}
        for score, correct in order
    ]


def test_curve_prefixes_grow_with_coverage():
    rows = structure_rows([(0.9, "1"), (0.8, "1"), (0.7, "0"), (0.6, "0"), (0.5, "0")])
    points = SEL.curve(rows, "loo_probability_full", "top1_correct_0p70")
    assert [point["coverage"] for point in points] == list(SEL.COVERAGE_GRID)
    sizes = [point["n"] for point in points]
    assert sizes == sorted(sizes)
    assert sizes[-1] == len(rows)


def test_curve_ranks_by_confidence_not_row_order():
    # The two correct rows are last in file order but most confident, so a curve that
    # respected the ranking finds them in the first prefix.
    rows = structure_rows([(0.1, "0"), (0.2, "0"), (0.3, "0"), (0.9, "1"), (0.95, "1")])
    points = SEL.curve(rows, "loo_probability_full", "top1_correct_0p70")
    assert points[0]["accuracy"] == 1.0
    assert points[-1]["accuracy"] == 0.4


def test_curve_skips_rows_without_a_ranking_value():
    rows = structure_rows([(0.9, "1"), (0.5, "0")])
    rows.append({"loo_probability_full": "", "top1_correct_0p70": "1"})
    assert SEL.curve(rows, "loo_probability_full", "top1_correct_0p70")[-1]["n"] == 2


def test_discrimination_is_zero_for_a_flat_curve():
    rows = structure_rows([(0.9, "1"), (0.7, "0"), (0.5, "1"), (0.3, "0")])
    flat = [{"accuracy": 0.5}, {"accuracy": 0.5}]
    assert SEL.discrimination(flat) == 0.0
    assert SEL.discrimination([]) is None
    assert SEL.discrimination(SEL.curve(rows, "loo_probability_full", "top1_correct_0p70")) is not None


def test_wilson_interval_brackets_the_observed_fraction():
    result = SEL.wilson(44, 92)
    assert result["accuracy"] == 0.4783
    assert result["ci95"][0] < result["accuracy"] < result["ci95"][1]
    assert SEL.wilson(0, 0)["accuracy"] is None


def pathway_rows(pairs: list[tuple[str, str]]) -> list[dict[str, str]]:
    return [
        {
            "true_super": true,
            "predicted_pathway": predicted,
            "true_sub": "Some Subpathway",
            "correct_sub": "false",
            "source": "neighbor_consensus",
            "in_reference_db": "false",
        }
        for true, predicted in pairs
    ]


def test_equivalence_group_merges_only_the_declared_classes():
    rows = pathway_rows([("Amino Acid", "Peptide"), ("Lipid", "Peptide")])
    raw = HARM.super_agreement(rows, ())
    harmonized = HARM.super_agreement(rows, HARM.SUPER_EQUIVALENCE)
    assert raw["correct"] == 0
    # The amino acid/peptide pair now agrees; lipid versus peptide still does not.
    assert harmonized["correct"] == 1


def test_canonical_leaves_unlisted_labels_untouched():
    assert HARM.canonical("lipid", HARM.SUPER_EQUIVALENCE) == "lipid"
    assert HARM.canonical("peptide", HARM.SUPER_EQUIVALENCE) == HARM.canonical(
        "amino acid", HARM.SUPER_EQUIVALENCE
    )


def test_unscoreable_subpathway_labels_are_excluded_only_when_requested():
    rows = [
        {"true_sub": "Ask Annie/Matt", "correct_sub": "false"},
        {"true_sub": "Ask DTS", "correct_sub": "false"},
        {"true_sub": "Unknown Steroid", "correct_sub": "false"},
        {"true_sub": "Tryptophan Metabolism", "correct_sub": "true"},
    ]
    assert HARM.sub_agreement(rows, drop_unscoreable=False)["n"] == 4
    scoreable = HARM.sub_agreement(rows, drop_unscoreable=True)
    assert scoreable["n"] == 1
    assert scoreable["accuracy"] == 1.0


def test_unscoreable_pattern_does_not_match_real_taxonomy_terms():
    for label in ("Tryptophan Metabolism", "Fatty Acid Metabolism", "Askorbate-like Pathway"):
        assert HARM.UNSCOREABLE_SUB.match(label) is None
    for label in ("Ask Annie/Matt", "ask dts", "Unknown", "Unknown Bile Acid"):
        assert HARM.UNSCOREABLE_SUB.match(label) is not None


def test_share_of_ceiling_is_the_ratio_of_observed_to_achievable():
    observed = {"accuracy": 0.5}
    ceiling = {"accuracy": 0.625}
    assert HARM.share_of_ceiling(observed, ceiling) == 0.8
    assert HARM.share_of_ceiling({"accuracy": None}, ceiling) is None
    assert HARM.share_of_ceiling(observed, {"accuracy": 0.0}) is None


def test_committed_analyses_match_a_fresh_computation():
    """The committed JSON must still follow from the committed per-control data."""

    import csv
    import json

    with (ROOT / "docs/paper/data/pathway/audit_pathway_accuracy.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    inferred = [row for row in rows if row.get("source") == "neighbor_consensus"]
    committed = json.loads((ROOT / "docs/paper/data/taxonomy_harmonization.json").read_text())

    assert HARM.super_agreement(inferred, ()) == committed["super_pathway"]["raw"]["inferred_only"]
    assert (
        HARM.sub_agreement(inferred, drop_unscoreable=True)
        == committed["subpathway"]["scoreable_only"]["inferred_only"]
    )
