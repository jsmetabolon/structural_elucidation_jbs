"""Unit tests for the L1-L5 identification-level figure builder."""
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


M = load_tool("build_msi_level_figure")


def test_a_specific_compound_name_counts_as_named():
    for name in (
        "glycocholate",
        "androsterone glucuronide",
        "2-octenoylcarnitine (C8)",
        "pyr-leu-phe",
        "LPI(20:2), 2-form",
        "sulfated resorcinol acetate OR catechol acetate",
        "Putative analog of succinate (Tanimoto 0.50)",
    ):
        assert M.is_named(name), name


def test_a_placeholder_does_not_count_as_named():
    for name in (
        "",
        None,
        "Tripeptide",
        "Small Peptide",
        "Dipeptide (His/Val composition)",
        "Tripeptide with Ala, Leu, Val",
        "Sulfate conjugate of C9H10O3",
        "glucuronide of C12H20O3 (5)",
        "1,3,7,9-tetramethylpurine-2,6,8-trione \u2014 IUPAC (generated) (C9H12N4O3)",
        "Amino Acid Derivative",
        "Sterol-like compound",
        "Unknown",
        "Confidence 0 - X-21785",
    ):
        assert not M.is_named(name), name


def _row(level, name, cohort="X"):
    return {"level": str(level), "name": name, "cohort": cohort}


def test_the_generous_rule_moves_only_named_l3_rows_to_l2():
    rows = [
        _row(3, "glycocholate"),
        _row(3, "Tripeptide"),
        _row(2, "taurine", cohort="Known"),
    ]
    delivered, _ = M.tally(rows, generous=False)
    generous, _ = M.tally(rows, generous=True)
    assert sum(delivered["2"].values()) == 1
    assert sum(delivered["3"].values()) == 2
    assert sum(generous["2"].values()) == 2
    assert sum(generous["3"].values()) == 1


def test_l1_l4_and_l5_are_untouched_by_the_generous_rule():
    rows = [_row(1, "taurine", cohort="Known"), _row(4, "isobutyl phthalate"), _row(5, "")]
    delivered, _ = M.tally(rows, generous=False)
    generous, _ = M.tally(rows, generous=True)
    for level in ("1", "4", "5"):
        assert delivered[level] == generous[level]


def test_rows_without_an_l1_l5_level_are_reported_rather_than_counted():
    rows = [_row("", "x"), _row("Exclude", "x"), _row("Internal Standard", "x"), _row(2, "taurine")]
    counts, unlevelled = M.tally(rows, generous=True)
    assert sum(sum(counter.values()) for counter in counts.values()) == 1
    assert unlevelled == {"(blank)": 1, "Exclude": 1, "Internal Standard": 1}


def test_a_caption_wraps_onto_lines_no_wider_than_asked_for():
    lines = M.wrap("one two three four five six seven eight nine ten", 14)
    assert lines == ["one two three", "four five six", "seven eight", "nine ten"]
