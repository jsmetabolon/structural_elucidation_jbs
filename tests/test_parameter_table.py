"""Tests for the exported evidence-score parameter table.

The table's whole claim is that it describes the model that actually runs. These
tests are what make that claim hold: they fail when a weight, a cutoff, or an
axis moves without the export being regenerated, which is the failure mode a
hand-maintained table has no way to notice.
"""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

from structural_elucidation.inference.mass_windows import (
    MASS_SCORE_PLATEAU_DA,
    MASS_SCORE_TOLERANCE_DA,
)
from structural_elucidation.reports.confidence import (
    AXIS_LABELS,
    MASS_WEAK,
    NEIGHBOR_SIMILARITY_FLOOR,
    RETENTION_IMPLAUSIBLE_Z,
    STRUCT_STRONG,
    WEIGHTS,
)

ROOT = Path(__file__).resolve().parent.parent
PARAMETERS = ROOT / "docs/paper/data/parameters"


def load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


T = load_tool("build_parameter_table")


def committed_json() -> dict:
    return json.loads((PARAMETERS / "evidence_parameters.json").read_text())


def committed_csv() -> list[dict[str, str]]:
    with (PARAMETERS / "evidence_parameters.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))


def test_committed_json_matches_a_fresh_computation():
    """A constant moving without a regeneration is exactly what this catches."""

    assert T.build_summary() == committed_json()


def test_committed_csv_matches_a_fresh_computation():
    fresh = [
        {key: str(value) for key, value in row.items()} for row in T.build_rows()
    ]
    assert fresh == committed_csv()


def test_every_scored_axis_appears_with_its_real_weight():
    rows = T.build_rows()
    scored = {row["axis_key"]: row["weight"] for row in rows if row["scored"]}
    assert scored == {key: round(weight, 2) for key, weight in WEIGHTS.items()}


def test_scored_weights_still_total_the_model_total():
    summary = T.build_summary()
    assert summary["weight_total"] == round(sum(WEIGHTS.values()), 2)
    assert summary["weight_total"] == 1.05


def test_every_axis_is_labelled_the_way_the_report_labels_it():
    """A parameter renamed in the report and not here would describe nothing."""

    rows = T.build_rows()
    for row in rows:
        if row["scored"]:
            assert row["parameter"] == AXIS_LABELS[row["axis_key"]]


def test_library_neighbors_is_carried_but_never_scored():
    rows = [row for row in T.build_rows() if row["parameter"] == "Library neighbors"]
    assert rows, "the neighbor count is documented even though it carries no weight"
    for row in rows:
        assert row["scored"] is False
        assert row["weight"] == 0.0
        assert row["axis_key"] == ""
    assert f"{NEIGHBOR_SIMILARITY_FLOOR:.2f}" in rows[0]["units"]


def test_every_band_carries_language_and_an_edge():
    for row in T.build_rows():
        assert row["report_language"].strip()
        assert row["low_end"].strip()
        assert row["high_end"].strip()
        assert row["units"].strip()


def test_mass_bands_follow_the_scoring_windows_rather_than_a_ppm_draft():
    """The drafted 0-1 / 1-3 / 3-5 ppm bands are not what the scorer implements."""

    rows = [row for row in T.build_rows() if row["axis_key"] == "mass"]
    strong, acceptable, borderline = rows
    assert strong["high_end"].startswith(f"{MASS_SCORE_PLATEAU_DA:.3f} Da")
    # The partial band ends where the linear arm crosses the weak-mass cutoff.
    assert acceptable["high_end"].startswith(
        f"{(1.0 - MASS_WEAK) * MASS_SCORE_TOLERANCE_DA:.4f} Da"
    )
    assert borderline["high_end"] == "no limit"


def test_ppm_equivalents_are_quoted_at_the_stated_reference_mass():
    summary = T.build_summary()
    reference = summary["ppm_reference_mass_da"]
    expected = MASS_SCORE_PLATEAU_DA / reference * 1e6
    strong = next(row for row in summary["parameters"] if row["axis_key"] == "mass")
    assert f"{expected:.1f} ppm" in strong["high_end"]


def test_retention_marks_the_band_that_demotes_on_its_own():
    rows = [row for row in T.build_rows() if row["axis_key"] == "retention"]
    borderline = rows[-1]
    assert borderline["low_end"] == f"{RETENTION_IMPLAUSIBLE_Z:.1f} sigma"
    assert "Tentative" in borderline["consequence"]


def test_structure_strong_band_starts_at_the_corroboration_cutoff():
    rows = [row for row in T.build_rows() if row["axis_key"] == "structure"]
    assert rows[0]["low_end"] == f"{STRUCT_STRONG:.2f}"


def test_supplement_renders_the_committed_table():
    supplement = load_tool("build_paper_supplement")
    lines = supplement.parameter_section(committed_json())
    body = "\n".join(lines)
    assert "**Table S10." in body
    for label in AXIS_LABELS.values():
        assert f"| {label} |" in body
    assert "| Library neighbors |" in body
