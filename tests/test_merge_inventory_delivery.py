"""The merged delivery keeps one source per eligible compound."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from merge_inventory_delivery import (  # noqa: E402
    ADDUCT_FALLBACK,
    EARLIER,
    NO_INPUTS,
    NO_PREDICTION,
    NO_SUPPORTED,
    STORED_MASS,
    assign_source,
    build_delivery,
    fill_missing_candidates,
)


def test_current_summary_wins_and_missing_compounds_are_listed() -> None:
    assert assign_source("X-1", current_ids={"X-1"}, settled_ids={"X-1"}, earlier_ids={"X-1"}, no_input_ids=set()) == (
        STORED_MASS,
        "",
    )
    assert assign_source("X-2", current_ids={"X-2"}, settled_ids=set(), earlier_ids={"X-2"}, no_input_ids=set()) == (
        ADDUCT_FALLBACK,
        "",
    )
    assert assign_source("X-3", current_ids=set(), settled_ids=set(), earlier_ids={"X-3"}, no_input_ids={"X-3"}) == (
        EARLIER,
        "",
    )
    assert assign_source("X-4", current_ids=set(), settled_ids=set(), earlier_ids=set(), no_input_ids={"X-4"})[1] == NO_INPUTS
    assert assign_source("X-5", current_ids=set(), settled_ids=set(), earlier_ids=set(), no_input_ids=set()) == (
        NO_PREDICTION,
        NO_SUPPORTED,
    )

    current = pd.DataFrame({"x_compound_id": ["X-1"], "reported_identification": ["Level 2: current"], "rank": ["1"]})
    earlier = pd.DataFrame(
        {"x_compound_id": ["X-1", "X-3"], "reported_identification": ["Level 4: old", "Level 4: kept"], "rank": ["1", "1"]}
    )
    known = pd.DataFrame({"Level": ["5"], "SubLevel": ["5"], "Tentative information": [""], "Cleaned External Info": [""]}, index=["X-4"])
    by_compound, top10, counts = build_delivery(
        eligible=["X-1", "X-3", "X-4"],
        current_by_compound=current.drop(columns="rank"),
        earlier_by_compound=earlier.drop(columns="rank"),
        current_top10=current,
        earlier_top10=earlier,
        settled_ids={"X-1"},
        no_input_ids={"X-4"},
        known=known,
    )
    assert list(by_compound["x_compound_id"]) == ["X-1", "X-3", "X-4"]
    assert by_compound.loc[by_compound.x_compound_id.eq("X-1"), "reported_identification"].item() == "Level 2: current"
    assert set(top10["x_compound_id"]) == {"X-1", "X-3"}
    assert by_compound.loc[by_compound.x_compound_id.eq("X-4"), "metabolon_level"].item() == "5"
    assert counts["compounds"] == 3
    assert counts[NO_PREDICTION] == 1


def test_gated_rows_fill_compounds_the_summary_exported_without_candidates() -> None:
    by_compound = pd.DataFrame({
        "x_compound_id": ["X-1", "X-2"],
        "prediction_source": ["stored_mass", "stored_mass"],
        "delivery_note": ["", ""],
        "reported_identification": ["Level 2: kept", "Level 4: filled"],
        "n_candidates": ["1", "0"],
        "rank1_candidate_smiles": ["C", ""],
    })
    top10 = pd.DataFrame({
        "x_compound_id": ["X-1"],
        "prediction_source": ["stored_mass"],
        "delivery_note": [""],
        "reported_identification": ["Level 2: kept"],
        "rank": ["1"],
        "candidate_smiles": ["C"],
        "candidate_formula": ["CH4"],
        "mass_match_score": ["1"],
        "abs_mass_error_da": ["0"],
        "abs_mass_error_ppm": ["1"],
        "evidence_confidence": ["0.5"],
        "evidence_confidence_basis": ["emitted"],
    })
    gated = pd.DataFrame({
        "control_id": ["X-1", "X-2"],
        "rank": ["9", "1"],
        "predicted_smiles": ["CCC", "CC"],
        "candidate_formula": ["C3", "C2"],
        "mass_match_score": ["0.1", "0.9"],
        "abs_mass_error_da": ["1", "0.001"],
        "abs_mass_error_ppm": ["50", "2"],
        "evidence_confidence": ["0.1", "0.4"],
        "evidence_confidence_basis": ["emitted", "emitted"],
        "pool": ["current", "current"],
        "lib": [400, 209],
    })
    filled = fill_missing_candidates(top10, by_compound, gated)
    assert list(filled["x_compound_id"]) == ["X-1", "X-2"]
    assert filled.loc[filled.x_compound_id.eq("X-1"), "candidate_smiles"].item() == "C"
    assert filled.loc[filled.x_compound_id.eq("X-2"), "candidate_smiles"].item() == "CC"
    assert filled.loc[filled.x_compound_id.eq("X-2"), "reported_identification"].item() == "Level 4: filled"
