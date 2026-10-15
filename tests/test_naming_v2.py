"""Identity gate regressions for catalog collisions and X annotations."""
from __future__ import annotations

import sqlite3

from structural_elucidation.reports.naming_v2 import CollisionIndex, decide
from structural_elucidation.utils.denovo import inchikey_from_smiles


def index(tmp_path):
    path = tmp_path / "library.sqlite"
    with sqlite3.connect(path) as connection:
        connection.execute(
            "CREATE TABLE metabolites (chemical_id TEXT, canonical_name TEXT, "
            "inchikeys TEXT, smiles TEXT, in_global_library TEXT)"
        )
        connection.executemany(
            "INSERT INTO metabolites VALUES (?, ?, ?, ?, ?)",
            [
                ("100006103", "3-hydroxyindolin-2-one sulfate",
                 inchikey_from_smiles("O=C1Nc2ccccc2C1OS(=O)(=O)O"),
                 "O=C1Nc2ccccc2C1OS(=O)(=O)O", "false | true"),
                ("2", "ethanol", inchikey_from_smiles("CCO"), "CCO", "true"),
            ],
        )
    return CollisionIndex(path)


OXINDOLE = "O=C1Nc2ccccc2C1OS(=O)(=O)O"


def candidate(rank, smiles, **extra):
    return {"rank": rank, "predicted_smiles": smiles, "mass_match_score": 1.0,
            "abs_mass_error_da": 0.0001, **extra}


def test_x12216_curated_name_and_collision(tmp_path):
    row = {"Formula": "C8H7NO5S", "Tentative information":
           "Possible 5,6-dihydroxyindole sulfate"}
    result = decide([candidate(1, OXINDOLE)], index(tmp_path), row)
    assert result["name_source"] == "known_x_tentative"
    assert result["curated_hypothesis"] == "5,6-dihydroxyindole sulfate"
    assert result["identity_candidate"] is None
    assert result["excluded_library_matches"][0]["matches"][0]["chemical_id"] == "100006103"


def test_x12729_and_x18935_share_collision_but_differ_in_name(tmp_path):
    library = index(tmp_path)
    x12729 = decide([candidate(1, OXINDOLE)], library, {
        "Tentative information": "Sulfate of C8H7NO2 Possible 5,6-dihydroxyindole sulfate",
        "Cleaned External Info": "5,6-dihydroxyindole-sulfate",
    })
    x18935 = decide([candidate(1, OXINDOLE)], library, {
        "Tentative information": "sulfate of C8H7NO2",
        "Cleaned External Info": "sulfate of C8H7NO2",
        "Formula": "C8H7NO5S",
    })
    assert x12729["curated_hypothesis"] == "5,6-dihydroxyindole-sulfate"
    assert x18935["name_source"] == "formula_or_mass"
    assert x12729["excluded_library_matches"][0]["matches"] == x18935["excluded_library_matches"][0]["matches"]


def test_selects_noncolliding_candidate_and_preserves_raw_rank(tmp_path):
    result = decide([candidate(1, "CCO"), candidate(2, "CCN")], index(tmp_path), {})
    assert result["raw_rank1_candidate"]["model_rank"] == 1
    assert result["identity_candidate_rank"] == 2
    assert result["identity_candidate"]["smiles"] == "CCN"


def test_only_authentic_standard_waives_collision(tmp_path):
    library = index(tmp_path)
    rows = [candidate(1, "CCO")]
    assert decide(rows, library, {}, assurance="curator_known")["identity_candidate"] is None
    confirmed = decide(rows, library, {}, assurance="standard_confirmed")
    assert confirmed["identity_candidate_rank"] == 1
    assert confirmed["exception_basis"] == "standard_confirmed"


def test_mismapped_library_name_is_quarantined(tmp_path):
    result = decide([candidate(1, "CCN", library_name="ethanol")], index(tmp_path), {})
    assert result["identity_candidate"] is None
    assert result["identity_exclusion_reason"] == "library_name_structure_conflict"


def test_formula_only_when_no_eligible_candidate(tmp_path):
    result = decide([candidate(1, "CCO")], index(tmp_path), {"Formula": "C2H6O"})
    assert result["name_source"] == "formula_or_mass"
    assert not result["identity_eligible"]


def test_standard_selects_verified_structure_not_raw_rank_one(tmp_path):
    result = decide(
        [candidate(1, "CCN"), candidate(2, "CCO")], index(tmp_path), {},
        assurance="standard_confirmed", confirmed_structure="CCO",
    )
    assert result["name_source"] == "authentic_standard"
    assert result["identity_candidate_rank"] == 2
    assert result["exception_basis"] == "standard_confirmed"


def test_known_library_synonym_is_quarantined(tmp_path):
    path = tmp_path / "library.sqlite"
    index(tmp_path)
    with sqlite3.connect(path) as connection:
        connection.execute(
            "CREATE TABLE synonym_index (chemical_id TEXT, entry_type TEXT, entry_value TEXT)"
        )
        connection.execute(
            "INSERT INTO synonym_index VALUES ('2', 'synonym', 'ethyl alcohol')"
        )
    library = CollisionIndex(path, supplied_names=["ethyl alcohol"])
    result = decide([candidate(1, "CCN", library_name="ethyl alcohol")], library, {})
    assert result["identity_exclusion_reason"] == "library_name_structure_conflict"


def test_curated_name_cannot_reassert_known_library_name(tmp_path):
    result = decide(
        [candidate(1, "CCO")], index(tmp_path),
        {"Tentative information": "Possible 3-hydroxyindolin-2-one sulfate"},
    )
    assert result["curated_hypothesis"] == ""
    assert result["name_source"] == "formula_or_mass"
