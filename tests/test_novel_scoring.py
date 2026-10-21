"""Deterministic score/report parity without Oracle, network, or model assets."""
from __future__ import annotations

import json
import sqlite3
import sys
from pathlib import Path

import pandas as pd
import pytest

from structural_elucidation.reports import compound, novel_scoring as scoring
from structural_elucidation.reports import confidence as confmod
from structural_elucidation.utils.denovo import inchikey_from_smiles
from tools import reported_identification as publication

CONTROL = "X-900000001"


def candidate(**overrides):
    return {"control_id": CONTROL, "rank": 1, "predicted_smiles": "CCO",
            "predicted_formula": "C2H6O", "target_formula": "C2H6O",
            "mass_match_score": 1.0, "abs_mass_error_da": 0.0,
            "evidence_confidence": 0.8, "unified_final_score": 0.5,
            "predicted_spectrum_score_v2": 0.6, "predicted_pathway": "Lipid",
            "pathway_best_reference_similarity": 0.8,
            "pathway_confidence_tier": "high", "pure_nce_score": 0.0,
            "nce_status": "zero", **overrides}


def metadata(**overrides):
    return {"control_id": CONTROL, "feature_id": 1, "parent_mz": 47.049141,
            "neutral_mass": 46.0418648, "lib_id": 400, "ion_mode": "positive",
            "feature_key": "feature_X_400_mz47.0491_ri1000.0_X-900000001",
            "ri": 1000.0, "peaks": [], **overrides}


@pytest.fixture
def reference_db(tmp_path):
    path = tmp_path / "library.sqlite"
    with sqlite3.connect(path) as connection:
        connection.execute("CREATE TABLE metabolites (chemical_id TEXT, canonical_name TEXT, smiles TEXT, inchikeys TEXT, in_global_library TEXT, super_pathway TEXT, pathway TEXT, metabolon_super TEXT, metabolon_sub TEXT, resolution_bucket TEXT, resolution_status TEXT, resolution_method TEXT, resolution_sources TEXT, resolution_query TEXT, resolution_match_name TEXT, resolution_score REAL, resolution_candidate_count INTEGER)")
        connection.execute("CREATE TABLE synonym_index (chemical_id TEXT, entry_type TEXT, entry_value TEXT, normalized_value TEXT, sources TEXT)")
    return path


def test_no_predictions_do_not_fabricate_scores():
    result = scoring.score_novel_feature([], metadata())
    assert result["status"] == "no_predictions"
    assert all(value is None for value in result["summary"].values())


def test_novel_handoff_cannot_inherit_standard_confirmation(reference_db):
    result = scoring.score_novel_feature([candidate()], metadata(
        compound_name="ethanol", prior_matched_name="ethanol", prior_match_level="MS2",
        metabolon_level="1", truth_match=True, assurance="standard_confirmed",
    ), db_path=reference_db)
    assert result["fields"]["assurance"] == "none"
    assert result["summary"]["se_msi_level"] != 1
    assert result["summary"]["se_reported_level"] != "1"
    assert result["fields"]["observed_neutral_mass"] == 46.0418648
    assert 0 <= result["summary"]["se_confidence"] <= 1
    assert 0 <= result["summary"]["se_evidence_score"] <= 1
    assert 0 <= result["summary"]["se_candidate_identity_confidence"] <= 100
    json.dumps(result, allow_nan=False)


def test_successful_zero_nce_is_distinct_from_failed_nce(reference_db):
    zero = scoring.score_novel_feature([candidate()], metadata(), db_path=reference_db)
    failed = scoring.score_novel_feature([candidate(nce_status="error")], metadata(), db_path=reference_db)
    assert zero["candidates"][0]["pure_nce_score"] == 0.0
    assert failed["candidates"][0]["pure_nce_score"] is None
    assert zero["summary"]["se_candidate_identity_confidence"] < failed["summary"]["se_candidate_identity_confidence"]
    assert failed["errors"][0]["stage"] == "pure_nce"
    assert "Pure NCE unavailable" in failed["summary"]["se_confidence_flags_json"]


def test_headless_runs_pure_nce_and_records_missing_evidence(monkeypatch, reference_db):
    from tools import rescore_pure_nce
    monkeypatch.setattr(rescore_pure_nce, "score_key", lambda key: (0.0, "no_observed_bins"))
    result = scoring.score_novel_feature([candidate(pure_nce_score=None, nce_status="")], metadata(), db_path=reference_db)
    assert result["candidates"][0]["pure_nce_score"] is None
    assert "no_observed_bins" in result["summary"]["se_confidence_flags_json"]
    assert not result["errors"]


def test_calibration_fallback_is_reported(monkeypatch, tmp_path, reference_db):
    path = tmp_path / "absent.json"
    monkeypatch.setattr(publication, "CALIBRATION_PATH", path)
    result = scoring.score_novel_feature([candidate()], metadata(), db_path=reference_db)
    assert result["calibration_checksum"] is None
    assert "using existing confidence bands" in result["summary"]["se_confidence_flags_json"]


def test_publication_level_2_star_stays_text(reference_db):
    result = scoring.score_novel_feature([candidate()], metadata(known_x_cleaned_external_info="glucuronide of C20H26O4"), db_path=reference_db)
    assert result["summary"]["se_reported_level"] == "2*"
    assert result["summary"]["se_candidate_identity_confidence"] >= 99


def test_identity_gate_reselects_candidate_preserving_model_rank(reference_db):
    with sqlite3.connect(reference_db) as connection:
        connection.execute("INSERT INTO metabolites (chemical_id, canonical_name, smiles, inchikeys, in_global_library) VALUES ('2','ethanol','CCO',?,'true')", (inchikey_from_smiles("CCO"),))
    result = scoring.score_novel_feature([candidate(), candidate(rank=2, predicted_smiles="COC")], metadata(), db_path=reference_db)
    assert result["naming_decision"]["identity_candidate_rank"] == 2
    assert result["candidates"][0]["model_rank"] == 2
    assert result["fields"]["smiles"] == "COC"
    assert result["raw_candidates"][0]["predicted_smiles"] == "CCO"


@pytest.mark.parametrize("level,changes", [
    (1, {"assurance": "standard_confirmed"}),
    (2, {"structure_similarity": 0.8, "ms2_spectrum_present": True, "ms2_score": 0.5, "mass_match_score": .98}),
    (3, {"structure_similarity": .6, "ms2_spectrum_present": False, "ms2_score": None, "mass_match_score": .95}),
    (4, {"structure_similarity": .6, "mass_match_score": 0., "abs_mass_error_da": 176., "observed_mz_available": True}),
    (5, {"structure_similarity": .2, "mass_match_score": None, "abs_mass_error_da": None, "observed_mz_available": False}),
])
def test_machine_summary_preserves_all_msi_levels(level, changes):
    fields = {"source": "predicted", "observed_mz_available": True,
              "pathway_tier": "high", "metabolon_level": "", **changes}
    result = scoring.summarize_scores(CONTROL, fields, [candidate()])
    assert result["se_msi_level"] == level
    assert result["se_msi_label"] == confmod.MSI_LABELS[level]
    assert result["se_confidence"] == confmod.compute_confidence(fields).overall


def report_args(monkeypatch, tmp_path, predictions, database, scoring_path):
    known = tmp_path / "known.csv"
    pd.DataFrame([{"Unknown": CONTROL, "Level": "", "KNOWN": ""}]).to_csv(known, index=False)
    captured = []
    original = compound.build_report
    monkeypatch.setattr(compound, "build_report", lambda args: captured.append(args))
    monkeypatch.setattr(sys, "argv", ["compound", "--control-id", CONTROL, "--predictions", str(predictions), "--db", str(database), "--known-x", str(known), "--known-x-all", str(tmp_path / "absent.csv"), "--curation-ledger", str(tmp_path / "absent-ledger.csv"), "--out-dir", str(tmp_path / "report"), "--scoring-json", str(scoring_path)])
    compound.main()
    monkeypatch.setattr(compound, "build_report", original)
    return captured[0]


@pytest.mark.parametrize("collision", [False, True])
def test_completed_scores_match_actual_rendered_report(monkeypatch, tmp_path, reference_db, collision):
    if collision:
        with sqlite3.connect(reference_db) as connection:
            connection.execute("INSERT INTO metabolites (chemical_id, canonical_name, smiles, inchikeys, in_global_library) VALUES ('2','ethanol','CCO',?,'true')", (inchikey_from_smiles("CCO"),))
    result = scoring.score_novel_feature([candidate()], metadata(), db_path=reference_db)
    path = tmp_path / "scores.json"
    path.write_text(json.dumps(result))
    predictions = tmp_path / "predictions.csv"
    pd.DataFrame([candidate()]).to_csv(predictions, index=False)
    args = report_args(monkeypatch, tmp_path, predictions, reference_db, path)
    manifest = compound.build_report(args)
    slug = compound.control_slug(CONTROL)
    saved = json.loads((Path(manifest["out_dir"]) / "data" / f"{slug}_confidence_inputs.json").read_text())
    assert saved == result["fields"]
    page = Path(manifest["html"]).read_text()
    assert str(result["summary"]["se_msi_level"]) in page
    assert result["summary"]["se_reported_identification"] in page
    assert "Candidate identity confidence" in page
    assert f'{round(result["summary"]["se_confidence"] * 100)}%' in page
    assert json.loads((Path(manifest["out_dir"]) / "data" / f"{slug}_scoring.json").read_text())["summary"] == result["summary"]


def test_precomputed_score_identity_must_match(tmp_path):
    path = tmp_path / "scores.json"
    path.write_text(json.dumps({"schema_version": 1, "control_id": "X-other", "status": "scored", "fields": {"msi_level": 2}}))
    with pytest.raises(ValueError, match="control_id"):
        scoring.read_scoring_payload(path, CONTROL)


def test_publication_accepts_ranked_csv_structure_and_formula_columns():
    fields = {"source": "predicted", "observed_mz_available": True, "mass_match_score": 1.0}
    result = scoring.summarize_scores(CONTROL, fields, [candidate()])
    assert result["se_reported_identification"] == "Level 4: C2H6O"
    assert result["se_candidate_identity_confidence"] is not None


@pytest.mark.parametrize("n_peaks,score,available", [(3, 0.0, True), (3, None, False), (0, 0.7, False)])
def test_observed_ms2_zero_remains_available(n_peaks, score, available):
    fields = scoring.prepare_confidence_fields(
        {"predicted_spectrum_score_v2": score, "parent_mz": 47.0}, [], {},
        pathway_vote_share=None, neighbors=[], n_peaks=n_peaks,
        target_smiles="CCO", target_exact_mass=46.0,
    )
    assert fields["ms2_spectrum_present"] is available
    assert fields["ms2_score"] == score


def test_unresolved_identity_does_not_score_rejected_structure(reference_db):
    with sqlite3.connect(reference_db) as connection:
        connection.execute("INSERT INTO metabolites (chemical_id, canonical_name, smiles, inchikeys, in_global_library) VALUES ('2','ethanol','CCO',?,'true')", (inchikey_from_smiles("CCO"),))
    result = scoring.score_novel_feature([candidate()], metadata(peaks=[[30.,1.]]), db_path=reference_db)
    assert result["naming_decision"]["identity_candidate"] is None
    assert result["fields"]["candidate_retention_proxy"] is None
    assert result["fields"]["ms2_score"] is None
    assert result["fields"]["mass_match_score"] is None
    assert result["summary"]["se_candidate_identity_confidence"] is None
    assert result["raw_candidates"][0]["predicted_smiles"] == "CCO"


def test_explicit_curation_can_confirm_a_novel_feature(tmp_path, reference_db):
    ledger = tmp_path / "ledger.csv"
    pd.DataFrame([{"control_id": CONTROL, "verdict": "standard_confirmed", "standard_run_id": "fixture-standard"}]).to_csv(ledger, index=False)
    result = scoring.score_novel_feature([candidate()], metadata(), db_path=reference_db, curation_ledger=ledger)
    assert result["summary"]["se_msi_level"] == 1
    assert result["summary"]["se_reported_level"] == "1"
    assert result["summary"]["se_candidate_identity_confidence"] == 100
    assert result["summary"]["se_confidence"] == 1
