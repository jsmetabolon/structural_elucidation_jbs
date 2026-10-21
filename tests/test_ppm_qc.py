"""Regenerate-and-compare tests for the ppm QC findings note."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load_tool(name: str):
    path = ROOT / "tools" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


Q = load_tool("audit_ppm_qc")
JSON_PATH = ROOT / "docs/paper/data/thresholds/audit_ppm_qc.json"
COMMITTED = json.loads(JSON_PATH.read_text())


def test_committed_counts_match_the_verification_list():
    slices = COMMITTED["slices"]
    assert slices["rank1_all"]["gt10"] == 721
    assert slices["all_top10"]["gt10"] == 13424
    assert slices["rank1_msi_1_3"]["gt10"] == 347
    assert slices["matched_truth"]["gt5_le10"] == 0
    pub = COMMITTED["truth_control_recovery"]["published"]
    assert (pub["n"], pub["rank1"], pub["top3"], pub["top10"]) == (
        92,
        44,
        59,
        68,
    )
    rec = COMMITTED["truth_control_recovery"]
    lost_10 = [x["x_compound_id"] for x in rec["match_le10_ppm"]["lost"]]
    lost_5 = [x["x_compound_id"] for x in rec["match_le5_ppm"]["lost"]]
    assert lost_10 == lost_5
    assert lost_10 == [
        "X-12811",
        "X-18345",
        "X-21470",
        "X-23080",
        "X-24137",
        "X-26094",
    ]


def test_regenerating_from_the_csvs_reproduces_the_committed_json():
    by_compound = Q.load_csv(Q.DEFAULT_BY_COMPOUND)
    top10 = Q.load_csv(Q.DEFAULT_TOP10)
    fresh = Q.audit(by_compound, top10)
    assert fresh == COMMITTED


def test_the_three_msi1_exceptions_are_named():
    items = COMMITTED["msi1_rank1_above_5_ppm"]
    ids = {item["x_compound_id"] for item in items}
    assert ids == {"X-18345", "X-21470", "X-24137"}
    note_path = ROOT / "docs/paper/data/thresholds/audit_ppm_qc.md"
    note = note_path.read_text()
    for xid in ids:
        assert f"`{xid}`" in note
    assert "10 ppm reporting ceiling was not applied" in note
    assert "no matches exceed 10 ppm" in note


def test_summary_csv_carries_the_headline_counts():
    by_compound = Q.load_csv(Q.DEFAULT_BY_COMPOUND)
    top10 = Q.load_csv(Q.DEFAULT_TOP10)
    rows = Q.csv_rows(Q.audit(by_compound, top10))
    by_key = {(r["section"], r["item"], r["metric"]): r["value"] for r in rows}
    assert by_key[("ppm_band", "rank1_all", "gt10")] == "721"
    assert by_key[("ppm_band", "all_top10", "gt10")] == "13424"
    assert by_key[("ppm_band", "rank1_msi_1_3", "gt10")] == "347"
    assert by_key[("ppm_band", "matched_truth", "gt5_le10")] == "0"
    assert by_key[("recovery", "published", "rank1")] == "44"
    committed = Q.load_csv(Q.DEFAULT_CSV)
    assert committed == rows


def test_filtered_candidates_are_within_10_ppm_and_keep_match_status():
    top10 = Q.load_csv(Q.DEFAULT_TOP10)
    rows = Q.filtered_candidate_rows(top10)
    assert len(rows) == 18211
    ppms = [float(r["abs_mass_error_ppm"]) for r in rows]
    assert max(ppms) <= 10
    assert min(ppms) >= 0
    match_ids = {
        r["x_compound_id"] for r in rows if r["matched_truth"] == "TRUE"
    }
    assert "X-18345" not in match_ids
    assert "X-21470" not in match_ids
    assert "X-24137" not in match_ids
    matches = [r for r in rows if r["matched_truth"] == "TRUE"]
    unknown = [r for r in rows if r["known_available"] == "FALSE"]
    known_miss = [r for r in rows if r["matched_truth"] == "FALSE"]
    assert len(matches) == 223
    assert unknown
    assert known_miss
    assert all(r["known_available"] == "TRUE" for r in matches)
    assert all(r["correct_result"] for r in matches)
    first = next(r for r in rows if r["x_compound_id"] == "X-01911")
    assert first["original_rank"] == "1"
    assert first["rank_after_10ppm_filter"] == "1"
    assert first["meets_qc_5ppm"] == "TRUE"
    assert first["matched_truth"] == ""
    assert first["correct_result"] == ""
    hit = next(
        r
        for r in rows
        if r["x_compound_id"] == "X-12127" and r["matched_truth"] == "TRUE"
    )
    assert "ddhc" in hit["correct_result"].lower()
    committed = Q.load_csv(Q.DEFAULT_FILTERED_CANDIDATES)
    assert committed == rows


def test_filtered_by_compound_records_surviving_known_matches():
    by_compound = Q.load_csv(Q.DEFAULT_BY_COMPOUND)
    top10 = Q.load_csv(Q.DEFAULT_TOP10)
    rows = Q.filtered_compound_rows(by_compound, top10)
    assert len(rows) == 3897
    by_id = {r["x_compound_id"]: r for r in rows}
    lost = by_id["X-18345"]
    assert lost["known_available"] == "TRUE"
    assert lost["correct_result"]
    assert lost["matched_truth"] == "FALSE"
    assert lost["match_survives_10ppm"] == "FALSE"
    assert lost["original_match_rank"] == "1"
    kept = by_id["X-12127"]
    assert "ddhc" in kept["correct_result"].lower()
    assert kept["matched_truth"] == "TRUE"
    assert kept["match_survives_10ppm"] == "TRUE"
    assert kept["match_survives_5ppm"] == "TRUE"
    assert kept["match_rank_after_10ppm_filter"] == "1"
    committed = Q.load_csv(Q.DEFAULT_FILTERED_BY_COMPOUND)
    assert committed == rows
