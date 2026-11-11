"""Unit tests for the x-compound-report agent skill helpers."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".cursor" / "skills" / "x-compound-report" / "scripts"))

import client_report as CR  # noqa: E402
import compare_reference as CMP  # noqa: E402
import resolve_latest as RL  # noqa: E402


def test_ids_are_canonicalized_from_messy_tokens():
    pairs = CR.ids_from_text("x 1911, X25385; X-23753 and X_00042")
    assert [cid for _, cid in pairs] == ["X-01911", "X-25385", "X-23753", "X-00042"]


def test_table_uses_named_column_and_ignores_other_columns():
    rows = [["Sample", "X Compound", "neighbor"], ["S1", "X-11677", "X-99998"], ["S2", "x 1911", ""]]
    assert [cid for _, cid in CR.ids_from_table(rows)] == ["X-11677", "X-01911"]


def test_table_without_header_picks_the_id_column():
    rows = [["a", "X-11677"], ["b", "X-12007"], ["c", "note"]]
    assert [cid for _, cid in CR.ids_from_table(rows)] == ["X-11677", "X-12007"]


def test_read_inputs_mixes_files_and_literals(tmp_path):
    listing = tmp_path / "ids.txt"
    listing.write_text("# client list\nX-11677\nnot an id\n")
    pairs, unparsed = CR.read_compound_inputs([str(listing), "X 2249"])
    assert [cid for _, cid in pairs] == ["X-11677", "X-02249"]
    assert unparsed == ["not an id"]


def test_screen_blocks_excluded_and_counts_duplicates():
    known = {
        "X-26724": [{"Unknown": "X-26724", "Level": "Exclude", "SubLevel": "Exclude", "LIB_ID": "209.0",
                     "Tentative information": "proprietary drug", "Cleaned External Info": "DO NOT INCLUDE"}],
        "X-11677": [{"Unknown": "X-11677", "Level": "1", "SubLevel": "KNOWN", "LIB_ID": "400.0",
                     "Tentative information": "", "Cleaned External Info": ""}],
    }
    rows = CR.screen([("X-11677", "X-11677"), ("X-26724", "X-26724"), ("x 11677", "X-11677"),
                      ("X-99999", "X-99999")], known, allow_excluded=False)
    by_id = {r["control_id"]: r for r in rows}
    assert by_id["X-26724"]["status"] == "excluded"
    assert "proprietary drug" in by_id["X-26724"]["exclusion_reason"]
    assert by_id["X-11677"]["duplicates"] == 1 and by_id["X-11677"]["known_x_libs"] == "400"
    assert by_id["X-99999"]["status"] == "included" and not by_id["X-99999"]["in_known_x"]
    kept = CR.screen([("X-26724", "X-26724")], known, allow_excluded=True)
    assert kept[0]["status"] == "included" and kept[0]["exclusion_reason"]


def test_shipped_known_x_blocks_client_proprietary_compounds():
    known = CR.load_known_x()
    rows = CR.screen([(c, c) for c in ("X-23358", "X-24933", "X-26724", "X-14536")], known, allow_excluded=False)
    assert {r["status"] for r in rows} == {"excluded"}


def test_package_copies_public_page_and_its_links_only(tmp_path):
    reports = tmp_path / "reports"
    bundle = reports / "X-00001"
    (bundle / "data").mkdir(parents=True)
    (bundle / "public.html").write_text('<a href="data/a.csv">a</a><a href="https://x.org">x</a><a href="#top">t</a>')
    (bundle / "internal.html").write_text("observed ions")
    (bundle / "data" / "a.csv").write_text("x\n")
    (bundle / "data" / "x-00001_scoring.json").write_text("{}")
    (reports / "index.html").write_text("<a href='X-00001/index.html'>X-00001</a>")
    result = CR.build_package(reports, tmp_path / "pkg", ["X-00001"], None)
    files = sorted(str(p.relative_to(tmp_path / "pkg")) for p in (tmp_path / "pkg").rglob("*") if p.is_file())
    assert files == ["X-00001/data/a.csv", "X-00001/index.html", "index.html"]
    assert result["problems"] == [] and result["n_reports"] == 1


def test_package_reports_missing_linked_file(tmp_path):
    bundle = tmp_path / "reports" / "X-00002"
    bundle.mkdir(parents=True)
    (bundle / "public.html").write_text('<a href="data/missing.csv">m</a>')
    result = CR.build_package(tmp_path / "reports", tmp_path / "pkg", ["X-00002"], None)
    assert result["problems"] == ["X-00002: public.html links missing file data/missing.csv"]


def _write_predictions(root: Path, name: str, rows: list[tuple[str, int, str]]) -> None:
    lib = root / "rerank_lib209"
    lib.mkdir(parents=True, exist_ok=True)
    with (lib / name).open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["control_id", "rank", "predicted_smiles"])
        writer.writerows(rows)


def test_resolver_prefers_ppm10_and_flags_conflicts_and_stale_scores(tmp_path):
    _write_predictions(tmp_path / "old", "unified_rerank_all_top10.csv", [("X-00001", 1, "CCO"), ("X-00001", 2, "CCC")])
    _write_predictions(tmp_path / "new", "unified_rerank_all_top10_ppm10.csv", [("X-00001", 1, "CCN")])
    current = hashlib.sha256(RL.CALIBRATION.read_bytes()).hexdigest()
    for name, checksum in (("fresh", current), ("stale", "0" * 64)):
        scores = tmp_path / name / "scores"
        scores.mkdir(parents=True)
        (scores / "X-00001.json").write_text(json.dumps(
            {"schema_version": 1, "control_id": "X-00001", "status": "scored", "calibration_checksum": checksum}))
    result = RL.resolve([tmp_path], {"X-00001", "X-00002"})
    assert result["chosen"]["predictions"]["root"].endswith("new")
    assert result["chosen"]["predictions_coverage"]["missing"] == ["X-00002"]
    assert result["chosen"]["scores"]["root"].endswith("fresh/scores")
    assert result["rank1_conflicts"][0]["control_id"] == "X-00001"
    assert any("stale" in w and "/stale/" in w for w in result["warnings"])


def test_compare_ignores_smiles_serialization_and_reports_real_changes():
    reference = {"X-1": {"x_compound_id": "X-1", "msi_level": "4", "rank1_candidate_smiles": "OCC"},
                 "X-2": {"x_compound_id": "X-2", "msi_level": "3", "rank1_candidate_smiles": "CCN"}}
    new = {"X-1": {"x_compound_id": "X-1", "msi_level": "4.0", "rank1_candidate_smiles": "CCO"},
           "X-2": {"x_compound_id": "X-2", "msi_level": "4", "rank1_candidate_smiles": "CCN"}}
    pytest.importorskip("rdkit")
    result = CMP.compare(new, reference)
    assert result["identical_compounds"] == 1
    assert result["diffs"] == [{"x_compound_id": "X-2", "field": "msi_level", "reference": "3", "new": "4"}]


def _write_v2_root(root: Path, xid: str, with_structure: bool) -> None:
    data = root / xid / "data"
    data.mkdir(parents=True)
    low = xid.lower()
    (data / f"{low}_raw_top_candidates.csv").write_text(f"feature_key\nfeature_X_209_mz128.0715_ri1157_{xid}\n")
    (data / f"{low}_pathway_prediction.csv").write_text(
        "control_id,predicted_pathway,predicted_subpathway,pathway_source\n"
        f"{xid},Lipid,Eicosanoid,known_x_annotation\n")
    (data / f"{low}_known_x_row.csv").write_text("Unknown,Formula,Level,SubLevel\n" f"{xid},C6H11NO2,3,3C\n")
    (data / f"{low}_top_candidates.csv").write_text("smiles\nCCO\n" if with_structure else "smiles\n")
    (data / "similarity_neighbors.csv").write_text("canonical_name,similarity_to_target\n")


def test_client_deck_sheet_rows_do_not_back_fill(tmp_path):
    import client_deck as CD

    _write_v2_root(tmp_path, "X-11111", with_structure=True)
    _write_v2_root(tmp_path, "X-22222", with_structure=False)
    summary = {"rank1_candidate_smiles": "CCO", "rank1_candidate_formula": "", "report_display_name": "name",
               "msi_level": "4.0", "msi_label": "Molecular formula"}
    with_structure = CD.build_row(tmp_path, "X-11111", summary)
    without = CD.build_row(tmp_path, "X-22222", summary)
    assert with_structure[0] == "X-11111" and with_structure[1] == 128.0715
    assert with_structure[2] == "C6H11NO2" and with_structure[5] == "CCO" and with_structure[9] == "4"
    assert with_structure[3:5] == ["Lipid", "Eicosanoid"] and with_structure[11] == "Curated X annotation"
    assert without[5] == "" and without[7] == "" and without[8] is None and without[12] == "3/3C"


def test_client_deck_refuses_stale_trees(tmp_path):
    import client_deck as CD

    for name in ("reports 4", "se_report_inputs"):
        stale = tmp_path / name
        stale.mkdir()
        with pytest.raises(SystemExit):
            CD.check_reports_root(stale)
    with pytest.raises(SystemExit):
        CD.check_reports_root(tmp_path)
