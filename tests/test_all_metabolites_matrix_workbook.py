"""Tests for the all-metabolites matrix workbook builder."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

from openpyxl import Workbook, load_workbook

ROOT = Path(__file__).resolve().parent.parent


def load_tool():
    name = "build_all_metabolites_matrix_workbook"
    tool_path = ROOT / "tools" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, tool_path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


M = load_tool()


def test_normalize_id_removes_spreadsheet_decimal_only():
    assert M.normalize_id(47607.0) == "47607"
    assert M.normalize_id("100000842.0") == "100000842"
    assert M.normalize_id("X-11677") == "X-11677"


def test_parse_rate_accepts_decimal_and_percent_scales():
    assert M.parse_rate(0.25) == 0.25
    assert M.parse_rate(25) == 0.25
    assert M.parse_rate("25%") == 0.25
    assert M.parse_rate(0) == 0


def test_parse_rate_keeps_unavailable_distinct_from_zero():
    assert M.parse_rate(None) is None
    assert M.parse_rate("Not on Plasma List") is None
    assert M.parse_rate("N/A") is None


def test_unique_x_collapses_compids_and_reports_metadata_conflicts():
    qa = []
    first = M.base_row(
        entity_type="X",
        comp_id="20",
        x_id="X-1",
        pathway="Lipid",
        name="candidate",
        msi_sublevel="2C",
        msi_sublevel_source="known_X.csv",
    )
    second = M.base_row(
        entity_type="X",
        comp_id="10",
        x_id="X-1",
        pathway="Amino Acid",
        name="candidate",
        msi_sublevel="Released",
        msi_sublevel_source="known_X.csv",
    )
    result = M.collapse_x_rows([first, second], qa)
    assert len(result) == 1
    assert result[0]["Comp ID"] == "10; 20"
    assert result[0]["Pathway"] == "Lipid"
    assert result[0]["MSI Sub-level"] == "2C"
    assert any(row["Check"] == "X metadata conflict" for row in qa)


def test_prediction_join_uses_canonical_confidence_and_report_msi():
    row = M.base_row(
        entity_type="X",
        x_id="X-1",
        msi="4",
        msi_source="fallback",
    )
    M.enrich_x_rows(
        [row],
        {
            "X-1": {
                "Top Match": "CCO",
                "Top Match Formula": "C2H6O",
                "Confidence Score": 0.91,
                "Confidence Basis": "backfilled",
                "msi": "3",
                "msi_sublevel": "2C",
                "feature_name": "ethanol-like",
            }
        },
        {"X-1": {"pathway": "Lipid", "subpathway": "Fatty Alcohol"}},
    )
    assert row["Top Match"] == "CCO"
    assert row["Confidence Score"] == 0.91
    assert row["MSI Level"] == "3"
    assert row["MSI Level Source"] == "delivered X report"
    assert row["MSI Sub-level"] == "2C"
    assert row["MSI Sub-level Source"] == "exec summary metabolon sub-level fallback"
    assert row["Pathway"] == "Lipid"


def test_curated_x_sublevel_is_not_overwritten_by_exec_summary():
    row = M.base_row(
        entity_type="X",
        x_id="X-1",
        msi_sublevel="3B",
        msi_sublevel_source="known_X.csv",
    )
    M.enrich_x_rows(
        [row],
        {"X-1": {"msi_sublevel": "2C"}},
        {},
    )
    assert row["MSI Sub-level"] == "3B"
    assert row["MSI Sub-level Source"] == "known_X.csv"


def test_x_known_tag_is_normalized_to_level_1():
    row = M.base_row(
        entity_type="X",
        x_id="X-1",
        msi="KNOWN",
        msi_sublevel="KNOWN",
    )
    assert row["MSI Level"] == "1"
    assert row["MSI Sub-level"] == "1"


def test_known_named_row_keeps_known_text_untouched():
    row = M.base_row(entity_type="Known", chem_id="5", msi="1", msi_sublevel="KNOWN")
    assert row["MSI Sub-level"] == "KNOWN"


def test_apply_known_comp_ids_backfills_only_blank_known_rows():
    known_blank = M.base_row(entity_type="Known", chem_id="35")
    known_filled = M.base_row(entity_type="Known", chem_id="36", comp_id="999")
    x_row = M.base_row(entity_type="X", x_id="X-1", chem_id="35")
    mapping = {"35": "42370", "36": "111"}
    M.apply_known_comp_ids([known_blank, known_filled, x_row], mapping)
    assert known_blank["Comp ID"] == "42370"
    assert known_filled["Comp ID"] == "999"
    assert x_row["Comp ID"] == ""


def test_named_msi_and_sublevel_are_joined_by_chemid():
    snapshot = M.base_row(
        entity_type="Known",
        chem_id="100",
        msi="1",
        msi_source="20250728 named library",
        msi_sublevel="KNOWN",
        msi_sublevel_source="20250728 named library",
    )
    latest = M.base_row(entity_type="Known", chem_id="100")
    M.add_named_msi([latest], [snapshot])
    assert latest["MSI Level"] == "1"
    assert latest["MSI Sub-level"] == "KNOWN"
    assert latest["MSI Sub-level Source"] == "20250728 named library"


def test_dedicated_feces_history_overrides_merged_detection(tmp_path):
    path = tmp_path / "feces.xlsx"
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "unk_feces_detection_rates"
    sheet.append(["CHEMICAL_NAME", "NUM_SETS_COMPLETE", "PCT"])
    sheet.append(["X-1", 57, 75])
    workbook.save(path)

    qa = []
    source = [
        {
            "Unknown": "X-1",
            "hPlasma Detection": 50,
            "hUrine Detection": "Not on Urine List",
            "hFecal Detection": 25,
        }
    ]
    result = M.build_x_detection(source, path, qa)
    assert result["X-1"]["Human Plasma"] == 0.5
    assert result["X-1"]["Human Feces"] == 0.75
    assert "Human Urine" not in result["X-1"]


def test_saliva_rate_counts_distinct_ssids(tmp_path):
    path = tmp_path / "saliva.xlsx"
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Select chro_curation_cpd"
    sheet.append(["SAMPLE_SET_ID", "CHEMICAL_ID"])
    sheet.append([1, 100])
    sheet.append([1, 100])
    sheet.append([2, 100])
    sheet.append([1, 200])
    workbook.save(path)

    rates, denominator = M.load_saliva_rates(path)
    assert denominator == 2
    assert rates["100"] == 1.0
    assert rates["200"] == 0.5


def test_membership_is_unique_by_logical_metabolite():
    known = M.base_row(entity_type="Known", chem_id="100", name="known")
    x_first = M.base_row(
        entity_type="X", comp_id="1", x_id="X-1", name="unknown"
    )
    x_second = M.base_row(
        entity_type="X", comp_id="2", x_id="X-1", name="unknown"
    )
    universes = {
        "Snapshot_20250728": [known, x_first],
        "Latest_Unique_X": [known, x_first],
        "Latest_PerCompID": [known, x_first, x_second],
    }
    membership = M.build_membership(universes)
    assert len(membership) == 2
    x_row = next(row for row in membership if row["Logical Key"] == "X:X-1")
    assert x_row["Membership Count"] == 3


def test_written_workbook_has_percent_format_and_required_sheets(tmp_path):
    row = M.base_row(entity_type="Known", chem_id="100", name="known")
    row["Human Plasma"] = 0.5
    universes = {
        "Snapshot_20250728": [dict(row)],
        "Latest_Unique_X": [dict(row)],
        "Latest_PerCompID": [dict(row)],
    }
    membership = M.build_membership(universes)
    output = tmp_path / "result.xlsx"
    qa = [M.qa_row("test", "PASS", 1, 1, "")]
    M.write_workbook(output, universes, membership, qa)

    workbook = load_workbook(output, read_only=False, data_only=True)
    assert workbook.sheetnames == [
        "Snapshot_20250728",
        "Latest_Unique_X",
        "Latest_PerCompID",
        "Universe_Summary",
        "Universe_Membership",
        "Data_Dictionary",
        "QA",
    ]
    sheet = workbook["Snapshot_20250728"]
    headers = {cell.value: cell.column for cell in sheet[1]}
    assert sheet.cell(2, headers["Human Plasma"]).value == 0.5
    assert sheet.cell(2, headers["Human Plasma"]).number_format == "0.0%"
    assert sheet.cell(2, headers["Mouse Serum"]).value is None
    assert "MSI Sub-level" in headers
    membership = workbook["Universe_Membership"]
    assert "MSI Sub-level" in {cell.value for cell in membership[1]}
    workbook.close()
