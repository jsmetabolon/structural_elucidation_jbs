#!/usr/bin/env python3
"""Build a joined all-metabolites workbook across three inventory universes.

The output is intentionally self-contained: source workbooks are read-only inputs,
all percentages are stored as decimal Excel percentages, and missing measurements
remain blank rather than being converted to measured zeroes.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MFP = Path.home() / "Downloads" / "MFP"
DEFAULT_DOWNLOADS = Path.home() / "Downloads"
DEFAULT_SYNONYM_DB = Path.home() / "Documents" / "Metabolite Synonym Database"

MATRIX_COLUMNS = [
    "Human Plasma",
    "Human Serum",
    "Human Feces",
    "Mouse Serum",
    "Human Feces (OMNI)",
    "Mouse Feces",
    "Human Urine",
    "Canine Feces",
    "Human CSF",
    "Mouse Tissue",
    "Mouse Urine",
    "Human Saliva",
]

UNAVAILABLE_MATRICES = {
    "Mouse Serum",
    "Canine Feces",
    "Human CSF",
    "Mouse Tissue",
    "Mouse Urine",
}

DATA_COLUMNS = [
    "Entity Type",
    "ChemID",
    "Comp ID",
    "X Compound ID",
    "Pathway",
    "Subpathway",
    "Name",
    "MSI Level",
    "MSI Level Source",
    "MSI Sub-level",
    "MSI Sub-level Source",
    "Top Match",
    "Top Match Formula",
    "Confidence Score",
    "Confidence Basis",
    "Pathway Source",
    *MATRIX_COLUMNS,
]

GLOBAL_RATE_HEADERS = {
    "Human Plasma": "hPlasma HD4",
    "Human Serum": "hSerum HD4",
    "Human Feces": "hFecal HD4",
    "Human Feces (OMNI)": "DNAG HD4",
    "Mouse Feces": "mouseFecal HD4",
    "Human Urine": "hUrine HD4",
}

X_RATE_HEADERS = {
    "Human Plasma": "hPlasma Detection",
    "Human Urine": "hUrine Detection",
    "Human Feces": "hFecal Detection",
}

NOTE_PATHWAYS = {
    "unknown",
    "do not include",
    "ask dts",
    "can add this",
    "none",
    "nan",
    "#ref!",
}

HEADER_FILL = PatternFill("solid", fgColor="1F4E79")
HEADER_FONT = Font(color="FFFFFF", bold=True)
PASS_FILL = PatternFill("solid", fgColor="C6E0B4")
WARN_FILL = PatternFill("solid", fgColor="FFE699")


def clean(value: Any) -> str:
    """Normalize a scalar to a stripped display string."""
    if value is None:
        return ""
    text = re.sub(r"\s+", " ", str(value).strip())
    return "" if text.lower() in {"nan", "none", "null"} else text


def normalize_id(value: Any) -> str:
    """Normalize spreadsheet numeric IDs without changing X identifiers."""
    text = clean(value)
    if not text:
        return ""
    try:
        number = float(text)
    except ValueError:
        return text
    if math.isfinite(number) and number.is_integer():
        return str(int(number))
    return text


def parse_rate(value: Any) -> float | None:
    """Parse either a decimal rate or a 0-100 percentage to a 0-1 rate."""
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, str):
        text = value.strip().replace("%", "")
        if not text or "not on" in text.lower() or text.lower() in {"n/a", "na", "nan", "none"}:
            return None
    else:
        text = value
    try:
        number = float(text)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number) or number < 0:
        return None
    if number > 1:
        number /= 100.0
    return number if number <= 1 else None


def first_present(values: Iterable[Any]) -> str:
    for value in values:
        text = clean(value)
        if text:
            return text
    return ""


def cleaned_name(*values: Any) -> str:
    text = first_present(values)
    return re.sub(r"^(KNOWN|Released)\s+", "", text, flags=re.IGNORECASE).strip()


def valid_pathway(value: Any) -> str:
    text = clean(value)
    return "" if text.lower() in NOTE_PATHWAYS else text


def sheet_rows(path: Path, sheet: str, header_row: int = 1) -> list[dict[str, Any]]:
    """Read one worksheet into dictionaries while preserving formulas as cached values."""
    workbook = load_workbook(path, read_only=True, data_only=True)
    worksheet = workbook[sheet]
    rows = worksheet.iter_rows(values_only=True)
    for _ in range(header_row - 1):
        next(rows)
    headers = [clean(value) for value in next(rows)]
    output = []
    for raw in rows:
        row = {header: raw[index] if index < len(raw) else None for index, header in enumerate(headers) if header}
        output.append(row)
    workbook.close()
    return output


def normalize_x_known(row: dict[str, Any]) -> None:
    """Represent the legacy X 'KNOWN' identity tag as MSI level 1."""
    if row.get("Entity Type") != "X":
        return
    for field in ("MSI Level", "MSI Sub-level"):
        if str(row.get(field) or "").strip().upper() == "KNOWN":
            row[field] = "1"


def base_row(
    *,
    entity_type: str,
    chem_id: Any = "",
    comp_id: Any = "",
    x_id: Any = "",
    pathway: Any = "",
    subpathway: Any = "",
    name: Any = "",
    msi: Any = "",
    msi_source: str = "",
    msi_sublevel: Any = "",
    msi_sublevel_source: str = "",
) -> dict[str, Any]:
    row = {column: None for column in DATA_COLUMNS}
    row.update(
        {
            "Entity Type": entity_type,
            "ChemID": normalize_id(chem_id),
            "Comp ID": normalize_id(comp_id),
            "X Compound ID": clean(x_id),
            "Pathway": valid_pathway(pathway),
            "Subpathway": valid_pathway(subpathway),
            "Name": clean(name),
            "MSI Level": normalize_id(msi),
            "MSI Level Source": msi_source,
            "MSI Sub-level": clean(msi_sublevel),
            "MSI Sub-level Source": msi_sublevel_source,
            "Pathway Source": "named library" if entity_type == "Known" else "curated X metadata",
        }
    )
    normalize_x_known(row)
    return row


def load_snapshot_named(path: Path) -> list[dict[str, Any]]:
    rows = sheet_rows(path, "MSI Levels for Library (5617)")
    output = []
    for source in rows:
        chem_id = normalize_id(source.get("CHEMICAL_ID"))
        if not chem_id:
            continue
        pathway = first_present(
            [source.get("New SuperPathway (only for partially characterized)"), source.get("SUPER_PATHWAY")]
        )
        subpathway = first_present(
            [source.get("New SubPathway (only for partially characterized)"), source.get("PATHWAY")]
        )
        output.append(
            base_row(
                entity_type="Known",
                chem_id=chem_id,
                pathway=pathway,
                subpathway=subpathway,
                name=source.get("CHEMICAL_NAME"),
                msi=source.get("MSI Level"),
                msi_source="20250728 named library",
                msi_sublevel=source.get("MSI SubLevel"),
                msi_sublevel_source="20250728 named library",
            )
        )
    return output


def snapshot_x_source_rows(path: Path) -> list[dict[str, Any]]:
    return sheet_rows(path, "JJ + MM Merged Info", header_row=2)


def load_snapshot_x(source_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output = []
    for source in source_rows:
        x_id = clean(source.get("Unknown"))
        if not x_id:
            continue
        output.append(
            base_row(
                entity_type="X",
                comp_id=source.get("CompID"),
                x_id=x_id,
                pathway=source.get("SuperPathway"),
                subpathway=first_present([source.get("Verified Subpathway"), source.get("Draft Subpathway")]),
                name=cleaned_name(
                    source.get("Cleaned External Info"),
                    source.get("Tentative information"),
                    source.get("KNOWN"),
                    source.get("External Information (Class or modification)"),
                ),
                msi=source.get("Level"),
                msi_source="20250728 X metadata",
                msi_sublevel=source.get("SubLevel"),
                msi_sublevel_source="20250728 X metadata",
            )
        )
    return output


def load_latest_named(global_path: Path) -> tuple[list[dict[str, Any]], dict[str, dict[str, float]]]:
    source_rows = sheet_rows(global_path, "GLOBAL 20260129", header_row=2)
    output = []
    detection: dict[str, dict[str, float]] = {}
    for source in source_rows:
        chem_id = normalize_id(source.get("CHEMICAL_ID"))
        if not chem_id:
            continue
        output.append(
            base_row(
                entity_type="Known",
                chem_id=chem_id,
                pathway=first_present([source.get("Animal Super"), source.get("Orig SuperPathway")]),
                subpathway=first_present([source.get("Animal Sub"), source.get("Orig SubPathway")]),
                name=source.get("Metabolite"),
                msi_source="",
            )
        )
        rates = {}
        for matrix, header in GLOBAL_RATE_HEADERS.items():
            rate = parse_rate(source.get(header))
            if rate is not None:
                rates[matrix] = rate
        detection[chem_id] = rates
    return output, detection


def load_current_x(path: Path) -> list[dict[str, Any]]:
    output = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for source in csv.DictReader(handle):
            x_id = clean(source.get("Unknown"))
            if not x_id:
                continue
            output.append(
                base_row(
                    entity_type="X",
                    comp_id=source.get("CompID"),
                    x_id=x_id,
                    pathway=source.get("SuperPathway"),
                    subpathway=source.get("Draft Subpathway"),
                    name=cleaned_name(
                        source.get("Cleaned External Info"),
                        source.get("Tentative information"),
                        source.get("KNOWN"),
                        source.get("External Information (Class or modification)"),
                    ),
                    msi=source.get("Level"),
                    msi_source="known_X.csv fallback",
                    msi_sublevel=source.get("SubLevel"),
                    msi_sublevel_source="known_X.csv",
                )
            )
    return output


def collapse_x_rows(rows: list[dict[str, Any]], qa: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[row["X Compound ID"]].append(row)
    output = []
    fields = ["Pathway", "Subpathway", "Name", "MSI Level", "MSI Sub-level"]
    for x_id, members in grouped.items():
        row = dict(members[0])
        row["Comp ID"] = "; ".join(sorted({clean(item["Comp ID"]) for item in members if clean(item["Comp ID"])}))
        for field in fields:
            values = [clean(item[field]) for item in members if clean(item[field])]
            distinct = list(dict.fromkeys(values))
            row[field] = distinct[0] if distinct else ""
            if len(distinct) > 1:
                qa.append(
                    qa_row(
                        "X metadata conflict",
                        "REVIEW",
                        x_id,
                        "",
                        f"{field}: {' | '.join(distinct)}; first non-empty retained",
                    )
                )
        output.append(row)
    return sorted(output, key=lambda row: row["X Compound ID"])


def load_predictions(path: Path) -> dict[str, dict[str, Any]]:
    output = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for source in csv.DictReader(handle):
            x_id = clean(source.get("x_compound_id"))
            if not x_id:
                continue
            output[x_id] = {
                "Top Match": clean(source.get("rank1_candidate_smiles")),
                "Top Match Formula": clean(source.get("rank1_candidate_formula")),
                "Confidence Score": parse_rate(source.get("rank1_evidence_confidence")),
                "Confidence Basis": clean(source.get("rank1_evidence_confidence_basis")),
                "msi": normalize_id(source.get("msi_level")),
                "msi_label": clean(source.get("msi_label")),
                "msi_sublevel": clean(source.get("metabolon_sublevel")),
                "feature_name": clean(source.get("feature_name")),
            }
    return output


def load_report_pathways(reports_root: Path | None) -> dict[str, dict[str, str]]:
    if reports_root is None or not reports_root.is_dir():
        return {}
    output = {}
    for directory in reports_root.iterdir():
        if not directory.is_dir() or not directory.name.startswith("X-"):
            continue
        x_id = directory.name
        path = directory / "data" / f"{x_id.lower()}_top_candidates.csv"
        if not path.exists():
            continue
        try:
            with path.open(newline="", encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle)
                first = next(reader, None)
        except (OSError, csv.Error):
            continue
        if not first:
            continue
        output[x_id] = {
            "pathway": valid_pathway(first.get("predicted_pathway")),
            "subpathway": valid_pathway(first.get("predicted_subpathway")),
        }
    return output


def enrich_x_rows(
    rows: list[dict[str, Any]],
    predictions: dict[str, dict[str, Any]],
    report_pathways: dict[str, dict[str, str]],
) -> None:
    for row in rows:
        if row["Entity Type"] != "X":
            continue
        x_id = row["X Compound ID"]
        prediction = predictions.get(x_id, {})
        for field in ["Top Match", "Top Match Formula", "Confidence Score", "Confidence Basis"]:
            row[field] = prediction.get(field)
        if prediction.get("msi"):
            row["MSI Level"] = prediction["msi"]
            row["MSI Level Source"] = "delivered X report"
        if not row["MSI Sub-level"] and prediction.get("msi_sublevel"):
            row["MSI Sub-level"] = prediction["msi_sublevel"]
            row["MSI Sub-level Source"] = "exec summary metabolon sub-level fallback"
        normalize_x_known(row)
        if not row["Name"] and prediction.get("feature_name"):
            row["Name"] = prediction["feature_name"]
        predicted_path = report_pathways.get(x_id, {})
        if not valid_pathway(row["Pathway"]) and predicted_path.get("pathway"):
            row["Pathway"] = predicted_path["pathway"]
            row["Pathway Source"] = "rank-1 prediction"
        if not valid_pathway(row["Subpathway"]) and predicted_path.get("subpathway"):
            row["Subpathway"] = predicted_path["subpathway"]
            if row["Pathway Source"] != "rank-1 prediction":
                row["Pathway Source"] = "curated pathway; predicted subpathway"


def add_named_msi(latest_named: list[dict[str, Any]], snapshot_named: list[dict[str, Any]]) -> None:
    by_id = {row["ChemID"]: row for row in snapshot_named}
    for row in latest_named:
        source = by_id.get(row["ChemID"], {})
        if source.get("MSI Level"):
            row["MSI Level"] = source["MSI Level"]
            row["MSI Level Source"] = "20250728 named library"
        if source.get("MSI Sub-level"):
            row["MSI Sub-level"] = source["MSI Sub-level"]
            row["MSI Sub-level Source"] = "20250728 named library"


def load_saliva_rates(path: Path) -> tuple[dict[str, float], int]:
    rows = sheet_rows(path, "Select chro_curation_cpd")
    sample_sets: set[str] = set()
    by_chem: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        chem_id = normalize_id(row.get("CHEMICAL_ID"))
        sample_set = normalize_id(row.get("SAMPLE_SET_ID"))
        if not chem_id or not sample_set:
            continue
        sample_sets.add(sample_set)
        by_chem[chem_id].add(sample_set)
    denominator = len(sample_sets)
    rates = {chem_id: len(found) / denominator for chem_id, found in by_chem.items()} if denominator else {}
    return rates, denominator


def build_x_detection(
    snapshot_sources: list[dict[str, Any]], feces_history: Path, qa: list[dict[str, Any]]
) -> dict[str, dict[str, float]]:
    values: dict[str, dict[str, set[float]]] = defaultdict(lambda: defaultdict(set))
    for row in snapshot_sources:
        x_id = clean(row.get("Unknown"))
        if not x_id:
            continue
        for matrix, header in X_RATE_HEADERS.items():
            rate = parse_rate(row.get(header))
            if rate is not None:
                values[x_id][matrix].add(rate)

    output: dict[str, dict[str, float]] = {}
    for x_id, matrices in values.items():
        output[x_id] = {}
        for matrix, rates in matrices.items():
            selected = max(rates)
            output[x_id][matrix] = selected
            if len(rates) > 1:
                qa.append(
                    qa_row(
                        "X detection conflict",
                        "REVIEW",
                        x_id,
                        matrix,
                        f"rates={sorted(rates)}; maximum retained",
                    )
                )

    history_rows = sheet_rows(feces_history, "unk_feces_detection_rates")
    for row in history_rows:
        x_id = clean(row.get("CHEMICAL_NAME"))
        rate = parse_rate(row.get("PCT"))
        if x_id and rate is not None:
            output.setdefault(x_id, {})["Human Feces"] = rate
    return output


def apply_detection(
    rows: list[dict[str, Any]],
    named_detection: dict[str, dict[str, float]],
    x_detection: dict[str, dict[str, float]],
) -> None:
    for row in rows:
        rates = named_detection.get(row["ChemID"], {}) if row["Entity Type"] == "Known" else x_detection.get(
            row["X Compound ID"], {}
        )
        for matrix in MATRIX_COLUMNS:
            row[matrix] = rates.get(matrix)


def logical_key(row: dict[str, Any]) -> str:
    return f"CHEM:{row['ChemID']}" if row["Entity Type"] == "Known" else f"X:{row['X Compound ID']}"


def qa_row(check: str, status: str, value: Any, expected: Any, details: str) -> dict[str, Any]:
    return {"Check": check, "Status": status, "Value": value, "Expected": expected, "Details": details}


def build_membership(universes: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    labels = list(universes)
    metadata: dict[str, dict[str, Any]] = {}
    memberships: dict[str, set[str]] = defaultdict(set)
    for label, rows in universes.items():
        for row in rows:
            key = logical_key(row)
            memberships[key].add(label)
            metadata.setdefault(key, row)
    output = []
    for key in sorted(metadata):
        row = metadata[key]
        output.append(
            {
                "Logical Key": key,
                "Entity Type": row["Entity Type"],
                "ChemID": row["ChemID"],
                "X Compound ID": row["X Compound ID"],
                "Name": row["Name"],
                "MSI Level": row["MSI Level"],
                "MSI Sub-level": row["MSI Sub-level"],
                **{label: label in memberships[key] for label in labels},
                "Membership Count": len(memberships[key]),
            }
        )
    return output


def summary_rows(universes: dict[str, list[dict[str, Any]]], membership: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output = []
    for label, rows in universes.items():
        keys = {logical_key(row) for row in rows}
        output.extend(
            [
                {"Scope": label, "Metric": "Rows", "Value": len(rows)},
                {"Scope": label, "Metric": "Unique logical metabolites", "Value": len(keys)},
                {
                    "Scope": label,
                    "Metric": "Known rows",
                    "Value": sum(row["Entity Type"] == "Known" for row in rows),
                },
                {"Scope": label, "Metric": "X rows", "Value": sum(row["Entity Type"] == "X" for row in rows)},
            ]
        )
    labels = list(universes)
    output.append({"Scope": "All universes", "Metric": "Union of logical metabolites", "Value": len(membership)})
    output.append(
        {
            "Scope": "All universes",
            "Metric": "Present in all three",
            "Value": sum(all(row[label] for label in labels) for row in membership),
        }
    )
    for label in labels:
        output.append(
            {
                "Scope": label,
                "Metric": "Exclusive logical metabolites",
                "Value": sum(row[label] and sum(bool(row[item]) for item in labels) == 1 for row in membership),
            }
        )
    return output


def data_dictionary_rows() -> list[dict[str, Any]]:
    rows = [
        {
            "Column": "ChemID",
            "Meaning": "Named-library CHEMICAL_ID. Blank for X compounds.",
            "Source": "20250728 MSI library or 20260130 Global",
        },
        {
            "Column": "Comp ID",
            "Meaning": (
                "Library compound ID; semicolon-joined in the unique-X universe. Knowns are "
                "backfilled from chem_annotations.csv where a library entry exists (partial coverage)."
            ),
            "Source": "JJ + MM Merged Info, data/known_X.csv, or chem_annotations.csv",
        },
        {
            "Column": "X Compound ID",
            "Meaning": "Stable X-##### identifier. Blank for named compounds.",
            "Source": "Unknown column / exec summary",
        },
        {
            "Column": "MSI Level",
            "Meaning": "Named-library MSI or delivered report MSI; curated X level is a labeled fallback.",
            "Source": "20250728 MSI library, exec summary, known_X.csv",
        },
        {
            "Column": "MSI Sub-level",
            "Meaning": "Defined named-library or X sub-level; text values such as KNOWN, Released, and 2C.",
            "Source": "20250728 MSI library, JJ + MM Merged Info, known_X.csv, or exec-summary fallback",
        },
        {
            "Column": "Top Match",
            "Meaning": "Rank-1 candidate SMILES for X compounds.",
            "Source": "exec_summary_by_compound.csv",
        },
        {
            "Column": "Confidence Score",
            "Meaning": "Canonical rank-1 evidence_confidence; blank when unavailable.",
            "Source": "exec_summary_by_compound.csv",
        },
    ]
    source_by_matrix = {
        matrix: f"20260130 Global {header}; X subset from 20250728 merged data where available"
        for matrix, header in GLOBAL_RATE_HEADERS.items()
    }
    source_by_matrix["Human Saliva"] = "Saliva Human (2018-2023) v2; distinct detected SSIDs / eligible SSIDs"
    source_by_matrix["Human Feces"] += "; dedicated X feces history takes precedence"
    source_by_matrix["Human Feces (OMNI)"] = "20260130 Global DNAG HD4 (DNA Genotek / OMNI)"
    for matrix in MATRIX_COLUMNS:
        rows.append(
            {
                "Column": matrix,
                "Meaning": "Detection percentage; blank means unavailable, numeric zero means measured zero.",
                "Source": "No source supplied; intentionally blank"
                if matrix in UNAVAILABLE_MATRICES
                else source_by_matrix.get(matrix, "Available source workbook"),
            }
        )
    return rows


def style_sheet(worksheet, *, percentage_columns: set[str] | None = None) -> None:
    percentage_columns = percentage_columns or set()
    worksheet.freeze_panes = "A2"
    worksheet.auto_filter.ref = worksheet.dimensions
    for cell in worksheet[1]:
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(wrap_text=True, vertical="top")
    headers = {cell.value: cell.column for cell in worksheet[1]}
    for header, column in headers.items():
        width = min(max(len(str(header or "")) + 2, 12), 42)
        sample_cells = worksheet.iter_cols(
            min_col=column,
            max_col=column,
            min_row=2,
            max_row=min(80, worksheet.max_row),
        )
        for cell in list(sample_cells)[0]:
            width = min(max(width, len(str(cell.value or "")) + 2), 42)
        worksheet.column_dimensions[get_column_letter(column)].width = width
        if header in percentage_columns:
            for cell in worksheet.iter_cols(min_col=column, max_col=column, min_row=2, max_row=worksheet.max_row):
                for item in cell:
                    item.number_format = "0.0%"
        if header == "Confidence Score":
            for cell in worksheet.iter_cols(min_col=column, max_col=column, min_row=2, max_row=worksheet.max_row):
                for item in cell:
                    item.number_format = "0.000000"
    for row in worksheet.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = Alignment(vertical="top", wrap_text=cell.column <= 14)


def add_table_sheet(workbook: Workbook, title: str, columns: list[str], rows: list[dict[str, Any]]) -> None:
    worksheet = workbook.create_sheet(title)
    worksheet.append(columns)
    for row in rows:
        worksheet.append([row.get(column) for column in columns])
    style_sheet(worksheet, percentage_columns=set(MATRIX_COLUMNS))


def write_workbook(
    output: Path,
    universes: dict[str, list[dict[str, Any]]],
    membership: list[dict[str, Any]],
    qa: list[dict[str, Any]],
) -> None:
    workbook = Workbook()
    workbook.remove(workbook.active)
    for title, rows in universes.items():
        add_table_sheet(workbook, title, DATA_COLUMNS, rows)

    summary = summary_rows(universes, membership)
    add_table_sheet(workbook, "Universe_Summary", ["Scope", "Metric", "Value"], summary)

    membership_columns = [
        "Logical Key",
        "Entity Type",
        "ChemID",
        "X Compound ID",
        "Name",
        "MSI Level",
        "MSI Sub-level",
        *universes.keys(),
        "Membership Count",
    ]
    add_table_sheet(workbook, "Universe_Membership", membership_columns, membership)
    add_table_sheet(workbook, "Data_Dictionary", ["Column", "Meaning", "Source"], data_dictionary_rows())
    add_table_sheet(workbook, "QA", ["Check", "Status", "Value", "Expected", "Details"], qa)

    qa_sheet = workbook["QA"]
    status_col = next(cell.column for cell in qa_sheet[1] if cell.value == "Status")
    for cell in list(qa_sheet.iter_cols(min_col=status_col, max_col=status_col, min_row=2))[0]:
        if cell.value == "PASS":
            cell.fill = PASS_FILL
        elif cell.value in {"REVIEW", "FAIL"}:
            cell.fill = WARN_FILL

    output.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(output)
    workbook.close()


def validate(
    universes: dict[str, list[dict[str, Any]]],
    named_detection: dict[str, dict[str, float]],
    membership: list[dict[str, Any]],
    qa: list[dict[str, Any]],
    saliva_denominator: int,
) -> None:
    for title, rows in universes.items():
        qa.append(qa_row("Universe row count", "PASS", len(rows), len(rows), title))
        keys = {logical_key(row) for row in rows}
        qa.append(qa_row("Universe unique logical count", "PASS", len(keys), len(keys), title))
        qa.append(
            qa_row(
                "MSI sub-level populated",
                "PASS",
                sum(bool(row["MSI Sub-level"]) for row in rows),
                "source-derived",
                title,
            )
        )
        if title == "Latest_Unique_X":
            x_ids = [row["X Compound ID"] for row in rows if row["Entity Type"] == "X"]
            status = "PASS" if len(x_ids) == len(set(x_ids)) else "FAIL"
            qa.append(qa_row("Latest unique X IDs", status, len(set(x_ids)), len(x_ids), title))
        for matrix in MATRIX_COLUMNS:
            populated = [row[matrix] for row in rows if row[matrix] is not None]
            valid = all(isinstance(value, (int, float)) and 0 <= value <= 1 for value in populated)
            qa.append(
                qa_row(
                    "Detection range",
                    "PASS" if valid else "FAIL",
                    len(populated),
                    "all values 0-1",
                    f"{title}: {matrix}",
                )
            )
    qa.append(
        qa_row(
            "Membership reconciliation",
            "PASS",
            len(membership),
            len({logical_key(row) for rows in universes.values() for row in rows}),
            "Union of logical keys",
        )
    )
    for matrix in sorted(UNAVAILABLE_MATRICES):
        populated = sum(row[matrix] is not None for rows in universes.values() for row in rows)
        qa.append(qa_row("Unavailable matrix blank", "PASS" if populated == 0 else "FAIL", populated, 0, matrix))
    qa.append(
        qa_row(
            "Saliva denominator",
            "PASS" if saliva_denominator else "FAIL",
            saliva_denominator,
            ">0",
            "SSID count",
        )
    )
    qa.append(qa_row("Named detection IDs", "PASS", len(named_detection), len(named_detection), "20260130 Global"))


def load_known_comp_ids(path: Path) -> dict[str, str]:
    """Map ChemID -> Comp ID from the synonym-database annotation export.

    Coverage is partial (only library-registered knowns), so callers must
    tolerate ChemIDs that are absent from this mapping.
    """
    mapping: dict[str, str] = {}
    if not path.exists():
        return mapping
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as handle:
        for source in csv.DictReader(handle):
            chem_id = normalize_id(source.get("CHEM_ID"))
            comp_id = normalize_id(source.get("COMP_ID"))
            if chem_id and comp_id:
                mapping.setdefault(chem_id, comp_id)
    return mapping


def apply_known_comp_ids(rows: Iterable[dict[str, Any]], mapping: dict[str, str]) -> None:
    """Backfill blank Comp IDs on Known rows from the ChemID->Comp ID mapping."""
    for row in rows:
        if row.get("Entity Type") != "Known":
            continue
        if row.get("Comp ID"):
            continue
        comp_id = mapping.get(str(row.get("ChemID") or ""))
        if comp_id:
            row["Comp ID"] = comp_id


def build(args: argparse.Namespace) -> dict[str, Any]:
    qa: list[dict[str, Any]] = []
    snapshot_named = load_snapshot_named(args.levels_workbook)
    snapshot_source_rows = snapshot_x_source_rows(args.levels_workbook)
    snapshot_x = load_snapshot_x(snapshot_source_rows)
    latest_named, named_detection = load_latest_named(args.global_workbook)
    current_x = load_current_x(args.known_x)
    latest_unique_x = collapse_x_rows(current_x, qa)
    add_named_msi(latest_named, snapshot_named)

    saliva_rates, saliva_denominator = load_saliva_rates(args.saliva_workbook)
    for chem_id, rate in saliva_rates.items():
        named_detection.setdefault(chem_id, {})["Human Saliva"] = rate
    x_detection = build_x_detection(snapshot_source_rows, args.x_feces_history, qa)

    predictions = load_predictions(args.exec_summary)
    report_pathways = load_report_pathways(args.reports_root)
    known_comp_ids = load_known_comp_ids(args.chem_annotations)
    universes = {
        "Snapshot_20250728": snapshot_named + snapshot_x,
        "Latest_Unique_X": [dict(row) for row in latest_named] + latest_unique_x,
        "Latest_PerCompID": [dict(row) for row in latest_named] + current_x,
    }
    for rows in universes.values():
        enrich_x_rows(rows, predictions, report_pathways)
        apply_known_comp_ids(rows, known_comp_ids)
        apply_detection(rows, named_detection, x_detection)

    membership = build_membership(universes)
    validate(universes, named_detection, membership, qa, saliva_denominator)
    write_workbook(args.output, universes, membership, qa)
    return {
        "output": args.output,
        "universes": {title: len(rows) for title, rows in universes.items()},
        "membership": len(membership),
        "qa_reviews": sum(row["Status"] == "REVIEW" for row in qa),
        "qa_failures": sum(row["Status"] == "FAIL" for row in qa),
    }


def parser() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument(
        "--global-workbook",
        type=Path,
        default=DEFAULT_MFP / "20260130 Microbiome vs. Global Detection (9).xlsx",
    )
    cli.add_argument(
        "--levels-workbook",
        type=Path,
        default=DEFAULT_DOWNLOADS / "20250728 Unknown Info with Levels & Detection.xlsx",
    )
    cli.add_argument(
        "--x-feces-history",
        type=Path,
        default=DEFAULT_DOWNLOADS / "20250728 hFeces X-compound (Unknown) Detection History.xlsx",
    )
    cli.add_argument("--saliva-workbook", type=Path, default=DEFAULT_MFP / "Saliva Human (2018-2023) v2.xlsx")
    cli.add_argument("--known-x", type=Path, default=PROJECT_ROOT / "data" / "known_X.csv")
    cli.add_argument("--chem-annotations", type=Path, default=DEFAULT_SYNONYM_DB / "chem_annotations.csv")
    cli.add_argument(
        "--exec-summary",
        type=Path,
        default=PROJECT_ROOT / "docs" / "paper" / "data" / "exec_summary" / "exec_summary_by_compound.csv",
    )
    cli.add_argument("--reports-root", type=Path, default=PROJECT_ROOT.parent / "reports 4")
    cli.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT / "docs" / "paper" / "data" / "all_metabolites_matrix_detection.xlsx",
    )
    return cli


def main() -> None:
    args = parser().parse_args()
    for name in [
        "global_workbook",
        "levels_workbook",
        "x_feces_history",
        "saliva_workbook",
        "known_x",
        "exec_summary",
    ]:
        path = getattr(args, name)
        if not path.exists():
            raise SystemExit(f"missing --{name.replace('_', '-')}: {path}")
    result = build(args)
    print(f"[all-metabolites] wrote {result['output']}")
    for title, count in result["universes"].items():
        print(f"[all-metabolites] {title}: {count:,} rows")
    print(
        f"[all-metabolites] membership: {result['membership']:,} logical metabolites | "
        f"QA reviews={result['qa_reviews']} failures={result['qa_failures']}"
    )


if __name__ == "__main__":
    main()
