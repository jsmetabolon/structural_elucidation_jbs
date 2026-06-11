#!/usr/bin/env python3
"""Legacy NCE workflow evidence helpers for the V2 integration branch.

This module ports the useful, portable parts of the older
NewCompoundElucidation/MS2SMILES workflow into the current spec2mol pipeline.
It intentionally avoids importing the old application runtime and avoids
optional spreadsheet dependencies such as openpyxl.
"""

from __future__ import annotations

import json
import math
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

from denovo_improvement_utils import (
    canonicalize_smiles,
    formula_counts,
    formula_from_smiles,
    pick_col,
    read_table,
)
from nce_evidence_utils import (
    formula_harmony_score,
    formula_mass,
    lib_ion_mode,
    mass_bin,
    motif_labels_for_legacy_fragment,
    neutral_mass_from_precursor,
    parse_json_list,
)
from x_domain_analog_utils import safe_float, safe_int, xid_from_feature_key


MS2SMILES_PATTERNS = [
    "ms2smiles-*.tsv",
    "ms2smiles-*.xlsx",
    "ms2smiles_predictions*.csv",
    "ms2-predictions-*.xlsx",
    "ms2-predictions-*.csv",
    "ms2-all-predictions*.xlsx",
    "ms2-all-predictions*.csv",
    "matches_*.xlsx",
    "double_*.xlsx",
    "smiles_table.csv",
    "unknown_list.csv",
]

REFERENCE_PATTERNS = [
    "fmla_finder_cache*",
    "harmony*",
    "ff_dbe*",
    "library_fragments*",
    "ms2smiles_models_*",
    "smiles2ms2_models_*",
    "ri_models_*",
    "chemont_models_*",
    "chemont_defs*",
    "mlims_tables*",
]

LEGACY_MOTIF_KEYWORDS = {
    "sulfate": ["sulfate", "sulfo", "s(=o)(=o)", "so3"],
    "glucuronide_like": ["glucuron", "hexuronic", "uronic"],
    "glycoside_like": ["glycos", "glucose", "hexose", "sugar", "ribose"],
    "phosphate": ["phosphate", "phospho", "p(=o)", "po3"],
    "peptide_glutathione": ["peptide", "glutathione", "cys", "cysteine", "gamma-glu", "glutam"],
    "taurine": ["taurine", "taur"],
    "carnitine": ["carnitine", "acylcarnitine"],
    "bile_acid_steroid": ["bile", "steroid", "chol", "sterol"],
    "glycine": ["glycine", "ncc(=o)"],
    "amide": ["amide", "c(=o)n"],
    "carboxyl": ["carboxyl", "c(=o)o"],
}


def normalize_column_name(value: Any) -> str:
    text = str(value or "").strip().strip('"').lstrip("\ufeff").lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


def _xlsx_column_index(cell_ref: str) -> int:
    letters = "".join(ch for ch in str(cell_ref) if ch.isalpha()).upper()
    out = 0
    for ch in letters:
        out = out * 26 + (ord(ch) - ord("A") + 1)
    return max(out - 1, 0)


def _xlsx_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    try:
        xml = zf.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    root = ET.fromstring(xml)
    ns = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    values: list[str] = []
    for si in root.findall("a:si", ns):
        text_parts = [node.text or "" for node in si.findall(".//a:t", ns)]
        values.append("".join(text_parts))
    return values


def _xlsx_sheet_names(zf: zipfile.ZipFile) -> list[str]:
    names = sorted(name for name in zf.namelist() if re.match(r"xl/worksheets/sheet\d+\.xml$", name))
    return names or [name for name in zf.namelist() if name.startswith("xl/worksheets/") and name.endswith(".xml")]


def read_xlsx_without_openpyxl(path: Path, max_sheets: int | None = 1) -> pd.DataFrame:
    """Read simple worksheet data from .xlsx OOXML using only stdlib XML.

    The legacy files use simple first-row headers and scalar cells, so this is
    enough and keeps the parser independent from openpyxl.
    """

    path = Path(path)
    with zipfile.ZipFile(path) as zf:
        shared_strings = _xlsx_shared_strings(zf)
        frames: list[pd.DataFrame] = []
        for sheet_idx, sheet_name in enumerate(_xlsx_sheet_names(zf)):
            if max_sheets is not None and sheet_idx >= int(max_sheets):
                break
            root = ET.fromstring(zf.read(sheet_name))
            ns = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
            rows: list[list[Any]] = []
            for row in root.findall(".//a:sheetData/a:row", ns):
                values: dict[int, Any] = {}
                for cell in row.findall("a:c", ns):
                    ref = cell.attrib.get("r", "")
                    idx = _xlsx_column_index(ref)
                    cell_type = cell.attrib.get("t", "")
                    value_node = cell.find("a:v", ns)
                    inline_node = cell.find("a:is/a:t", ns)
                    raw = value_node.text if value_node is not None else inline_node.text if inline_node is not None else ""
                    if cell_type == "s":
                        value = shared_strings[safe_int(raw, -1)] if 0 <= safe_int(raw, -1) < len(shared_strings) else ""
                    else:
                        value = raw or ""
                    values[idx] = value
                if values:
                    width = max(values) + 1
                    rows.append([values.get(i, "") for i in range(width)])
            if not rows:
                continue
            header = [normalize_column_name(value) for value in rows[0]]
            data = rows[1:]
            width = len(header)
            normalized = [row + [""] * max(width - len(row), 0) for row in data]
            frame = pd.DataFrame([row[:width] for row in normalized], columns=header)
            frame = frame.loc[:, [col for col in frame.columns if col]]
            frame["_legacy_sheet"] = Path(sheet_name).stem
            frames.append(frame)
        return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def read_legacy_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".xlsx":
        return read_xlsx_without_openpyxl(path)
    df = read_table(path)
    df.columns = [normalize_column_name(col) for col in df.columns]
    return df


def classify_legacy_source(path: Path, df: pd.DataFrame) -> str:
    name = path.name.lower()
    cols = set(df.columns)
    if name.startswith("matches_"):
        return "known_match"
    if name.startswith("double_"):
        return "double_match"
    if name == "smiles_table.csv":
        return "smiles_table"
    if name == "unknown_list.csv":
        return "unknown_metadata"
    if "fragment" in cols or "frag_chem_name" in cols or "prediction" in cols:
        return "fragment_prediction"
    if name.startswith("ms2"):
        return "fragment_prediction"
    return "other"


def first_present(row: pd.Series, names: Iterable[str], default: Any = "") -> Any:
    for name in names:
        if name in row.index:
            value = row.get(name)
            if value is not None and not (isinstance(value, float) and np.isnan(value)) and str(value) != "":
                return value
    return default


def motif_labels_for_legacy_text(*values: Any) -> tuple[str, ...]:
    text = " ".join(str(value or "").lower() for value in values)
    labels = set(motif_labels_for_legacy_fragment(text))
    for label, needles in LEGACY_MOTIF_KEYWORDS.items():
        if any(needle in text for needle in needles):
            labels.add(label)
    return tuple(sorted(labels))


def _float_or_blank(value: Any) -> float | str:
    out = safe_float(value, math.nan)
    return out if np.isfinite(out) else ""


def normalize_legacy_rows(path: Path, df: pd.DataFrame) -> pd.DataFrame:
    kind = classify_legacy_source(path, df)
    rows: list[dict[str, Any]] = []
    for _, row in df.iterrows():
        comp_id = safe_int(first_present(row, ["comp_id", "unknown_comp_id", "chemical_id"]), 0)
        known_comp_id = safe_int(first_present(row, ["known_comp_id", "best_comp_id", "second_comp_id"]), 0)
        lib_id = safe_int(first_present(row, ["lib_id", "library"]), 0)
        parent_mass = _float_or_blank(first_present(row, ["parent_mass", "unknown_parent_mass", "parent_mz"]))
        fragment = first_present(row, ["fragment", "frag_smiles", "prediction", "subsmiles", "predicted_fragment"])
        frag_name = first_present(row, ["frag_chem_name", "fragment_name"])
        details = first_present(row, ["details", "historical_notes", "notes"])
        unknown_name = first_present(row, ["unknown_name", "name", "unknown"])
        matched_known = first_present(row, ["known_name", "best_name", "second_name", "name"])
        repeated_header_values = {
            str(value or "").strip().lower()
            for value in [unknown_name, details, matched_known, fragment, frag_name]
            if str(value or "").strip()
        }
        if repeated_header_values and repeated_header_values.issubset(
            {"unknown_name", "details", "best_name", "second_name", "fragment", "frag_chem_name", "name"}
        ):
            continue
        motif_labels = motif_labels_for_legacy_text(fragment, frag_name, details, unknown_name, matched_known)
        rows.append(
            {
                "legacy_source_path": str(path),
                "legacy_source_file": path.name,
                "legacy_source_kind": kind,
                "comp_id": comp_id if comp_id else "",
                "lib_id": lib_id if lib_id else "",
                "ion_mode": lib_ion_mode(lib_id) if lib_id else "",
                "parent_mass": parent_mass,
                "unknown_name": unknown_name,
                "details": details,
                "predicted_fragment": fragment,
                "frag_chem_name": frag_name,
                "frag_len": first_present(row, ["frag_len"]),
                "fragment_mass": _float_or_blank(first_present(row, ["frag_exact_mass", "fragment_mass"])),
                "complementary_mass": _float_or_blank(first_present(row, ["complementary_mass"]),
                ),
                "probability": _float_or_blank(first_present(row, ["probability", "prediction_probability", "score", "confidence"])),
                "precision": _float_or_blank(first_present(row, ["precision"])),
                "recall": _float_or_blank(first_present(row, ["recall"])),
                "f1": _float_or_blank(first_present(row, ["f1"])),
                "matched_known_comp_id": known_comp_id if known_comp_id else "",
                "matched_known_name": matched_known,
                "match_score": _float_or_blank(first_present(row, ["score", "best_score"])),
                "second_score": _float_or_blank(first_present(row, ["second_score"])),
                "legacy_motif_labels": "|".join(motif_labels),
                "historical_notes": " | ".join(str(v) for v in [details, first_present(row, ["parent_name"])] if str(v or "").strip()),
            }
        )
    out = pd.DataFrame(rows)
    return out.dropna(how="all").reset_index(drop=True)


def discover_legacy_files(ms2smiles_dir: Path) -> list[Path]:
    root = Path(ms2smiles_dir).expanduser()
    if not root.exists():
        return []
    files: set[Path] = set()
    for pattern in MS2SMILES_PATTERNS:
        files.update(root.glob(pattern))
    return sorted(path for path in files if path.is_file())


def build_legacy_catalog(ms2smiles_dir: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in discover_legacy_files(ms2smiles_dir):
        try:
            df = read_legacy_table(path)
            if df.empty:
                continue
            normalized = normalize_legacy_rows(path, df)
            if not normalized.empty:
                frames.append(normalized)
        except Exception as exc:
            frames.append(
                pd.DataFrame(
                    [
                        {
                            "legacy_source_path": str(path),
                            "legacy_source_file": path.name,
                            "legacy_source_kind": "parse_error",
                            "parse_error": str(exc),
                        }
                    ]
                )
            )
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True, sort=False)
    if "comp_id" in out.columns:
        out["_sort_comp"] = pd.to_numeric(out["comp_id"], errors="coerce").fillna(0).astype(int)
        out = out.sort_values(["_sort_comp", "legacy_source_file"], kind="mergesort").drop(columns=["_sort_comp"])
    return out.reset_index(drop=True)


def _json_top_counter(counter: Counter, n: int = 20) -> str:
    return json.dumps(dict(counter.most_common(n)), sort_keys=True)


def summarize_legacy_catalog(catalog: pd.DataFrame) -> pd.DataFrame:
    if catalog.empty or "comp_id" not in catalog.columns:
        return pd.DataFrame(columns=["comp_id"])
    work = catalog.copy()
    work["_comp_id"] = pd.to_numeric(work["comp_id"], errors="coerce")
    work = work[work["_comp_id"].notna()].copy()
    if work.empty:
        return pd.DataFrame(columns=["comp_id"])
    rows: list[dict[str, Any]] = []
    for comp_id, group in work.groupby("_comp_id", dropna=True):
        motifs = Counter()
        for text in group.get("legacy_motif_labels", pd.Series(dtype=str)).fillna("").astype(str):
            for label in text.split("|"):
                if label:
                    motifs[label] += 1
        frag_masses = [
            round(value, 6)
            for value in pd.to_numeric(group.get("fragment_mass"), errors="coerce").dropna().astype(float).tolist()
            if value > 0
        ]
        comp_masses = [
            round(value, 6)
            for value in pd.to_numeric(group.get("complementary_mass"), errors="coerce").dropna().astype(float).tolist()
            if value > 0
        ]
        probs = pd.to_numeric(group.get("probability"), errors="coerce").dropna()
        match_scores = pd.to_numeric(group.get("match_score"), errors="coerce").dropna()
        match_names = [str(v) for v in group.get("matched_known_name", pd.Series(dtype=str)).fillna("").astype(str) if str(v)]
        rows.append(
            {
                "comp_id": int(comp_id),
                "legacy_record_count_v2": int(group.shape[0]),
                "legacy_fragment_count_v2": int(group["predicted_fragment"].fillna("").astype(str).ne("").sum())
                if "predicted_fragment" in group.columns
                else 0,
                "legacy_match_count_v2": int(group["matched_known_comp_id"].fillna("").astype(str).ne("").sum())
                if "matched_known_comp_id" in group.columns
                else 0,
                "legacy_top_motif_v2": motifs.most_common(1)[0][0] if motifs else "",
                "legacy_motifs_json_v2": _json_top_counter(motifs),
                "legacy_max_probability_v2": float(probs.max()) if not probs.empty else 0.0,
                "legacy_mean_probability_v2": float(probs.mean()) if not probs.empty else 0.0,
                "legacy_fragment_masses_json_v2": json.dumps(sorted(set(frag_masses))[:80]),
                "legacy_fragment_bins_json_v2": json.dumps(sorted({mass_bin(v) for v in frag_masses if v > 0})[:80]),
                "legacy_complementary_masses_json_v2": json.dumps(sorted(set(comp_masses))[:80]),
                "legacy_loss_bins_json_v2": json.dumps(sorted({mass_bin(v) for v in comp_masses if v > 0})[:80]),
                "legacy_top_match_name_v2": Counter(match_names).most_common(1)[0][0] if match_names else "",
                "legacy_top_match_score_v2": float(match_scores.max()) if not match_scores.empty else 0.0,
                "legacy_source_files_v2": "|".join(sorted(set(group["legacy_source_file"].fillna("").astype(str)))),
                "legacy_historical_notes_v2": " | ".join(
                    sorted({str(v) for v in group.get("historical_notes", pd.Series(dtype=str)).fillna("").astype(str) if str(v).strip()})[:8]
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("comp_id").reset_index(drop=True)


def load_known_x_comp_id_map(known_x: Path) -> dict[str, int]:
    path = Path(known_x)
    if not path.exists():
        return {}
    df = read_table(path)
    unknown_col = pick_col(df, ["Unknown", "unknown", "x_compound_id", "XID"])
    comp_col = pick_col(df, ["CompID", "comp_id", "compound_id"])
    if unknown_col is None or comp_col is None:
        return {}
    out: dict[str, int] = {}
    for _, row in df.iterrows():
        xid = str(row.get(unknown_col, "") or "").strip()
        comp_id = safe_int(row.get(comp_col), 0)
        if xid and comp_id:
            out[xid] = comp_id
    return out


def assign_missing_comp_ids_from_xids(feature_df: pd.DataFrame, known_x: Path) -> pd.DataFrame:
    if feature_df.empty or "feature_key" not in feature_df.columns:
        return feature_df.copy()
    x_map = load_known_x_comp_id_map(known_x)
    if not x_map:
        return feature_df.copy()
    out = feature_df.copy()
    if "comp_id" not in out.columns:
        out["comp_id"] = 0
    comp = pd.to_numeric(out["comp_id"], errors="coerce").fillna(0).astype(int)
    xids = out["feature_key"].astype(str).map(xid_from_feature_key)
    replacement = xids.map(lambda value: x_map.get(str(value), 0) if value else 0)
    out["comp_id"] = np.where(comp.gt(0), comp, replacement)
    return out


def summarize_buddy_df(buddy_path: Path, feature_df: pd.DataFrame, tolerance_da: float = 0.01) -> pd.DataFrame:
    path = Path(buddy_path)
    if not path.exists() or feature_df.empty or "parent_mz" not in feature_df.columns:
        return pd.DataFrame(columns=["feature_key", "lib_id"])
    try:
        buddy = read_table(path)
    except Exception:
        return pd.DataFrame(columns=["feature_key", "lib_id"])
    if "mz_feature" not in buddy.columns:
        return pd.DataFrame(columns=["feature_key", "lib_id"])
    buddy = buddy.copy()
    buddy["_mz_feature"] = pd.to_numeric(buddy["mz_feature"], errors="coerce")
    buddy = buddy[buddy["_mz_feature"].notna()].copy()
    if buddy.empty:
        return pd.DataFrame(columns=["feature_key", "lib_id"])
    rows: list[dict[str, Any]] = []
    for _, row in feature_df.iterrows():
        parent_mz = safe_float(row.get("parent_mz"), math.nan)
        if not np.isfinite(parent_mz):
            continue
        hits = buddy[(buddy["_mz_feature"] - float(parent_mz)).abs().le(float(tolerance_da))].copy()
        if hits.empty:
            continue
        hits["_abs_delta"] = (hits["_mz_feature"] - float(parent_mz)).abs()
        hits = hits.sort_values(["_abs_delta"], kind="mergesort")
        rows.append(
            {
                "feature_key": str(row.get("feature_key", "")),
                "lib_id": safe_int(row.get("lib_id"), 0),
                "buddy_match_count_v2": int(hits.shape[0]),
                "buddy_top_delta_v2": safe_float(hits.iloc[0].get("delta"), math.nan),
                "buddy_top_theor_delta_v2": safe_float(hits.iloc[0].get("theor_delta"), math.nan),
                "buddy_top_ratio_v2": safe_float(hits.iloc[0].get("ratio"), 0.0),
                "buddy_top_cov_v2": safe_float(hits.iloc[0].get("c_o_v"), 0.0),
                "buddy_top_tally_v2": safe_float(hits.iloc[0].get("tally"), 0.0),
                "buddy_inferred_m_v2": str(hits.iloc[0].get("inferred_m", "") or ""),
                "ms1_buddy_status_v2": "matched_legacy_buddy_df_by_parent_mz",
            }
        )
    return pd.DataFrame(rows)


def audit_missing_nce_reference_files(nce_dir: Path) -> pd.DataFrame:
    root = Path(nce_dir).expanduser()
    rows: list[dict[str, Any]] = []
    search_roots = [root, root / "elucidation", root / "reference", root / "data"]
    for pattern in REFERENCE_PATTERNS:
        matches: list[Path] = []
        for search_root in search_roots:
            if search_root.exists():
                matches.extend(search_root.rglob(pattern))
        matches = sorted({path for path in matches if path.is_file()})
        rows.append(
            {
                "reference_pattern": pattern,
                "status": "found" if matches else "missing",
                "n_found": int(len(matches)),
                "example_paths": "|".join(str(path) for path in matches[:8]),
                "fallback_used": "rebuilt_proxy_from_current_data" if not matches else "portable_file_available",
            }
        )
    return pd.DataFrame(rows)


def merge_v1_and_legacy_feature_evidence(
    *,
    v1_feature_df: pd.DataFrame,
    legacy_summary: pd.DataFrame,
    buddy_summary: pd.DataFrame,
    known_x: Path,
) -> pd.DataFrame:
    base = v1_feature_df.copy() if not v1_feature_df.empty else pd.DataFrame()
    if base.empty:
        base = pd.DataFrame(columns=["feature_key", "lib_id", "comp_id", "parent_mz", "primary_neutral_mass"])
    base = assign_missing_comp_ids_from_xids(base, known_x)
    if not legacy_summary.empty and "comp_id" in base.columns:
        base["_comp_id_int"] = pd.to_numeric(base["comp_id"], errors="coerce").fillna(0).astype(int)
        base = base.merge(legacy_summary, left_on="_comp_id_int", right_on="comp_id", how="left", suffixes=("", "_legacy_v2"))
        base = base.drop(columns=["_comp_id_int", "comp_id_legacy_v2"], errors="ignore")
    if not buddy_summary.empty and {"feature_key", "lib_id"}.issubset(base.columns):
        buddy = buddy_summary.copy()
        buddy["feature_key"] = buddy["feature_key"].astype(str)
        buddy["lib_id"] = pd.to_numeric(buddy["lib_id"], errors="coerce").fillna(0).astype(int)
        base["feature_key"] = base["feature_key"].astype(str)
        base["lib_id"] = pd.to_numeric(base["lib_id"], errors="coerce").fillna(0).astype(int)
        base = base.merge(buddy, on=["feature_key", "lib_id"], how="left")
    defaults: dict[str, Any] = {
        "legacy_record_count_v2": 0,
        "legacy_fragment_count_v2": 0,
        "legacy_match_count_v2": 0,
        "legacy_top_motif_v2": "",
        "legacy_motifs_json_v2": "{}",
        "legacy_max_probability_v2": 0.0,
        "legacy_mean_probability_v2": 0.0,
        "legacy_fragment_masses_json_v2": "[]",
        "legacy_fragment_bins_json_v2": "[]",
        "legacy_complementary_masses_json_v2": "[]",
        "legacy_loss_bins_json_v2": "[]",
        "legacy_top_match_name_v2": "",
        "legacy_top_match_score_v2": 0.0,
        "legacy_historical_notes_v2": "",
        "buddy_match_count_v2": 0,
        "buddy_top_ratio_v2": 0.0,
        "buddy_top_cov_v2": 0.0,
        "buddy_top_tally_v2": 0.0,
        "ms1_buddy_status_v2": "not_available_or_no_parent_mz_match",
    }
    for col, default in defaults.items():
        if col not in base.columns:
            base[col] = default
        base[col] = base[col].fillna(default)
    if "primary_neutral_mass" not in base.columns and {"parent_mz", "lib_id"}.issubset(base.columns):
        base["primary_neutral_mass"] = [
            neutral_mass_from_precursor(mz, safe_int(lib_id, 0))
            for mz, lib_id in zip(base["parent_mz"], base["lib_id"])
        ]
    base["nce_workflow_v2_evidence_status"] = np.where(
        pd.to_numeric(base["legacy_record_count_v2"], errors="coerce").fillna(0).gt(0),
        "legacy_ms2smiles_or_match_evidence",
        "current_feature_evidence_only",
    )
    return base.sort_values(["lib_id", "feature_key"], kind="mergesort").reset_index(drop=True)


def formula_fragment_consistency_score(candidate_formula: Any, fragment_masses_json: Any, complementary_masses_json: Any) -> float:
    cand_mass = formula_mass(candidate_formula)
    if not np.isfinite(cand_mass) or cand_mass <= 0:
        return 0.5
    frag_masses = [safe_float(value, math.nan) for value in parse_json_list(fragment_masses_json)]
    comp_masses = [safe_float(value, math.nan) for value in parse_json_list(complementary_masses_json)]
    frag_masses = [value for value in frag_masses if np.isfinite(value) and value > 0]
    comp_masses = [value for value in comp_masses if np.isfinite(value) and value > 0]
    if not frag_masses and not comp_masses:
        return 0.5
    scores: list[float] = []
    for value in frag_masses:
        scores.append(1.0 if value < cand_mass - 1.0 else 0.1)
    for value in comp_masses:
        scores.append(1.0 if value < cand_mass - 1.0 else 0.1)
    for frag in frag_masses:
        for comp in comp_masses:
            err = abs((frag + comp) - cand_mass)
            if err <= 0.01:
                scores.append(1.0)
            elif err <= 0.1:
                scores.append(0.9)
            elif err <= 1.0:
                scores.append(0.7)
            elif err <= 5.0:
                scores.append(0.35)
            else:
                scores.append(0.0)
    return float(np.clip(np.mean(scores), 0.0, 1.0)) if scores else 0.5


def element_formula_overlap(left_formula: Any, right_formula: Any) -> float:
    left = formula_counts(left_formula)
    right = formula_counts(right_formula)
    if not left or not right:
        return 0.5
    keys = set(left) | set(right)
    distance = sum(abs(int(left.get(k, 0)) - int(right.get(k, 0))) for k in keys)
    size = sum(max(int(left.get(k, 0)), int(right.get(k, 0))) for k in keys)
    return float(np.clip(1.0 - distance / max(size, 1), 0.0, 1.0))


def legacy_motif_candidate_score(candidate_labels: Iterable[str], motifs_json: Any, top_motif: Any = "") -> float:
    labels = {str(label) for label in candidate_labels if str(label)}
    motifs: set[str] = set()
    try:
        parsed = json.loads(str(motifs_json or "{}"))
        if isinstance(parsed, dict):
            motifs.update(str(key) for key, value in parsed.items() if safe_float(value, 0.0) > 0)
    except Exception:
        pass
    top = str(top_motif or "").strip()
    if top:
        motifs.add(top)
    if not labels and not motifs:
        return 0.5
    if not labels or not motifs:
        return 0.25
    return float(len(labels & motifs) / max(len(labels | motifs), 1))


def likely_next_operator_from_row(row: pd.Series) -> str:
    bucket = str(row.get("failure_label", "") or "")
    motif_text = " ".join(
        str(row.get(col, "") or "")
        for col in ["legacy_top_motif_v2", "legacy_motifs_json_v2", "legacy_historical_notes_v2"]
    ).lower()
    top10 = safe_float(row.get("top10_best_tanimoto"), 0.0)
    oracle = safe_float(row.get("full_pool_oracle_best_tanimoto"), 0.0)
    if top10 >= 0.7:
        return "validate_current_top10"
    if oracle >= 0.7 and top10 < 0.7:
        return "rerank_only_oracle_candidate"
    if "sulfate" in motif_text or "phosphate" in motif_text:
        return "targeted_conjugation_sulfate_phosphate"
    if "glucuron" in motif_text or "glycoside" in motif_text:
        return "glycoside_glucuronide_local_assembly"
    if "taurine" in motif_text or "carnitine" in motif_text:
        return "acyl_conjugate_operator"
    if "glycine" in motif_text or "glutathione" in motif_text or "peptide" in motif_text:
        return "peptide_glutathione_fragment_assembly"
    if bucket in {"no_candidate", "weak_oracle", "no_mass_compatible_novel"}:
        return "operator_redesign_for_candidate_pool"
    return "collect_ms1_ms2_reference_data"
