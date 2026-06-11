#!/usr/bin/env python3
"""Unified readers for the local WGCNA cohort exports."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd

from structural_elucidation.paths import REPO_ROOT

CACHE_ROOT = REPO_ROOT / "data" / "wgcna_cache"

COHORT_INPUTS: dict[str, dict[str, Any]] = {
    "HELIOS": {
        "module": CACHE_ROOT / "HELIOS" / "WGCNA_HELIOS_all_6_0.25_3.csv",
        "corr": CACHE_ROOT / "HELIOS" / "full_correlation_matrix.csv",
        "putative": CACHE_ROOT / "HELIOS" / "unk_info.csv",
        "module_col": "MODULE",
    },
    "Twins_UK": {
        "module": CACHE_ROOT / "Twins_UK" / "WGCNA_Twins_UK_2420_5_0.25_3.csv",
        "corr": CACHE_ROOT / "Twins_UK" / "corrMatrix.csv",
        "putative": CACHE_ROOT / "Twins_UK" / "network_modules_with_details.csv",
        "module_col": "MODULE",
    },
    "PRIMM_feces": {
        "module": CACHE_ROOT / "PRIMM_feces" / "WGCNA_PRIMMfull_feces_23920_5_0.25_3.xlsx",
        "corr": CACHE_ROOT / "PRIMM_feces" / "corrTN_PRIMMfull_feces_23920_0.3_0.03_0.01.xlsx",
        "putative": CACHE_ROOT / "PRIMM_feces" / "WGCNA_PRIMMfull_feces_23920_5_0.25_3.xlsx",
        "module_col": "MODULE",
    },
    "PRIMM_serum": {
        "module": CACHE_ROOT / "PRIMM_serum" / "WGCNA_PRIMMfull_serum_23629_5_0.25_3_withUnkInfo.xlsx",
        "corr": CACHE_ROOT / "PRIMM_serum" / "corrTN_PRIMMfull_serum_23629_0.3_0.03_0.01_withUnkInfo.xlsx",
        "putative": CACHE_ROOT / "PRIMM_serum" / "WGCNA_PRIMMfull_serum_23629_5_0.25_3_withUnkInfo.xlsx",
        "module_col": "MODULE",
    },
}

NORMALIZED_MODULE_COLUMNS = [
    "cohort",
    "chem_id",
    "chem_id2",
    "module",
    "type",
    "super_pathway",
    "sub_pathway",
    "inchikey",
    "smiles",
    "plot_name",
    "lib_id",
    "comp_id",
]


def _read_table(path: Path, **kwargs: Any) -> pd.DataFrame:
    path = Path(path).expanduser()
    if not path.exists():
        return pd.DataFrame()
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path, **kwargs)
    return pd.read_csv(path, low_memory=False, **kwargs)


def _first_present(row: pd.Series, names: list[str], default: Any = "") -> Any:
    for name in names:
        if name in row.index and pd.notna(row.get(name)):
            return row.get(name)
    return default


def normalize_control_id(value: Any) -> str:
    text = str(value or "").strip()
    match = re.search(r"\bX[-_ ]?0*(\d+)\b", text, flags=re.IGNORECASE)
    if match:
        return f"X-{int(match.group(1))}"
    return ""


def normalize_corr_label(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if text.startswith("X-"):
        return f"X{int(text.split('-', 1)[1])}" if text.split("-", 1)[1].isdigit() else text.replace("-", "")
    if text.upper().startswith("X") and text[1:].isdigit():
        return "X" + str(int(text[1:]))
    try:
        return "X" + str(int(float(text)))
    except Exception:
        return text


def _normalize_type(row: pd.Series) -> str:
    raw = str(row.get("TYPE", "") or "").strip().upper()
    if raw:
        return raw
    control_id = normalize_control_id(
        _first_present(row, ["PLOT_NAME", "CHEMICAL_NAME", "COMMENTS", "PUTATIVE_DETAILS"])
    )
    return "UNNAMED" if control_id else "NAMED"


def load_module_assignments(cohort: str) -> pd.DataFrame:
    """Load one WGCNA cohort module table with stable normalized columns."""
    if cohort not in COHORT_INPUTS:
        raise KeyError(f"Unknown WGCNA cohort: {cohort}")
    spec = COHORT_INPUTS[cohort]
    raw = _read_table(spec["module"])
    rows: list[dict[str, Any]] = []
    module_col = str(spec.get("module_col", "MODULE"))
    for _, row in raw.iterrows():
        module = _first_present(row, [module_col, "MODULE", "MODULE1", "MODULE_ID"], 0)
        try:
            if int(float(module)) == 0:
                continue
        except Exception:
            pass
        plot_name = _first_present(row, ["PLOT_NAME", "CHEMICAL_NAME", "COMMENTS"], "")
        control_id = normalize_control_id(plot_name)
        if control_id:
            plot_name = control_id
        chem_id = _first_present(row, ["CHEM_ID"], "")
        rows.append(
            {
                "cohort": cohort,
                "chem_id": str(chem_id or ""),
                "chem_id2": str(_first_present(row, ["CHEM_ID2"], chem_id) or ""),
                "module": str(module),
                "type": _normalize_type(row),
                "super_pathway": str(_first_present(row, ["SUPER_PATHWAY"], "") or ""),
                "sub_pathway": str(_first_present(row, ["SUB_PATHWAY"], "") or ""),
                "inchikey": str(_first_present(row, ["INCHIKEY"], "") or ""),
                "smiles": str(_first_present(row, ["SMILES"], "") or ""),
                "plot_name": str(plot_name or ""),
                "lib_id": _first_present(row, ["LIB_ID"], ""),
                "comp_id": _first_present(row, ["COMP_ID"], ""),
            }
        )
    out = pd.DataFrame(rows, columns=NORMALIZED_MODULE_COLUMNS)
    if not out.empty:
        out["lib_id"] = pd.to_numeric(out["lib_id"], errors="coerce").astype("Int64")
        out["comp_id"] = pd.to_numeric(out["comp_id"], errors="coerce").astype("Int64")
    return out


def _write_cache(df: pd.DataFrame, cache_path: Path) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        df.to_parquet(cache_path)
    except Exception:
        df.to_pickle(cache_path)


def _read_cache(cache_path: Path) -> pd.DataFrame:
    try:
        return pd.read_parquet(cache_path)
    except Exception:
        return pd.read_pickle(cache_path)


def _load_raw_corr(cohort: str) -> pd.DataFrame:
    raw = _read_table(COHORT_INPUTS[cohort]["corr"])
    if raw.empty:
        return raw
    first_col = str(raw.columns[0])
    if first_col.lower().startswith("unnamed") or first_col in {"", "row", "index"}:
        raw = raw.set_index(first_col)
    elif raw.shape[0] == raw.shape[1] - 1:
        raw = raw.set_index(first_col)
    if raw.index.dtype == object or raw.index.name is not None:
        raw.index = [normalize_corr_label(value) for value in raw.index]
    raw.columns = [normalize_corr_label(value) for value in raw.columns]
    return raw


def load_correlation_matrix(cohort: str) -> pd.DataFrame:
    """Load and cache a cohort correlation matrix/table."""
    if cohort not in COHORT_INPUTS:
        raise KeyError(f"Unknown WGCNA cohort: {cohort}")
    cache_path = CACHE_ROOT / f"{cohort}_corr.parquet"
    if cache_path.exists():
        return _read_cache(cache_path)
    corr = _load_raw_corr(cohort)
    _write_cache(corr, cache_path)
    return corr


def load_putative_details(cohort: str) -> pd.DataFrame:
    """Return putative free-text rows with normalized identifiers."""
    if cohort not in COHORT_INPUTS:
        raise KeyError(f"Unknown WGCNA cohort: {cohort}")
    raw = _read_table(COHORT_INPUTS[cohort]["putative"])
    rows: list[dict[str, Any]] = []
    for _, row in raw.iterrows():
        text = _first_present(row, ["PUTATIVE_DETAILS", "putative_details"], "")
        if pd.isna(text) or not str(text).strip():
            continue
        plot_name = _first_present(row, ["PLOT_NAME", "CHEMICAL_NAME", "COMMENTS"], "")
        control_id = normalize_control_id(plot_name) or normalize_control_id(text)
        if control_id:
            plot_name = control_id
        rows.append(
            {
                "cohort": cohort,
                "chro_lib_entry_id": _first_present(row, ["CHRO_LIB_ENTRY_ID"], ""),
                "chem_id": str(_first_present(row, ["CHEM_ID"], "") or ""),
                "plot_name": str(plot_name or ""),
                "control_id": control_id,
                "putative_details": str(text),
            }
        )
    return pd.DataFrame(
        rows, columns=["cohort", "chro_lib_entry_id", "chem_id", "plot_name", "control_id", "putative_details"]
    )


def all_cohorts() -> list[str]:
    return list(COHORT_INPUTS)


if __name__ == "__main__":
    for cohort_name in all_cohorts():
        modules = load_module_assignments(cohort_name)
        corr = load_correlation_matrix(cohort_name)
        putative = load_putative_details(cohort_name)
        print(f"{cohort_name}: modules={modules.shape} corr={corr.shape} putative={putative.shape}")
