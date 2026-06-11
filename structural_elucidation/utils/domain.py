#!/usr/bin/env python3
"""Shared helpers for X-domain hard holdout and analog diagnostics."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from structural_elucidation.utils.chemistry import canonicalize_smiles, inchikey_from_smiles
from structural_elucidation.utils.io import dataframe_to_markdown, read_table, write_csv, write_json

NEGATIVE_LIBS = {209, 305}
POSITIVE_LIBS = {400, 402}
ION_LABELS = {"negative": 0.0, "positive": 1.0}
DEFAULT_HARD_LIMITS = {209: 58, 305: 12, 400: 34, 402: 7}


def lib_mode(lib_id: int) -> str:
    return "negative" if int(lib_id) in NEGATIVE_LIBS else "positive"


def parse_libs(text: str) -> list[int]:
    return [int(part.strip()) for part in str(text).split(",") if part.strip()]


def parse_limits(text: str) -> dict[int, int]:
    out = dict(DEFAULT_HARD_LIMITS)
    if not str(text or "").strip():
        return out
    for item in str(text).split(","):
        if not item.strip():
            continue
        lib_text, value_text = item.split(":", 1)
        out[int(lib_text)] = int(value_text)
    return out


def safe_float(value: Any, default: float = math.nan) -> float:
    try:
        out = float(value)
        return out if np.isfinite(out) else default
    except Exception:
        return default


def safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(value))
    except Exception:
        return default


def xid_from_feature_key(value: Any) -> str:
    match = re.search(r"(X-\d+)", str(value or ""))
    return match.group(1) if match else ""


UNKNOWN_LABEL_TO_XNUM = re.compile(r"Unknown\s*-\s*(\d+)", re.I)
UNKNOWN_LABEL_SUFFIX = re.compile(r"Unknown\s*-\s*\d+\s*[-_]", re.I)
NMD_NAME_PREFIX = re.compile(r"^NMD", re.I)


def unknown_label_to_xid(value: Any) -> str:
    match = UNKNOWN_LABEL_TO_XNUM.search(str(value or "").strip())
    return f"X-{match.group(1)}" if match else ""


def is_nmd_compound_name(compound_name: Any, report_name: Any = "") -> bool:
    for raw in (compound_name, report_name):
        text = str(raw or "").strip()
        if text and NMD_NAME_PREFIX.match(text):
            return True
    return False


def resolve_canonical_control_id(compound_name: Any, report_name: Any = "") -> tuple[str, str, str]:
    """Return (canonical_control_id, canonical_compound_name, name_normalization)."""
    original = str(compound_name or "").strip()
    for raw in (compound_name, report_name):
        xid = unknown_label_to_xid(raw)
        if xid:
            norm = "unknown_to_x_suffix" if UNKNOWN_LABEL_SUFFIX.search(str(raw or "")) else "unknown_to_x"
            return xid, xid, norm
    for raw in (compound_name, report_name):
        xid = xid_from_feature_key(raw)
        if xid:
            return xid, xid, "already_x"
    if original:
        return "", original, "named"
    return "", "", "empty"


def normalize_ion_mode(value: Any) -> str:
    text = str(value or "").strip().lower()
    if text in {"0", "0.0", "negative", "neg"}:
        return "negative"
    if text in {"1", "1.0", "positive", "pos"}:
        return "positive"
    return text


def has_lib(row: pd.Series, lib_id: int) -> bool:
    raw = str(row.get("lib_ids", "") or "")
    parts = {part.strip() for part in raw.replace(",", ";").split(";") if part.strip()}
    if str(int(lib_id)) in parts:
        return True
    return safe_int(row.get("parent_mz_source_lib"), -1) == int(lib_id)


def feature_summary(df: pd.DataFrame, *, default_scope: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    work = df.copy()
    if "feature_key" not in work.columns:
        raise ValueError("feature table must include feature_key")
    if "comp_id" not in work.columns:
        work["comp_id"] = np.nan
    if "lib_id" not in work.columns:
        work["lib_id"] = np.nan
    if "ion_mode_label" in work.columns:
        work["_ion_mode"] = work["ion_mode_label"].map(normalize_ion_mode)
    elif "ion_mode" in work.columns:
        work["_ion_mode"] = work["ion_mode"].map(normalize_ion_mode)
    else:
        work["_ion_mode"] = ""
    if "ri" not in work.columns and "rt" in work.columns:
        work["ri"] = work["rt"]
    if "parent_mz" not in work.columns:
        work["parent_mz"] = np.nan
    grouped = (
        work.groupby("feature_key", dropna=False)
        .agg(
            comp_id=("comp_id", "first"),
            lib_id=("lib_id", "first"),
            ion_mode=("_ion_mode", "first"),
            parent_mz=("parent_mz", "first"),
            ri=("ri", "first") if "ri" in work.columns else ("parent_mz", "first"),
            rt=("rt", "first") if "rt" in work.columns else ("parent_mz", "first"),
            peak_count=("mass", "count") if "mass" in work.columns else ("feature_key", "count"),
        )
        .reset_index()
    )
    grouped["feature_key"] = grouped["feature_key"].astype(str)
    grouped["xid"] = grouped["feature_key"].map(xid_from_feature_key)
    grouped["holdout_scope"] = str(default_scope)
    grouped["comp_id"] = pd.to_numeric(grouped["comp_id"], errors="coerce")
    grouped["lib_id"] = pd.to_numeric(grouped["lib_id"], errors="coerce")
    grouped["parent_mz"] = pd.to_numeric(grouped["parent_mz"], errors="coerce")
    grouped["ri"] = pd.to_numeric(grouped["ri"], errors="coerce")
    grouped["rt"] = pd.to_numeric(grouped["rt"], errors="coerce")
    grouped["peak_count"] = pd.to_numeric(grouped["peak_count"], errors="coerce")
    return grouped


def load_known_x_comp_ids(known_x: Path) -> set[int]:
    if not known_x.exists():
        return set()
    df = read_table(known_x)
    if "CompID" not in df.columns:
        return set()
    return {int(value) for value in pd.to_numeric(df["CompID"], errors="coerce").dropna().astype(int)}


def load_x_structure_keys(known_x_all: Path) -> tuple[set[str], set[str]]:
    if not known_x_all.exists():
        return set(), set()
    df = read_table(known_x_all)
    smiles_col = pick_col(df, ["smiles", "SMILES", "canonical_smiles"])
    if smiles_col is None:
        return set(), set()
    keyed = add_structure_keys(df, smiles_col)
    return (
        {value for value in keyed["canonical_smiles_key"].astype(str) if value},
        {value for value in keyed["inchikey_key"].astype(str) if value},
    )


def load_existing_holdout_keys(holdout_dir: Path) -> tuple[set[int], set[str], set[str]]:
    comp_ids: set[int] = set()
    canonicals: set[str] = set()
    inchikeys: set[str] = set()
    if not holdout_dir.exists():
        return comp_ids, canonicals, inchikeys
    for path in sorted(holdout_dir.glob("*.csv")):
        try:
            df = read_table(path)
        except Exception:
            continue
        if "comp_id" in df.columns:
            comp_ids.update(int(value) for value in pd.to_numeric(df["comp_id"], errors="coerce").dropna().astype(int))
        smiles_col = pick_col(df, ["canonical_smiles", "smiles", "SMILES"])
        if smiles_col is not None:
            keyed = add_structure_keys(df, smiles_col)
            canonicals.update(value for value in keyed["canonical_smiles_key"].astype(str) if value)
            inchikeys.update(value for value in keyed["inchikey_key"].astype(str) if value)
        if "inchikey" in df.columns:
            inchikeys.update(value for value in df["inchikey"].astype(str) if value and value != "nan")
    return comp_ids, canonicals, inchikeys


def load_train_exclude_comp_ids(holdout_dir: Path) -> set[int]:
    out: set[int] = set()
    if not holdout_dir.exists():
        return out
    for path in sorted(holdout_dir.glob("train_exclude_*.csv")):
        try:
            df = read_table(path)
        except Exception:
            continue
        if "comp_id" in df.columns:
            out.update(int(value) for value in pd.to_numeric(df["comp_id"], errors="coerce").dropna().astype(int))
    return out


def load_training_structures(repo_root: Path, libs: Iterable[int]) -> pd.DataFrame:
    paths = {
        "negative": repo_root / "merged_encoder/data/decoder_train_neg_209_305_smiles.csv",
        "positive": repo_root / "merged_encoder/data/decoder_train_pos_400_402_smiles.csv",
    }
    rows: list[pd.DataFrame] = []
    for mode, path in paths.items():
        if not path.exists():
            continue
        df = read_table(path)
        comp_col = pick_col(df, ["comp_id", "CompID", "id"])
        smiles_col = pick_col(df, ["canonical_smiles", "smiles", "SMILES"])
        if comp_col is None or smiles_col is None:
            continue
        out = df.copy()
        out["comp_id"] = pd.to_numeric(out[comp_col], errors="coerce").astype("Int64")
        out["canonical_smiles"] = out[smiles_col].astype(str)
        out = out.dropna(subset=["comp_id"]).copy()
        out = add_structure_keys(out, "canonical_smiles")
        out["formula"] = out["canonical_smiles_key"].map(formula_from_smiles)
        out["asset_mode"] = mode
        out["asset_path"] = str(path)
        rows.append(out)
    if not rows:
        return pd.DataFrame()
    combined = pd.concat(rows, ignore_index=True)
    keep_rows: list[pd.DataFrame] = []
    for lib_id in libs:
        mode = lib_mode(int(lib_id))
        scoped = combined[
            combined["asset_mode"].eq(mode) & combined.apply(lambda row: has_lib(row, int(lib_id)), axis=1)
        ].copy()
        scoped["lib_id"] = int(lib_id)
        scoped["ion_mode"] = mode
        keep_rows.append(scoped)
    return pd.concat(keep_rows, ignore_index=True) if keep_rows else pd.DataFrame()


def domain_distance(
    *,
    target_parent_mz: float,
    candidate_parent_mz: float,
    target_ri: float,
    candidate_ri: float,
    target_peak_count: float,
    candidate_peak_count: float,
    target_latent_distance: float = math.nan,
    candidate_latent_distance: float = math.nan,
) -> dict[str, float]:
    mz_delta = abs(safe_float(candidate_parent_mz, 0.0) - safe_float(target_parent_mz, 0.0)) / 50.0
    ri_delta = abs(safe_float(candidate_ri, 0.0) - safe_float(target_ri, 0.0)) / 1000.0
    peak_delta = abs(
        math.log1p(max(safe_float(candidate_peak_count, 0.0), 0.0))
        - math.log1p(max(safe_float(target_peak_count, 0.0), 0.0))
    )
    latent_delta = 0.0
    if np.isfinite(safe_float(target_latent_distance)) and np.isfinite(safe_float(candidate_latent_distance)):
        latent_delta = abs(float(candidate_latent_distance) - float(target_latent_distance))
    score = mz_delta + ri_delta + 0.5 * peak_delta + 0.25 * latent_delta
    return {
        "domain_distance": float(score),
        "mz_distance_component": float(mz_delta),
        "ri_distance_component": float(ri_delta),
        "peak_count_distance_component": float(peak_delta),
        "latent_distance_component": float(latent_delta),
    }


def canonical_key(smiles: Any) -> str:
    return canonicalize_smiles(smiles) or ""


def inchikey_key(smiles: Any) -> str:
    canonical = canonical_key(smiles)
    return inchikey_from_smiles(canonical) if canonical else ""
