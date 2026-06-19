#!/usr/bin/env python3
"""Utilities for de novo generation improvement experiments.

Pipeline context:
    Builds formula feature maps and structure keys for de novo candidate
    tables. Re-exports core chemistry helpers so experiment scripts import
    from one module.

Key I/O:
    ``load_formula_feature_map``: comp_id-keyed formula vectors from CSV.
    ``add_structure_keys``: canonical SMILES and InChIKey columns for joins.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from structural_elucidation.utils.chemistry import (
    FORMULA_ELEMENTS,
    FORMULA_FEATURE_COLUMNS,
    FORMULA_FEATURE_SCALE,
    HAS_RDKIT,
    canonicalize_smiles,
    exact_mass_from_smiles,
    formula_counts,
    formula_from_smiles,
    formula_vector_from_counts,
    formula_vector_from_formula,
    formula_vector_from_smiles,
    inchikey_from_smiles,
    mol_from_smiles,
    zero_formula_vector,
)
from structural_elucidation.utils.io import pick_col, read_table

# Re-export chemistry symbols for backward-compatible denovo imports.
__all__ = [
    "FORMULA_ELEMENTS",
    "FORMULA_FEATURE_COLUMNS",
    "FORMULA_FEATURE_SCALE",
    "HAS_RDKIT",
    "add_structure_keys",
    "adduct_vector",
    "canonicalize_smiles",
    "exact_mass_from_smiles",
    "formula_counts",
    "formula_from_smiles",
    "formula_vector_from_counts",
    "formula_vector_from_formula",
    "formula_vector_from_smiles",
    "inchikey_from_smiles",
    "load_formula_feature_map",
    "mol_from_smiles",
    "row_formula",
    "zero_formula_vector",
]


# ── Ion-mode encoding ────────────────────────────────────────────────
def adduct_vector(ion_mode: str) -> list[float]:
    """One-hot encode ion mode for model features: ``[neg, pos]``.

    Args:
        ion_mode: ``negative``/``neg`` or ``positive``/``pos`` prefix.

    Returns:
        Two-element float list; ``[0, 0]`` when mode is unrecognized.
    """
    mode = str(ion_mode or "").strip().lower()
    if mode.startswith("neg"):
        return [1.0, 0.0]
    if mode.startswith("pos"):
        return [0.0, 1.0]
    return [0.0, 0.0]


# ── Formula resolution from heterogeneous tables ─────────────────────
def row_formula(row: pd.Series) -> str:
    """Resolve a molecular formula from common column names or SMILES fallback.

    Args:
        row: Candidate or feature metadata row.

    Returns:
        Formula string, or ``""`` when neither explicit formula nor SMILES resolves.
    """
    for col in ["formula", "Formula", "predicted_formula", "candidate_formula"]:
        if col in row.index:
            value = str(row.get(col, "") or "").strip()
            if value and value.lower() != "nan":
                return value
    for col in ["canonical_smiles", "smiles", "SMILES", "predicted_smiles"]:
        if col in row.index:
            formula = formula_from_smiles(row.get(col))
            if formula:
                return formula
    return ""


def load_formula_feature_map(path: Path, key_candidates: list[str]) -> dict[str, dict[str, Any]]:
    """Load per-key formula counts and ML vectors from a CSV export.

    Args:
        path: Formula feature CSV path.
        key_candidates: Column names to try as join keys (first match wins).

    Returns:
        Mapping from key string to ``{formula, counts, vector, source}``.

    Raises:
        KeyError: When no ``key_candidates`` column exists in the file.
    """
    if path is None or not Path(path).exists():
        return {}
    df = read_table(Path(path))
    key_col = pick_col(df, key_candidates)
    if key_col is None:
        raise KeyError(f"Could not find any key column {key_candidates} in formula CSV {path}")
    out: dict[str, dict[str, Any]] = {}
    for _, row in df.iterrows():
        key = str(row.get(key_col, "")).strip()
        if not key or key.lower() == "nan":
            continue
        formula = row_formula(row)
        if not formula:
            continue
        counts = formula_counts(formula)
        out[key] = {
            "formula": formula,
            "counts": counts,
            "vector": formula_vector_from_counts(counts),
            "source": str(path),
        }
    return out


# ── Structure keys for deduplication joins ───────────────────────────
def add_structure_keys(df: pd.DataFrame, smiles_col: str) -> pd.DataFrame:
    """Append ``canonical_smiles_key`` and ``inchikey_key`` columns.

    Args:
        df: Input table containing ``smiles_col``.
        smiles_col: Column with raw SMILES strings.

    Returns:
        Copy of ``df`` with structure key columns (empty strings when column missing).
    """
    out = df.copy()
    if smiles_col not in out.columns:
        out["canonical_smiles_key"] = ""
        out["inchikey_key"] = ""
        return out
    out["canonical_smiles_key"] = out[smiles_col].map(lambda value: canonicalize_smiles(value) or "")
    out["inchikey_key"] = out["canonical_smiles_key"].map(inchikey_from_smiles)
    return out
