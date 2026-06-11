#!/usr/bin/env python3
"""Utilities for de novo generation improvement experiments."""

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


def adduct_vector(ion_mode: str) -> list[float]:
    mode = str(ion_mode or "").strip().lower()
    if mode.startswith("neg"):
        return [1.0, 0.0]
    if mode.startswith("pos"):
        return [0.0, 1.0]
    return [0.0, 0.0]


def row_formula(row: pd.Series) -> str:
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


def add_structure_keys(df: pd.DataFrame, smiles_col: str) -> pd.DataFrame:
    out = df.copy()
    if smiles_col not in out.columns:
        out["canonical_smiles_key"] = ""
        out["inchikey_key"] = ""
        return out
    out["canonical_smiles_key"] = out[smiles_col].map(lambda value: canonicalize_smiles(value) or "")
    out["inchikey_key"] = out["canonical_smiles_key"].map(inchikey_from_smiles)
    return out
