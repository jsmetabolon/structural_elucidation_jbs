#!/usr/bin/env python3
"""Utilities for de novo generation improvement experiments."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors

    HAS_RDKIT = True
except ImportError:  # pragma: no cover - runtime dependency
    Chem = None
    Descriptors = None
    rdMolDescriptors = None
    HAS_RDKIT = False


FORMULA_ELEMENTS = ["C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"]
FORMULA_FEATURE_COLUMNS = [f"formula_count_{element}" for element in FORMULA_ELEMENTS]
FORMULA_FEATURE_SCALE = 64.0


def read_table(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception:
        df = pd.read_csv(path)
    df.columns = [str(col).strip().strip('"').lstrip("\ufeff") for col in df.columns]
    return df


def pick_col(df: pd.DataFrame, options: list[str]) -> str | None:
    lookup = {str(col).strip().lower(): str(col) for col in df.columns}
    for option in options:
        found = lookup.get(str(option).strip().lower())
        if found is not None:
            return found
    return None


def canonicalize_smiles(smiles: Any) -> str | None:
    if not HAS_RDKIT:
        return None
    text = str(smiles or "").strip()
    if not text or text.lower() == "nan":
        return None
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return None


def mol_from_smiles(smiles: Any):
    if not HAS_RDKIT:
        return None
    canonical = canonicalize_smiles(smiles)
    if canonical is None:
        return None
    return Chem.MolFromSmiles(canonical)


def inchikey_from_smiles(smiles: Any) -> str:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return ""
    try:
        return str(Chem.MolToInchiKey(mol))
    except Exception:
        return ""


def formula_from_smiles(smiles: Any) -> str:
    mol = mol_from_smiles(smiles)
    if mol is None or rdMolDescriptors is None:
        return ""
    try:
        return str(rdMolDescriptors.CalcMolFormula(mol))
    except Exception:
        return ""


def exact_mass_from_smiles(smiles: Any) -> float:
    mol = mol_from_smiles(smiles)
    if mol is None or Descriptors is None:
        return float("nan")
    try:
        return float(Descriptors.ExactMolWt(mol))
    except Exception:
        return float("nan")


_FORMULA_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def formula_counts(formula: Any) -> dict[str, int]:
    text = str(formula or "").strip()
    if not text or text.lower() == "nan":
        return {}
    counts: dict[str, int] = {}
    for element, raw_count in _FORMULA_RE.findall(text):
        if element not in FORMULA_ELEMENTS:
            continue
        count = int(raw_count) if raw_count else 1
        counts[element] = counts.get(element, 0) + count
    return counts


def formula_vector_from_counts(counts: dict[str, int]) -> list[float]:
    return [
        float(np.clip(float(counts.get(element, 0)) / FORMULA_FEATURE_SCALE, 0.0, 1.0))
        for element in FORMULA_ELEMENTS
    ]


def formula_vector_from_formula(formula: Any) -> list[float]:
    return formula_vector_from_counts(formula_counts(formula))


def formula_vector_from_smiles(smiles: Any) -> list[float]:
    return formula_vector_from_formula(formula_from_smiles(smiles))


def zero_formula_vector() -> list[float]:
    return [0.0 for _ in FORMULA_ELEMENTS]


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
        raise KeyError(
            f"Could not find any key column {key_candidates} in formula CSV {path}"
        )
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
