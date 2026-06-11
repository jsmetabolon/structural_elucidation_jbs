"""Shared chemistry helpers."""

from __future__ import annotations

import re
from typing import Any

import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors

    HAS_RDKIT = True
except ImportError:  # pragma: no cover
    Chem = None
    Descriptors = None
    rdMolDescriptors = None
    HAS_RDKIT = False


FORMULA_ELEMENTS = ["C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"]
FORMULA_FEATURE_COLUMNS = [f"formula_count_{element}" for element in FORMULA_ELEMENTS]
FORMULA_FEATURE_SCALE = 64.0


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
