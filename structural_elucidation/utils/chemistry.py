"""Shared chemistry helpers for SMILES normalization and formula features.

Pipeline context:
    Canonical SMILES, InChIKey, and elemental formula vectors are reused by
    de novo scoring, domain holdout keys, merge policies, and enrichment
    joins. RDKit is optional; callers receive empty or NaN sentinels when
    cheminformatics are unavailable.

Key I/O:
    ``FORMULA_FEATURE_COLUMNS``: fixed-order count features for ML models.
    ``canonicalize_smiles`` / ``formula_vector_from_smiles``: structure keys.
"""

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

# ── Formula feature schema ───────────────────────────────────────────
FORMULA_ELEMENTS = ["C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"]
FORMULA_FEATURE_COLUMNS = [f"formula_count_{element}" for element in FORMULA_ELEMENTS]
FORMULA_FEATURE_SCALE = 64.0


# ── SMILES normalization ─────────────────────────────────────────────
def canonicalize_smiles(smiles: Any) -> str | None:
    """Return isomeric canonical SMILES, or ``None`` when invalid or RDKit is absent.

    Args:
        smiles: Raw SMILES string or scalar coercible to string.

    Returns:
        Canonical SMILES, or ``None`` on empty, NaN, parse, or sanitize failure.
    """
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


def mol_from_smiles(smiles: Any) -> Any | None:
    """Build a sanitized RDKit molecule from ``smiles`` via canonicalization.

    Args:
        smiles: Raw SMILES string.

    Returns:
        RDKit ``Mol``, or ``None`` when RDKit is unavailable or SMILES is invalid.
    """
    if not HAS_RDKIT:
        return None
    canonical = canonicalize_smiles(smiles)
    if canonical is None:
        return None
    return Chem.MolFromSmiles(canonical)


def inchikey_from_smiles(smiles: Any) -> str:
    """Derive InChIKey from ``smiles``; returns ``""`` on failure.

    Args:
        smiles: Raw or canonical SMILES.

    Returns:
        InChIKey string, or empty when structure cannot be resolved.
    """
    mol = mol_from_smiles(smiles)
    if mol is None:
        return ""
    try:
        return str(Chem.MolToInchiKey(mol))
    except Exception:
        return ""


def formula_from_smiles(smiles: Any) -> str:
    """Return molecular formula for ``smiles``; empty string on failure.

    Args:
        smiles: Raw or canonical SMILES.

    Returns:
        Hill-system formula string, or ``""`` when unavailable.
    """
    mol = mol_from_smiles(smiles)
    if mol is None or rdMolDescriptors is None:
        return ""
    try:
        return str(rdMolDescriptors.CalcMolFormula(mol))
    except Exception:
        return ""


def exact_mass_from_smiles(smiles: Any) -> float:
    """Return monoisotopic exact mass for ``smiles``.

    Args:
        smiles: Raw or canonical SMILES.

    Returns:
        Exact molecular weight in Da, or ``nan`` when unavailable.
    """
    mol = mol_from_smiles(smiles)
    if mol is None or Descriptors is None:
        return float("nan")
    try:
        return float(Descriptors.ExactMolWt(mol))
    except Exception:
        return float("nan")


# ── Formula parsing and ML feature vectors ───────────────────────────
_FORMULA_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def formula_counts(formula: Any) -> dict[str, int]:
    """Parse elemental counts from a formula string for supported elements.

    Only ``FORMULA_ELEMENTS`` are retained so feature vectors stay aligned
    with training schema.

    Args:
        formula: Hill-system or similar formula text.

    Returns:
        Element symbol to integer count mapping (may be empty).
    """
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
    """Scale elemental counts into clipped [0, 1] features for ML models.

    Args:
        counts: Element symbol to count mapping.

    Returns:
        Fixed-order vector aligned with ``FORMULA_ELEMENTS``.
    """
    return [
        float(np.clip(float(counts.get(element, 0)) / FORMULA_FEATURE_SCALE, 0.0, 1.0)) for element in FORMULA_ELEMENTS
    ]


def formula_vector_from_formula(formula: Any) -> list[float]:
    """Build a formula feature vector directly from a formula string.

    Args:
        formula: Molecular formula text.

    Returns:
        Fixed-order scaled count vector.
    """
    return formula_vector_from_counts(formula_counts(formula))


def formula_vector_from_smiles(smiles: Any) -> list[float]:
    """Build a formula feature vector from ``smiles`` via RDKit formula derivation.

    Args:
        smiles: Raw or canonical SMILES.

    Returns:
        Fixed-order scaled count vector (zeros when formula cannot be resolved).
    """
    return formula_vector_from_formula(formula_from_smiles(smiles))


def zero_formula_vector() -> list[float]:
    """Return an all-zero formula feature vector matching ``FORMULA_ELEMENTS`` order."""
    return [0.0 for _ in FORMULA_ELEMENTS]
