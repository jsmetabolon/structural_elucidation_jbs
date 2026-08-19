#!/usr/bin/env python3
"""Canonical ``evidence_confidence`` definition, and backfill for pooled rows.

``evidence_confidence`` is a reranker feature (``pipeline/rerank.py``) and an F2
promotion-gate input (``chemistry/diagnostics.py``), so rows from every
generator have to carry it on one scale. Generators that never emit the column
reach the pool as 0.0 through the pool's missing-column fill, which is
indistinguishable from a measured zero.

The definition below is the one in ``pipeline/candidates.py``: the strongest
single piece of evidence available for the row, plus a small MS1-support bonus.
``chemistry/conjugates.py`` computes its own ``0.35 + 0.35 * prior`` blend, so
values from that track are not on this scale.

A formula prior is only recomputed when a finite target neutral mass is
available. ``formula_prior_score`` falls back to element-ratio plausibility
alone when the mass is missing and returns 0.65 for any well-formed formula,
which would manufacture evidence for rows that carry no mass at all.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd

from structural_elucidation.chemistry.formula import formula_prior_score
from structural_elucidation.utils.domain import safe_float

MS1_SUPPORT_BONUS: float = 0.05

# Candidate columns are searched in order; the first non-empty one wins.
FORMULA_COLUMNS: tuple[str, ...] = ("predicted_formula", "candidate_formula", "peptide_formula")
NEUTRAL_MASS_COLUMNS: tuple[str, ...] = (
    "target_neutral_mass",
    "neutral_mass",
    "mass_hypothesis_neutral_mass",
    "theoretical_neutral_mass",
)
SPECTRAL_COLUMNS: tuple[str, ...] = ("forward_spectrum_score", "ms2_match_score", "nce_pipeline_score")

BASIS_COLUMN: str = "evidence_confidence_basis"
BASIS_EMITTED: str = "emitted"
BASIS_BACKFILLED: str = "backfilled"
BASIS_UNAVAILABLE: str = "unavailable"

_TRUE_STRINGS: frozenset[str] = frozenset({"true", "t", "yes", "y", "1"})


def evidence_confidence(
    *,
    formula_prior: Any = 0.0,
    forward_spectrum: Any = 0.0,
    ms2_match: Any = 0.0,
    nce_pipeline: Any = 0.0,
    ms1_supported: Any = False,
) -> float:
    """Strongest available evidence term for one candidate, plus MS1 support."""
    best = max(
        _finite(formula_prior),
        _finite(forward_spectrum),
        _finite(ms2_match),
        _finite(nce_pipeline),
    )
    bonus = MS1_SUPPORT_BONUS if _truthy(ms1_supported) else 0.0
    return float(np.clip(best + bonus, 0.0, 1.0))


def backfill_evidence_confidence(
    df: pd.DataFrame, *, tolerance_da: float = 0.05
) -> tuple[pd.Series, pd.Series]:
    """Fill ``evidence_confidence`` for rows whose generator did not emit it.

    Returns the completed score and a basis label per row. Positive emitted
    values are never overwritten; rows with no usable input are labelled
    ``unavailable`` and left at 0.0 so downstream layers can distinguish them
    from a computed zero.
    """
    if df.empty:
        return pd.Series(dtype=float), pd.Series(dtype=str)

    emitted = _numeric(df, "evidence_confidence")
    prior = _numeric(df, "formula_prior_score")
    formula = _first_present(df, FORMULA_COLUMNS)
    mass = _first_present_numeric(df, NEUTRAL_MASS_COLUMNS)

    recompute = prior.le(0.0) & formula.ne("") & np.isfinite(mass)
    for index in df.index[recompute]:
        prior.at[index] = formula_prior_score(formula.at[index], mass.at[index], tolerance_da=tolerance_da)

    spectral = [_numeric(df, column) for column in SPECTRAL_COLUMNS]
    strongest = pd.concat([prior, *spectral], axis=1).max(axis=1)
    bonus = _truthy_series(df, "ms1_buddy_supported").map({True: MS1_SUPPORT_BONUS, False: 0.0})
    computed = (strongest + bonus).clip(0.0, 1.0)

    missing = emitted.le(0.0)
    score = emitted.where(~missing, computed).astype(float)
    basis = pd.Series(BASIS_EMITTED, index=df.index, dtype=object)
    basis[missing & computed.gt(0.0)] = BASIS_BACKFILLED
    basis[missing & computed.le(0.0)] = BASIS_UNAVAILABLE
    return score, basis


def _finite(value: Any) -> float:
    """Coerce to a finite float, treating unparseable and NaN inputs as 0.0."""
    number = safe_float(value, math.nan)
    return float(number) if np.isfinite(number) else 0.0


def _truthy(value: Any) -> bool:
    """Parse mixed bool/numeric/string truth values from candidate CSVs."""
    if isinstance(value, bool):
        return value
    number = safe_float(value, math.nan)
    if np.isfinite(number):
        return number > 0.0
    return str(value).strip().lower() in _TRUE_STRINGS


def _numeric(df: pd.DataFrame, column: str) -> pd.Series:
    if column not in df.columns:
        return pd.Series(0.0, index=df.index, dtype=float)
    return pd.to_numeric(df[column], errors="coerce").fillna(0.0).astype(float)


def _truthy_series(df: pd.DataFrame, column: str) -> pd.Series:
    if column not in df.columns:
        return pd.Series(False, index=df.index, dtype=bool)
    return df[column].map(_truthy).fillna(False).astype(bool)


def _first_present(df: pd.DataFrame, columns: tuple[str, ...]) -> pd.Series:
    """First non-empty value across *columns*, as stripped strings."""
    out = pd.Series("", index=df.index, dtype=object)
    for column in columns:
        if column not in df.columns:
            continue
        values = df[column].fillna("").astype(str).str.strip()
        out = out.where(out.ne(""), values)
    return out.replace({"nan": ""})


def _first_present_numeric(df: pd.DataFrame, columns: tuple[str, ...]) -> pd.Series:
    """First finite numeric value across *columns*; NaN where none is present."""
    out = pd.Series(np.nan, index=df.index, dtype=float)
    for column in columns:
        if column not in df.columns:
            continue
        values = pd.to_numeric(df[column], errors="coerce")
        out = out.where(out.notna(), values)
    return out
