#!/usr/bin/env python3
"""Portable NCE-inspired formula prior features."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd

from structural_elucidation.nce.evidence import formula_harmony_score, formula_mass
from structural_elucidation.utils.domain import safe_float


def formula_mass_fit_score(formula: Any, neutral_mass: Any, tolerance_da: float = 0.05) -> float:
    """Score how closely a formula's monoisotopic mass matches the target neutral mass.

    Returns a clipped linear decay from 1.0 at exact match to 0.0 at
    ``tolerance_da`` or beyond, giving a cheap mass-consistency prior.
    """
    mass = formula_mass(formula)
    neutral = safe_float(neutral_mass, math.nan)
    if not np.isfinite(mass) or not np.isfinite(neutral):
        return 0.0
    return float(np.clip(1.0 - abs(mass - neutral) / max(float(tolerance_da), 1e-9), 0.0, 1.0))


def formula_prior_score(formula: Any, neutral_mass: Any, tolerance_da: float = 0.05) -> float:
    """Combine formula harmony and mass-fit into a single [0, 1] prior score.

    Harmony captures element-ratio plausibility; mass fit anchors the score to
    the observed neutral mass so implausible formulas are down-weighted early.
    """
    harmony = formula_harmony_score(formula, neutral_mass)
    fit = formula_mass_fit_score(formula, neutral_mass, tolerance_da)
    return float(np.clip(0.65 * harmony + 0.35 * fit, 0.0, 1.0))


def annotate_formula_prior(
    df: pd.DataFrame,
    *,
    formula_col: str = "target_formula",
    neutral_mass_col: str = "neutral_mass",
    output_col: str = "formula_prior_score",
    tolerance_da: float = 0.05,
) -> pd.DataFrame:
    """Add a ``formula_prior_score`` column by row-wise prior evaluation.

    Keeps the input frame immutable aside from the new column so downstream
    enrichment steps can chain multiple annotation helpers safely.
    """
    out = df.copy()
    scores = []
    for _, row in out.iterrows():
        scores.append(formula_prior_score(row.get(formula_col, ""), row.get(neutral_mass_col), tolerance_da))
    out[output_col] = scores
    return out
