#!/usr/bin/env python3
# flake8: noqa
"""Portable NCE-inspired formula prior features."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd

from nce_evidence_utils import formula_harmony_score, formula_mass
from x_domain_analog_utils import safe_float


def formula_mass_fit_score(formula: Any, neutral_mass: Any, tolerance_da: float = 0.05) -> float:
    mass = formula_mass(formula)
    neutral = safe_float(neutral_mass, math.nan)
    if not np.isfinite(mass) or not np.isfinite(neutral):
        return 0.0
    return float(np.clip(1.0 - abs(mass - neutral) / max(float(tolerance_da), 1e-9), 0.0, 1.0))


def formula_prior_score(formula: Any, neutral_mass: Any, tolerance_da: float = 0.05) -> float:
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
    out = df.copy()
    scores = []
    for _, row in out.iterrows():
        scores.append(formula_prior_score(row.get(formula_col, ""), row.get(neutral_mass_col), tolerance_da))
    out[output_col] = scores
    return out
