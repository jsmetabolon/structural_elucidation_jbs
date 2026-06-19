#!/usr/bin/env python3
# pylint: disable=broad-exception-caught
"""Portable predicted-spectrum verifier for de novo candidates.

This mirrors the usable core of NewCompoundElucidation's
``PredictedSpectrumMatcher``: predict fragment/loss bins from a candidate
SMILES and score their overlap with the observed MS2 bins. The private
smiles2ms2 model loader is intentionally not imported here.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from structural_elucidation.nce.evidence import candidate_fragment_bins, overlap_score, parse_json_list
from structural_elucidation.utils.domain import safe_int


def _parse_bins(value: Any) -> list[int]:
    """Parse JSON or list bin columns into sorted unique positive integers."""
    out: list[int] = []
    for item in parse_json_list(value):
        try:
            intval = int(item)
        except Exception:
            continue
        if intval:
            out.append(intval)
    return sorted(set(out))


def predicted_spectrum_score(
    smiles: Any,
    observed_ion_bins: Any,
    observed_loss_bins: Any,
    *,
    lib_id: int = 209,
) -> float:
    """Score predicted fragment/loss bins against observed MS2 bins."""
    ion_bins = _parse_bins(observed_ion_bins)
    loss_bins = _parse_bins(observed_loss_bins)
    if not ion_bins and not loss_bins:
        return 0.5
    predicted_ions, predicted_losses = candidate_fragment_bins(str(smiles or ""), int(lib_id))
    ion_score = overlap_score(predicted_ions, ion_bins)
    loss_score = overlap_score(predicted_losses, loss_bins)
    return float(np.clip(0.55 * ion_score + 0.45 * loss_score, 0.0, 1.0))


def annotate_predicted_spectrum_scores(
    candidates: pd.DataFrame,
    *,
    smiles_col: str = "predicted_smiles",
    lib_id_col: str = "lib_id",
    output_col: str = "predicted_spectrum_verifier_score",
) -> pd.DataFrame:
    """Add a verifier score column to a candidate frame."""
    out = candidates.copy()
    if out.empty:
        out[output_col] = []
        return out
    scores: list[float] = []
    for _, row in out.iterrows():
        scores.append(
            predicted_spectrum_score(
                row.get(smiles_col, ""),
                row.get("observed_ion_bins_json", "[]"),
                row.get("observed_loss_bins_json", "[]"),
                lib_id=safe_int(row.get(lib_id_col), 209),
            )
        )
    out[output_col] = scores
    return out
