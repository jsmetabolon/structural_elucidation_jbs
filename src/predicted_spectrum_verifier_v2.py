#!/usr/bin/env python3
# flake8: noqa
"""Real NCE smiles2ms2-backed predicted-spectrum verifier."""

from __future__ import annotations

from typing import Any

import pandas as pd

import json

import numpy as np

from nce_evidence_utils import PROTON_MASS, candidate_fragment_bins, mass_bin, overlap_score, parse_json_list
from nce_local_runtime_utils import score_smiles_against_observed_bins
from predicted_spectrum_verifier_v1 import predicted_spectrum_score as fallback_predicted_spectrum_score
from x_domain_analog_utils import safe_int


def _parse_bins(value: Any) -> list[Any]:
    return parse_json_list(value)


def _parse_ms2_pairs(value: Any) -> list[tuple[float, float]]:
    try:
        parsed = json.loads(str(value or "[]"))
    except Exception:
        return []
    out: list[tuple[float, float]] = []
    for item in parsed if isinstance(parsed, list) else []:
        if not isinstance(item, (list, tuple)) or not item:
            continue
        try:
            mz = float(item[0])
            intensity = float(item[1]) if len(item) > 1 else 1.0
        except Exception:
            continue
        if np.isfinite(mz) and mz > 0:
            out.append((mz, intensity if np.isfinite(intensity) else 1.0))
    return out


def _positive_adduct_bins_from_pairs(parent_mz: Any, ms2_pairs_json: Any) -> tuple[list[int], list[int]]:
    try:
        parent = float(parent_mz)
    except Exception:
        return ([], [])
    if not np.isfinite(parent) or parent <= PROTON_MASS:
        return ([], [])
    parent_neutral = parent - PROTON_MASS
    ions: set[int] = set()
    losses: set[int] = set()
    for mz, _intensity in _parse_ms2_pairs(ms2_pairs_json):
        fragment_neutral = mz - PROTON_MASS
        if fragment_neutral <= 0:
            continue
        ions.add(mass_bin(fragment_neutral))
        losses.add(mass_bin(max(parent_neutral - fragment_neutral, 1e-9)))
    return (sorted(ions), sorted(losses))


def _numeric_score(value: Any, default: float = 0.0) -> float:
    try:
        out = float(value)
    except Exception:
        return default
    return out if np.isfinite(out) else default


def _mass_or_adduct_score(row: pd.Series) -> float:
    mass_match = _numeric_score(row.get("mass_match_score"), default=float("nan"))
    if np.isfinite(mass_match):
        return float(np.clip(mass_match, 0.0, 1.0))
    abs_error = _numeric_score(row.get("abs_mass_error_da"), default=float("nan"))
    if np.isfinite(abs_error):
        return float(np.clip(np.exp(-abs(abs_error) / 0.05), 0.0, 1.0))
    return 0.0


def _strict_fragment_overlap_score(smiles: Any, ion_bins: list[Any], loss_bins: list[Any], lib_id: int) -> float:
    predicted_ions, predicted_losses = candidate_fragment_bins(str(smiles or ""), int(lib_id))
    ion_values = [int(v) for v in ion_bins if str(v) not in {"", "0", "nan"}]
    loss_values = [int(v) for v in loss_bins if str(v) not in {"", "0", "nan"}]
    if (not predicted_ions and not predicted_losses) or (not ion_values and not loss_values):
        return 0.0
    ion_score = overlap_score(predicted_ions, ion_values) if predicted_ions and ion_values else 0.0
    loss_score = overlap_score(predicted_losses, loss_values) if predicted_losses and loss_values else 0.0
    return float(np.clip(0.55 * ion_score + 0.45 * loss_score, 0.0, 1.0))


def positive_mode_adduct_aware_rescue_score(smiles: Any, row: pd.Series, observed: dict[str, Any]) -> float:
    """Positive-mode rescue score for precursor-dominated spectra.

    This is deliberately labelled as a rescue component: it combines strict
    fragment-bin overlap with precursor/adduct consistency. It should not be
    interpreted as pure smiles2ms2 fragment evidence.
    """

    lib_id = safe_int(row.get("lib_id", observed.get("lib_id", 400)), 400)
    if int(lib_id) not in {400, 402}:
        return 0.0
    current_ions = _parse_bins(observed.get("observed_ion_bins_json", "[]"))
    current_losses = _parse_bins(observed.get("observed_loss_bins_json", "[]"))
    raw_ions, raw_losses = _positive_adduct_bins_from_pairs(
        observed.get("parent_mz", row.get("parent_mz", 0.0)),
        observed.get("ms2_pairs_json", "[]"),
    )
    fragment_score = max(
        _strict_fragment_overlap_score(smiles, current_ions, current_losses, int(lib_id)),
        _strict_fragment_overlap_score(smiles, raw_ions, raw_losses, int(lib_id)),
    )
    mass_score = _mass_or_adduct_score(row)
    quality_score = 1.0 if len(_parse_ms2_pairs(observed.get("ms2_pairs_json", "[]"))) >= 3 else 0.0
    return float(np.clip(0.45 * fragment_score + 0.45 * mass_score + 0.10 * quality_score, 0.0, 1.0))


def predicted_spectrum_score_v2(
    smiles: Any,
    observed_ion_bins: Any,
    observed_loss_bins: Any,
    *,
    lib_id: int = 209,
    fallback_on_error: bool = True,
) -> float:
    try:
        score = score_smiles_against_observed_bins(
            smiles,
            _parse_bins(observed_ion_bins),
            _parse_bins(observed_loss_bins),
            lib_id=int(lib_id),
        )
    except Exception:
        score = 0.0
    if score <= 0 and fallback_on_error:
        return float(fallback_predicted_spectrum_score(smiles, observed_ion_bins, observed_loss_bins, lib_id=int(lib_id)))
    return float(score)


def annotate_predicted_spectrum_scores_v2(
    candidates: pd.DataFrame,
    observed_lookup: pd.DataFrame,
    *,
    smiles_col: str = "predicted_smiles",
    output_col: str = "predicted_spectrum_score_v2",
) -> pd.DataFrame:
    out = candidates.copy()
    if out.empty:
        out[output_col] = []
        return out
    if not observed_lookup.empty and "control_id" in observed_lookup.columns:
        observed_lookup = observed_lookup.copy()
        observed_lookup["control_id"] = observed_lookup["control_id"].fillna("").astype(str)
        observed_lookup = observed_lookup[observed_lookup["control_id"].ne("")].drop_duplicates("control_id", keep="first")
    lookup = observed_lookup.set_index("control_id").to_dict(orient="index") if not observed_lookup.empty and "control_id" in observed_lookup.columns else {}
    cache: dict[tuple[str, str, str, int], float] = {}
    scores: list[float] = []
    sources: list[str] = []
    raw_scores: list[float] = []
    adduct_scores: list[float] = []
    for _, row in out.iterrows():
        control_id = str(row.get("control_id", "") or "")
        observed = lookup.get(control_id, {})
        if not observed:
            scores.append(0.0)
            sources.append("missing_observed_ms2")
            raw_scores.append(0.0)
            adduct_scores.append(0.0)
            continue
        smiles = str(row.get(smiles_col, "") or "")
        ion_json = str(observed.get("observed_ion_bins_json", "[]") or "[]")
        loss_json = str(observed.get("observed_loss_bins_json", "[]") or "[]")
        lib_id = safe_int(row.get("lib_id", observed.get("lib_id", 209)), 209)
        key = (smiles, ion_json, loss_json, int(lib_id))
        if key not in cache:
            cache[key] = predicted_spectrum_score_v2(smiles, ion_json, loss_json, lib_id=int(lib_id))
        raw_score = float(cache[key])
        adduct_score = positive_mode_adduct_aware_rescue_score(smiles, row, observed)
        final_score = max(raw_score, adduct_score)
        raw_scores.append(raw_score)
        adduct_scores.append(adduct_score)
        scores.append(final_score)
        if adduct_score > raw_score:
            sources.append(f"{observed.get('spectrum_source', '')}|v19_positive_adduct_mass_fragment_rescue")
        else:
            sources.append(str(observed.get("spectrum_source", "")))
    out[output_col] = scores
    out[f"{output_col}_raw_nce_or_portable"] = raw_scores
    out[f"{output_col}_v19_positive_adduct_aware"] = adduct_scores
    out["predicted_spectrum_score_v2_source"] = sources
    return out
