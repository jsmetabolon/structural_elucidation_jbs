#!/usr/bin/env python3
"""Supervised peptide-spectrum calibration utilities.

This module keeps peptide spectral evidence separate from the existing
permutation heuristic.  It learns only from non-X known peptide spectra, then
applies the calibrated event weights to X candidates without using X truth.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from peptide_ranker_v1_utils import (
    IMMONIUM_MASS,
    NEUTRAL_LOSS_MASS,
    PROTON_MASS,
    RESIDUE_MASS,
    WATER_MASS,
    _bin_neighbors,
    _mass_matches_bins,
    mass_bin,
    parse_bins,
    prefix_suffix_masses,
    sequence_tokens,
)
from peptide_permutation_v1_utils import normalize_ordered_peptide_sequence
from x_domain_analog_utils import safe_float, write_json


DEFAULT_OUTPUT_ROOT = Path("models/peptide_spectrum_supervision_v1")
POLICY = "protected_champion_peptide_spectrum_calibrated_guarded"

EVENT_TYPES = [
    "b_ion",
    "y_ion",
    "b_loss",
    "y_loss",
    "immonium",
    "neutral_loss",
    "n_terminal_immonium",
    "c_terminal_immonium",
    "pyr_specific",
]


def load_calibration(path: Path) -> dict[str, Any]:
    if path.is_dir():
        path = path / "peptide_spectrum_calibration.json"
    if not path.exists():
        raise FileNotFoundError(f"Missing peptide spectrum calibration: {path}")
    return json.loads(path.read_text())


def expected_fragment_events(tokens: list[str]) -> list[dict[str, Any]]:
    events: list[dict[str, Any]] = []
    b_masses, y_masses = prefix_suffix_masses(tokens)
    for idx, mass in enumerate(b_masses, start=1):
        events.append({"event_type": "b_ion", "event_label": f"b{idx}", "mass": float(mass), "bin_kind": "ion"})
        events.append({"event_type": "b_loss", "event_label": f"b{idx}_loss", "mass": float(mass), "bin_kind": "loss"})
    for idx, mass in enumerate(y_masses, start=1):
        events.append({"event_type": "y_ion", "event_label": f"y{idx}", "mass": float(mass), "bin_kind": "ion"})
        events.append({"event_type": "y_loss", "event_label": f"y{idx}_loss", "mass": float(mass), "bin_kind": "loss"})
    for token in tokens:
        if token in IMMONIUM_MASS:
            events.append({"event_type": "immonium", "event_label": f"{token}_immonium", "mass": float(IMMONIUM_MASS[token]), "bin_kind": "ion"})
    if tokens:
        if tokens[0] in IMMONIUM_MASS:
            events.append({"event_type": "n_terminal_immonium", "event_label": f"N_{tokens[0]}_immonium", "mass": float(IMMONIUM_MASS[tokens[0]]), "bin_kind": "ion"})
        if tokens[-1] in IMMONIUM_MASS:
            events.append({"event_type": "c_terminal_immonium", "event_label": f"C_{tokens[-1]}_immonium", "mass": float(IMMONIUM_MASS[tokens[-1]]), "bin_kind": "ion"})
    for name, mass in NEUTRAL_LOSS_MASS.items():
        events.append({"event_type": "neutral_loss", "event_label": name, "mass": float(mass), "bin_kind": "loss"})
    if "Pyr" in tokens:
        events.extend(
            [
                {"event_type": "pyr_specific", "event_label": "Pyr_immonium", "mass": float(IMMONIUM_MASS["Pyr"]), "bin_kind": "ion"},
                {"event_type": "pyr_specific", "event_label": "Pyr_residue", "mass": float(RESIDUE_MASS["Pyr"] + PROTON_MASS), "bin_kind": "ion"},
                {"event_type": "pyr_specific", "event_label": "Pyr_water_loss", "mass": float(WATER_MASS), "bin_kind": "loss"},
            ]
        )
    return events


def event_match(event: dict[str, Any], ion_bins: set[int], loss_bins: set[int], window: int = 1) -> bool:
    mass = safe_float(event.get("mass"), math.nan)
    if not np.isfinite(mass) or mass <= 0:
        return False
    bins = ion_bins if str(event.get("bin_kind")) == "ion" else loss_bins
    return bool(_bin_neighbors(mass_bin(mass), int(window)) & bins)


def sequence_event_rows(
    *,
    sequence: Any,
    ion_bins: set[int],
    loss_bins: set[int],
    example_id: str,
    group_id: str,
    label: int,
    example_type: str,
    source_feature_key: Any,
) -> list[dict[str, Any]]:
    tokens = sequence_tokens(sequence)
    rows: list[dict[str, Any]] = []
    for event in expected_fragment_events(tokens):
        rows.append(
            {
                "example_id": example_id,
                "group_id": group_id,
                "label": int(label),
                "example_type": example_type,
                "source_feature_key": str(source_feature_key),
                "ordered_sequence": normalize_ordered_peptide_sequence(sequence),
                "event_type": event["event_type"],
                "event_label": event["event_label"],
                "event_mass": round(float(event["mass"]), 6),
                "bin_kind": event["bin_kind"],
                "observed_match": int(event_match(event, ion_bins, loss_bins)),
            }
        )
    return rows


def event_type_scores(sequence: Any, ion_bins: set[int], loss_bins: set[int]) -> dict[str, float]:
    tokens = sequence_tokens(sequence)
    events = expected_fragment_events(tokens)
    out: dict[str, float] = {}
    for event_type in EVENT_TYPES:
        scoped = [event for event in events if event["event_type"] == event_type]
        if not scoped:
            out[f"spectrum_{event_type}_score"] = 0.0
            out[f"spectrum_{event_type}_count"] = 0
            continue
        hits = sum(event_match(event, ion_bins, loss_bins) for event in scoped)
        out[f"spectrum_{event_type}_score"] = float(hits / max(len(scoped), 1))
        out[f"spectrum_{event_type}_count"] = int(hits)
    return out


def calibrate_event_weights(events: pd.DataFrame, alpha: float = 0.5) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for event_type in EVENT_TYPES:
        scoped = events[events.get("event_type", pd.Series(dtype=str)).astype(str).eq(event_type)].copy()
        if scoped.empty:
            rows.append({"event_type": event_type, "n_true": 0, "n_decoy": 0, "true_hit_rate": 0.0, "decoy_hit_rate": 0.0, "weight": 0.0})
            continue
        true = scoped[scoped["label"].astype(int).eq(1)]
        decoy = scoped[scoped["label"].astype(int).eq(0)]
        true_hits = float(pd.to_numeric(true.get("observed_match"), errors="coerce").fillna(0.0).sum())
        decoy_hits = float(pd.to_numeric(decoy.get("observed_match"), errors="coerce").fillna(0.0).sum())
        true_rate = (true_hits + alpha) / (max(true.shape[0], 0) + 2 * alpha)
        decoy_rate = (decoy_hits + alpha) / (max(decoy.shape[0], 0) + 2 * alpha)
        effect = max(float(true_rate - decoy_rate), 0.0)
        odds = math.log((true_rate * (1.0 - decoy_rate) + 1e-9) / (decoy_rate * (1.0 - true_rate) + 1e-9))
        weight = float(max(0.0, min(4.0, odds)) * max(effect, 0.05 if effect > 0 else 0.0))
        rows.append(
            {
                "event_type": event_type,
                "n_true": int(true.shape[0]),
                "n_decoy": int(decoy.shape[0]),
                "true_hit_rate": round(float(true_rate), 6),
                "decoy_hit_rate": round(float(decoy_rate), 6),
                "hit_rate_delta": round(float(effect), 6),
                "log_odds_lift": round(float(odds), 6),
                "weight": round(float(weight), 6),
            }
        )
    calib = pd.DataFrame(rows)
    weights = {str(row["event_type"]): float(row["weight"]) for _, row in calib.iterrows()}
    if sum(weights.values()) <= 1e-9:
        weights = {event_type: 1.0 for event_type in EVENT_TYPES}
    payload = {
        "event_weights": weights,
        "event_types": EVENT_TYPES,
        "calibration_alpha": float(alpha),
    }
    return calib, payload


def calibrated_spectrum_score(sequence: Any, ion_bins: set[int], loss_bins: set[int], calibration: dict[str, Any]) -> tuple[float, dict[str, float]]:
    features = event_type_scores(sequence, ion_bins, loss_bins)
    weights = {str(k): safe_float(v, 0.0) for k, v in dict(calibration.get("event_weights", {})).items()}
    numerator = 0.0
    denom = 0.0
    for event_type in EVENT_TYPES:
        weight = max(float(weights.get(event_type, 0.0)), 0.0)
        value = safe_float(features.get(f"spectrum_{event_type}_score"), 0.0)
        numerator += weight * value
        denom += weight
    score = numerator / denom if denom > 1e-12 else 0.0
    return float(np.clip(score, 0.0, 1.0)), features


def add_calibrated_spectrum_scores(df: pd.DataFrame, calibration: dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    score_rows: list[dict[str, Any]] = []
    for _, row in out.iterrows():
        sequence = row.get("candidate_ordered_sequence", row.get("peptide_composition", ""))
        ions = set()
        losses = set()
        for col in ["dual_channel__ion_bins_json", "light_top50__ion_bins_json", "raw_all__ion_bins_json", "observed_ion_bins_json"]:
            if col in row.index:
                ions |= parse_bins(row.get(col))
        for col in ["dual_channel__loss_bins_json", "light_top50__loss_bins_json", "raw_all__loss_bins_json", "observed_loss_bins_json"]:
            if col in row.index:
                losses |= parse_bins(row.get(col))
        score, features = calibrated_spectrum_score(sequence, ions, losses, calibration)
        score_rows.append({"peptide_spectrum_calibrated_score": score, **features})
    scores = pd.DataFrame(score_rows, index=out.index)
    for col in scores.columns:
        out[col] = pd.to_numeric(scores[col], errors="coerce").fillna(0.0)
    return out


def write_calibration(root: Path, calibration: dict[str, Any], event_calibration: pd.DataFrame) -> None:
    root.mkdir(parents=True, exist_ok=True)
    event_calibration.to_csv(root / "peptide_spectrum_event_calibration.csv", index=False)
    write_json(root / "peptide_spectrum_calibration.json", calibration)
