#!/usr/bin/env python3
"""Shared feature code for peptide-specific ranking V1."""

from __future__ import annotations

import json
import math
import hashlib
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from denovo_improvement_utils import canonicalize_smiles, exact_mass_from_smiles, formula_from_smiles, inchikey_from_smiles
from nce_evidence_utils import bin_base, ion_and_loss_bins, mass_bin, parse_json_list
from peptide_composition_v1_utils import normalize_peptide_composition, numeric_series
from peptide_permutation_v1_utils import add_peptide_permutation_features, infer_sequence_from_smiles, normalize_ordered_peptide_sequence
from x_domain_analog_utils import safe_float, safe_int


POLICY = "protected_champion_peptide_model_guarded"
DEFAULT_OUTPUT_ROOT = Path("models/peptide_specific_ranker_v1")
REGRESSION_PROTECTED_209 = [
    "X-17919",
    "X-18345",
    "X-21410",
    "X-21821",
    "X-24137",
    "X-24972",
    "X-25459",
    "X-25523",
    "X-25744",
]

RESIDUE_MASS = {
    "Ala": 71.037113805,
    "Arg": 156.10111105,
    "Asn": 114.04292747,
    "Asp": 115.026943065,
    "Cys": 103.009184505,
    "Gln": 128.05857754,
    "Glu": 129.042593135,
    "Gly": 57.021463735,
    "His": 137.058911875,
    "Lys": 128.09496305,
    "Met": 131.040484645,
    "Phe": 147.068413945,
    "Pro": 97.052763875,
    "Pyr": 111.032028405,
    "Ser": 87.032028435,
    "Thr": 101.047678505,
    "Trp": 186.07931298,
    "Tyr": 163.063328575,
    "Val": 99.068413945,
    "Xle": 113.084064015,
}
WATER_MASS = 18.010564684
PROTON_MASS = 1.007276467

IMMONIUM_MASS = {
    "Ala": 44.0500,
    "Arg": 129.1135,
    "Asn": 87.0553,
    "Asp": 88.0393,
    "Cys": 76.0221,
    "Gln": 101.0709,
    "Glu": 102.0549,
    "Gly": 30.0344,
    "His": 110.0713,
    "Lys": 101.1073,
    "Met": 104.0534,
    "Phe": 120.0813,
    "Pro": 70.0657,
    "Pyr": 84.0444,
    "Ser": 60.0449,
    "Thr": 74.0606,
    "Trp": 159.0922,
    "Tyr": 136.0762,
    "Val": 72.0813,
    "Xle": 86.0969,
}
NEUTRAL_LOSS_MASS = {
    "water": 18.0106,
    "ammonia": 17.0265,
    "co2": 43.9898,
}

FEATURE_COLUMNS = [
    "feat_composition_match",
    "feat_same_length_as_feature",
    "feat_candidate_length",
    "feat_mass_match",
    "feat_existing_rank",
    "feat_seed_priority",
    "feat_generation",
    "feat_evidence_confidence",
    "feat_connectivity_confidence",
    "feat_pyr_terminal",
    "feat_branch_penalty_inverse",
    "feat_b_ion_score",
    "feat_y_ion_score",
    "feat_b_loss_score",
    "feat_y_loss_score",
    "feat_by_ion_score",
    "feat_by_loss_score",
    "feat_sequence_spectrum_score",
    "feat_immonium_score",
    "feat_neutral_loss_score",
    "feat_terminal_residue_score",
    "feat_n_terminal_immonium_score",
    "feat_c_terminal_immonium_score",
    "feat_pyr_specific_score",
    "feat_ms2_motif_score",
    "feat_nce_fragment_score",
    "feat_recursive_fragment_score",
    "feat_ordered_sequence_support",
    "feat_composition_support",
    "feat_non_peptide_seed",
]
FEATURE_FAMILIES = {
    "composition": ["feat_composition_match", "feat_same_length_as_feature", "feat_candidate_length"],
    "mass_rank_seed": ["feat_mass_match", "feat_existing_rank", "feat_seed_priority", "feat_generation", "feat_evidence_confidence"],
    "connectivity": ["feat_connectivity_confidence", "feat_pyr_terminal", "feat_branch_penalty_inverse", "feat_ordered_sequence_support", "feat_composition_support"],
    "fragment": [
        "feat_b_ion_score",
        "feat_y_ion_score",
        "feat_b_loss_score",
        "feat_y_loss_score",
        "feat_by_ion_score",
        "feat_by_loss_score",
        "feat_sequence_spectrum_score",
        "feat_immonium_score",
        "feat_neutral_loss_score",
        "feat_terminal_residue_score",
        "feat_n_terminal_immonium_score",
        "feat_c_terminal_immonium_score",
        "feat_pyr_specific_score",
    ],
    "external": ["feat_ms2_motif_score", "feat_nce_fragment_score", "feat_recursive_fragment_score", "feat_non_peptide_seed"],
}


def parse_csv(text: Any) -> list[str]:
    return [part.strip() for part in str(text or "").split(",") if part.strip()]


def parse_libs(text: Any) -> list[int]:
    return [int(part) for part in parse_csv(text)]


def resolve(path: Path, repo_root: Path) -> Path:
    return path if path.is_absolute() else repo_root / path


def sequence_tokens(sequence: Any) -> list[str]:
    text = normalize_ordered_peptide_sequence(sequence)
    return [part for part in text.split("-") if part]


def peptide_neutral_mass(tokens: list[str]) -> float:
    if not tokens or any(token not in RESIDUE_MASS for token in tokens):
        return float("nan")
    return float(sum(RESIDUE_MASS[token] for token in tokens) + WATER_MASS)


def prefix_suffix_masses(tokens: list[str]) -> tuple[list[float], list[float]]:
    if len(tokens) < 2:
        return [], []
    b: list[float] = []
    y: list[float] = []
    for idx in range(1, len(tokens)):
        prefix = tokens[:idx]
        suffix = tokens[idx:]
        if all(token in RESIDUE_MASS for token in prefix):
            b.append(float(sum(RESIDUE_MASS[token] for token in prefix) + PROTON_MASS))
        if all(token in RESIDUE_MASS for token in suffix):
            y.append(float(sum(RESIDUE_MASS[token] for token in suffix) + WATER_MASS + PROTON_MASS))
    return b, y


def _bin_neighbors(value: int, window: int = 1) -> set[int]:
    if not value:
        return set()
    return {int(value) + delta for delta in range(-int(window), int(window) + 1)}


def _mass_matches_bins(masses: list[float], bins: set[int]) -> int:
    hits = 0
    for mass in masses:
        if not np.isfinite(mass) or mass <= 0:
            continue
        b = mass_bin(mass)
        if _bin_neighbors(b) & bins:
            hits += 1
    return hits


def spectrum_scores(tokens: list[str], ion_bins: set[int], loss_bins: set[int]) -> dict[str, float]:
    b_masses, y_masses = prefix_suffix_masses(tokens)
    b_total = max(len(b_masses), 1)
    y_total = max(len(y_masses), 1)
    by_total = max(len(b_masses) + len(y_masses), 1)
    b_ion_hits = _mass_matches_bins(b_masses, ion_bins)
    y_ion_hits = _mass_matches_bins(y_masses, ion_bins)
    b_loss_hits = _mass_matches_bins(b_masses, loss_bins)
    y_loss_hits = _mass_matches_bins(y_masses, loss_bins)
    by_ion_hits = b_ion_hits + y_ion_hits
    by_loss_hits = b_loss_hits + y_loss_hits
    immonium = [IMMONIUM_MASS[token] for token in tokens if token in IMMONIUM_MASS]
    immonium_score = _mass_matches_bins(immonium, ion_bins) / max(len(immonium), 1)
    losses = list(NEUTRAL_LOSS_MASS.values())
    neutral_loss_score = _mass_matches_bins(losses, loss_bins) / max(len(losses), 1)
    terminal = []
    n_terminal = []
    c_terminal = []
    if tokens:
        n_terminal = [IMMONIUM_MASS.get(tokens[0], 0.0)]
        c_terminal = [IMMONIUM_MASS.get(tokens[-1], 0.0)]
        terminal.extend(n_terminal + c_terminal)
    terminal = [value for value in terminal if value > 0]
    n_terminal = [value for value in n_terminal if value > 0]
    c_terminal = [value for value in c_terminal if value > 0]
    terminal_score = _mass_matches_bins(terminal, ion_bins) / max(len(terminal), 1)
    n_terminal_score = _mass_matches_bins(n_terminal, ion_bins) / max(len(n_terminal), 1)
    c_terminal_score = _mass_matches_bins(c_terminal, ion_bins) / max(len(c_terminal), 1)
    pyr_specific = 0.0
    if "Pyr" in tokens:
        pyr_specific = float(
            max(
                _mass_matches_bins([IMMONIUM_MASS["Pyr"]], ion_bins),
                _mass_matches_bins([RESIDUE_MASS["Pyr"] + PROTON_MASS], ion_bins),
                _mass_matches_bins([WATER_MASS], loss_bins),
            )
            > 0
        )
    return {
        "feat_b_ion_score": float(b_ion_hits / b_total),
        "feat_y_ion_score": float(y_ion_hits / y_total),
        "feat_b_loss_score": float(b_loss_hits / b_total),
        "feat_y_loss_score": float(y_loss_hits / y_total),
        "feat_by_ion_score": float(by_ion_hits / by_total),
        "feat_by_loss_score": float(by_loss_hits / by_total),
        "feat_sequence_spectrum_score": float((by_ion_hits + by_loss_hits) / max(2 * by_total, 1)),
        "feat_immonium_score": float(immonium_score),
        "feat_neutral_loss_score": float(neutral_loss_score),
        "feat_terminal_residue_score": float(terminal_score),
        "feat_n_terminal_immonium_score": float(n_terminal_score),
        "feat_c_terminal_immonium_score": float(c_terminal_score),
        "feat_pyr_specific_score": float(pyr_specific),
    }


def parse_bins(value: Any) -> set[int]:
    out: set[int] = set()
    for item in parse_json_list(value):
        try:
            out.add(int(float(item)))
        except Exception:
            continue
    return out


def feature_observed_bins(row: pd.Series) -> tuple[set[int], set[int]]:
    ion_cols = [
        "dual_channel__ion_bins_json",
        "light_top50__ion_bins_json",
        "raw_all__ion_bins_json",
        "observed_ion_bins_json",
    ]
    loss_cols = [
        "dual_channel__loss_bins_json",
        "light_top50__loss_bins_json",
        "raw_all__loss_bins_json",
        "observed_loss_bins_json",
    ]
    ions: set[int] = set()
    losses: set[int] = set()
    for col in ion_cols:
        if col in row.index:
            ions |= parse_bins(row.get(col))
    for col in loss_cols:
        if col in row.index:
            losses |= parse_bins(row.get(col))
    return ions, losses


def feature_summary_to_bins(summary: pd.DataFrame) -> dict[str, tuple[set[int], set[int]]]:
    out: dict[str, tuple[set[int], set[int]]] = {}
    if summary.empty or "feature_key" not in summary.columns:
        return out
    for _, row in summary.iterrows():
        ions, losses = feature_observed_bins(row)
        out[str(row.get("feature_key", ""))] = (ions, losses)
        comp_id = str(safe_int(row.get("comp_id"), -1))
        if comp_id != "-1":
            out[comp_id] = (ions, losses)
    return out


def load_spectrum_evidence(path: Path) -> pd.DataFrame:
    if path.is_dir():
        candidate = path / "feature_spectrum_evidence.csv"
        if candidate.exists():
            path = candidate
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, low_memory=False)


def compute_feature_row(
    *,
    feature_key: Any,
    candidate_smiles: Any,
    candidate_sequence: Any = "",
    feature_composition: Any = "",
    feature_sequence: Any = "",
    ion_bins: set[int] | None = None,
    loss_bins: set[int] | None = None,
    source_row: pd.Series | None = None,
) -> dict[str, Any]:
    source = source_row if source_row is not None else pd.Series(dtype=object)
    parsed = infer_sequence_from_smiles(candidate_smiles)
    seq = str(parsed.get("ordered_peptide_sequence", "") or "") or normalize_ordered_peptide_sequence(candidate_sequence)
    tokens = sequence_tokens(seq)
    candidate_comp = normalize_peptide_composition(seq or candidate_sequence or source.get("peptide_composition", ""))
    feature_comp = normalize_peptide_composition(feature_composition or feature_sequence)
    ions = ion_bins or set()
    losses = loss_bins or set()
    spectrum = spectrum_scores(tokens, ions, losses)
    rank = safe_float(source.get("rank"), 999999.0)
    existing_rank = 1.0 / max(rank, 1.0) if np.isfinite(rank) else 0.0
    score_mass = safe_float(source.get("mass_match_score"), math.nan)
    if not np.isfinite(score_mass):
        score_mass = safe_float(source.get("score_mass_base"), 0.0)
    composition_len = len([part for part in feature_comp.split("-") if part])
    candidate_len = len(tokens) if tokens else len([part for part in candidate_comp.split("-") if part])
    out = {
        "feature_key": str(feature_key),
        "candidate_sequence": seq,
        "candidate_composition": candidate_comp,
        "feature_composition": feature_comp,
        "candidate_length": int(candidate_len),
        "feature_length": int(composition_len),
        "sequence_parse_status": parsed.get("sequence_parse_status", ""),
        "feat_composition_match": float(bool(candidate_comp and feature_comp and candidate_comp == feature_comp)),
        "feat_same_length_as_feature": float(bool(composition_len and candidate_len == composition_len)),
        "feat_candidate_length": float(min(candidate_len, 4) / 4.0),
        "feat_mass_match": float(np.clip(score_mass, 0.0, 1.0)),
        "feat_existing_rank": float(np.clip(existing_rank, 0.0, 1.0)),
        "feat_seed_priority": float(np.clip(safe_float(source.get("seed_priority_score"), safe_float(source.get("seed_source_score"), 0.0)), 0.0, 1.0)),
        "feat_generation": float(np.clip(safe_float(source.get("generation_score"), safe_float(source.get("assembly_score"), 0.0)), 0.0, 1.0)),
        "feat_evidence_confidence": float(np.clip(safe_float(source.get("evidence_confidence"), 0.0), 0.0, 1.0)),
        "feat_connectivity_confidence": float(np.clip(safe_float(parsed.get("peptide_connectivity_confidence"), 0.0), 0.0, 1.0)),
        "feat_pyr_terminal": float(np.clip(safe_float(parsed.get("pyr_terminal_score"), 0.0), 0.0, 1.0)),
        "feat_branch_penalty_inverse": float(1.0 - np.clip(safe_float(parsed.get("peptide_branch_or_cycle_penalty"), 0.0), 0.0, 1.0)),
        "feat_ms2_motif_score": float(np.clip(safe_float(source.get("ms2__ms2smiles_motif_score"), safe_float(source.get("ms2smiles_motif_score"), 0.0)), 0.0, 1.0)),
        "feat_nce_fragment_score": float(np.clip(safe_float(source.get("nce__score_fragment_gluing_v4"), safe_float(source.get("score_nce_fragment_base"), 0.0)), 0.0, 1.0)),
        "feat_recursive_fragment_score": float(np.clip(safe_float(source.get("recursive__recursive_fragment_group_score"), safe_float(source.get("score_recursive_base"), 0.0)), 0.0, 1.0)),
        "feat_ordered_sequence_support": float(np.clip(safe_float(source.get("ordered_sequence_support_score"), 0.0), 0.0, 1.0)),
        "feat_composition_support": float(np.clip(safe_float(source.get("composition_support_score"), 0.0), 0.0, 1.0)),
        "feat_non_peptide_seed": float(np.clip(safe_float(source.get("score_non_peptide_seed_base"), 0.0), 0.0, 1.0)),
    }
    out.update(spectrum)
    return out


def add_feature_columns(frame: pd.DataFrame, *, feature_composition_col: str = "feature_composition") -> pd.DataFrame:
    out = add_peptide_permutation_features(frame, top_n=30)
    rows: list[dict[str, Any]] = []
    for _, row in out.iterrows():
        ions, losses = feature_observed_bins(row)
        rows.append(
            compute_feature_row(
                feature_key=row.get("feature_key", ""),
                candidate_smiles=row.get("predicted_smiles", row.get("candidate_smiles", "")),
                candidate_sequence=row.get("candidate_ordered_sequence", row.get("peptide_composition", "")),
                feature_composition=row.get(feature_composition_col, row.get("normalized_peptide_composition", "")),
                ion_bins=ions,
                loss_bins=losses,
                source_row=row,
            )
        )
    feat = pd.DataFrame(rows, index=out.index)
    for col in FEATURE_COLUMNS:
        out[col] = pd.to_numeric(feat.get(col, pd.Series(0.0, index=out.index)), errors="coerce").fillna(0.0).astype(float)
    for col in ["candidate_sequence", "candidate_composition", "feature_composition", "sequence_parse_status"]:
        if col in feat.columns:
            out[col] = feat[col]
    return out


def model_input_matrix(df: pd.DataFrame, feature_columns: list[str] | None = None) -> np.ndarray:
    cols = feature_columns or FEATURE_COLUMNS
    work = df.copy()
    for col in cols:
        if col not in work.columns:
            work[col] = 0.0
    return work[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=np.float32)


def split_for_key(key: Any, val_fraction: float = 0.2) -> str:
    return "val" if stable_bucket(key, 10_000) / 10_000.0 < float(val_fraction) else "train"


def stable_bucket(value: Any, modulo: int) -> int:
    digest = hashlib.sha1(str(value).encode("utf-8")).hexdigest()
    return int(digest[:12], 16) % max(int(modulo), 1)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def candidate_exact_mass(row: pd.Series) -> float:
    mass = safe_float(row.get("predicted_exact_mass"), math.nan)
    if np.isfinite(mass):
        return mass
    return exact_mass_from_smiles(row.get("predicted_smiles", row.get("candidate_smiles", "")))


def structure_keys(smiles: Any) -> tuple[str, str]:
    canonical = canonicalize_smiles(smiles) or ""
    return canonical, inchikey_from_smiles(canonical) if canonical else ""
