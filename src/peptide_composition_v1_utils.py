#!/usr/bin/env python3
"""Utilities for peptide composition normalization and control-blind scoring."""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from denovo_improvement_utils import canonicalize_smiles
from evaluate_denovo_controls import tanimoto
from evaluate_denovo_failure_modes import load_x_controls
from run_retrieval_anchored_denovo_assembly_v4 import resolve
from x_domain_analog_utils import safe_float, safe_int


GOOD_TANIMOTO = 0.7
RESIDUES = {
    "ALA": "Ala",
    "ARG": "Arg",
    "ASN": "Asn",
    "ASP": "Asp",
    "CYS": "Cys",
    "GLN": "Gln",
    "GLU": "Glu",
    "GLY": "Gly",
    "HIS": "His",
    "ILE": "Xle",
    "LEU": "Xle",
    "LYS": "Lys",
    "MET": "Met",
    "PHE": "Phe",
    "PRO": "Pro",
    "PYR": "Pyr",
    "SER": "Ser",
    "THR": "Thr",
    "TRP": "Trp",
    "TYR": "Tyr",
    "VAL": "Val",
    "XLE": "Xle",
}
FULL_RESIDUES = {
    "ALANINE": "Ala",
    "ARGININE": "Arg",
    "ASPARAGINE": "Asn",
    "ASPARTIC": "Asp",
    "CYSTEINE": "Cys",
    "GLUTAMINE": "Gln",
    "GLUTAMIC": "Glu",
    "GLYCINE": "Gly",
    "HISTIDINE": "His",
    "ISOLEUCINE": "Xle",
    "LEUCINE": "Xle",
    "LYSINE": "Lys",
    "METHIONINE": "Met",
    "PHENYLALANINE": "Phe",
    "PROLINE": "Pro",
    "PYROGLUTAMATE": "Pyr",
    "PYROGLUTAMIC": "Pyr",
    "SERINE": "Ser",
    "THREONINE": "Thr",
    "TRYPTOPHAN": "Trp",
    "TYROSINE": "Tyr",
    "VALINE": "Val",
}
RESIDUE_ORDER = {
    residue: idx
    for idx, residue in enumerate(
        ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Lys", "Met", "Phe", "Pro", "Pyr", "Ser", "Thr", "Trp", "Tyr", "Val", "Xle"]
    )
}


def parse_residue_tokens(value: Any) -> list[str]:
    text = str(value or "")
    if not text.strip() or text.strip().lower() == "nan":
        return []
    upper = text.upper()
    upper = re.sub(r"\bKNOWN\b|\bRELEASED\b|\bPUTATIVE\b|\bTENTATIVE\b", " ", upper)
    upper = re.sub(r"\b[DL]-", "", upper)
    tokens: list[str] = []
    pattern = "|".join(sorted([*FULL_RESIDUES.keys(), *RESIDUES.keys()], key=len, reverse=True))
    for match in re.finditer(rf"(?<![A-Z])({pattern})(?![A-Z])", upper):
        raw = match.group(1)
        token = FULL_RESIDUES.get(raw, RESIDUES.get(raw, ""))
        if token:
            tokens.append(token)
    return tokens


def normalize_peptide_composition(value: Any) -> str:
    tokens = parse_residue_tokens(value)
    if not tokens:
        return ""
    return "-".join(sorted(tokens, key=lambda token: (RESIDUE_ORDER.get(token, 999), token)))


def composition_length(value: Any) -> int:
    normalized = normalize_peptide_composition(value) if "-" in str(value or "") or str(value or "") in RESIDUE_ORDER else str(value or "")
    return len([part for part in normalized.split("-") if part])


def composition_match(left: Any, right: Any) -> bool:
    left_norm = normalize_peptide_composition(left)
    right_norm = normalize_peptide_composition(right)
    if not left_norm or not right_norm:
        return False
    return left_norm == right_norm


def numeric_series(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default).astype(float)


def peptide_truth_from_known_x(known_x: Path) -> pd.DataFrame:
    if not known_x.exists():
        return pd.DataFrame(columns=["control_id", "truth_peptide_composition", "truth_source", "truth_text"])
    df = pd.read_csv(known_x, low_memory=False)
    rows: list[dict[str, str]] = []
    text_cols = [
        col
        for col in ["KNOWN", "Tentative information", "External Information (Class or modification)", "Cleaned External Info", "SuperPathway", "Draft Subpathway"]
        if col in df.columns
    ]
    for _, row in df.iterrows():
        xid = str(row.get("Unknown", "") or "").strip()
        if not xid:
            continue
        text = ""
        normalized = ""
        for col in text_cols:
            candidate_text = str(row.get(col, "") or "")
            candidate = normalize_peptide_composition(candidate_text)
            if candidate and len(candidate.split("-")) >= 2:
                text = candidate_text
                normalized = candidate
                break
        if normalized and len(normalized.split("-")) >= 2:
            rows.append({"control_id": xid, "truth_peptide_composition": normalized, "truth_source": "known_x_annotation", "truth_text": text})
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["control_id", "truth_peptide_composition", "truth_source", "truth_text"])
    return out.drop_duplicates(["control_id", "truth_peptide_composition"]).reset_index(drop=True)


def s0_truth_map(repo_root: Path, known_x: Path, known_x_all: Path, lib_id: int = 209) -> dict[str, list[str]]:
    controls = load_x_controls(known_x=resolve(known_x, repo_root), known_x_all=resolve(known_x_all, repo_root), lib_id=int(lib_id), ion_mode="negative")
    return {control.control_id: list(control.smiles) for control in controls}


def best_tanimoto_to_truth(smiles: Any, truth_smiles: list[str]) -> float:
    canonical = canonicalize_smiles(smiles) or ""
    if not canonical:
        return 0.0
    return max([tanimoto(canonical, truth) for truth in truth_smiles] or [0.0])


def add_peptide_composition_features(df: pd.DataFrame, *, top_n: int = 30) -> pd.DataFrame:
    out = df.copy()
    if out.empty:
        return out
    if "control_id" not in out.columns:
        out["control_id"] = ""
    out["normalized_peptide_composition"] = out.get("peptide_composition", pd.Series("", index=out.index)).map(normalize_peptide_composition)
    out["composition_length"] = out["normalized_peptide_composition"].map(lambda value: len([part for part in str(value).split("-") if part]))
    out["is_peptide_candidate"] = out["normalized_peptide_composition"].astype(str).ne("")
    if "peptide_mass_error_da" in out.columns:
        mass_error_values = out["peptide_mass_error_da"]
    elif "seed_mass_error_da" in out.columns:
        mass_error_values = out["seed_mass_error_da"]
    else:
        mass_error_values = pd.Series(math.nan, index=out.index)
    out["composition_mass_error_da"] = pd.to_numeric(mass_error_values, errors="coerce").abs()
    seed_priority = numeric_series(out, "seed_priority_score", math.nan).fillna(numeric_series(out, "seed_source_score", 0.0)).fillna(0.0)
    generation = numeric_series(out, "generation_score", math.nan).fillna(numeric_series(out, "assembly_score", 0.0)).fillna(0.0)
    evidence = numeric_series(out, "evidence_confidence", 0.0)
    mass = numeric_series(out, "mass_match_score", math.nan).fillna(numeric_series(out, "score_mass_base", 0.0)).fillna(0.0)
    source = numeric_series(out, "seed_source_score", 0.0)
    out["composition_candidate_score"] = (0.30 * seed_priority + 0.22 * generation + 0.18 * evidence + 0.18 * mass + 0.12 * source).clip(0.0, 1.0)
    out["_rank_for_composition"] = pd.to_numeric(out.get("rank"), errors="coerce").fillna(999999).astype(float)
    out["_topn_for_composition"] = out.groupby("feature_key", dropna=False)["_rank_for_composition"].rank(method="first", ascending=True).le(int(top_n))
    peptide_top = out[out["is_peptide_candidate"] & out["_topn_for_composition"]].copy()
    if peptide_top.empty:
        out["composition_support_count"] = 0
        out["composition_support_score"] = 0.0
        out["composition_within_group_rank"] = 999999
        out["composition_length_plausibility"] = 0.0
        return out
    group_cols = ["feature_key", "normalized_peptide_composition"]
    grouped = (
        peptide_top.groupby(group_cols, dropna=False)
        .agg(
            composition_support_count=("predicted_smiles", "count"),
            composition_support_score=("composition_candidate_score", "sum"),
            composition_best_rank=("_rank_for_composition", "min"),
            composition_length=("composition_length", "first"),
        )
        .reset_index()
    )
    max_score = grouped.groupby("feature_key")["composition_support_score"].transform("max").replace(0, np.nan)
    grouped["composition_support_score_norm"] = (grouped["composition_support_score"] / max_score).fillna(0.0).clip(0.0, 1.0)
    length_support = grouped.groupby(["feature_key", "composition_length"], dropna=False)["composition_support_score"].sum().reset_index()
    length_support = length_support.sort_values(["feature_key", "composition_support_score", "composition_length"], ascending=[True, False, True], kind="mergesort")
    dominant_length = length_support.drop_duplicates("feature_key", keep="first")[["feature_key", "composition_length"]].rename(columns={"composition_length": "dominant_composition_length"})
    grouped = grouped.merge(dominant_length, on="feature_key", how="left")
    grouped["composition_length_plausibility"] = np.where(
        grouped["composition_length"].eq(grouped["dominant_composition_length"]),
        1.0,
        1.0 / (1.0 + (grouped["composition_length"] - grouped["dominant_composition_length"]).abs().clip(lower=0)),
    )
    out = out.merge(
        grouped[
            [
                "feature_key",
                "normalized_peptide_composition",
                "composition_support_count",
                "composition_support_score_norm",
                "composition_best_rank",
                "composition_length_plausibility",
            ]
        ],
        on=["feature_key", "normalized_peptide_composition"],
        how="left",
    )
    out["composition_support_count"] = numeric_series(out, "composition_support_count", 0).astype(int)
    out["composition_support_score"] = numeric_series(out, "composition_support_score_norm", 0.0)
    out["composition_within_group_rank"] = (
        out[out["is_peptide_candidate"]]
        .groupby(["feature_key", "normalized_peptide_composition"], dropna=False)["_rank_for_composition"]
        .rank(method="first", ascending=True)
        .reindex(out.index)
        .fillna(999999)
    )
    out["composition_length_plausibility"] = numeric_series(out, "composition_length_plausibility", 0.0)
    return out.drop(columns=["composition_support_score_norm"], errors="ignore")


def discover_protected_candidate_paths(candidate_root: Path, profiles: list[str] | None = None) -> list[Path]:
    base = candidate_root / "runs/protected_early_fusion/lib_209_negative/selfies_profiles"
    if not base.exists():
        return []
    wanted = set(profiles or [])
    out: list[Path] = []
    for path in sorted(base.glob("*/candidate_metrics_all.csv")):
        profile = path.parent.name
        if wanted and profile not in wanted:
            continue
        out.append(path)
    return out


def profile_from_candidate_path(path: Path) -> str:
    parts = list(path.parts)
    if "selfies_profiles" in parts:
        idx = parts.index("selfies_profiles")
        if idx + 1 < len(parts):
            return parts[idx + 1]
    return path.parent.name


def top_compositions_for_feature(scoped: pd.DataFrame, top_n: int) -> pd.DataFrame:
    empty_cols = [
        "normalized_peptide_composition",
        "composition_support_count",
        "composition_support_score",
        "best_permutation_rank",
        "top_permutation_smiles",
        "composition_length",
    ]
    work = scoped[scoped["is_peptide_candidate"]].copy()
    if work.empty:
        return pd.DataFrame(columns=empty_cols)
    work = work.sort_values(["_rank_for_composition", "composition_candidate_score", "predicted_smiles"], ascending=[True, False, True], kind="mergesort").head(int(top_n))
    grouped = (
        work.groupby("normalized_peptide_composition", dropna=False)
        .agg(
            composition_support_count=("predicted_smiles", "count"),
            composition_support_score=("composition_candidate_score", "sum"),
            best_permutation_rank=("_rank_for_composition", "min"),
            top_permutation_smiles=("predicted_smiles", "first"),
            composition_length=("composition_length", "first"),
        )
        .reset_index()
    )
    return grouped.sort_values(["composition_support_count", "composition_support_score", "best_permutation_rank"], ascending=[False, False, True], kind="mergesort")


def peptide_score_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = add_peptide_composition_features(df, top_n=30)
    seed_priority = numeric_series(out, "seed_priority_score", math.nan).fillna(numeric_series(out, "seed_source_score", 0.0)).fillna(0.0)
    generation = numeric_series(out, "generation_score", math.nan).fillna(numeric_series(out, "assembly_score", 0.0)).fillna(0.0)
    evidence = numeric_series(out, "evidence_confidence", 0.0)
    mass = numeric_series(out, "mass_match_score", math.nan).fillna(numeric_series(out, "score_mass_base", 0.0)).fillna(0.0)
    champion = numeric_series(out, "protected_champion", 0.0)
    rank_score = numeric_series(out, "score_rank", 0.0)
    peptide_gate = out["is_peptide_candidate"].astype(float)
    length_penalty = np.where(out["composition_length"].gt(3), 0.08 * (out["composition_length"] - 3), 0.0)
    out["score_peptide_composition_base"] = (
        peptide_gate
        * (
            0.30 * out["composition_support_score"]
            + 0.20 * seed_priority
            + 0.16 * mass
            + 0.12 * generation
            + 0.10 * evidence
            + 0.07 * out["composition_length_plausibility"]
            + 0.05 * (1.0 / out["composition_within_group_rank"].replace(0, 1).clip(lower=1))
            - length_penalty
        )
    ).clip(0.0, 1.0)
    out["score_protected_champion_base"] = champion.clip(0.0, 1.0)
    if "score_mass_base" not in out.columns:
        out["score_mass_base"] = mass.clip(0.0, 1.0)
    if "seed_priority_score" not in out.columns:
        out["seed_priority_score"] = seed_priority.clip(0.0, 1.0)
    out["score_peptide_protected_guarded"] = (10.0 * champion + 0.55 * out["score_peptide_composition_base"] + 0.18 * mass + 0.13 * seed_priority + 0.08 * rank_score + 0.06 * generation)
    return out
