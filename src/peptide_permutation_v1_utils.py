#!/usr/bin/env python3
"""Utilities for peptide permutation/connectivity ranking."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from peptide_composition_v1_utils import (
    RESIDUE_ORDER,
    add_peptide_composition_features,
    normalize_peptide_composition,
    numeric_series,
    parse_residue_tokens,
    profile_from_candidate_path,
)
from x_domain_analog_utils import safe_float


try:
    from rdkit import Chem
except Exception:  # pragma: no cover - exercised on systems without RDKit
    Chem = None  # type: ignore[assignment]


POLICY = "protected_champion_peptide_permutation_guarded"


def normalize_ordered_peptide_sequence(value: Any) -> str:
    tokens = parse_residue_tokens(value)
    return "-".join(tokens)


def _heavy_counts(mol: Any, start_idx: int, blocked: set[int]) -> dict[str, int]:
    stack = [start_idx]
    seen = set(blocked)
    counts: dict[str, int] = {}
    while stack:
        idx = stack.pop()
        if idx in seen:
            continue
        seen.add(idx)
        atom = mol.GetAtomWithIdx(int(idx))
        symbol = atom.GetSymbol()
        if symbol != "H":
            counts[symbol] = counts.get(symbol, 0) + 1
        for nbr in atom.GetNeighbors():
            nidx = nbr.GetIdx()
            if nidx not in seen:
                stack.append(nidx)
    return counts


def _has_double_oxygen(atom: Any) -> bool:
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetSymbol() == "O" and str(bond.GetBondType()) == "DOUBLE":
            return True
    return False


def _is_carbonyl(atom: Any) -> bool:
    return atom.GetSymbol() == "C" and _has_double_oxygen(atom)


def _sidechain_label(counts: dict[str, int]) -> str:
    c = counts.get("C", 0)
    n = counts.get("N", 0)
    o = counts.get("O", 0)
    s = counts.get("S", 0)
    if c == 0 and n == 0 and o == 0 and s == 0:
        return "Gly"
    if c == 1 and n == 0 and o == 0 and s == 0:
        return "Ala"
    if c == 3 and n == 1 and o == 1:
        return "Gln"
    if c == 3 and n == 0 and o >= 2:
        return "Glu"
    if c == 4 and n == 1:
        return "Lys"
    if c == 3 and n == 0 and o == 0 and s == 0:
        return "Val"
    if c == 4 and n == 0 and o == 0 and s == 0:
        return "Xle"
    if c == 3 and s == 1:
        return "Met"
    return ""


def _detect_pyr_residues(mol: Any) -> list[dict[str, Any]]:
    residues: list[dict[str, Any]] = []
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        ring_set = set(int(idx) for idx in ring)
        if len(ring_set) != 5:
            continue
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring_set]
        n_atoms = [atom for atom in ring_atoms if atom.GetSymbol() == "N"]
        carbonyls = [atom for atom in ring_atoms if _is_carbonyl(atom)]
        if len(n_atoms) != 1 or not carbonyls:
            continue
        n_atom = n_atoms[0]
        alpha_candidates = []
        for atom in ring_atoms:
            if atom.GetSymbol() != "C" or _is_carbonyl(atom):
                continue
            external_carbonyls = [
                nbr
                for nbr in atom.GetNeighbors()
                if nbr.GetIdx() not in ring_set and _is_carbonyl(nbr)
            ]
            if external_carbonyls:
                alpha_candidates.append((atom, external_carbonyls[0]))
        if not alpha_candidates:
            continue
        alpha, carboxyl = alpha_candidates[0]
        residues.append(
            {
                "label": "Pyr",
                "alpha_idx": alpha.GetIdx(),
                "n_idx": n_atom.GetIdx(),
                "carbonyl_idx": carboxyl.GetIdx(),
                "is_pyr": True,
            }
        )
    return residues


def _detect_regular_residues(mol: Any, pyr_alpha: set[int], pyr_carbonyl: set[int]) -> list[dict[str, Any]]:
    residues: list[dict[str, Any]] = []
    for carbonyl in mol.GetAtoms():
        if not _is_carbonyl(carbonyl) or carbonyl.GetIdx() in pyr_carbonyl:
            continue
        alpha_candidates = [
            nbr
            for nbr in carbonyl.GetNeighbors()
            if nbr.GetSymbol() == "C" and nbr.GetIdx() not in pyr_alpha and not _is_carbonyl(nbr)
        ]
        for alpha in alpha_candidates:
            n_candidates = [nbr for nbr in alpha.GetNeighbors() if nbr.GetSymbol() == "N"]
            if not n_candidates:
                continue
            n_atom = n_candidates[0]
            side_neighbors = [
                nbr
                for nbr in alpha.GetNeighbors()
                if nbr.GetIdx() not in {carbonyl.GetIdx(), n_atom.GetIdx()}
                and nbr.GetSymbol() not in {"H"}
            ]
            counts: dict[str, int] = {}
            for side in side_neighbors:
                sub_counts = _heavy_counts(mol, side.GetIdx(), {alpha.GetIdx()})
                for key, value in sub_counts.items():
                    counts[key] = counts.get(key, 0) + int(value)
            label = _sidechain_label(counts)
            if not label:
                continue
            residues.append(
                {
                    "label": label,
                    "alpha_idx": alpha.GetIdx(),
                    "n_idx": n_atom.GetIdx(),
                    "carbonyl_idx": carbonyl.GetIdx(),
                    "is_pyr": False,
                }
            )
    return residues


def infer_sequence_from_smiles(smiles: Any) -> dict[str, Any]:
    if Chem is None:
        return {
            "ordered_peptide_sequence": "",
            "sequence_parse_status": "no_rdkit",
            "peptide_bond_count": 0,
            "peptide_connectivity_confidence": 0.0,
            "pyr_terminal_score": 0.0,
            "peptide_branch_or_cycle_penalty": 0.0,
        }
    text = str(smiles or "").strip()
    if not text or text.lower() == "nan":
        return {
            "ordered_peptide_sequence": "",
            "sequence_parse_status": "empty_smiles",
            "peptide_bond_count": 0,
            "peptide_connectivity_confidence": 0.0,
            "pyr_terminal_score": 0.0,
            "peptide_branch_or_cycle_penalty": 0.0,
        }
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return {
            "ordered_peptide_sequence": "",
            "sequence_parse_status": "invalid_smiles",
            "peptide_bond_count": 0,
            "peptide_connectivity_confidence": 0.0,
            "pyr_terminal_score": 0.0,
            "peptide_branch_or_cycle_penalty": 0.0,
        }
    pyr = _detect_pyr_residues(mol)
    pyr_alpha = {int(row["alpha_idx"]) for row in pyr}
    pyr_carbonyl = {int(row["carbonyl_idx"]) for row in pyr}
    residues = pyr + _detect_regular_residues(mol, pyr_alpha, pyr_carbonyl)
    if not residues:
        return {
            "ordered_peptide_sequence": "",
            "sequence_parse_status": "no_residue_units",
            "peptide_bond_count": 0,
            "peptide_connectivity_confidence": 0.0,
            "pyr_terminal_score": 0.0,
            "peptide_branch_or_cycle_penalty": 0.0,
        }
    residue_by_n = {int(row["n_idx"]): idx for idx, row in enumerate(residues)}
    edges: dict[int, int] = {}
    incoming: dict[int, int] = {}
    peptide_bonds = 0
    for idx, row in enumerate(residues):
        carbonyl = mol.GetAtomWithIdx(int(row["carbonyl_idx"]))
        for nbr in carbonyl.GetNeighbors():
            nidx = nbr.GetIdx()
            if nbr.GetSymbol() == "N" and nidx in residue_by_n and residue_by_n[nidx] != idx:
                target = residue_by_n[nidx]
                edges[idx] = target
                incoming[target] = idx
                peptide_bonds += 1
                break
    starts = [idx for idx in range(len(residues)) if idx not in incoming]
    if len(starts) != 1:
        labels = [str(row["label"]) for row in residues]
        sequence = "-".join(labels)
        status = "ambiguous_connectivity"
    else:
        ordered: list[str] = []
        seen: set[int] = set()
        cur = starts[0]
        while cur not in seen:
            seen.add(cur)
            ordered.append(str(residues[cur]["label"]))
            if cur not in edges:
                break
            cur = edges[cur]
        sequence = "-".join(ordered)
        status = "parsed" if len(ordered) == len(residues) else "partial_connectivity"
    labels = [part for part in sequence.split("-") if part]
    expected_bonds = max(0, len(labels) - 1)
    bond_score = 1.0 / (1.0 + abs(float(peptide_bonds - expected_bonds)))
    if "Pyr" in labels:
        pyr_terminal_score = 1.0 if labels and labels[0] == "Pyr" else 0.15
    else:
        pyr_terminal_score = 0.5
    branch_penalty = 0.0
    if len(edges) > len(set(edges.values())):
        branch_penalty = 0.5
    if status in {"ambiguous_connectivity", "partial_connectivity"}:
        branch_penalty = max(branch_penalty, 0.25)
    return {
        "ordered_peptide_sequence": sequence,
        "sequence_parse_status": status,
        "peptide_bond_count": int(peptide_bonds),
        "peptide_connectivity_confidence": float(bond_score),
        "pyr_terminal_score": float(pyr_terminal_score),
        "peptide_branch_or_cycle_penalty": float(branch_penalty),
    }


def _topn_rank(df: pd.DataFrame) -> pd.Series:
    return pd.to_numeric(df.get("rank"), errors="coerce").fillna(999999).astype(float)


def add_peptide_permutation_features(df: pd.DataFrame, *, top_n: int = 30) -> pd.DataFrame:
    out = add_peptide_composition_features(df, top_n=top_n)
    if out.empty:
        return out
    sequence_from_label = out.get("peptide_composition", pd.Series("", index=out.index)).map(normalize_ordered_peptide_sequence)
    parsed_rows = out.get("predicted_smiles", pd.Series("", index=out.index)).map(infer_sequence_from_smiles)
    parsed = pd.DataFrame(list(parsed_rows), index=out.index)
    for col in parsed.columns:
        out[col] = parsed[col]
    out["candidate_ordered_sequence"] = np.where(
        out["ordered_peptide_sequence"].astype(str).ne(""),
        out["ordered_peptide_sequence"],
        sequence_from_label,
    )
    out["candidate_sequence_composition"] = out["candidate_ordered_sequence"].map(normalize_peptide_composition)
    out["sequence_matches_declared_composition"] = (
        out["candidate_sequence_composition"].astype(str).eq(out["normalized_peptide_composition"].astype(str))
        & out["normalized_peptide_composition"].astype(str).ne("")
    )
    out["permutation_rank_source"] = np.where(out["ordered_peptide_sequence"].astype(str).ne(""), "smiles_connectivity", "peptide_composition_field")
    out["_rank_for_permutation"] = _topn_rank(out)
    peptide_top = out[out["is_peptide_candidate"] & out["_rank_for_permutation"].le(int(top_n))].copy()
    if peptide_top.empty:
        out["ordered_sequence_support_count"] = 0
        out["ordered_sequence_support_score"] = 0.0
        out["ordered_sequence_within_composition_rank"] = 999999
    else:
        peptide_top["_sequence_candidate_score"] = (
            0.34 * numeric_series(peptide_top, "seed_priority_score", 0.0)
            + 0.22 * numeric_series(peptide_top, "mass_match_score", 0.0)
            + 0.18 * numeric_series(peptide_top, "evidence_confidence", 0.0)
            + 0.14 * numeric_series(peptide_top, "generation_score", 0.0)
            + 0.12 * numeric_series(peptide_top, "peptide_connectivity_confidence", 0.0)
        ).clip(0.0, 1.0)
        grouped = (
            peptide_top.groupby(["feature_key", "normalized_peptide_composition", "candidate_ordered_sequence"], dropna=False)
            .agg(
                ordered_sequence_support_count=("predicted_smiles", "count"),
                ordered_sequence_support_score=("_sequence_candidate_score", "sum"),
                ordered_sequence_best_rank=("_rank_for_permutation", "min"),
            )
            .reset_index()
        )
        max_score = grouped.groupby(["feature_key", "normalized_peptide_composition"], dropna=False)["ordered_sequence_support_score"].transform("max").replace(0, np.nan)
        grouped["ordered_sequence_support_score_norm"] = (grouped["ordered_sequence_support_score"] / max_score).fillna(0.0).clip(0.0, 1.0)
        grouped = grouped.sort_values(
            ["feature_key", "normalized_peptide_composition", "ordered_sequence_support_score", "ordered_sequence_best_rank"],
            ascending=[True, True, False, True],
            kind="mergesort",
        )
        grouped["ordered_sequence_within_composition_rank"] = grouped.groupby(["feature_key", "normalized_peptide_composition"], dropna=False).cumcount() + 1
        out = out.merge(
            grouped[
                [
                    "feature_key",
                    "normalized_peptide_composition",
                    "candidate_ordered_sequence",
                    "ordered_sequence_support_count",
                    "ordered_sequence_support_score_norm",
                    "ordered_sequence_best_rank",
                    "ordered_sequence_within_composition_rank",
                ]
            ],
            on=["feature_key", "normalized_peptide_composition", "candidate_ordered_sequence"],
            how="left",
        )
        out["ordered_sequence_support_count"] = numeric_series(out, "ordered_sequence_support_count", 0).astype(int)
        out["ordered_sequence_support_score"] = numeric_series(out, "ordered_sequence_support_score_norm", 0.0)
        out["ordered_sequence_within_composition_rank"] = numeric_series(out, "ordered_sequence_within_composition_rank", 999999).astype(int)
        out = out.drop(columns=["ordered_sequence_support_score_norm"], errors="ignore")
    out["permutation_length_penalty"] = np.where(out["composition_length"].gt(3), 0.16 * (out["composition_length"] - 3), 0.0)
    out["permutation_extra_residue_penalty"] = np.where(
        out["normalized_peptide_composition"].astype(str).str.contains("Ala-Gly", regex=False)
        & out["composition_length"].gt(3),
        0.12,
        0.0,
    )
    out["score_peptide_permutation_base"] = (
        out["is_peptide_candidate"].astype(float)
        * (
            0.24 * numeric_series(out, "composition_support_score", 0.0)
            + 0.22 * numeric_series(out, "pyr_terminal_score", 0.0)
            + 0.16 * numeric_series(out, "peptide_connectivity_confidence", 0.0)
            + 0.13 * numeric_series(out, "ordered_sequence_support_score", 0.0)
            + 0.10 * numeric_series(out, "seed_priority_score", 0.0)
            + 0.08 * numeric_series(out, "mass_match_score", 0.0)
            + 0.07 * numeric_series(out, "evidence_confidence", 0.0)
            - numeric_series(out, "permutation_length_penalty", 0.0)
            - numeric_series(out, "permutation_extra_residue_penalty", 0.0)
            - 0.20 * numeric_series(out, "peptide_branch_or_cycle_penalty", 0.0)
        )
    ).clip(0.0, 1.0)
    non_peptide_gate = (~out["is_peptide_candidate"].astype(bool)).astype(float)
    out["score_non_peptide_seed_base"] = (
        non_peptide_gate
        * (
            0.30 * numeric_series(out, "seed_priority_score", 0.0)
            + 0.22 * numeric_series(out, "generation_score", 0.0)
            + 0.20 * numeric_series(out, "mass_match_score", 0.0)
            + 0.15 * numeric_series(out, "evidence_confidence", 0.0)
            + 0.08 * numeric_series(out, "fusion_source_priority", 0.0)
            + 0.05 * numeric_series(out, "score_rank", 0.0)
        )
    ).clip(0.0, 1.0)
    protected = numeric_series(out, "protected_champion", 0.0).clip(0.0, 1.0)
    out["score_protected_champion_base"] = protected
    rank_score = numeric_series(out, "score_rank", 0.0)
    out["score_peptide_permutation_guarded"] = (
        10.0 * protected
        + 0.70 * out["score_peptide_permutation_base"]
        + 0.62 * out["score_non_peptide_seed_base"]
        + 0.15 * numeric_series(out, "score_peptide_composition_base", 0.0)
        + 0.08 * rank_score
        + 0.05 * (numeric_series(out, "score_mass_base", np.nan).fillna(numeric_series(out, "mass_match_score", 0.0)))
        + 0.02 * numeric_series(out, "fusion_source_priority", 0.0)
    )
    return out


def discover_permutation_feature_paths(feature_root: Path, profiles: list[str] | None = None) -> list[Path]:
    base = feature_root / "runs/peptide_permutation_features/lib_209_negative/selfies_profiles"
    if not base.exists():
        return []
    wanted = set(profiles or [])
    paths: list[Path] = []
    for path in sorted(base.glob("*/candidate_metrics_all.csv")):
        profile = profile_from_candidate_path(path)
        if wanted and profile not in wanted:
            continue
        paths.append(path)
    return paths


def rank_permutation_pool(scored: pd.DataFrame, *, top_k: int | None = None) -> pd.DataFrame:
    if scored.empty:
        return scored.copy()
    out = scored.copy()
    out["feature_key"] = out["feature_key"].astype(str)
    out["predicted_smiles"] = out["predicted_smiles"].astype(str)
    sort_cols = [
        "feature_key",
        "score_peptide_permutation_guarded",
        "score_protected_champion_base",
        "score_peptide_permutation_base",
        "pyr_terminal_score",
        "peptide_connectivity_confidence",
        "ordered_sequence_support_score",
        "composition_support_score",
        "seed_priority_score",
        "rank",
        "predicted_smiles",
    ]
    for col in sort_cols:
        if col not in out.columns:
            out[col] = "" if col in {"feature_key", "predicted_smiles"} else 0.0
    out = out.sort_values(
        sort_cols,
        ascending=[True, False, False, False, False, False, False, False, False, True, True],
        kind="mergesort",
    ).copy()
    out["rank"] = out.groupby("feature_key", dropna=False).cumcount() + 1
    out["rerank_policy"] = POLICY
    out["rerank_score"] = out["score_peptide_permutation_guarded"]
    if top_k is not None:
        return out[out["rank"].le(int(top_k))].reset_index(drop=True)
    return out.reset_index(drop=True)


def top_ordered_sequences_for_control(scoped: pd.DataFrame, top_n: int) -> pd.DataFrame:
    work = scoped[scoped.get("is_peptide_candidate", pd.Series(False, index=scoped.index)).astype(bool)].copy()
    if work.empty:
        return pd.DataFrame()
    work["_rank_tmp"] = pd.to_numeric(work.get("rank"), errors="coerce").fillna(999999)
    work = work[work["_rank_tmp"].le(int(top_n))]
    if work.empty:
        return pd.DataFrame()
    grouped = (
        work.groupby(["normalized_peptide_composition", "candidate_ordered_sequence"], dropna=False)
        .agg(
            support_count=("predicted_smiles", "count"),
            support_score=("score_peptide_permutation_base", "sum"),
            best_rank=("_rank_tmp", "min"),
            best_smiles=("predicted_smiles", "first"),
            best_connectivity_confidence=("peptide_connectivity_confidence", "max"),
            best_pyr_terminal_score=("pyr_terminal_score", "max"),
        )
        .reset_index()
    )
    return grouped.sort_values(["support_count", "support_score", "best_rank"], ascending=[False, False, True], kind="mergesort")


def dataframe_report_table(df: pd.DataFrame, max_rows: int = 40) -> str:
    if df.empty:
        return "_No rows._"
    from x_domain_analog_utils import dataframe_to_markdown

    return dataframe_to_markdown(df.head(int(max_rows)))


def metric_bool(value: Any) -> bool:
    return str(value).lower() in {"1", "true", "yes", "y"}


def as_float(value: Any, default: float = 0.0) -> float:
    try:
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return default
        return float(value)
    except Exception:
        return default
