#!/usr/bin/env python3
"""Build HMDB/NCE-snapshot anchored non-peptide candidates for missed S0 controls.

Generates library-mass-matched structures (plus conjugate expansions) for X
controls that failed primary recovery, using diagnosis CSV conjugate hints.

Pipeline context:
    Upstream of ``verify`` and ``pool``. Targets controls where
    ``primary_stage != recovered`` in the diagnosis table. Writes under
    ``models/library_anchored_candidates_v1/runs/library_anchored_v1/``.

Key I/O:
    Inputs: ``control_stage_diagnosis.csv``, feature rows, neutral-mass
    hypotheses, HMDB/vendor libraries, NCE full-smiles snapshot.
    Outputs: ``candidate_metrics_all.csv``, ``denovo_report.csv``, per-control
    metrics, source summary, and markdown report.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from structural_elucidation.chemistry.conjugates import CONJUGATE_REACTIONS, apply_reaction, mass_match_score
from structural_elucidation.utils.denovo import canonicalize_smiles, exact_mass_from_smiles, formula_from_smiles
from structural_elucidation.evaluation.failure_modes import load_x_controls
from structural_elucidation.chemistry.formula import formula_prior_score
from structural_elucidation.evaluation.rerank_assembly import eval_reranked, write_denovo_report
from structural_elucidation.paths import REPO_ROOT, resolve
from structural_elucidation.utils.domain import (
    dataframe_to_markdown,
    safe_float,
    safe_int,
    write_csv,
    write_json,
    xid_from_feature_key,
)

# ── Default library paths and run label ──────────────────────────────
DEFAULT_HMDB: Path = REPO_ROOT / "vendor/nce/elucidation/table_builders/hmdb_minimal.csv"
DEFAULT_HMDB_EXTENDED: Path = (
    REPO_ROOT / "vendor/nce/elucidation/table_builders/HMDB_36_classyfire_21_annotations.csv"
)
DEFAULT_NCE_FULL: Path = REPO_ROOT / "data/nce_full_smiles_candidates.csv"
RUN_NAME: str = "library_anchored_v1"  # output subdirectory under runs/


# ── Column coercion helpers ──────────────────────────────────────────
def numeric(df: pd.DataFrame, col: str, default: float = math.nan) -> pd.Series:
    """Coerce a DataFrame column to float, filling missing with *default*."""
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def bool_value(value: Any) -> bool:
    """Parse common truthy string representations."""
    return str(value or "").strip().lower() in {"true", "1", "yes", "y"}


def ppm_hit_mask(masses: pd.Series, target: float, ppm: float) -> pd.Series:
    """Return a boolean mask for masses within *ppm* of *target*.

    Uses absolute ppm tolerance (not asymmetric) so small targets are not
    over-matched relative to large ones.
    """
    values = pd.to_numeric(masses, errors="coerce")
    tol = abs(float(target)) * float(ppm) / 1_000_000.0
    return values.sub(float(target)).abs().le(tol)


# ── Structure library loaders ────────────────────────────────────────
def prepare_hmdb(path: Path) -> pd.DataFrame:
    """Load the minimal HMDB table into the unified library schema."""
    if not path.is_file():
        return pd.DataFrame()
    raw = pd.read_csv(path, low_memory=False)
    smiles_col = "canonical_smiles" if "canonical_smiles" in raw.columns else "smiles"
    out = raw.copy()
    out["library_source"] = "hmdb_minimal"
    out["library_name"] = out.get("name", pd.Series("", index=out.index))
    out["library_smiles"] = out[smiles_col].astype(str)
    out["library_formula"] = out.get("chemical_formula", pd.Series("", index=out.index))
    out["library_exact_mass"] = numeric(out, "mass")
    return out[
        ["library_source", "library_name", "library_smiles", "library_formula", "library_exact_mass"]
        + [col for col in ["accession", "chemical_id"] if col in out.columns]
    ].dropna(subset=["library_exact_mass"])


def first_existing_col(df: pd.DataFrame, names: list[str]) -> str | None:
    """Return the first column name from *names* present in *df* (case-insensitive)."""
    lower = {col.lower(): col for col in df.columns}
    for name in names:
        if name in df.columns:
            return name
        if name.lower() in lower:
            return lower[name.lower()]
    return None


def exact_mass_series(smiles: pd.Series) -> pd.Series:
    """Compute monoisotopic mass per SMILES row when the CSV lacks a mass column."""
    return pd.Series([exact_mass_from_smiles(value) for value in smiles], index=smiles.index, dtype=float)


def prepare_structure_library(path: Path, source_label: str) -> pd.DataFrame:
    """Normalize an arbitrary structure CSV into the unified library schema."""
    if not path.is_file():
        return pd.DataFrame()
    raw = pd.read_csv(path, low_memory=False)
    if raw.empty:
        return pd.DataFrame()
    smiles_col = first_existing_col(
        raw, ["canonical_smiles", "smiles", "SMILES", "Smiles", "structure_smiles", "library_smiles"]
    )
    if smiles_col is None:
        return pd.DataFrame()
    name_col = first_existing_col(raw, ["name", "Name", "compound_name", "library_name", "Common_Name", "COMMON_NAME"])
    formula_col = first_existing_col(
        raw, ["chemical_formula", "formula", "Formula", "molecular_formula", "candidate_formula"]
    )
    mass_col = first_existing_col(
        raw,
        ["exact_mass", "mass", "monoisotopic_mass", "MonoisotopicMass", "library_exact_mass", "candidate_exact_mass"],
    )
    out = raw.copy()
    out["library_source"] = str(source_label)
    out["library_name"] = out[name_col].astype(str) if name_col else ""
    out["library_smiles"] = out[smiles_col].astype(str)
    out["library_formula"] = out[formula_col].astype(str) if formula_col else ""
    out["library_exact_mass"] = numeric(out, mass_col) if mass_col else exact_mass_series(out["library_smiles"])
    keep = ["library_source", "library_name", "library_smiles", "library_formula", "library_exact_mass"]
    for col in [
        "accession",
        "chemical_id",
        "hmdb_id",
        "CompoundID",
        "chebi_id",
        "lipidmaps_id",
        "drugbank_id",
        "public_spectral_prevalence",
    ]:
        if col in out.columns:
            keep.append(col)
    return out[keep].dropna(subset=["library_exact_mass"])


def prepare_hmdb_extended(path: Path) -> pd.DataFrame:
    """Load HMDB extended annotations, joining spreadsheet names when needed."""
    if not path.is_file():
        return pd.DataFrame()
    raw = pd.read_csv(path, low_memory=False)
    if {"CompoundID", "Smiles"}.issubset(raw.columns) and "exact_mass" not in raw.columns:
        spreadsheet = path.parent / "hmdb_spreadsheet.csv"
        if spreadsheet.exists():
            base = pd.read_csv(spreadsheet, low_memory=False)
            base_key = first_existing_col(base, ["hmdb_id", "secondary_ids"])
            if base_key:
                lookup = base.copy()
                lookup["_hmdb_join_key"] = lookup[base_key].astype(str).str.replace("HMDB", "HMDB", regex=False)
                annot = raw[["CompoundID", "ParentName"]].drop_duplicates("CompoundID").copy()
                annot["_hmdb_join_key"] = annot["CompoundID"].astype(str)
                merged = lookup.merge(annot, on="_hmdb_join_key", how="left")
                if "ParentName" in merged.columns:
                    merged["library_name"] = (
                        merged.get("name", pd.Series("", index=merged.index)).astype(str)
                        + " | "
                        + merged["ParentName"].fillna("").astype(str)
                    )
                return prepare_structure_library_from_frame(merged, "hmdb_extended")
    return prepare_structure_library(path, "hmdb_extended")


def prepare_structure_library_from_frame(raw: pd.DataFrame, source_label: str) -> pd.DataFrame:
    """Normalize an in-memory frame (post-join) into the library schema."""
    if raw.empty:
        return pd.DataFrame()
    smiles_col = first_existing_col(raw, ["canonical_smiles", "smiles", "SMILES", "Smiles", "library_smiles"])
    if smiles_col is None:
        return pd.DataFrame()
    name_col = first_existing_col(raw, ["library_name", "name", "Name", "compound_name"])
    formula_col = first_existing_col(raw, ["chemical_formula", "formula", "Formula", "molecular_formula"])
    mass_col = first_existing_col(raw, ["exact_mass", "mass", "monoisotopic_mass", "library_exact_mass"])
    out = raw.copy()
    out["library_source"] = str(source_label)
    out["library_name"] = out[name_col].astype(str) if name_col else ""
    out["library_smiles"] = out[smiles_col].astype(str)
    out["library_formula"] = out[formula_col].astype(str) if formula_col else ""
    out["library_exact_mass"] = numeric(out, mass_col) if mass_col else exact_mass_series(out["library_smiles"])
    keep = ["library_source", "library_name", "library_smiles", "library_formula", "library_exact_mass"]
    for col in ["accession", "chemical_id", "hmdb_id", "CompoundID"]:
        if col in out.columns:
            keep.append(col)
    return out[keep].dropna(subset=["library_exact_mass"])


def prepare_external_libraries(args: argparse.Namespace) -> pd.DataFrame:
    """Concatenate all configured vendor libraries, deduplicating by canonical SMILES."""
    frames = [
        prepare_hmdb(Path(str(args.hmdb_csv)).expanduser()),
        prepare_hmdb_extended(Path(str(args.hmdb_extended_csv)).expanduser()),
        prepare_structure_library(Path(str(args.chebi_csv)).expanduser(), "chebi"),
        prepare_structure_library(Path(str(args.lipidmaps_csv)).expanduser(), "lipidmaps"),
        prepare_structure_library(Path(str(args.drugbank_csv)).expanduser(), "drugbank_open"),
        prepare_structure_library(Path(str(args.masst_structures_csv)).expanduser(), "masst_public_structures"),
    ]
    frames = [frame for frame in frames if not frame.empty]
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True, sort=False)
    out["library_smiles"] = out["library_smiles"].astype(str)
    out["_canonical_library_smiles"] = out["library_smiles"].map(lambda value: canonicalize_smiles(value) or "")
    out = out[out["_canonical_library_smiles"].astype(bool)].copy()
    # Dedup per source so the same structure from HMDB and ChEBI stays distinct.
    out = out.drop_duplicates(["library_source", "_canonical_library_smiles"], keep="first")
    return out.drop(columns=["_canonical_library_smiles"]).reset_index(drop=True)


def prepare_nce_snapshot(path: Path, lib_id: int) -> pd.DataFrame:
    """Load NCE full-smiles snapshot rows scoped to *lib_id*."""
    if not path.is_file():
        return pd.DataFrame()
    raw = pd.read_csv(path, low_memory=False)
    if "lib_id" in raw.columns:
        raw = raw[numeric(raw, "lib_id").eq(int(lib_id))].copy()
    if raw.empty:
        return pd.DataFrame()
    out = raw.copy()
    out["library_source"] = "nce_snapshot_full_smiles"
    out["library_name"] = out.get("candidate_name", pd.Series("", index=out.index))
    out["library_smiles"] = out.get("smiles", pd.Series("", index=out.index)).astype(str)
    out["library_formula"] = out.get("candidate_formula", pd.Series("", index=out.index))
    out["library_exact_mass"] = numeric(out, "candidate_exact_mass")
    out["control_id"] = out.get("x_compound_id", pd.Series("", index=out.index)).astype(str)
    keep = [
        "library_source",
        "library_name",
        "library_smiles",
        "library_formula",
        "library_exact_mass",
        "control_id",
        "feature_key",
        "comp_id",
        "source_track",
        "modification",
        "forward_spectrum_score",
        "ms2_match_score",
        "nce_pipeline_score",
        "neutral_mass",
        "candidate_comp_id",
    ]
    return out[[col for col in keep if col in out.columns]].dropna(subset=["library_exact_mass"])


# ── Control selection and mass hypotheses ────────────────────────────
def target_controls(diagnosis: pd.DataFrame) -> list[str]:
    """Return control IDs that did not reach ``recovered`` in the primary stage."""
    primary = diagnosis.get("primary_stage", pd.Series("", index=diagnosis.index)).astype(str)
    return sorted(set(diagnosis.loc[~primary.eq("recovered"), "control_id"].astype(str)))


def parse_controls(text: Any) -> list[str]:
    """Parse a comma-separated control ID list, preserving order and uniqueness."""
    controls = [part.strip() for part in str(text or "").split(",") if part.strip()]
    return sorted(dict.fromkeys(controls))


def mass_hypotheses(features: pd.DataFrame, adducts: pd.DataFrame) -> pd.DataFrame:
    """Union neutral-mass targets from feature fallbacks and adduct hypotheses."""
    feature_map = features.copy()
    if "control_id" not in feature_map.columns:
        feature_map["control_id"] = feature_map.get("feature_key", pd.Series("", index=feature_map.index)).map(
            xid_from_feature_key
        )
    adduct_work = adducts.copy()
    if not adduct_work.empty:
        adduct_work["control_id"] = adduct_work.get("feature_key", pd.Series("", index=adduct_work.index)).map(
            xid_from_feature_key
        )
    rows: list[dict[str, Any]] = []
    for _, row in feature_map.iterrows():
        rows.append(
            {
                "control_id": str(row.get("control_id", "")),
                "feature_key": row.get("feature_key", ""),
                "comp_id": row.get("comp_id", ""),
                "lib_id": row.get("lib_id", 209),
                "ion_mode": row.get("ion_mode", "negative"),
                "parent_mz": row.get("parent_mz", ""),
                "target_neutral_mass": row.get("neutral_mass_fallback", ""),
                "adduct_label": "neutral_mass_fallback",
                "formula_prior_score": 0.0,
                "ms1_buddy_supported": False,
            }
        )
    for _, row in adduct_work.iterrows():
        rows.append(
            {
                "control_id": str(row.get("control_id", "")),
                "feature_key": row.get("feature_key", ""),
                "comp_id": row.get("comp_id", ""),
                "lib_id": row.get("lib_id", 209),
                "ion_mode": row.get("ion_mode", "negative"),
                "parent_mz": row.get("parent_mz", ""),
                "target_neutral_mass": row.get("neutral_mass", ""),
                "adduct_label": row.get("adduct_label", ""),
                "formula_prior_score": row.get("formula_prior_score", 0.0),
                "ms1_buddy_supported": row.get("ms1_buddy_supported", False),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["target_neutral_mass"] = pd.to_numeric(out["target_neutral_mass"], errors="coerce")
    return (
        out[out["target_neutral_mass"].notna()]
        .drop_duplicates(["control_id", "feature_key", "target_neutral_mass"])
        .reset_index(drop=True)
    )


# ── Candidate row construction ───────────────────────────────────────
def append_candidate(
    rows: list[dict[str, Any]],
    *,
    target: pd.Series,
    hit: pd.Series,
    predicted_smiles: str,
    operator_family: str,
    library_match_type: str,
    args: argparse.Namespace,
) -> None:
    """Append one scored candidate row if SMILES canonicalizes and mass is finite."""
    canonical = canonicalize_smiles(predicted_smiles)
    if not canonical:
        return
    exact_mass = exact_mass_from_smiles(canonical)
    if not np.isfinite(exact_mass):
        exact_mass = safe_float(hit.get("library_exact_mass"), math.nan)
    target_mass = safe_float(target.get("target_neutral_mass"), math.nan)
    abs_error = (
        abs(float(exact_mass) - float(target_mass))
        if np.isfinite(exact_mass) and np.isfinite(target_mass)
        else math.nan
    )
    formula = formula_from_smiles(canonical) or str(hit.get("library_formula", ""))
    formula_prior = formula_prior_score(formula, target_mass, tolerance_da=float(args.formula_tolerance_da))
    mass_score = mass_match_score(abs_error, float(args.mass_score_tolerance_da))
    nce_forward = safe_float(hit.get("forward_spectrum_score"), 0.0)
    nce_ms2 = safe_float(hit.get("ms2_match_score"), 0.0)
    nce_pipeline = safe_float(hit.get("nce_pipeline_score"), 0.0)
    ms1_supported = bool_value(target.get("ms1_buddy_supported"))
    # evidence_confidence feeds downstream rerank features and F2 promotion gates.
    evidence = float(
        np.clip(max(formula_prior, nce_forward, nce_ms2, nce_pipeline) + (0.05 if ms1_supported else 0.0), 0.0, 1.0)
    )
    # Mass-heavy blend: library anchor is primarily a mass/formula rescue track.
    generation = float(
        np.clip(
            0.40 * mass_score + 0.25 * formula_prior + 0.20 * nce_pipeline + 0.10 * nce_forward + 0.05 * nce_ms2,
            0.0,
            1.0,
        )
    )
    feature_key = str(
        target.get("feature_key", "") or hit.get("feature_key", "") or f"library_anchor_{target.get('control_id', '')}"
    )
    rows.append(
        {
            "feature_key": feature_key,
            "control_id": str(target.get("control_id", "")),
            "comp_id": safe_int(target.get("comp_id", hit.get("comp_id", -1)), -1),
            "lib_id": safe_int(target.get("lib_id"), 209),
            "ion_mode": str(target.get("ion_mode", "negative") or "negative"),
            "profile": str(args.profile),
            "track": str(args.track),
            "predicted_smiles": canonical,
            "predicted_exact_mass": round(float(exact_mass), 7) if np.isfinite(exact_mass) else "",
            "predicted_formula": formula,
            "target_neutral_mass": round(float(target_mass), 7) if np.isfinite(target_mass) else "",
            "abs_mass_error_da": round(float(abs_error), 7) if np.isfinite(abs_error) else "",
            "mass_match_score": round(float(mass_score), 6),
            "evidence_confidence": round(evidence, 6),
            "forward_spectrum_score": nce_forward,
            "ms2_match_score": nce_ms2,
            "nce_pipeline_score": nce_pipeline,
            "operator_family": operator_family,
            "generation_score": round(generation, 6),
            "novelty_score": 1.0,
            "formula_prior_score": round(float(formula_prior), 6),
            "ms1_buddy_supported": ms1_supported,
            "library_source": hit.get("library_source", ""),
            "library_name": hit.get("library_name", ""),
            "library_match_type": library_match_type,
            "library_smiles_present": True,
            "source_track": hit.get("source_track", ""),
            "modification": hit.get("modification", ""),
            "candidate_source": "library_anchored_candidates_v1",
        }
    )


def build_candidates(
    controls: list[str],
    targets: pd.DataFrame,
    hmdb: pd.DataFrame,
    nce: pd.DataFrame,
    diagnosis: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    """Mass-match library hits per control and expand conjugates from diagnosis hints."""
    rows: list[dict[str, Any]] = []
    hint_map = diagnosis.set_index("control_id").get("conjugate_hints", pd.Series(dtype=str)).astype(str).to_dict()
    for control_id in controls:
        scoped_targets = targets[targets["control_id"].astype(str).eq(str(control_id))].copy()
        direct_nce = nce[nce.get("control_id", pd.Series(dtype=str)).astype(str).eq(str(control_id))].head(
            int(args.max_nce_direct_per_control)
        )
        if scoped_targets.empty and not direct_nce.empty:
            scoped_targets = pd.DataFrame(
                [
                    {
                        "control_id": control_id,
                        "feature_key": direct_nce.iloc[0].get("feature_key", f"library_anchor_{control_id}"),
                        "comp_id": direct_nce.iloc[0].get("comp_id", -1),
                        "lib_id": int(args.lib_id),
                        "ion_mode": str(args.ion_mode),
                        "target_neutral_mass": direct_nce.iloc[0].get(
                            "neutral_mass", direct_nce.iloc[0].get("library_exact_mass", math.nan)
                        ),
                        "adduct_label": "nce_snapshot_control",
                        "formula_prior_score": 0.0,
                        "ms1_buddy_supported": False,
                    }
                ]
            )
        library_hits: list[tuple[pd.Series, str]] = []
        for _, target in scoped_targets.iterrows():
            target_mass = safe_float(target.get("target_neutral_mass"), math.nan)
            if not np.isfinite(target_mass):
                continue
            if not hmdb.empty:
                for _, hit in (
                    hmdb[ppm_hit_mask(hmdb["library_exact_mass"], target_mass, float(args.ppm_tolerance))]
                    .head(int(args.max_hits_per_control))
                    .iterrows()
                ):
                    library_hits.append((hit, "hmdb_mass_match"))
            if not nce.empty:
                for _, hit in (
                    nce[ppm_hit_mask(nce["library_exact_mass"], target_mass, float(args.ppm_tolerance))]
                    .head(int(args.max_hits_per_control))
                    .iterrows()
                ):
                    library_hits.append((hit, "nce_snapshot_mass_match"))
        for _, hit in direct_nce.iterrows():
            library_hits.append((hit, "nce_snapshot_control_match"))
        seen_hits: set[tuple[str, str, str]] = set()
        for _, target in scoped_targets.iterrows():
            for hit, match_type in library_hits:
                base_smiles = str(hit.get("library_smiles", "") or "")
                key = (str(target.get("feature_key", "")), canonicalize_smiles(base_smiles) or base_smiles, match_type)
                if key in seen_hits:
                    continue
                seen_hits.add(key)
                append_candidate(
                    rows,
                    target=target,
                    hit=hit,
                    predicted_smiles=base_smiles,
                    operator_family=f"library_anchor:{match_type}",
                    library_match_type=match_type,
                    args=args,
                )
                hints = [part.strip().lower() for part in str(hint_map.get(control_id, "")).split(",") if part.strip()]
                for hint in hints:
                    for reaction_name in CONJUGATE_REACTIONS.get(hint, tuple()):
                        for product in apply_reaction(base_smiles, reaction_name, int(args.max_products_per_reaction)):
                            append_candidate(
                                rows,
                                target=target,
                                hit=hit,
                                predicted_smiles=product,
                                operator_family=f"library_anchor_conjugate:{reaction_name}",
                                library_match_type=f"{match_type}+{hint}",
                                args=args,
                            )
    candidates = pd.DataFrame(rows)
    if candidates.empty:
        return candidates
    candidates = candidates.drop_duplicates(["feature_key", "predicted_smiles", "operator_family"], keep="first").copy()
    candidates = candidates.sort_values(
        [
            "feature_key",
            "mass_match_score",
            "generation_score",
            "formula_prior_score",
            "nce_pipeline_score",
            "abs_mass_error_da",
            "predicted_smiles",
        ],
        ascending=[True, False, False, False, False, True, True],
        kind="mergesort",
    )
    candidates["rank"] = candidates.groupby("feature_key", dropna=False).cumcount() + 1
    return candidates[candidates["rank"].le(int(args.max_candidates_per_feature))].reset_index(drop=True)


# ── CLI entry point ──────────────────────────────────────────────────
def run(args: argparse.Namespace) -> dict[str, Any]:
    """Execute library-anchored candidate generation and evaluation."""
    repo_root = REPO_ROOT
    output_root = resolve(args.output_root, repo_root)
    output_root.mkdir(parents=True, exist_ok=True)
    diagnosis = pd.read_csv(resolve(args.diagnosis_csv, repo_root), low_memory=False)
    features = pd.read_csv(resolve(args.feature_rows, repo_root), low_memory=False)
    adducts = pd.read_csv(resolve(args.neutral_mass_hypotheses, repo_root), low_memory=False)
    controls = parse_controls(getattr(args, "controls", "")) or target_controls(diagnosis)
    hmdb = prepare_external_libraries(args)
    nce = prepare_nce_snapshot(Path(str(args.nce_full_smiles)).expanduser(), int(args.lib_id))
    targets = mass_hypotheses(features, adducts)
    candidates = build_candidates(controls, targets, hmdb, nce, diagnosis, args)
    out_dir = (
        output_root
        / "runs"
        / RUN_NAME
        / f"track_{args.track}"
        / f"lib_{int(args.lib_id)}_{args.ion_mode}"
        / "selfies_profiles"
        / str(args.profile)
    )
    candidate_path = out_dir / "candidate_metrics_all.csv"
    denovo_path = out_dir / "denovo_report.csv"
    write_csv(candidate_path, candidates)
    ranked = (
        candidates[candidates.get("rank", pd.Series(dtype=int)).le(int(args.top_k))].copy()
        if not candidates.empty
        else candidates
    )
    if not ranked.empty and "_rank" not in ranked.columns:
        ranked["_rank"] = pd.to_numeric(ranked.get("rank"), errors="coerce").fillna(999999).astype(int)
    write_denovo_report(ranked, denovo_path)
    controls_eval = load_x_controls(
        known_x=repo_root / "known_X.csv",
        known_x_all=repo_root / "known_X_all_available_structures_simple.csv",
        lib_id=int(args.lib_id),
        ion_mode=str(args.ion_mode),
    )
    if ranked.empty:
        metrics, per_control = [], []
    else:
        metrics, per_control = eval_reranked(
            candidate_df=candidates,
            ranked_df=ranked,
            controls=controls_eval,
            lib_id=int(args.lib_id),
            ion_mode=str(args.ion_mode),
            profile=str(args.profile),
            policy=RUN_NAME,
            top_k=int(args.top_k),
            candidate_path=candidate_path,
            denovo_path=denovo_path,
        )
    metrics_df = pd.DataFrame(metrics)
    per_df = pd.DataFrame(per_control)
    write_csv(output_root / "library_anchored_metrics.csv", metrics_df)
    write_csv(output_root / "library_anchored_per_control.csv", per_df)
    summary = pd.DataFrame(
        [
            {
                "n_controls_attempted": int(len(controls)),
                "n_controls_with_candidates": int(candidates["control_id"].nunique()) if not candidates.empty else 0,
                "n_candidate_rows": int(candidates.shape[0]),
                "n_library_rows": int(hmdb.shape[0]),
                "n_nce_snapshot_rows": int(nce.shape[0]),
                "candidate_metrics_all": str(candidate_path),
            }
        ]
    )
    by_source = (
        candidates.groupby(["library_source", "library_match_type"], dropna=False)
        .agg(n_controls=("control_id", "nunique"), n_rows=("predicted_smiles", "count"))
        .reset_index()
        if not candidates.empty
        else pd.DataFrame()
    )
    write_csv(output_root / "library_anchored_candidate_summary.csv", summary)
    write_csv(output_root / "library_anchored_source_summary.csv", by_source)
    report = [
        "# Library-Anchored Candidates V1",
        "",
        "## Summary",
        dataframe_to_markdown(summary),
        "",
        "## Source Summary",
        dataframe_to_markdown(by_source),
        "",
        "## Metrics",
        dataframe_to_markdown(metrics_df),
        "",
    ]
    (output_root / "library_anchored_candidate_report.md").write_text("\n".join(report))
    payload = summary.iloc[0].to_dict()
    write_json(output_root / "library_anchored_candidate_summary.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    """Parse CLI arguments and invoke :func:`run`."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=Path("models/library_anchored_candidates_v1"))
    parser.add_argument(
        "--diagnosis-csv", type=Path, default=Path("models/denovo_dashboard_v1/diagnosis/control_stage_diagnosis.csv")
    )
    parser.add_argument(
        "--feature-rows", type=Path, default=Path("models/peptide_first_denovo_v1_1/feature_peptide_rows.csv")
    )
    parser.add_argument(
        "--neutral-mass-hypotheses",
        type=Path,
        default=Path("models/upstream_adduct_evidence_v1/neutral_mass_hypotheses.csv"),
    )
    parser.add_argument("--hmdb-csv", type=Path, default=DEFAULT_HMDB)
    parser.add_argument("--hmdb-extended-csv", type=Path, default=DEFAULT_HMDB_EXTENDED)
    parser.add_argument("--chebi-csv", type=Path, default=Path(""))
    parser.add_argument("--lipidmaps-csv", type=Path, default=Path(""))
    parser.add_argument("--drugbank-csv", type=Path, default=Path(""))
    parser.add_argument("--masst-structures-csv", type=Path, default=Path(""))
    parser.add_argument("--nce-full-smiles", type=Path, default=DEFAULT_NCE_FULL)
    parser.add_argument("--lib-id", type=int, default=209)
    parser.add_argument("--ion-mode", default="negative")
    parser.add_argument("--controls", default="", help="Optional comma-separated X control IDs to process.")
    parser.add_argument("--profile", default="strict_1da")
    parser.add_argument("--track", default="control_blind")
    parser.add_argument("--ppm-tolerance", type=float, default=5.0)
    parser.add_argument("--mass-score-tolerance-da", type=float, default=1.0)
    parser.add_argument("--formula-tolerance-da", type=float, default=0.05)
    parser.add_argument("--max-hits-per-control", type=int, default=120)
    parser.add_argument("--max-nce-direct-per-control", type=int, default=80)
    parser.add_argument("--max-products-per-reaction", type=int, default=3)
    parser.add_argument("--max-candidates-per-feature", type=int, default=250)
    parser.add_argument("--top-k", type=int, default=10)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
