#!/usr/bin/env python3
"""Score existing candidate pools with the real NCE predicted-spectrum verifier.

Runs ``annotate_predicted_spectrum_scores_v2`` across all enumerator outputs,
writes per-source ``*_with_verifier.csv`` files, and concatenates a combined
scored pool for downstream merging.

Pipeline context:
    Bridges enumerator ``candidate_metrics_all.csv`` artifacts to ``pool``.
    Default inputs cover lib 209 negative; other scopes use
    ``default_inputs_for_scope``.

Key I/O:
    Inputs: observed MS2 lookup, comma-separated candidate_metrics CSVs.
    Outputs: ``combined_scored_candidates.csv``, per-input verifier CSVs,
    manifest, winner-rank diagnostic, summary JSON/markdown.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import pandas as pd

from structural_elucidation.utils.denovo import canonicalize_smiles
from structural_elucidation.nce.score_cache import VerifierScoreCache
from structural_elucidation.nce.verifier_v2 import annotate_predicted_spectrum_scores_v2
from structural_elucidation.paths import REPO_ROOT, resolve
from structural_elucidation.utils.domain import dataframe_to_markdown, write_csv, write_json
from structural_elucidation.utils.timing import stage_timer

# ── Default verifier input paths (lib 209 negative) ──────────────────
DEFAULT_INPUTS: tuple[str, ...] = (
    "models/peptide_first_denovo_v1_1/rerank/runs/peptide_first_rerank/lib_209_negative/selfies_profiles/protected_union_control_blind/candidate_metrics_all.csv",
    "models/peptide_first_denovo_v1_1/rerank/runs/peptide_first_rerank/lib_209_negative/selfies_profiles/protected_union_annotation_assisted/candidate_metrics_all.csv",
    "models/peptide_first_denovo_v1_1/rerank_verifier/runs/peptide_first_rerank/lib_209_negative/selfies_profiles/protected_union_control_blind/candidate_metrics_all.csv",
    "models/peptide_first_denovo_v1_1/rerank_verifier/runs/peptide_first_rerank/lib_209_negative/selfies_profiles/protected_union_annotation_assisted/candidate_metrics_all.csv",
    "models/library_anchored_candidates_v1/runs/library_anchored_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/nonpeptide_gated_v1/runs/nonpeptide_gated_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/multi_conjugate_candidates_v1/runs/multi_conjugate_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/biotransformation_gated_v2/runs/biotransformation_gated_v2/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/biotransformation_gated_v2/runs/biotransformation_gated_v2/track_annotation_assisted/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/nce_direct_runs_v1/nce_direct_top10.csv",
    "models/mlims_anchored_candidates_v1/runs/mlims_anchored_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/mlims_strict_overlay_v1/candidate_metrics_all.csv",
    "models/targeted_combinatorial_v1/candidate_metrics_all.csv",
    "models/mass_aware_combinatorial_v1/candidate_metrics_all.csv",
    "models/lactoyl_amino_acid_candidates_v1/runs/lactoyl_amino_acid_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/sulfate_position_candidates_v1/runs/sulfate_position_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/taurine_variant_candidates_v1/runs/taurine_variant_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/n_acetyl_candidates_v1/runs/n_acetyl_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/glucuronide_candidates_v1/runs/glucuronide_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/glycine_candidates_v1/runs/glycine_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/glutathione_candidates_v1/runs/glutathione_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/cysteinyl_candidates_v1/runs/cysteinyl_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/hydroxyl_candidates_v1/runs/hydroxyl_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/dehydro_candidates_v1/runs/dehydro_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/demethyl_candidates_v1/runs/demethyl_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/methylation_candidates_v1/runs/methylation_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/wgcna_seeded_candidates_v1/lib_209_negative/candidate_metrics_all.csv",
)


def default_inputs_for_scope(lib_id: int, ion_mode: str) -> tuple[str, ...]:
    """Return default candidate_metrics paths for non-209-negative scopes."""
    if (int(lib_id), str(ion_mode)) == (209, "negative"):
        return DEFAULT_INPUTS
    scope = f"lib_{int(lib_id)}_{str(ion_mode)}"
    return (
        f"models/ec2_heavy_library_anchored_v1_lib{int(lib_id)}/runs/library_anchored_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/ec2_heavy_mlims_anchored_v1_lib{int(lib_id)}/runs/mlims_anchored_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/ec2_heavy_mlims_strict_overlay_v1_lib{int(lib_id)}/candidate_metrics_all.csv",
        f"models/ec2_heavy_targeted_combinatorial_v1_lib{int(lib_id)}/candidate_metrics_all.csv",
        f"models/ec2_heavy_mass_aware_combinatorial_v1_lib{int(lib_id)}/candidate_metrics_all.csv",
        f"models/ec2_heavy_multi_conjugate_v1_lib{int(lib_id)}/runs/multi_conjugate_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/ec2_heavy_nonpeptide_gated_v1_lib{int(lib_id)}/runs/nonpeptide_gated_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/ec2_heavy_biotransformation_gated_v2_lib{int(lib_id)}/runs/biotransformation_gated_v2/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/ec2_heavy_biotransformation_gated_v2_lib{int(lib_id)}/runs/biotransformation_gated_v2/track_annotation_assisted/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/integrated_denovo_audit_v1/ec2_bundle_cross_lib_v1/{scope}/candidate_metrics_all.csv",
        f"models/wgcna_seeded_candidates_v1/{scope}/candidate_metrics_all.csv",
        f"models/lactoyl_amino_acid_candidates_v1/{scope}/runs/lactoyl_amino_acid_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/sulfate_position_candidates_v1/{scope}/runs/sulfate_position_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/taurine_variant_candidates_v1/{scope}/runs/taurine_variant_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/ec2_heavy_peptide_first_v1_lib{int(lib_id)}/runs/peptide_first_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/nucleoside_candidates_v1/{scope}/runs/nucleoside_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/carnitine_acyl_candidates_v1/{scope}/runs/carnitine_acyl_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/glutathione_candidates_v1/{scope}/direct_conjugate_v18/runs/glutathione_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/n_acetyl_candidates_v1/{scope}/direct_conjugate_v18/runs/n_acetyl_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
        f"models/glycine_candidates_v1/{scope}/direct_conjugate_v18/runs/glycine_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    )


def parse_inputs(text: str, lib_id: int = 209, ion_mode: str = "negative") -> list[Path]:
    """Parse comma-separated input paths, deduplicating while preserving order."""
    if not str(text or "").strip():
        values = default_inputs_for_scope(int(lib_id), str(ion_mode))
    else:
        values = [part.strip() for part in str(text).split(",") if part.strip()]
    paths: list[Path] = []
    seen: set[str] = set()
    for value in values:
        path = Path(value)
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        paths.append(path)
    return paths


def s0_controls_for_scope(path: Path, lib_id: int) -> set[str]:
    """Return S0-tier X control IDs for splitting combined scored output."""
    if not path.exists():
        return set()
    df = pd.read_csv(path, low_memory=False)
    if "lib_id" in df.columns:
        df = df[pd.to_numeric(df["lib_id"], errors="coerce").eq(int(lib_id))].copy()
    if "evidence_tier" not in df.columns or "x_compound_id" not in df.columns:
        return set()
    return set(df[df["evidence_tier"].astype(str).str.contains("S0", na=False)]["x_compound_id"].astype(str))


def merge_all_x_observed_lookup(
    observed: pd.DataFrame, extra_path: Path, lib_id: int, ion_mode: str
) -> tuple[pd.DataFrame, int]:
    """Extend an S0-only observed lookup with all-X spectra for include-all-x runs.

    Returns:
        Merged observed frame and count of newly added controls.
    """
    if observed.empty or not extra_path.exists():
        return observed, 0
    try:
        extra = pd.read_csv(extra_path, low_memory=False)
    except (pd.errors.EmptyDataError, OSError):
        return observed, 0
    if extra.empty or "control_id" not in extra.columns:
        return observed, 0
    scoped = extra.copy()
    if "lib_id" in scoped.columns:
        scoped = scoped[pd.to_numeric(scoped["lib_id"], errors="coerce").eq(int(lib_id))].copy()
    if "ion_mode" in scoped.columns:
        scoped = scoped[scoped["ion_mode"].astype(str).eq(str(ion_mode))].copy()
    if scoped.empty:
        return observed, 0
    existing = set(observed.get("control_id", pd.Series(dtype=str)).dropna().astype(str))
    scoped = scoped[~scoped["control_id"].astype(str).isin(existing)].copy()
    if scoped.empty:
        return observed, 0
    merged = pd.concat([observed, scoped], ignore_index=True, sort=False)
    return merged, int(scoped["control_id"].nunique())


def _boolish_series(series: pd.Series) -> pd.Series:
    """Parse common truthy string values in observed-lookup columns."""
    return series.astype(str).str.lower().isin({"true", "1", "yes", "y"})


def prepare_observed_lookup_for_scope(
    observed: pd.DataFrame, lib_id: int, ion_mode: str
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Scope, deduplicate, and rank observed MS2 rows per control for verifier join.

    When multiple observed rows exist per control, prefer richer MS2 evidence.
    """
    stats = {
        "n_observed_input_rows": int(observed.shape[0]),
        "n_observed_scoped_rows": 0,
        "n_observed_duplicate_controls": 0,
        "n_observed_duplicate_rows_dropped": 0,
    }
    if observed.empty or "control_id" not in observed.columns:
        return observed.copy(), stats

    work = observed.copy()
    work["control_id"] = work["control_id"].fillna("").astype(str)
    work = work[work["control_id"].ne("")].copy()
    if "lib_id" in work.columns:
        work = work[pd.to_numeric(work["lib_id"], errors="coerce").eq(int(lib_id))].copy()
    if "ion_mode" in work.columns:
        work = work[work["ion_mode"].fillna("").astype(str).eq(str(ion_mode))].copy()
    if "ms2_present" in work.columns:
        work = work[_boolish_series(work["ms2_present"])].copy()

    stats["n_observed_scoped_rows"] = int(work.shape[0])
    if work.empty:
        return work, stats

    counts = work["control_id"].value_counts()
    stats["n_observed_duplicate_controls"] = int(counts.gt(1).sum())
    work["_observed_input_order"] = range(len(work))
    work["_observed_n_ms2_pairs"] = pd.to_numeric(
        work.get("n_ms2_pairs", pd.Series(-1, index=work.index)),
        errors="coerce",
    ).fillna(-1)
    ion_bins = pd.to_numeric(work.get("n_observed_ion_bins", pd.Series(0, index=work.index)), errors="coerce").fillna(0)
    loss_bins = pd.to_numeric(work.get("n_observed_loss_bins", pd.Series(0, index=work.index)), errors="coerce").fillna(
        0
    )
    work["_observed_n_bins"] = ion_bins + loss_bins
    if "spectrum_source" in work.columns:
        source_order = {
            str(value): idx
            for idx, value in enumerate(work["spectrum_source"].fillna("").astype(str).drop_duplicates())
        }
        work["_observed_source_order"] = (
            work["spectrum_source"].fillna("").astype(str).map(source_order).fillna(len(source_order))
        )
    else:
        work["_observed_source_order"] = 0
    work = (
        work.sort_values(
            [
                "control_id",
                "_observed_n_ms2_pairs",
                "_observed_n_bins",
                "_observed_source_order",
                "_observed_input_order",
            ],
            ascending=[True, False, False, True, True],
            kind="mergesort",
        )
        .drop_duplicates("control_id", keep="first")
        .drop(columns=["_observed_input_order", "_observed_n_ms2_pairs", "_observed_n_bins", "_observed_source_order"])
        .reset_index(drop=True)
    )
    stats["n_observed_duplicate_rows_dropped"] = int(stats["n_observed_scoped_rows"] - work.shape[0])
    return work, stats


def source_label_for_path(path: Path) -> str:
    """Map an input CSV path to a stable reranker source_label string."""
    text = str(path)
    if "merged_structuremasst_bundle" in text and "selfies_profiles" in text:
        return "merged_structuremasst_selfies_v1"
    if "merged_structuremasst_bundle" in text and "/retrieval/" in text:
        return "merged_structuremasst_retrieval_v1"
    if "ec2_bundle_cross_lib_v1" in text:
        return "ec2_bundle_retrieval_selfies_v1"
    if "nce_direct_runs_v1" in text:
        return "nce_direct_v1"
    if "library_anchored_candidates_v1" in text or "ec2_heavy_library_anchored_v1" in text:
        return "library_anchored_v1"
    if "nonpeptide_gated_v1" in text:
        return "nonpeptide_gated_v1"
    if "multi_conjugate_candidates_v1" in text or "ec2_heavy_multi_conjugate_v1" in text:
        return "multi_conjugate_v2_smarts"
    if "biotransformation_gated_v2_multi_conjugate" in text:
        return "biotransformation_gated_v2_multi_conjugate"
    if "biotransformation_gated_v2" in text:
        return "biotransformation_gated_v2"
    if "mlims_anchored_candidates_v1" in text or "ec2_heavy_mlims_anchored_v1" in text:
        return "mlims_anchored_v1"
    if "mlims_strict_overlay_v1" in text:
        return "mlims_strict_overlay_v1"
    if "targeted_combinatorial_v1" in text:
        return "targeted_combinatorial_v1"
    if "mass_aware_combinatorial_v1" in text:
        return "mass_aware_combinatorial_v1"
    if "wgcna_seeded_candidates_v1" in text:
        return "wgcna_seeded_v1"
    if "sulfate_position_candidates_v" in text:
        return "sulfate_position_v1"
    if "taurine_variant_candidates_v" in text:
        return "taurine_variant_v1"
    if "lactoyl_amino_acid_candidates_v" in text:
        return "lactoyl_amino_acid_v1"
    if "n_acetyl_candidates_v" in text:
        return "n_acetyl_v1"
    if "glucuronide_candidates_v" in text:
        return "glucuronide_v1"
    if "glycine_candidates_v" in text:
        return "glycine_v1"
    if "glutathione_candidates_v" in text:
        return "glutathione_v1"
    if "cysteinyl_candidates_v" in text:
        return "cysteinyl_v1"
    if "hydroxyl_candidates_v" in text:
        return "hydroxyl_v1"
    if "dehydro_candidates_v" in text:
        return "dehydro_v1"
    if "demethyl_candidates_v" in text:
        return "demethyl_v1"
    if "methylation_candidates_v" in text:
        return "methylation_v1"
    if "rerank_verifier" in text:
        return "peptide_first_v1_1_verifier"
    if "peptide_first_denovo_v1_1" in text:
        return "peptide_first_v1_1"
    return path.parent.name


def normalize_frame(df: pd.DataFrame, source_label: str, lib_id: int = 209, ion_mode: str = "negative") -> pd.DataFrame:
    """Normalize candidate_metrics schema before verifier scoring."""
    out = df.copy()
    is_unified_pool = "source_labels" in out.columns and (
        "unified_pool_rank" in out.columns or source_label.startswith("unified_candidate_pool")
    )
    if is_unified_pool and "predicted_spectrum_score_v2" in out.columns:
        out["_input_predicted_spectrum_score_v2"] = pd.to_numeric(out["predicted_spectrum_score_v2"], errors="coerce")
    if "predicted_smiles" not in out.columns and "smiles" in out.columns:
        out["predicted_smiles"] = out["smiles"]
    if "control_id" not in out.columns:
        out["control_id"] = (
            out.get("feature_key", pd.Series("", index=out.index))
            .astype(str)
            .str.extract(r"(X-\d+)", expand=False)
            .fillna("")
        )
    if "feature_key" not in out.columns:
        out["feature_key"] = out["control_id"].astype(str).map(lambda value: f"{source_label}_{value}")
    if "rank" not in out.columns:
        out["rank"] = out.groupby("control_id", dropna=False).cumcount() + 1
    if "lib_id" not in out.columns:
        out["lib_id"] = int(lib_id)
    if "ion_mode" not in out.columns:
        out["ion_mode"] = str(ion_mode)
    else:
        out["ion_mode"] = out["ion_mode"].fillna(str(ion_mode)).astype(str)
    if is_unified_pool:
        out["source_label"] = (
            out["source_labels"].fillna(out.get("source_label", pd.Series("", index=out.index))).astype(str)
        )
        out["candidate_source"] = out.get("candidate_source", out["source_label"]).fillna(out["source_label"])
    else:
        out["source_label"] = source_label
        out["candidate_source"] = out.get("candidate_source", pd.Series(source_label, index=out.index)).fillna(
            source_label
        )
    return out


def trim_for_scoring(df: pd.DataFrame, max_rank: int) -> pd.DataFrame:
    """Limit rows per source to top *max_rank* candidates before verifier cost."""
    if df.empty or int(max_rank) <= 0:
        return df
    out = df.copy()
    out["_rank_numeric"] = pd.to_numeric(out.get("rank"), errors="coerce").fillna(999999)
    return out[out["_rank_numeric"].le(int(max_rank))].drop(columns=["_rank_numeric"]).copy()


# ── CLI entry point ──────────────────────────────────────────────────
def run(args: argparse.Namespace) -> dict[str, Any]:
    """Score all configured inputs and write combined verifier outputs."""
    repo_root = REPO_ROOT
    output_root = resolve(args.output_root, repo_root)
    output_root.mkdir(parents=True, exist_ok=True)
    observed = pd.read_csv(resolve(args.observed_lookup, repo_root), low_memory=False)
    n_all_x_observed_added = 0
    if bool(args.include_all_x):
        observed, n_all_x_observed_added = merge_all_x_observed_lookup(
            observed,
            resolve(args.all_x_observed_lookup, repo_root),
            int(args.lib_id),
            str(args.ion_mode),
        )
        if bool(getattr(args, "verbose", False)):
            print(
                json.dumps(
                    {
                        "event": "observed_lookup_loaded",
                        "n_observed_controls": (
                            int(observed["control_id"].nunique()) if "control_id" in observed.columns else 0
                        ),
                        "n_all_x_observed_added": int(n_all_x_observed_added),
                    },
                    sort_keys=True,
                )
            )
    observed, observed_stats = prepare_observed_lookup_for_scope(observed, int(args.lib_id), str(args.ion_mode))
    workers = max(1, int(getattr(args, "workers", 1) or 1))
    score_cache: VerifierScoreCache | None = None
    if not bool(getattr(args, "no_score_cache", False)):
        cache_path = getattr(args, "score_cache", None)
        if cache_path is None:
            cache_path = Path("cache/verifier_scores") / f"lib{int(args.lib_id)}_{args.ion_mode}_scores.csv"
        score_cache = VerifierScoreCache(resolve(Path(cache_path), repo_root))
    n_scores_preloaded = len(score_cache) if score_cache is not None else 0
    manifest_rows: list[dict[str, Any]] = []
    scored_frames: list[pd.DataFrame] = []
    for rel in parse_inputs(args.input_candidate_metrics, args.lib_id, args.ion_mode):
        path = resolve(rel, repo_root)
        if not path.exists():
            manifest_rows.append(
                {"input_csv": str(path), "status": "missing", "n_rows_in": 0, "n_rows_scored": 0, "output_csv": ""}
            )
            continue
        source_label = source_label_for_path(path)
        try:
            raw = pd.read_csv(path, low_memory=False)
        except pd.errors.EmptyDataError:
            manifest_rows.append(
                {
                    "input_csv": str(path),
                    "source_label": source_label,
                    "status": "empty",
                    "n_rows_in": 0,
                    "n_rows_scored": 0,
                    "output_csv": "",
                }
            )
            continue
        if raw.empty:
            manifest_rows.append(
                {
                    "input_csv": str(path),
                    "source_label": source_label,
                    "status": "empty",
                    "n_rows_in": 0,
                    "n_rows_scored": 0,
                    "output_csv": "",
                }
            )
            continue
        norm = normalize_frame(raw, source_label, lib_id=int(args.lib_id), ion_mode=str(args.ion_mode))
        work = trim_for_scoring(norm, int(args.max_rank_per_source))
        if work.empty or "predicted_smiles" not in work.columns:
            manifest_rows.append(
                {
                    "input_csv": str(path),
                    "source_label": source_label,
                    "status": "empty",
                    "n_rows_in": int(raw.shape[0]),
                    "n_rows_scored": 0,
                    "output_csv": "",
                }
            )
            continue
        scored = annotate_predicted_spectrum_scores_v2(work, observed, workers=workers, score_cache=score_cache)
        if score_cache is not None:
            score_cache.save()  # persist incrementally so interrupted runs keep progress
        if scored.empty or "predicted_smiles" not in scored.columns:
            manifest_rows.append(
                {
                    "input_csv": str(path),
                    "source_label": source_label,
                    "status": "empty",
                    "n_rows_in": int(raw.shape[0]),
                    "n_rows_scored": 0,
                    "output_csv": "",
                }
            )
            continue
        if "_input_predicted_spectrum_score_v2" in scored.columns:
            rescored = pd.to_numeric(scored["predicted_spectrum_score_v2"], errors="coerce").fillna(0.0)
            prior = pd.to_numeric(scored["_input_predicted_spectrum_score_v2"], errors="coerce").fillna(0.0)
            scored["predicted_spectrum_score_v2_rescored"] = rescored
            # Keep the better of pre-verifier and rescored scores for unified pools.
            scored["predicted_spectrum_score_v2"] = pd.concat([prior, rescored], axis=1).max(axis=1)
        scored["canonical_2d_smiles"] = scored["predicted_smiles"].map(lambda value: canonicalize_smiles(value) or "")
        out_path = path.with_name(path.stem + "_with_verifier.csv")
        write_csv(out_path, scored)
        scored["verifier_scored_csv"] = str(out_path)
        scored_frames.append(scored)
        manifest_rows.append(
            {
                "input_csv": str(path),
                "source_label": source_label,
                "status": "scored",
                "n_rows_in": int(raw.shape[0]),
                "n_rows_scored": int(scored.shape[0]),
                "output_csv": str(out_path),
            }
        )
    manifest = pd.DataFrame(manifest_rows)
    write_csv(output_root / "verifier_scoring_manifest.csv", manifest)
    combined = pd.concat(scored_frames, ignore_index=True, sort=False) if scored_frames else pd.DataFrame()
    write_csv(output_root / "combined_scored_candidates.csv", combined)
    write_csv(output_root / "scored_candidates.csv", combined)
    s0_controls = s0_controls_for_scope(resolve(args.all_x_catalog, repo_root), int(args.lib_id))
    if not combined.empty and s0_controls:
        write_csv(
            output_root / "combined_scored_candidates_s0.csv",
            combined[combined["control_id"].astype(str).isin(s0_controls)].copy(),
        )
        write_csv(
            output_root / "combined_scored_candidates_non_s0.csv",
            combined[~combined["control_id"].astype(str).isin(s0_controls)].copy(),
        )
    else:
        write_csv(output_root / "combined_scored_candidates_s0.csv", pd.DataFrame(columns=combined.columns))
        write_csv(output_root / "combined_scored_candidates_non_s0.csv", pd.DataFrame(columns=combined.columns))

    current = pd.read_csv(resolve(args.current_integrated_best, repo_root), low_memory=False)
    winners = current[pd.to_numeric(current.get("top10_best_tanimoto"), errors="coerce").fillna(0.0).ge(0.7)].copy()
    rank_rows: list[dict[str, Any]] = []
    if not combined.empty:
        combined["_verifier_rank"] = (
            combined.sort_values(
                ["control_id", "predicted_spectrum_score_v2"], ascending=[True, False], kind="mergesort"
            )
            .groupby("control_id", dropna=False)
            .cumcount()
            + 1
        )
    for _, row in winners.iterrows():
        control_id = str(row.get("control_id", ""))
        winner_smiles = canonicalize_smiles(row.get("top10_smiles", "")) or ""
        scoped = combined[
            (combined.get("control_id", pd.Series(dtype=str)).astype(str).eq(control_id))
            & (combined.get("canonical_2d_smiles", pd.Series(dtype=str)).astype(str).eq(winner_smiles))
        ]
        rank_rows.append(
            {
                "control_id": control_id,
                "source_label": row.get("source_label", ""),
                "winner_smiles": row.get("top10_smiles", ""),
                "winner_found_in_scored_pool": bool(not scoped.empty),
                "best_verifier_score": (
                    float(pd.to_numeric(scoped.get("predicted_spectrum_score_v2"), errors="coerce").fillna(0.0).max())
                    if not scoped.empty
                    else 0.0
                ),
                "best_verifier_rank": (
                    int(pd.to_numeric(scoped.get("_verifier_rank"), errors="coerce").fillna(999999).min())
                    if not scoped.empty
                    else 999999
                ),
            }
        )
    rank_diag = pd.DataFrame(rank_rows)
    write_csv(output_root / "current_winner_verifier_rank_diagnostic.csv", rank_diag)
    if not rank_diag.empty:
        rank_values = pd.to_numeric(rank_diag.get("best_verifier_rank"), errors="coerce")
        rank_values = rank_values[rank_values.lt(999999)]
        median_rank = float(rank_values.median()) if not rank_values.empty else 999999.0
    else:
        median_rank = 999999.0
    summary = {
        "n_inputs": int(manifest.shape[0]),
        "n_scored_inputs": int(manifest["status"].eq("scored").sum()) if not manifest.empty else 0,
        "n_rows_scored": int(sum(row.get("n_rows_scored", 0) for row in manifest_rows)),
        "verifier_workers": int(workers),
        "score_cache_path": str(score_cache.path) if score_cache is not None else "",
        "n_score_cache_preloaded": int(n_scores_preloaded),
        "n_score_cache_entries": int(len(score_cache)) if score_cache is not None else 0,
        "include_all_x": bool(args.include_all_x),
        "n_observed_controls": int(observed["control_id"].nunique()) if "control_id" in observed.columns else 0,
        "n_all_x_observed_added": int(n_all_x_observed_added),
        **observed_stats,
        "n_current_winners": int(rank_diag.shape[0]),
        "n_current_winners_found": int(rank_diag["winner_found_in_scored_pool"].sum()) if not rank_diag.empty else 0,
        "median_current_winner_verifier_rank": median_rank,
        "acceptance_pass": (
            bool(median_rank <= float(args.max_acceptance_median_rank)) if median_rank != 999999.0 else False
        ),  # median verifier rank of current winners must stay low
    }
    write_json(output_root / "verifier_scoring_summary.json", summary)
    report = [
        "# Candidate Verifier Scoring V1",
        "",
        "## Summary",
        dataframe_to_markdown(pd.DataFrame([summary])),
        "",
        "## Manifest",
        dataframe_to_markdown(manifest),
        "",
        "## Current Winner Verifier Ranks",
        dataframe_to_markdown(
            rank_diag.sort_values("best_verifier_rank").head(40) if not rank_diag.empty else rank_diag
        ),
        "",
    ]
    (output_root / "verifier_scoring_report.md").write_text("\n".join(report))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def main() -> None:
    """Parse CLI arguments and invoke :func:`run`; exit 2 on acceptance failure."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=Path("models/verifier_scoring_v1"))
    parser.add_argument("--observed-lookup", type=Path, default=Path("models/observed_ms2_lookup_v1/per_control.csv"))
    parser.add_argument(
        "--current-integrated-best",
        type=Path,
        default=Path("models/integrated_denovo_audit_v1/integrated_per_control_best.csv"),
    )
    parser.add_argument("--input-candidate-metrics", default="")
    parser.add_argument("--lib-id", type=int, default=209)
    parser.add_argument("--ion-mode", default="negative")
    parser.add_argument("--max-rank-per-source", type=int, default=50)
    parser.add_argument("--max-acceptance-median-rank", type=float, default=5.0)
    parser.add_argument("--include-all-x", action="store_true")
    parser.add_argument(
        "--all-x-catalog", type=Path, default=Path("models/x_known_curation_v1/known_x_evidence_catalog.csv")
    )
    parser.add_argument(
        "--all-x-observed-lookup", type=Path, default=Path("models/observed_ms2_lookup_v8/per_control_all_x.csv")
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("SE_VERIFIER_WORKERS", "1") or 1),
        help="Process count for NCE spectrum scoring (default: SE_VERIFIER_WORKERS env or 1; results are identical at any count).",
    )
    parser.add_argument(
        "--score-cache",
        type=Path,
        default=None,
        help="Persistent verifier score cache CSV (default: cache/verifier_scores/lib{lib}_{ion}_scores.csv).",
    )
    parser.add_argument("--no-score-cache", action="store_true", help="Disable the persistent score cache.")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    with stage_timer("verify", lib_id=int(args.lib_id), ion_mode=str(args.ion_mode)):
        payload = run(args)
    raise SystemExit(0 if payload["acceptance_pass"] else 2)


if __name__ == "__main__":
    main()
