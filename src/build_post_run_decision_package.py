#!/usr/bin/env python3
"""Build a post-run quality decision package for merged retrieval/SELFIES runs.

The report separates selected-output quality from full candidate-pool oracle
quality so we can tell whether SELFIES is failing to generate useful structures
or mostly failing to rank/select them.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from benchmark_selfies_merge_policies import rank_selfies_pool
from evaluate_denovo_controls import (
    HAS_RDKIT,
    canonicalize,
    evaluate_controls,
    extract_xid,
    load_controls,
    load_predictions,
    tanimoto,
)


PROFILE_ORDER = ["strict_1da", "balanced", "coverage_explore"]
THRESHOLDS = [0.5, 0.7, 0.9]
METRIC_COLUMNS = [
    "lib_id",
    "ion_mode",
    "source_type",
    "source_name",
    "lane",
    "policy",
    "metric_scope",
    "n_controls",
    "n_controls_covered",
    "n_exact_match",
    "exact_match_rate",
    "n_tanimoto_ge_0p5",
    "n_tanimoto_ge_0p7",
    "n_tanimoto_ge_0p9",
    "tanimoto_ge_0p7_rate",
    "mean_best_tanimoto",
    "median_best_tanimoto",
    "path",
]


@dataclass(frozen=True)
class LibContext:
    lib_dir: Path
    lib_id: int
    ion_mode: str
    manifest: Dict[str, Any]


def read_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    """Render a compact Markdown table without pandas' optional tabulate dependency."""
    if df.empty:
        return ""
    cols = [str(col) for col in df.columns]

    def clean(value: Any) -> str:
        if pd.isna(value):
            return ""
        return str(value).replace("\n", " ").replace("|", "\\|")

    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for row in df.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(clean(value) for value in row) + " |")
    return "\n".join(lines)


def existing_path(path_like: Any, run_root: Path) -> Optional[Path]:
    if not path_like:
        return None
    raw = str(path_like)
    path = Path(raw)
    if path.exists():
        return path
    if not path.is_absolute():
        candidate = run_root / path
        if candidate.exists():
            return candidate
    marker = run_root.name
    parts = path.parts
    if marker in parts:
        rel = Path(*parts[parts.index(marker) + 1 :])
        candidate = run_root / rel
        if candidate.exists():
            return candidate
    return None


def discover_libs(run_root: Path) -> List[LibContext]:
    contexts: List[LibContext] = []
    for lib_dir in sorted(run_root.glob("lib_*_*")):
        manifest_path = lib_dir / "merged_dual_lane_manifest.json"
        if not manifest_path.exists():
            continue
        manifest = read_json(manifest_path)
        lib_id = int(manifest.get("lib_id") or str(lib_dir.name).split("_")[1])
        ion_mode = str(manifest.get("ion_mode") or str(lib_dir.name).split("_", 2)[2])
        contexts.append(LibContext(lib_dir=lib_dir, lib_id=lib_id, ion_mode=ion_mode, manifest=manifest))
    return contexts


def metric_row(
    *,
    ctx: LibContext,
    source_type: str,
    source_name: str,
    metric_scope: str,
    summary: Dict[str, Any],
    path: Optional[Path] = None,
    lane: str = "",
    policy: str = "",
) -> Dict[str, Any]:
    return {
        "lib_id": int(ctx.lib_id),
        "ion_mode": ctx.ion_mode,
        "source_type": source_type,
        "source_name": source_name,
        "lane": lane,
        "policy": policy,
        "metric_scope": metric_scope,
        "n_controls": int(summary.get("n_controls", 0) or 0),
        "n_controls_covered": int(summary.get("n_controls_covered", 0) or 0),
        "n_exact_match": int(summary.get("n_exact_match", 0) or 0),
        "exact_match_rate": float(summary.get("exact_match_rate", 0.0) or 0.0),
        "n_tanimoto_ge_0p5": int(
            summary.get("n_topk_tanimoto_ge_0.5", summary.get("n_oracle_tanimoto_ge_0.5", 0)) or 0
        ),
        "n_tanimoto_ge_0p7": int(
            summary.get("n_topk_tanimoto_ge_0.7", summary.get("n_oracle_tanimoto_ge_0.7", 0)) or 0
        ),
        "n_tanimoto_ge_0p9": int(
            summary.get("n_topk_tanimoto_ge_0.9", summary.get("n_oracle_tanimoto_ge_0.9", 0)) or 0
        ),
        "tanimoto_ge_0p7_rate": float(
            summary.get("topk_tanimoto_ge_0.7_rate", summary.get("oracle_tanimoto_ge_0.7_rate", 0.0)) or 0.0
        ),
        "mean_best_tanimoto": float(summary.get("mean_best_tanimoto", 0.0) or 0.0),
        "median_best_tanimoto": float(summary.get("median_best_tanimoto", 0.0) or 0.0),
        "path": str(path) if path is not None else "",
    }


def evaluate_prediction_path(
    *,
    ctx: LibContext,
    controls: Dict[str, List[str]],
    pred_path: Path,
    source_type: str,
    source_name: str,
    top_k: int,
    lane: str = "",
    policy: str = "",
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    preds = load_predictions(pred_path)
    rows: List[Dict[str, Any]] = []
    per_control_rows: List[Dict[str, Any]] = []
    for scope_name, k in [("rank1", 1), (f"top{top_k}", top_k)]:
        per_control, summary = evaluate_controls(
            controls,
            preds,
            track_name=f"{source_name}_{scope_name}",
            top_k=int(k),
        )
        rows.append(
            metric_row(
                ctx=ctx,
                source_type=source_type,
                source_name=source_name,
                metric_scope=scope_name,
                summary=summary,
                path=pred_path,
                lane=lane,
                policy=policy,
            )
        )
        for item in per_control:
            item = dict(item)
            item.update(
                {
                    "lib_id": int(ctx.lib_id),
                    "ion_mode": ctx.ion_mode,
                    "source_type": source_type,
                    "source_name": source_name,
                    "lane": lane,
                    "policy": policy,
                    "metric_scope": scope_name,
                    "path": str(pred_path),
                }
            )
            per_control_rows.append(item)
    return rows, per_control_rows


def numeric_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(df[col], errors="coerce")


def mass_coverage_row(ctx: LibContext, profile: str, profile_dir: Path, candidate_df: pd.DataFrame) -> Dict[str, Any]:
    abs_err = numeric_series(candidate_df, "abs_mass_error_da")
    feature_keys = candidate_df["feature_key"].astype(str) if "feature_key" in candidate_df.columns else pd.Series(dtype=str)
    row: Dict[str, Any] = {
        "lib_id": int(ctx.lib_id),
        "ion_mode": ctx.ion_mode,
        "profile": profile,
        "candidate_rows": int(candidate_df.shape[0]),
        "candidate_features": int(feature_keys.nunique()) if not feature_keys.empty else 0,
        "path": str(profile_dir / "candidate_metrics_all.csv"),
    }
    denom = max(int(candidate_df.shape[0]), 1)
    for window in [0.5, 1.0, 2.0, 5.0]:
        key = str(window).replace(".", "p")
        mask = abs_err.le(float(window)) if not abs_err.empty else pd.Series(dtype=bool)
        row[f"within_{key}_da"] = int(mask.sum()) if not mask.empty else 0
        row[f"within_{key}_da_rate"] = float(row[f"within_{key}_da"] / denom)
        row[f"features_within_{key}_da"] = (
            int(candidate_df.loc[mask, "feature_key"].astype(str).nunique())
            if not mask.empty and "feature_key" in candidate_df.columns
            else 0
        )
    return row


def make_rank_maps(
    candidate_df: pd.DataFrame,
    denovo_path: Optional[Path],
    mass_tolerance_da: float,
) -> Tuple[Dict[Tuple[str, str], int], Dict[Tuple[str, str], int], Dict[Tuple[str, str], int]]:
    report_rank_map: Dict[Tuple[str, str], int] = {}
    if denovo_path is not None and denovo_path.exists():
        report_df = pd.read_csv(denovo_path)
        if {"feature_key", "predicted_smiles", "rank"}.issubset(report_df.columns):
            for row in report_df[["feature_key", "predicted_smiles", "rank"]].itertuples(index=False):
                key = (str(row.feature_key), str(row.predicted_smiles))
                rank = int(pd.to_numeric(pd.Series([row.rank]), errors="coerce").fillna(999999).iloc[0])
                report_rank_map[key] = min(rank, report_rank_map.get(key, rank))

    window_rank_map: Dict[Tuple[str, str], int] = {}
    all_score_rank_map: Dict[Tuple[str, str], int] = {}
    required = {"feature_key", "predicted_smiles", "abs_mass_error_da", "mass_match_score", "generation_support_fraction", "mean_log_prob"}
    if required.issubset(candidate_df.columns):
        ranked_window = rank_selfies_pool(candidate_df, tolerance_da=float(mass_tolerance_da))
        if not ranked_window.empty and "selfies_rank" in ranked_window.columns:
            for row in ranked_window[["feature_key", "predicted_smiles", "selfies_rank"]].itertuples(index=False):
                key = (str(row.feature_key), str(row.predicted_smiles))
                rank = int(row.selfies_rank)
                window_rank_map[key] = min(rank, window_rank_map.get(key, rank))

        ranked_all = rank_selfies_pool(candidate_df, tolerance_da=float("inf"))
        if not ranked_all.empty and "selfies_rank" in ranked_all.columns:
            for row in ranked_all[["feature_key", "predicted_smiles", "selfies_rank"]].itertuples(index=False):
                key = (str(row.feature_key), str(row.predicted_smiles))
                rank = int(row.selfies_rank)
                all_score_rank_map[key] = min(rank, all_score_rank_map.get(key, rank))
    return report_rank_map, window_rank_map, all_score_rank_map


def best_control_match(pred_can: Optional[str], control_smiles: Iterable[str]) -> Tuple[bool, float]:
    if pred_can is None:
        return False, 0.0
    exact = False
    best = 0.0
    for ctrl_smi in control_smiles:
        if pred_can == ctrl_smi:
            exact = True
        best = max(best, tanimoto(pred_can, ctrl_smi))
    return exact, best


def evaluate_full_pool_oracle(
    *,
    ctx: LibContext,
    controls: Dict[str, List[str]],
    profile: str,
    profile_dir: Path,
    candidate_df: pd.DataFrame,
    denovo_path: Optional[Path],
    mass_tolerance_da: float,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    if candidate_df.empty or "feature_key" not in candidate_df.columns or "predicted_smiles" not in candidate_df.columns:
        rows = [
            {
                "lib_id": int(ctx.lib_id),
                "ion_mode": ctx.ion_mode,
                "profile": profile,
                "x_compound_id": xid,
                "n_candidates": 0,
                "exact_match": False,
                "best_tanimoto": 0.0,
                "failure_bucket": "no_predictions",
            }
            for xid in sorted(controls)
        ]
        return oracle_summary(controls, rows), rows

    pool = candidate_df.copy()
    pool["feature_key"] = pool["feature_key"].astype(str)
    pool["predicted_smiles"] = pool["predicted_smiles"].astype(str)
    pool["_xid"] = pool["feature_key"].apply(extract_xid)
    pool = pool[pool["_xid"].isin(set(controls))].copy()

    report_rank_map, window_rank_map, all_score_rank_map = make_rank_maps(pool, denovo_path, mass_tolerance_da)
    canonical_cache: Dict[str, Optional[str]] = {}

    rows: List[Dict[str, Any]] = []
    for xid, control_smiles in sorted(controls.items()):
        xid_pool = pool[pool["_xid"] == xid].copy()
        if xid_pool.empty:
            rows.append(
                {
                    "lib_id": int(ctx.lib_id),
                    "ion_mode": ctx.ion_mode,
                    "profile": profile,
                    "x_compound_id": xid,
                    "n_candidates": 0,
                    "exact_match": False,
                    "best_tanimoto": 0.0,
                    "failure_bucket": "no_predictions",
                    "oracle_report_rank": None,
                    "oracle_mass_window_rank": None,
                    "oracle_all_score_rank": None,
                    "oracle_in_profile_mass_window": False,
                }
            )
            continue

        candidate_rows: List[Dict[str, Any]] = []
        for source_index, row in xid_pool.reset_index(drop=True).iterrows():
            smi = str(row.get("predicted_smiles", "")).strip()
            if smi not in canonical_cache:
                canonical_cache[smi] = canonicalize(smi)
            exact, best_tan = best_control_match(canonical_cache[smi], control_smiles)
            feature_key = str(row.get("feature_key", ""))
            key = (feature_key, smi)
            abs_mass_error = pd.to_numeric(pd.Series([row.get("abs_mass_error_da")]), errors="coerce").iloc[0]
            candidate_rows.append(
                {
                    "feature_key": feature_key,
                    "predicted_smiles": smi,
                    "canonical_smiles": canonical_cache[smi],
                    "exact_match": bool(exact),
                    "best_tanimoto": float(best_tan),
                    "abs_mass_error_da": None if pd.isna(abs_mass_error) else float(abs_mass_error),
                    "predicted_formula": row.get("predicted_formula", row.get("candidate_formula", "")),
                    "mass_error_hypothesis": row.get("mass_error_hypothesis", ""),
                    "mass_hypothesis_adduct": row.get("mass_hypothesis_adduct", ""),
                    "mass_match_score": row.get("mass_match_score", ""),
                    "generation_support_fraction": row.get("generation_support_fraction", ""),
                    "mean_log_prob": row.get("mean_log_prob", ""),
                    "oracle_report_rank": report_rank_map.get(key),
                    "oracle_mass_window_rank": window_rank_map.get(key),
                    "oracle_all_score_rank": all_score_rank_map.get(key),
                    "source_order": int(source_index + 1),
                }
            )

        def sort_key(item: Dict[str, Any]) -> Tuple[Any, ...]:
            abs_mass = item["abs_mass_error_da"]
            report_rank = item["oracle_report_rank"]
            window_rank = item["oracle_mass_window_rank"]
            return (
                -float(item["best_tanimoto"]),
                not bool(item["exact_match"]),
                math.inf if abs_mass is None else float(abs_mass),
                math.inf if report_rank is None else int(report_rank),
                math.inf if window_rank is None else int(window_rank),
                int(item["source_order"]),
            )

        best = sorted(candidate_rows, key=sort_key)[0]
        best_tan = float(best["best_tanimoto"])
        if bool(best["exact_match"]):
            bucket = "exact_match"
        elif best_tan >= 0.7:
            bucket = "structural_hit"
        elif best_tan >= 0.3:
            bucket = "low_similarity"
        else:
            bucket = "structural_miss"

        rows.append(
            {
                "lib_id": int(ctx.lib_id),
                "ion_mode": ctx.ion_mode,
                "profile": profile,
                "x_compound_id": xid,
                "n_candidates": int(len(xid_pool)),
                "exact_match": bool(best["exact_match"]),
                "best_tanimoto": round(best_tan, 4),
                "failure_bucket": bucket,
                "oracle_feature_key": best["feature_key"],
                "oracle_smiles": best["predicted_smiles"],
                "oracle_canonical_smiles": best["canonical_smiles"],
                "oracle_predicted_formula": best["predicted_formula"],
                "oracle_abs_mass_error_da": best["abs_mass_error_da"],
                "oracle_mass_error_hypothesis": best["mass_error_hypothesis"],
                "oracle_mass_hypothesis_adduct": best["mass_hypothesis_adduct"],
                "oracle_mass_match_score": best["mass_match_score"],
                "oracle_generation_support_fraction": best["generation_support_fraction"],
                "oracle_mean_log_prob": best["mean_log_prob"],
                "oracle_report_rank": best["oracle_report_rank"],
                "oracle_mass_window_rank": best["oracle_mass_window_rank"],
                "oracle_all_score_rank": best["oracle_all_score_rank"],
                "oracle_in_profile_mass_window": bool(best["oracle_mass_window_rank"] is not None),
                "candidate_metrics_csv": str(profile_dir / "candidate_metrics_all.csv"),
                "denovo_csv": str(denovo_path) if denovo_path is not None else "",
            }
        )

    return oracle_summary(controls, rows), rows


def oracle_summary(controls: Dict[str, List[str]], rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    covered = [row for row in rows if int(row.get("n_candidates", 0) or 0) > 0]
    best_tans = [float(row.get("best_tanimoto", 0.0) or 0.0) for row in covered]
    exact_count = sum(1 for row in covered if bool(row.get("exact_match")))
    summary: Dict[str, Any] = {
        "n_controls": int(len(controls)),
        "n_controls_covered": int(len(covered)),
        "n_exact_match": int(exact_count),
        "exact_match_rate": round(exact_count / max(len(controls), 1), 4),
        "mean_best_tanimoto": round(float(np.mean(best_tans)) if best_tans else 0.0, 4),
        "median_best_tanimoto": round(float(np.median(best_tans)) if best_tans else 0.0, 4),
    }
    for threshold in THRESHOLDS:
        count = sum(1 for value in best_tans if value >= threshold)
        suffix = str(threshold)
        summary[f"n_oracle_tanimoto_ge_{suffix}"] = int(count)
        summary[f"oracle_tanimoto_ge_{suffix}_rate"] = round(count / max(len(controls), 1), 4)
    return summary


def profile_mass_tolerance(profile_dir: Path, fallback: float = 5.0) -> float:
    summary = read_json(profile_dir / "selfies_run_summary.json")
    value = summary.get("mass_tolerance_da", fallback)
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(fallback)


def collect_profile_metrics(
    *,
    ctx: LibContext,
    controls: Dict[str, List[str]],
    top_k: int,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    metric_rows: List[Dict[str, Any]] = []
    per_control_rows: List[Dict[str, Any]] = []
    oracle_rows: List[Dict[str, Any]] = []
    mass_rows: List[Dict[str, Any]] = []

    profiles_root = ctx.lib_dir / "selfies_profiles"
    profile_names = [name for name in PROFILE_ORDER if (profiles_root / name).exists()]
    profile_names.extend(sorted({p.name for p in profiles_root.iterdir() if p.is_dir()} - set(profile_names)) if profiles_root.exists() else [])

    for profile in profile_names:
        profile_dir = profiles_root / profile
        denovo_path = profile_dir / "denovo_report.csv"
        candidate_path = profile_dir / "candidate_metrics_all.csv"
        if denovo_path.exists():
            rows, per_rows = evaluate_prediction_path(
                ctx=ctx,
                controls=controls,
                pred_path=denovo_path,
                source_type="selfies_profile",
                source_name=profile,
                top_k=top_k,
            )
            metric_rows.extend(rows)
            per_control_rows.extend(per_rows)
        if candidate_path.exists():
            candidate_df = pd.read_csv(candidate_path)
            mass_rows.append(mass_coverage_row(ctx, profile, profile_dir, candidate_df))
            tolerance = profile_mass_tolerance(profile_dir, fallback=5.0)
            oracle_summary_payload, profile_oracle_rows = evaluate_full_pool_oracle(
                ctx=ctx,
                controls=controls,
                profile=profile,
                profile_dir=profile_dir,
                candidate_df=candidate_df,
                denovo_path=denovo_path if denovo_path.exists() else None,
                mass_tolerance_da=tolerance,
            )
            metric_rows.append(
                metric_row(
                    ctx=ctx,
                    source_type="selfies_profile",
                    source_name=profile,
                    metric_scope="full_pool_oracle",
                    summary=oracle_summary_payload,
                    path=candidate_path,
                )
            )
            oracle_rows.extend(profile_oracle_rows)
    return metric_rows, per_control_rows, oracle_rows, mass_rows


def collect_retrieval_metrics(
    *,
    ctx: LibContext,
    controls: Dict[str, List[str]],
    top_k: int,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    retrieval = ctx.manifest.get("retrieval", {}) if isinstance(ctx.manifest.get("retrieval"), dict) else {}
    retrieval_path = existing_path(retrieval.get("merged_predictions_csv"), ctx.lib_dir.parent)
    if retrieval_path is None:
        candidates = sorted((ctx.lib_dir / "retrieval").glob("*merged_predictions.csv"))
        retrieval_path = candidates[0] if candidates else ctx.lib_dir / "retrieval" / "merged_predictions.csv"
    if not retrieval_path.exists():
        return [], []
    return evaluate_prediction_path(
        ctx=ctx,
        controls=controls,
        pred_path=retrieval_path,
        source_type="retrieval",
        source_name="retrieval",
        top_k=top_k,
    )


def policy_metric_lookup(comparison_path: Path) -> Dict[str, Dict[str, Any]]:
    payload = read_json(comparison_path)
    out: Dict[str, Dict[str, Any]] = {}
    for item in payload.get("policy_metrics", []) or []:
        policy = str(item.get("policy", ""))
        if policy:
            out[policy] = item
    return out


def collect_policy_metrics(
    *,
    ctx: LibContext,
    controls: Dict[str, List[str]],
    top_k: int,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    metric_rows: List[Dict[str, Any]] = []
    per_control_rows: List[Dict[str, Any]] = []
    lane_decision_rows: List[Dict[str, Any]] = []

    for lane_dir in sorted(ctx.lib_dir.glob("*_lane")):
        lane = lane_dir.name.replace("_lane", "")
        comparison_path = lane_dir / "comparison_summary.json"
        lookup = policy_metric_lookup(comparison_path)
        comparison = read_json(comparison_path)
        decision = comparison.get("comparison_decision", {}) or comparison.get("decision", {}) or {}
        winner = decision.get("winning_policy_overall", {}) or {}
        blend = decision.get("winning_blend_policy", {}) or {}
        lane_decision_rows.append(
            {
                "lib_id": int(ctx.lib_id),
                "ion_mode": ctx.ion_mode,
                "lane": lane,
                "keep_selfies_track": bool(decision.get("keep_selfies_track", False)),
                "strong_success": bool(decision.get("strong_success", False)),
                "winning_policy_overall": winner.get("policy", ""),
                "winning_policy_mean_best_tanimoto": winner.get("mean_best_tanimoto", ""),
                "winning_policy_n_tanimoto_ge_0p7": winner.get("n_topk_tanimoto_ge_0.7", ""),
                "winning_blend_policy": blend.get("policy", ""),
                "winning_blend_mean_best_tanimoto": blend.get("mean_best_tanimoto", ""),
                "winning_blend_n_tanimoto_ge_0p7": blend.get("n_topk_tanimoto_ge_0.7", ""),
                "close_reason": decision.get("close_reason", ""),
                "comparison_summary_json": str(comparison_path) if comparison_path.exists() else "",
            }
        )
        policies_root = lane_dir / "policies"
        for policy_dir in sorted(p for p in policies_root.glob("*") if p.is_dir()) if policies_root.exists() else []:
            pred_path = policy_dir / "merged_predictions.csv"
            if not pred_path.exists():
                candidates = sorted(policy_dir.glob("*merged_predictions.csv"))
                pred_path = candidates[0] if candidates else pred_path
            if not pred_path.exists():
                continue
            rows, per_rows = evaluate_prediction_path(
                ctx=ctx,
                controls=controls,
                pred_path=pred_path,
                source_type="merge_policy",
                source_name=policy_dir.name,
                top_k=top_k,
                lane=lane,
                policy=policy_dir.name,
            )
            for row in rows:
                prior = lookup.get(policy_dir.name, {})
                row["diagnostic"] = bool(prior.get("diagnostic", False))
                row["is_blend_policy"] = bool(prior.get("is_blend_policy", False))
            metric_rows.extend(rows)
            per_control_rows.extend(per_rows)
    return metric_rows, per_control_rows, lane_decision_rows


def selected_vs_augmented_row(ctx: LibContext, run_root: Path) -> Dict[str, Any]:
    selected = ctx.manifest.get("selected_delivery", {}) if isinstance(ctx.manifest.get("selected_delivery"), dict) else {}
    augmented = ctx.manifest.get("denovo_augmented_delivery", {}) if isinstance(ctx.manifest.get("denovo_augmented_delivery"), dict) else {}
    selected_path = existing_path(selected.get("merged_predictions_csv"), run_root)
    augmented_path = existing_path(augmented.get("csv"), run_root)

    row: Dict[str, Any] = {
        "lib_id": int(ctx.lib_id),
        "ion_mode": ctx.ion_mode,
        "selected_lane": selected.get("lane", ""),
        "selected_policy": selected.get("policy", ""),
        "selected_reason": selected.get("reason", ""),
        "selected_path": str(selected_path) if selected_path is not None else "",
        "augmented_path": str(augmented_path) if augmented_path is not None else "",
        "augmented_manifest_rows": int(augmented.get("rows", 0) or 0),
        "augmented_manifest_features_with_predictions": int(augmented.get("features_with_predictions", 0) or 0),
        "selected_is_retrieval_only": str(selected.get("lane", "")) == "retrieval_only",
    }
    if selected_path is not None and selected_path.exists():
        selected_df = pd.read_csv(selected_path)
        row["selected_rows"] = int(selected_df.shape[0])
        row["selected_features"] = int(selected_df["feature_key"].astype(str).nunique()) if "feature_key" in selected_df.columns else 0
        row["selected_track_counts"] = json.dumps(selected_df.get("track", pd.Series(dtype=str)).astype(str).value_counts().to_dict(), sort_keys=True)
    else:
        row["selected_rows"] = 0
        row["selected_features"] = 0
        row["selected_track_counts"] = "{}"
    if augmented_path is not None and augmented_path.exists():
        aug_df = pd.read_csv(augmented_path)
        row["augmented_rows"] = int(aug_df.shape[0])
        row["augmented_features"] = int(aug_df["feature_key"].astype(str).nunique()) if "feature_key" in aug_df.columns else 0
        row["augmented_track_counts"] = json.dumps(aug_df.get("track", pd.Series(dtype=str)).astype(str).value_counts().to_dict(), sort_keys=True)
        row["augmented_profile_counts"] = json.dumps(aug_df.get("source_profile", pd.Series(dtype=str)).astype(str).value_counts().to_dict(), sort_keys=True)
    else:
        row["augmented_rows"] = 0
        row["augmented_features"] = 0
        row["augmented_track_counts"] = "{}"
        row["augmented_profile_counts"] = "{}"
    return row


def feature_scope_row(ctx: LibContext) -> Dict[str, Any]:
    scope = ctx.manifest.get("feature_scope", {}) if isinstance(ctx.manifest.get("feature_scope"), dict) else {}
    expected = scope.get("expected_features")
    output = scope.get("output_features")
    return {
        "lib_id": int(ctx.lib_id),
        "ion_mode": ctx.ion_mode,
        "expected_count_enforced": bool(scope.get("expected_count_enforced", False)),
        "expected_features": expected,
        "output_features": output,
        "input_features": scope.get("input_features", ""),
        "scoped_features_before_limit": scope.get("scoped_features_before_limit", ""),
        "feature_scope_ok": bool(expected == output) if expected is not None and output is not None else False,
        "requested_lib_id": scope.get("requested_lib_id", ""),
        "requested_ion_mode": scope.get("requested_ion_mode", ""),
        "scoped_feature_table": scope.get("scoped_feature_table", ""),
    }


def latent_alignment_row(ctx: LibContext) -> Dict[str, Any]:
    selfies = ctx.manifest.get("selfies", {}) if isinstance(ctx.manifest.get("selfies"), dict) else {}
    audit = selfies.get("latent_alignment_audit", {}) if isinstance(selfies.get("latent_alignment_audit"), dict) else {}
    return {
        "lib_id": int(ctx.lib_id),
        "ion_mode": ctx.ion_mode,
        "status": audit.get("status", ""),
        "feature_count": audit.get("feature_count", ""),
        "matched_retrieval_neighbor_count": audit.get("matched_retrieval_neighbor_count", ""),
        "feature_norm_mean": audit.get("feature_norm_mean", ""),
        "rank1_embedding_norm_mean": audit.get("rank1_embedding_norm_mean", ""),
        "rank1_feature_cosine_mean": audit.get("rank1_feature_cosine_mean", ""),
        "calibrated_transform_required_for_delivery": audit.get("calibrated_transform_required_for_delivery", ""),
        "calibrated_transform_reason": audit.get("calibrated_transform_reason", ""),
        "json": audit.get("json", ""),
        "csv": audit.get("csv", ""),
    }


def global_augmented_audit(run_root: Path) -> Dict[str, Any]:
    path = run_root / "metabolon_denovo_augmented_predictions.csv"
    if not path.exists():
        return {"path": str(path), "exists": False}
    df = pd.read_csv(path)
    payload: Dict[str, Any] = {
        "path": str(path),
        "exists": True,
        "rows": int(df.shape[0]),
        "columns": list(df.columns),
    }
    for col in ["lib_id", "ion_mode", "track", "source_profile", "selection_reason", "top_candidate_track"]:
        if col in df.columns:
            payload[f"{col}_counts"] = {str(k): int(v) for k, v in df[col].astype(str).value_counts().to_dict().items()}
    return payload


def choose_next_branch(metrics_df: pd.DataFrame) -> Tuple[str, str, Dict[str, Any]]:
    profile_df = metrics_df[metrics_df["source_type"].eq("selfies_profile")].copy()
    if profile_df.empty:
        return "insufficient_data", "No SELFIES profile metrics were available.", {}

    pivot = profile_df.pivot_table(
        index=["lib_id", "ion_mode", "source_name"],
        columns="metric_scope",
        values=["mean_best_tanimoto", "n_tanimoto_ge_0p7", "n_controls_covered"],
        aggfunc="max",
    )
    improvements: List[Dict[str, Any]] = []
    weak_oracles: List[Dict[str, Any]] = []
    for index, row in pivot.iterrows():
        lib_id, ion_mode, profile = index
        oracle_mean = float(row.get(("mean_best_tanimoto", "full_pool_oracle"), 0.0) or 0.0)
        top10_mean = float(row.get(("mean_best_tanimoto", "top10"), 0.0) or 0.0)
        oracle_hits = int(row.get(("n_tanimoto_ge_0p7", "full_pool_oracle"), 0) or 0)
        top10_hits = int(row.get(("n_tanimoto_ge_0p7", "top10"), 0) or 0)
        item = {
            "lib_id": int(lib_id),
            "ion_mode": ion_mode,
            "profile": profile,
            "oracle_mean_best_tanimoto": round(oracle_mean, 4),
            "top10_mean_best_tanimoto": round(top10_mean, 4),
            "oracle_n_tanimoto_ge_0p7": oracle_hits,
            "top10_n_tanimoto_ge_0p7": top10_hits,
            "oracle_minus_top10_mean": round(oracle_mean - top10_mean, 4),
            "oracle_minus_top10_hits": int(oracle_hits - top10_hits),
        }
        if oracle_hits > top10_hits or (oracle_mean - top10_mean) >= 0.05:
            improvements.append(item)
        else:
            weak_oracles.append(item)

    details = {
        "profiles_with_oracle_headroom": improvements,
        "profiles_without_oracle_headroom": weak_oracles,
    }
    if improvements:
        reason = (
            "Full-pool oracle has structural headroom in at least one profile, "
            "so ranking/selection should be improved before launching expanded generation."
        )
        return "reranking_first", reason, details
    reason = (
        "Full-pool oracle does not materially improve over current top-10 metrics, "
        "so generation constraints should be improved before expanded runs."
    )
    return "generation_constraints_first", reason, details


def write_markdown_report(
    *,
    output_path: Path,
    run_root: Path,
    metrics_df: pd.DataFrame,
    selected_df: pd.DataFrame,
    lane_df: pd.DataFrame,
    feature_scope_df: pd.DataFrame,
    latent_df: pd.DataFrame,
    branch: str,
    branch_reason: str,
    has_rdkit: bool,
) -> None:
    lines: List[str] = []
    lines.append("# Post-Run Decision Package")
    lines.append("")
    lines.append(f"- Run root: `{run_root}`")
    lines.append(f"- RDKit available: `{has_rdkit}`")
    lines.append(f"- Recommendation: `{branch}`")
    lines.append(f"- Rationale: {branch_reason}")
    lines.append("- Expanded run decision: `do_not_launch_expanded_yet`")
    lines.append("")

    lines.append("## Selected Delivery")
    if selected_df.empty:
        lines.append("No selected delivery rows found.")
    else:
        cols = ["lib_id", "ion_mode", "selected_lane", "selected_policy", "selected_reason", "selected_rows", "augmented_rows"]
        lines.append(dataframe_to_markdown(selected_df[cols]))
    lines.append("")

    lines.append("## Retrieval vs SELFIES Quality")
    if metrics_df.empty:
        lines.append("No metric rows found.")
    else:
        display = metrics_df[
            metrics_df["source_type"].isin(["retrieval", "selfies_profile"])
            & metrics_df["metric_scope"].isin(["rank1", "top10", "full_pool_oracle"])
        ].copy()
        display = display[
            [
                "lib_id",
                "ion_mode",
                "source_type",
                "source_name",
                "metric_scope",
                "n_controls_covered",
                "n_exact_match",
                "n_tanimoto_ge_0p7",
                "mean_best_tanimoto",
            ]
        ].sort_values(["lib_id", "source_type", "source_name", "metric_scope"])
        lines.append(dataframe_to_markdown(display))
    lines.append("")

    lines.append("## Lane Decisions")
    if lane_df.empty:
        lines.append("No lane decision rows found.")
    else:
        cols = [
            "lib_id",
            "ion_mode",
            "lane",
            "keep_selfies_track",
            "strong_success",
            "winning_policy_overall",
            "winning_policy_mean_best_tanimoto",
            "winning_policy_n_tanimoto_ge_0p7",
            "close_reason",
        ]
        lines.append(dataframe_to_markdown(lane_df[cols]))
    lines.append("")

    lines.append("## Scope And Latent Checks")
    if not feature_scope_df.empty:
        lines.append(dataframe_to_markdown(feature_scope_df[["lib_id", "ion_mode", "expected_features", "output_features", "feature_scope_ok"]]))
    if not latent_df.empty:
        lines.append("")
        lines.append(dataframe_to_markdown(latent_df[["lib_id", "ion_mode", "status", "matched_retrieval_neighbor_count", "rank1_feature_cosine_mean", "calibrated_transform_required_for_delivery"]]))
    lines.append("")

    lines.append("## Next Work")
    if branch == "reranking_first":
        lines.extend(
            [
                "1. Implement formula/adduct/subformula-aware reranking on frozen candidates.",
                "2. Add retrieval-neighbor anchoring and spectrum-consistency features to the ranker.",
                "3. Re-run this package and require rank-1/top-10 gains without retrieval regression before promotion.",
            ]
        )
    elif branch == "generation_constraints_first":
        lines.extend(
            [
                "1. Improve generation constraints before expanding: formula/adduct-constrained sampling and retrieval-anchored generation.",
                "2. Keep retrieval as selected delivery while using SELFIES only as augmented exploratory candidates.",
                "3. Re-run oracle analysis after constrained generation to check whether the pool quality improves.",
            ]
        )
    else:
        lines.append("1. Fill missing artifacts and rerun this package.")
    lines.append("")
    output_path.write_text("\n".join(lines) + "\n")


def build_decision_package(
    *,
    run_root: Path,
    known_x: Path,
    known_x_all: Path,
    output_dir: Path,
    top_k: int,
) -> Dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    libs = discover_libs(run_root)
    if not libs:
        raise RuntimeError(f"No lib manifests found under {run_root}")

    all_metric_rows: List[Dict[str, Any]] = []
    all_per_control_rows: List[Dict[str, Any]] = []
    all_oracle_rows: List[Dict[str, Any]] = []
    all_mass_rows: List[Dict[str, Any]] = []
    all_policy_rows: List[Dict[str, Any]] = []
    all_policy_per_control_rows: List[Dict[str, Any]] = []
    lane_rows: List[Dict[str, Any]] = []
    selected_rows: List[Dict[str, Any]] = []
    feature_scope_rows: List[Dict[str, Any]] = []
    latent_rows: List[Dict[str, Any]] = []

    for ctx in libs:
        controls, control_scope = load_controls(
            known_x,
            known_x_all,
            control_lib_id=int(ctx.lib_id),
            fail_on_mixed_lib_id=False,
        )
        retrieval_rows, retrieval_per_control = collect_retrieval_metrics(ctx=ctx, controls=controls, top_k=top_k)
        profile_rows, profile_per_control, oracle_rows, mass_rows = collect_profile_metrics(
            ctx=ctx,
            controls=controls,
            top_k=top_k,
        )
        policy_rows, policy_per_control, lib_lane_rows = collect_policy_metrics(ctx=ctx, controls=controls, top_k=top_k)

        for row in retrieval_rows + profile_rows + policy_rows:
            row["control_scope_xids"] = control_scope.get("control_scope_xids", "")
            row["controls_loaded"] = control_scope.get("controls_loaded", "")

        all_metric_rows.extend(retrieval_rows)
        all_metric_rows.extend(profile_rows)
        all_policy_rows.extend(policy_rows)
        all_per_control_rows.extend(retrieval_per_control)
        all_per_control_rows.extend(profile_per_control)
        all_policy_per_control_rows.extend(policy_per_control)
        all_oracle_rows.extend(oracle_rows)
        all_mass_rows.extend(mass_rows)
        lane_rows.extend(lib_lane_rows)
        selected_rows.append(selected_vs_augmented_row(ctx, run_root))
        feature_scope_rows.append(feature_scope_row(ctx))
        latent_rows.append(latent_alignment_row(ctx))

    metrics_df = pd.DataFrame(all_metric_rows)
    policy_df = pd.DataFrame(all_policy_rows)
    per_control_df = pd.DataFrame(all_per_control_rows)
    policy_per_control_df = pd.DataFrame(all_policy_per_control_rows)
    oracle_df = pd.DataFrame(all_oracle_rows)
    mass_df = pd.DataFrame(all_mass_rows)
    lane_df = pd.DataFrame(lane_rows)
    selected_df = pd.DataFrame(selected_rows)
    feature_scope_df = pd.DataFrame(feature_scope_rows)
    latent_df = pd.DataFrame(latent_rows)

    if not metrics_df.empty:
        metrics_df = metrics_df.reindex(columns=[*METRIC_COLUMNS, *[c for c in metrics_df.columns if c not in METRIC_COLUMNS]])
    metrics_df.to_csv(output_dir / "quality_metrics_by_scope.csv", index=False)
    policy_df.to_csv(output_dir / "policy_metrics_by_scope.csv", index=False)
    per_control_df.to_csv(output_dir / "selected_output_per_control_metrics.csv", index=False)
    policy_per_control_df.to_csv(output_dir / "policy_per_control_metrics.csv", index=False)
    oracle_df.to_csv(output_dir / "candidate_pool_oracle_per_control.csv", index=False)
    mass_df.to_csv(output_dir / "profile_mass_coverage.csv", index=False)
    lane_df.to_csv(output_dir / "lane_decisions.csv", index=False)
    selected_df.to_csv(output_dir / "selected_vs_augmented_audit.csv", index=False)
    feature_scope_df.to_csv(output_dir / "feature_scope_audit.csv", index=False)
    latent_df.to_csv(output_dir / "latent_alignment_audit.csv", index=False)

    branch, branch_reason, branch_details = choose_next_branch(metrics_df)
    global_aug = global_augmented_audit(run_root)
    summary = {
        "run_root": str(run_root),
        "output_dir": str(output_dir),
        "top_k": int(top_k),
        "rdkit_available": bool(HAS_RDKIT),
        "libs": [{"lib_id": int(ctx.lib_id), "ion_mode": ctx.ion_mode, "lib_dir": str(ctx.lib_dir)} for ctx in libs],
        "recommendation": {
            "next_branch": branch,
            "reason": branch_reason,
            "expanded_run_decision": "do_not_launch_expanded_yet",
            "details": branch_details,
        },
        "global_augmented_audit": global_aug,
        "outputs": {
            "quality_metrics_by_scope_csv": str(output_dir / "quality_metrics_by_scope.csv"),
            "candidate_pool_oracle_per_control_csv": str(output_dir / "candidate_pool_oracle_per_control.csv"),
            "profile_mass_coverage_csv": str(output_dir / "profile_mass_coverage.csv"),
            "policy_metrics_by_scope_csv": str(output_dir / "policy_metrics_by_scope.csv"),
            "selected_vs_augmented_audit_csv": str(output_dir / "selected_vs_augmented_audit.csv"),
            "lane_decisions_csv": str(output_dir / "lane_decisions.csv"),
            "feature_scope_audit_csv": str(output_dir / "feature_scope_audit.csv"),
            "latent_alignment_audit_csv": str(output_dir / "latent_alignment_audit.csv"),
            "markdown_report": str(output_dir / "post_run_decision_report.md"),
        },
    }
    write_json(output_dir / "post_run_decision_summary.json", summary)
    write_markdown_report(
        output_path=output_dir / "post_run_decision_report.md",
        run_root=run_root,
        metrics_df=metrics_df,
        selected_df=selected_df,
        lane_df=lane_df,
        feature_scope_df=feature_scope_df,
        latent_df=latent_df,
        branch=branch,
        branch_reason=branch_reason,
        has_rdkit=bool(HAS_RDKIT),
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", required=True, type=Path, help="Completed merged run artifact root.")
    parser.add_argument("--known-x", required=True, type=Path, help="known_X.csv path.")
    parser.add_argument("--known-x-all", required=True, type=Path, help="known_X_all_available_structures_simple.csv path.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to <run-root>/post_run_decision_package.",
    )
    parser.add_argument("--top-k", type=int, default=10, help="Selected-output top-k scope to compare with rank-1 and oracle.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_root = args.run_root.resolve()
    output_dir = args.output_dir.resolve() if args.output_dir is not None else run_root / "post_run_decision_package"
    summary = build_decision_package(
        run_root=run_root,
        known_x=args.known_x.resolve(),
        known_x_all=args.known_x_all.resolve(),
        output_dir=output_dir,
        top_k=int(args.top_k),
    )
    print(json.dumps(summary["recommendation"], indent=2, sort_keys=True))
    print(f"Wrote decision package: {output_dir}")


if __name__ == "__main__":
    main()
