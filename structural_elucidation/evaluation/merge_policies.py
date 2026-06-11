#!/usr/bin/env python3
"""Benchmark fixed retrieval/SELFIES merge policies on frozen X209 artifacts."""

from __future__ import annotations

from structural_elucidation.paths import REPO_ROOT

import argparse
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List

import numpy as np
import pandas as pd
import structural_elucidation.utils.merge as merge_policy


CORE_POLICY_ORDER = [
    "retrieval_only",
    "selfies_only_5da",
    "pin_high_rank1_joint_rest",
    "top5_selfies_then_joint_rest",
    "high_protected_medium_low_joint",
    "retrieval_rank1_plus_best_selfies",
]
DIAGNOSTIC_POLICY_ORDER = [
    "selfies_only_1da_diag",
    "selfies_only_0p5da_diag",
]
BLEND_POLICY_NAMES = {
    "pin_high_rank1_joint_rest",
    "top5_selfies_then_joint_rest",
    "high_protected_medium_low_joint",
    "retrieval_rank1_plus_best_selfies",
}
HITMAX_REFERENCE_JSON = Path("models/x209_known_x_hitmax_reference.json")
HITMAX_REFERENCE_CEILING_POLICY = "selfies_only_5da"
HITMAX_REFERENCE_COVERAGE_POLICY = "pin_high_rank1_joint_rest"
HITMAX_PROMOTION_GATES = {
    "min_n_controls_covered": 48,
    "ceiling_policy": HITMAX_REFERENCE_CEILING_POLICY,
    "coverage_preserving_policy": HITMAX_REFERENCE_COVERAGE_POLICY,
}
PROMOTION_GATES = {
    "min_mean_best_tanimoto": 0.3080,
    "min_n_topk_tanimoto_ge_0p7": 3,
    "min_n_controls_covered": 50,
}
STRONG_SUCCESS_GATES = {
    "selfies_mean_best_tanimoto_baseline": 0.3568,
    "selfies_n_topk_tanimoto_ge_0p7_baseline": 4,
}
PLAN_HARD_WINDOW_REFERENCE = {
    "selfies_only_1da_diag": {
        "n_controls_covered": 32,
        "mean_best_tanimoto": 0.3376,
    },
    "selfies_only_0p5da_diag": {
        "n_controls_covered": 21,
        "mean_best_tanimoto": 0.2809,
    },
}


@dataclass
class PolicyResult:
    name: str
    diagnostic: bool
    ranked_csv: Path
    merged_csv: Path
    control_eval_summary_json: Path
    control_eval_per_xid_csv: Path
    metric: dict


def artifact_filename(prefix: str, stem: str, suffix: str) -> str:
    normalized = str(prefix).strip()
    if not normalized:
        return f"{stem}{suffix}"
    return f"{normalized}_{stem}{suffix}"


def resolve_existing_artifact(base_dir: Path, candidates: Iterable[str]) -> Path:
    for candidate in candidates:
        path = base_dir / candidate
        if path.exists():
            return path
    return base_dir / next(iter(candidates))


def sigmoid(x: float | pd.Series) -> float | pd.Series:
    return 1.0 / (1.0 + np.exp(-x))


def compute_selfies_score(
    mass_match_score: float | pd.Series,
    generation_support_fraction: float | pd.Series,
    mean_log_prob: float | pd.Series,
) -> float | pd.Series:
    """Composite score requested in PLAN.md."""
    is_series = isinstance(mass_match_score, pd.Series)
    index = mass_match_score.index if is_series else None
    mm = pd.to_numeric(mass_match_score, errors="coerce")
    support = pd.to_numeric(generation_support_fraction, errors="coerce")
    log_prob = pd.to_numeric(mean_log_prob, errors="coerce")
    lp_term = sigmoid(np.nan_to_num(log_prob, nan=-99.0) + 1.0)
    out = 0.55 * np.nan_to_num(mm, nan=0.0) + 0.25 * np.nan_to_num(support, nan=0.0) + 0.20 * lp_term
    if is_series:
        return pd.Series(np.asarray(out, dtype=float), index=index)
    return float(np.asarray(out).item())


def normalize_score(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(0.0).clip(lower=0.0, upper=1.0)


def rank_selfies_pool(candidate_metrics_df: pd.DataFrame, tolerance_da: float) -> pd.DataFrame:
    """Rank SELFIES pool for a given mass window and enforce per-feature SMILES uniqueness."""
    pool = candidate_metrics_df.copy()
    pool = pool[pd.to_numeric(pool["abs_mass_error_da"], errors="coerce") <= float(tolerance_da)].copy()
    if pool.empty:
        return pool

    pool["selfies_score"] = compute_selfies_score(
        pool["mass_match_score"],
        pool["generation_support_fraction"],
        pool["mean_log_prob"],
    )
    pool["normalized_merge_score"] = normalize_score(pool["selfies_score"])
    pool["abs_mass_error_da"] = pd.to_numeric(pool["abs_mass_error_da"], errors="coerce").fillna(np.inf)
    pool["generation_support_fraction"] = (
        pd.to_numeric(pool["generation_support_fraction"], errors="coerce").fillna(0.0)
    )
    pool["predicted_smiles"] = pool["predicted_smiles"].astype(str)

    pool = (
        pool.sort_values(
            [
                "feature_key",
                "selfies_score",
                "abs_mass_error_da",
                "generation_support_fraction",
                "predicted_smiles",
            ],
            ascending=[True, False, True, False, True],
            kind="mergesort",
        )
        .drop_duplicates(subset=["feature_key", "predicted_smiles"], keep="first")
        .copy()
    )
    pool["selfies_rank"] = pool.groupby("feature_key").cumcount() + 1
    return pool.reset_index(drop=True)


def compare_selfies_pool_to_report(selfies_pool_df: pd.DataFrame, selfies_report_df: pd.DataFrame) -> dict:
    report_pairs = {
        (str(row.feature_key), str(row.predicted_smiles))
        for row in selfies_report_df[["feature_key", "predicted_smiles"]].itertuples(index=False)
    }
    pool_pairs = {
        (str(row.feature_key), str(row.predicted_smiles))
        for row in selfies_pool_df[["feature_key", "predicted_smiles"]].itertuples(index=False)
    }
    return {
        "pool_rows": int(len(pool_pairs)),
        "report_rows": int(len(report_pairs)),
        "pool_minus_report": int(len(pool_pairs - report_pairs)),
        "report_minus_pool": int(len(report_pairs - pool_pairs)),
        "exact_match": bool(pool_pairs == report_pairs),
        "selfies_only_duplicates_removed": int(
            selfies_pool_df.duplicated(subset=["feature_key", "predicted_smiles"]).sum()
        ),
    }


def load_track_metric(summary_json_path: Path, track_name: str) -> dict:
    payload = json.loads(summary_json_path.read_text())
    for track in payload.get("tracks", []):
        if str(track.get("track")) == str(track_name):
            return track
    raise KeyError(f"Track '{track_name}' not found in {summary_json_path}")


def score_tuple(metric: dict) -> tuple[float, int]:
    return (
        float(metric.get("mean_best_tanimoto", 0.0)),
        int(metric.get("n_topk_tanimoto_ge_0.7", 0)),
    )


def metric_from_policy_eval(summary_json_path: Path) -> dict:
    payload = json.loads(summary_json_path.read_text())
    merged_tracks = [t for t in payload.get("tracks", []) if str(t.get("track")) == "merged"]
    if merged_tracks:
        return dict(merged_tracks[0])
    if payload.get("tracks"):
        return dict(payload["tracks"][0])
    return {}


def map_retrieval_rows_for_merge(retrieval_df: pd.DataFrame) -> pd.DataFrame:
    out = retrieval_df.copy()
    out["track"] = "retrieval"
    out["smiles"] = out["smiles"].astype(str)
    out["normalized_merge_score"] = normalize_score(out["retrieval_score"])
    out["score_for_sort"] = out["normalized_merge_score"]
    out["abs_error_for_sort"] = pd.to_numeric(out["abs_mass_diff"], errors="coerce").fillna(np.inf)
    out["support_for_sort"] = 0.0
    out["source_rank"] = pd.to_numeric(out["rank"], errors="coerce").fillna(9999).astype(int)
    out["selfies_score"] = np.nan
    out["mass_match_score"] = np.nan
    out["mean_log_prob"] = np.nan
    out["generation_support_fraction"] = np.nan
    out["predicted_exact_mass"] = np.nan
    out["candidate_formula"] = out.get("candidate_formula", np.nan)
    return out


def map_selfies_rows_for_merge(selfies_pool_df: pd.DataFrame, feature_meta_df: pd.DataFrame) -> pd.DataFrame:
    out = selfies_pool_df.copy()
    out["track"] = "selfies"
    out["smiles"] = out["predicted_smiles"].astype(str)
    out["source_rank"] = pd.to_numeric(out["selfies_rank"], errors="coerce").fillna(9999).astype(int)
    out["score_for_sort"] = out["normalized_merge_score"]
    out["abs_error_for_sort"] = pd.to_numeric(out["abs_mass_error_da"], errors="coerce").fillna(np.inf)
    out["support_for_sort"] = pd.to_numeric(out["generation_support_fraction"], errors="coerce").fillna(0.0)
    out["retrieval_score"] = np.nan
    out["confidence_tier"] = None
    out["recommended"] = None
    out["latent_cosine"] = np.nan
    out["gate_used"] = np.nan
    out["abs_mass_diff"] = np.nan
    out["fragment_match_count"] = np.nan
    out["fragment_matched_peak_fraction"] = np.nan
    out["fragment_explained_intensity_fraction"] = np.nan
    out["candidate_latent_cosine"] = np.nan
    out["name"] = None
    out["status"] = "ok"
    out["candidate_formula"] = out.get("predicted_formula", np.nan)
    if not feature_meta_df.empty and "feature_key" in feature_meta_df.columns:
        merge_cols = ["feature_key"]
        for col in ["comp_id"]:
            if col in feature_meta_df.columns:
                merge_cols.append(col)
        out = out.merge(feature_meta_df[merge_cols].drop_duplicates("feature_key"), on="feature_key", how="left")
    if "comp_id" not in out.columns:
        out["comp_id"] = np.nan
    return out


def joint_sort_rows(pool_df: pd.DataFrame) -> pd.DataFrame:
    if pool_df.empty:
        return pool_df
    sort_df = pool_df.copy()
    sort_df["track_priority"] = sort_df["track"].map({"retrieval": 0, "selfies": 1}).fillna(9).astype(int)
    sort_df = sort_df.sort_values(
        [
            "score_for_sort",
            "abs_error_for_sort",
            "support_for_sort",
            "track_priority",
            "source_rank",
            "smiles",
        ],
        ascending=[False, True, False, True, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    sort_df = sort_df.drop(columns=["track_priority"], errors="ignore")
    return sort_df


def finalize_policy_rows(ranked_rows_df: pd.DataFrame, topk: int = 10) -> pd.DataFrame:
    if ranked_rows_df.empty:
        return ranked_rows_df
    out = ranked_rows_df.copy()
    out["policy_rank"] = out.groupby("feature_key").cumcount() + 1
    out["include_topk"] = out["policy_rank"] <= int(topk)
    return out


def to_merged_predictions_csv_schema(policy_rows_df: pd.DataFrame, topk: int = 10) -> pd.DataFrame:
    if policy_rows_df.empty:
        return pd.DataFrame(
            columns=[
                "feature_key",
                "track",
                "rank",
                "authoritative_track_default",
                "comp_id",
                "smiles",
                "name",
                "latent_cosine",
                "gate_used",
                "abs_mass_diff",
                "confidence_tier",
                "recommended",
                "candidate_formula",
                "retrieval_score",
                "fragment_match_count",
                "fragment_matched_peak_fraction",
                "fragment_explained_intensity_fraction",
                "predicted_exact_mass",
                "mass_match_score",
                "mean_log_prob",
                "generation_support_fraction",
                "candidate_latent_cosine",
                "merged_rank",
                "is_authoritative_track",
                "normalized_merge_score",
                "selfies_score",
                "policy_rank",
            ]
        )

    top_df = policy_rows_df[policy_rows_df["policy_rank"] <= int(topk)].copy()
    top_df["merged_rank"] = top_df["policy_rank"].astype(int)
    top_df["rank"] = top_df["source_rank"].astype(int)
    top_df["authoritative_track_default"] = "retrieval"
    top_df["is_authoritative_track"] = top_df["track"].eq("retrieval")
    top_df["candidate_formula"] = top_df.get("candidate_formula", np.nan)
    top_df["predicted_exact_mass"] = top_df.get("predicted_exact_mass", np.nan)
    top_df["mass_match_score"] = top_df.get("mass_match_score", np.nan)
    top_df["mean_log_prob"] = top_df.get("mean_log_prob", np.nan)
    top_df["generation_support_fraction"] = top_df.get("generation_support_fraction", np.nan)
    top_df["candidate_latent_cosine"] = top_df.get("candidate_latent_cosine", np.nan)

    columns = [
        "feature_key",
        "track",
        "rank",
        "authoritative_track_default",
        "comp_id",
        "smiles",
        "name",
        "latent_cosine",
        "gate_used",
        "abs_mass_diff",
        "confidence_tier",
        "recommended",
        "candidate_formula",
        "retrieval_score",
        "fragment_match_count",
        "fragment_matched_peak_fraction",
        "fragment_explained_intensity_fraction",
        "predicted_exact_mass",
        "mass_match_score",
        "mean_log_prob",
        "generation_support_fraction",
        "candidate_latent_cosine",
        "merged_rank",
        "is_authoritative_track",
        "normalized_merge_score",
        "selfies_score",
        "policy_rank",
    ]
    for col in columns:
        if col not in top_df.columns:
            top_df[col] = np.nan
    return top_df[columns].sort_values(["feature_key", "merged_rank"], kind="mergesort").reset_index(drop=True)


def select_selfies_rows_per_feature(selfies_merge_df: pd.DataFrame, max_per_feature: int) -> pd.DataFrame:
    if selfies_merge_df.empty:
        return selfies_merge_df
    return (
        selfies_merge_df.sort_values(["feature_key", "source_rank"], kind="mergesort")
        .groupby("feature_key", as_index=False, group_keys=False)
        .head(int(max_per_feature))
        .reset_index(drop=True)
    )


def build_policy_retrieval_only(
    feature_keys: Iterable[str],
    retrieval_merge_df: pd.DataFrame,
) -> pd.DataFrame:
    chunks: List[pd.DataFrame] = []
    grouped = retrieval_merge_df.sort_values(["feature_key", "source_rank", "smiles"], kind="mergesort").groupby(
        "feature_key", sort=False
    )
    for feature_key in feature_keys:
        part = grouped.get_group(feature_key).copy() if feature_key in grouped.groups else pd.DataFrame()
        if not part.empty:
            chunks.append(part)
    if not chunks:
        return pd.DataFrame()
    return pd.concat(chunks, ignore_index=True)


def build_policy_selfies_only(
    feature_keys: Iterable[str],
    selfies_merge_df: pd.DataFrame,
) -> pd.DataFrame:
    chunks: List[pd.DataFrame] = []
    grouped = selfies_merge_df.sort_values(
        ["feature_key", "source_rank", "smiles"], kind="mergesort"
    ).groupby("feature_key", sort=False)
    for feature_key in feature_keys:
        part = grouped.get_group(feature_key).copy() if feature_key in grouped.groups else pd.DataFrame()
        if not part.empty:
            chunks.append(part)
    if not chunks:
        return pd.DataFrame()
    return pd.concat(chunks, ignore_index=True)


def build_policy_pin_high_rank1_joint_rest(
    feature_keys: Iterable[str],
    retrieval_merge_df: pd.DataFrame,
    selfies_cap_df: pd.DataFrame,
    rank1_tier_map: Dict[str, str],
) -> pd.DataFrame:
    chunks: List[pd.DataFrame] = []
    retr_group = retrieval_merge_df.groupby("feature_key", sort=False)
    self_group = selfies_cap_df.groupby("feature_key", sort=False)
    empty_retr = retrieval_merge_df.iloc[0:0].copy()
    empty_self = selfies_cap_df.iloc[0:0].copy()

    for feature_key in feature_keys:
        retr = retr_group.get_group(feature_key).copy() if feature_key in retr_group.groups else empty_retr.copy()
        selfs = self_group.get_group(feature_key).copy() if feature_key in self_group.groups else empty_self.copy()
        tier = str(rank1_tier_map.get(feature_key, "")).lower()

        if tier == "high" and not retr.empty:
            pinned = retr[retr["source_rank"] == 1].head(1).copy()
            remainder = pd.concat([retr[retr["source_rank"] != 1], selfs], ignore_index=True)
            sorted_rest = joint_sort_rows(remainder)
            ranked = pd.concat([pinned, sorted_rest], ignore_index=True)
        else:
            ranked = joint_sort_rows(pd.concat([retr, selfs], ignore_index=True))
        if not ranked.empty:
            chunks.append(ranked)

    if not chunks:
        return pd.DataFrame()
    return pd.concat(chunks, ignore_index=True)


def build_policy_high_protected_medium_low_joint(
    feature_keys: Iterable[str],
    retrieval_merge_df: pd.DataFrame,
    selfies_cap_df: pd.DataFrame,
    rank1_tier_map: Dict[str, str],
) -> pd.DataFrame:
    chunks: List[pd.DataFrame] = []
    retr_group = retrieval_merge_df.groupby("feature_key", sort=False)
    self_group = selfies_cap_df.groupby("feature_key", sort=False)
    empty_retr = retrieval_merge_df.iloc[0:0].copy()
    empty_self = selfies_cap_df.iloc[0:0].copy()

    for feature_key in feature_keys:
        retr = retr_group.get_group(feature_key).copy() if feature_key in retr_group.groups else empty_retr.copy()
        selfs = self_group.get_group(feature_key).copy() if feature_key in self_group.groups else empty_self.copy()
        tier = str(rank1_tier_map.get(feature_key, "")).lower()
        if tier == "high":
            ranked = retr.sort_values(["source_rank", "smiles"], kind="mergesort")
        else:
            ranked = joint_sort_rows(pd.concat([retr, selfs], ignore_index=True))
        if not ranked.empty:
            chunks.append(ranked)

    if not chunks:
        return pd.DataFrame()
    return pd.concat(chunks, ignore_index=True)


def build_policy_retrieval_rank1_plus_best_selfies(
    feature_keys: Iterable[str],
    retrieval_merge_df: pd.DataFrame,
    selfies_merge_df: pd.DataFrame,
) -> pd.DataFrame:
    chunks: List[pd.DataFrame] = []
    retr_group = retrieval_merge_df.groupby("feature_key", sort=False)
    self_group = selfies_merge_df.groupby("feature_key", sort=False)
    empty_retr = retrieval_merge_df.iloc[0:0].copy()
    empty_self = selfies_merge_df.iloc[0:0].copy()

    for feature_key in feature_keys:
        retr = retr_group.get_group(feature_key).copy() if feature_key in retr_group.groups else empty_retr.copy()
        selfs = self_group.get_group(feature_key).copy() if feature_key in self_group.groups else empty_self.copy()
        ranked_parts: List[pd.DataFrame] = []

        retr = retr.sort_values(["source_rank", "smiles"], kind="mergesort")
        selfs = selfs.sort_values(["source_rank", "smiles"], kind="mergesort")

        if not retr.empty:
            rank1 = retr[retr["source_rank"] == 1].head(1).copy()
            if not rank1.empty:
                ranked_parts.append(rank1)

        if not selfs.empty:
            selfie_eligible = selfs[
                (pd.to_numeric(selfs["abs_mass_error_da"], errors="coerce") <= 1.0)
                | (pd.to_numeric(selfs["selfies_score"], errors="coerce") >= 0.75)
            ].copy()
            if not selfie_eligible.empty:
                seen_smiles = set()
                if ranked_parts:
                    seen_smiles = set(str(x) for x in ranked_parts[0]["smiles"].astype(str).tolist())
                for row in selfie_eligible.itertuples(index=False):
                    if str(row.smiles) not in seen_smiles:
                        ranked_parts.append(pd.DataFrame([row._asdict()]))
                        break

        if not retr.empty:
            remaining_retr = retr[retr["source_rank"] != 1].copy()
            ranked_parts.append(remaining_retr)

        if ranked_parts:
            ranked = pd.concat(ranked_parts, ignore_index=True)
            chunks.append(ranked)

    if not chunks:
        return pd.DataFrame()
    return pd.concat(chunks, ignore_index=True)


# Keep benchmark/predictor merge behavior in lockstep by delegating shared helpers.
compute_selfies_score = merge_policy.compute_selfies_score
normalize_score = merge_policy.normalize_score
rank_selfies_pool = merge_policy.rank_selfies_pool
compare_selfies_pool_to_report = merge_policy.compare_selfies_pool_to_report
map_retrieval_rows_for_merge = merge_policy.map_retrieval_rows_for_merge
map_selfies_rows_for_merge = merge_policy.map_selfies_rows_for_merge
joint_sort_rows = merge_policy.joint_sort_rows
finalize_policy_rows = merge_policy.finalize_policy_rows
to_merged_predictions_csv_schema = merge_policy.to_merged_predictions_csv_schema
select_selfies_rows_per_feature = merge_policy.select_selfies_rows_per_feature
build_policy_retrieval_only = merge_policy.build_policy_retrieval_only
build_policy_selfies_only = merge_policy.build_policy_selfies_only
build_policy_pin_high_rank1_joint_rest = merge_policy.build_policy_pin_high_rank1_joint_rest
build_policy_top5_selfies_then_joint_rest = merge_policy.build_policy_top5_selfies_then_joint_rest
build_policy_high_protected_medium_low_joint = merge_policy.build_policy_high_protected_medium_low_joint
build_policy_retrieval_rank1_plus_best_selfies = merge_policy.build_policy_retrieval_rank1_plus_best_selfies


def resolve_subprocess_python() -> Path:
    candidates = []
    virtual_env = os.environ.get("VIRTUAL_ENV", "").strip()
    if virtual_env:
        candidates.append(Path(virtual_env) / "bin" / "python")
    repo_python = REPO_ROOT / ".venv" / "bin" / "python"
    candidates.append(repo_python)
    candidates.append(Path(sys.executable))
    candidates.append(Path(sys.prefix) / "bin" / "python")

    for candidate in candidates:
        if candidate.exists():
            return candidate.absolute()
    return Path(sys.executable).absolute()


def run_control_eval(
    eval_script: Path,
    policy_dir: Path,
    merged_csv: Path,
    known_x: Path,
    known_x_all: Path,
    control_lib_id: int,
) -> tuple[Path, Path]:
    cmd = [
        str(resolve_subprocess_python()),
        str(eval_script),
        "--prediction_dir",
        str(policy_dir),
        "--merged_csv",
        str(merged_csv),
        "--output_dir",
        str(policy_dir),
        "--known_x",
        str(known_x),
        "--known_x_all",
        str(known_x_all),
        "--control_lib_id",
        str(int(control_lib_id)),
    ]
    print("Running control eval:", " ".join(cmd))
    completed = subprocess.run(cmd, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"Control eval failed for {merged_csv} with code {completed.returncode}")
    return policy_dir / "control_eval_summary.json", policy_dir / "control_eval_per_xid.csv"


def ensure_track_metric_summary(
    eval_script: Path,
    prediction_dir: Path,
    summary_json_path: Path,
    known_x: Path,
    known_x_all: Path,
    control_lib_id: int,
    denovo_csv: Path | None = None,
    retrieval_csv: Path | None = None,
    merged_csv: Path | None = None,
) -> Path:
    if summary_json_path.exists():
        return summary_json_path

    cmd = [
        str(resolve_subprocess_python()),
        str(eval_script),
        "--prediction_dir",
        str(prediction_dir),
        "--output_dir",
        str(prediction_dir),
        "--known_x",
        str(known_x),
        "--known_x_all",
        str(known_x_all),
        "--control_lib_id",
        str(int(control_lib_id)),
    ]
    if denovo_csv is not None and denovo_csv.exists():
        cmd.extend(["--denovo_csv", str(denovo_csv)])
    if retrieval_csv is not None and retrieval_csv.exists():
        cmd.extend(["--retrieval_csv", str(retrieval_csv)])
    if merged_csv is not None and merged_csv.exists():
        cmd.extend(["--merged_csv", str(merged_csv)])

    print("Generating missing baseline control eval:", " ".join(cmd))
    completed = subprocess.run(cmd, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"Baseline control eval failed for {prediction_dir} with code {completed.returncode}"
        )
    if not summary_json_path.exists():
        raise FileNotFoundError(f"Expected baseline control eval summary at {summary_json_path}")
    return summary_json_path


def build_comparison_payload(
    policy_results: List[PolicyResult],
    retrieval_baseline_metric: dict,
    selfies_baseline_metric: dict,
    pool_alignment: dict,
) -> dict:
    all_metrics = []
    for result in policy_results:
        metric = dict(result.metric)
        metric["policy"] = result.name
        metric["diagnostic"] = bool(result.diagnostic)
        metric["is_blend_policy"] = result.name in BLEND_POLICY_NAMES
        all_metrics.append(metric)

    core_metrics = [m for m in all_metrics if not m["diagnostic"]]
    blend_metrics = [m for m in core_metrics if m["is_blend_policy"]]

    winner_any = max(core_metrics, key=score_tuple) if core_metrics else None
    winner_blend = max(blend_metrics, key=score_tuple) if blend_metrics else None

    def promotable(metric: dict | None) -> bool:
        if metric is None:
            return False
        return int(metric.get("n_controls_covered", 0)) >= int(PROMOTION_GATES["min_n_controls_covered"])

    def passes_keep_gates(metric: dict | None) -> bool:
        if metric is None:
            return False
        return (
            int(metric.get("n_controls_covered", 0)) >= int(PROMOTION_GATES["min_n_controls_covered"])
            and float(metric.get("mean_best_tanimoto", 0.0)) >= float(PROMOTION_GATES["min_mean_best_tanimoto"])
            and int(metric.get("n_topk_tanimoto_ge_0.7", 0))
            >= int(PROMOTION_GATES["min_n_topk_tanimoto_ge_0p7"])
        )

    keep_selfies_track = passes_keep_gates(winner_blend)
    strong_success = bool(
        winner_blend
        and (
            float(winner_blend.get("mean_best_tanimoto", 0.0))
            > float(STRONG_SUCCESS_GATES["selfies_mean_best_tanimoto_baseline"])
            or int(winner_blend.get("n_topk_tanimoto_ge_0.7", 0))
            > int(STRONG_SUCCESS_GATES["selfies_n_topk_tanimoto_ge_0p7_baseline"])
        )
    )

    diagnostics = {
        metric["policy"]: {
            "actual": {
                "n_controls_covered": int(metric.get("n_controls_covered", 0)),
                "mean_best_tanimoto": float(metric.get("mean_best_tanimoto", 0.0)),
            },
            "plan_reference": PLAN_HARD_WINDOW_REFERENCE.get(metric["policy"], {}),
        }
        for metric in all_metrics
        if metric["diagnostic"]
    }

    decision = {
        "winning_policy_overall": winner_any,
        "winning_blend_policy": winner_blend,
        "winning_blend_promotable": promotable(winner_blend),
        "keep_selfies_track": keep_selfies_track,
        "strong_success": strong_success,
        "close_reason": None if keep_selfies_track else "no_merged_policy_cleared_followup_gates",
    }
    return {
        "policy_metrics": all_metrics,
        "retrieval_baseline": retrieval_baseline_metric,
        "selfies_baseline": selfies_baseline_metric,
        "pool_alignment_check": pool_alignment,
        "promotion_gates": PROMOTION_GATES,
        "strong_success_gate_reference": STRONG_SUCCESS_GATES,
        "diagnostic_hard_windows": diagnostics,
        "decision": decision,
    }


def policy_metric_by_name(comparison_payload: dict, policy_name: str) -> dict:
    for metric in comparison_payload.get("policy_metrics", []):
        if str(metric.get("policy")) == str(policy_name):
            return dict(metric)
    return {}


def build_known_x_hitmax_reference_payload(
    comparison_payload: dict,
    coverage_reference_payload: dict,
    coverage_reference_path: Path,
    output_dir: Path,
    selfies_dir: Path,
    retrieval_dir: Path,
    selfies_max_per_feature: int,
    lane_policy_name: str | None = None,
) -> dict:
    hitmax_ceiling = policy_metric_by_name(comparison_payload, HITMAX_REFERENCE_CEILING_POLICY)
    winning_blend = policy_metric_by_name(comparison_payload, lane_policy_name) if lane_policy_name else {}
    if not winning_blend:
        winning_blend = dict((comparison_payload.get("decision") or {}).get("winning_blend_policy") or {})
    coverage_reference = dict(
        ((coverage_reference_payload.get("current_winner") or {}).get("merged")) or {}
    )
    pin_high_cap_candidate = policy_metric_by_name(comparison_payload, HITMAX_REFERENCE_COVERAGE_POLICY)

    hitmax_gate_pass = bool(
        winning_blend
        and int(winning_blend.get("n_controls_covered", 0)) >= int(HITMAX_PROMOTION_GATES["min_n_controls_covered"])
        and (
            float(winning_blend.get("mean_best_tanimoto", 0.0)) > float(hitmax_ceiling.get("mean_best_tanimoto", 0.0))
            or int(winning_blend.get("n_topk_tanimoto_ge_0.7", 0))
            > int(hitmax_ceiling.get("n_topk_tanimoto_ge_0.7", 0))
        )
    )
    coverage_followup_pass = bool(
        pin_high_cap_candidate
        and coverage_reference
        and int(pin_high_cap_candidate.get("n_controls_covered", 0))
        >= int(coverage_reference.get("n_controls_covered", 0))
        and (
            float(pin_high_cap_candidate.get("mean_best_tanimoto", 0.0))
            > float(coverage_reference.get("mean_best_tanimoto", 0.0))
            or int(pin_high_cap_candidate.get("n_topk_tanimoto_ge_0.7", 0))
            > int(coverage_reference.get("n_topk_tanimoto_ge_0.7", 0))
            or int(pin_high_cap_candidate.get("n_exact_match", 0))
            > int(coverage_reference.get("n_exact_match", 0))
        )
    )

    if hitmax_gate_pass:
        recommendation = "promote_merged_hitmax_candidate"
    elif coverage_followup_pass:
        recommendation = "keep_pin_high_rank1_joint_rest_cap10_as_coverage_preserving_followup_candidate"
    else:
        recommendation = "leave_selfies_only_5da_as_hitmax_ceiling"

    return {
        "reference_name": "x209_known_x_hitmax_reference",
        "status": "research_hitmax_branch",
        "artifacts": {
            "output_dir": str(output_dir),
            "selfies_dir": str(selfies_dir),
            "retrieval_dir": str(retrieval_dir),
            "comparison_summary_json": str(output_dir / "comparison_summary.json"),
            "coverage_reference_json": str(coverage_reference_path),
        },
        "benchmark_config": {
            "selfies_max_per_feature": int(selfies_max_per_feature),
            "top_k": 10,
        },
        "hitmax_ceiling": {
            "policy": HITMAX_REFERENCE_CEILING_POLICY,
            "metric": hitmax_ceiling,
        },
        "coverage_preserving_reference": {
            "policy": coverage_reference.get("policy", HITMAX_REFERENCE_COVERAGE_POLICY),
            "metric": coverage_reference,
        },
        "winning_blend_candidate": winning_blend,
        "coverage_preserving_followup_candidate": pin_high_cap_candidate,
        "promotion_rule": {
            "promote_merged_hitmax_candidate_if": (
                "mean_best_tanimoto improves above the hit-max ceiling or "
                "n_topk_tanimoto_ge_0.7 improves above the hit-max ceiling, "
                "while n_controls_covered >= 48"
            ),
            "coverage_preserving_followup_rule": (
                "If no merged candidate clears the hit-max ceiling, keep pin_high_rank1_joint_rest "
                "only when it improves over the seeded coverage-preserving reference."
            ),
        },
        "decision": {
            "hitmax_gate_pass": hitmax_gate_pass,
            "coverage_preserving_followup_pass": coverage_followup_pass,
            "recommended_action": recommendation,
        },
    }


def comparison_table_rows(comparison_payload: dict) -> list[dict]:
    rows = [
        {
            "policy": "baseline_retrieval",
            "diagnostic": False,
            **comparison_payload["retrieval_baseline"],
        },
        {
            "policy": "baseline_selfies",
            "diagnostic": False,
            **comparison_payload["selfies_baseline"],
        },
    ]
    rows.extend(comparison_payload.get("policy_metrics", []))
    df = pd.DataFrame(rows)
    if df.empty:
        return []
    if "mean_best_tanimoto" in df.columns:
        df["mean_best_tanimoto"] = pd.to_numeric(df["mean_best_tanimoto"], errors="coerce")
    if "n_topk_tanimoto_ge_0.7" in df.columns:
        df["n_topk_tanimoto_ge_0.7"] = pd.to_numeric(df["n_topk_tanimoto_ge_0.7"], errors="coerce")
    if "n_controls_covered" in df.columns:
        df["n_controls_covered"] = pd.to_numeric(df["n_controls_covered"], errors="coerce")
    df = df.sort_values(
        ["diagnostic", "mean_best_tanimoto", "n_topk_tanimoto_ge_0.7", "policy"],
        ascending=[True, False, False, True],
        kind="mergesort",
    )
    return df.to_dict(orient="records")


def write_followup_memo(markdown_path: Path, comparison_payload: dict) -> None:
    decision = comparison_payload["decision"]
    winner = decision.get("winning_blend_policy") or {}
    keep_close = "KEEP" if decision.get("keep_selfies_track") else "CLOSE"
    lines = [
        "# X209 SELFIES Follow-Up Merge Benchmark",
        "",
        "## Outcome",
        f"- Decision: **{keep_close}** SELFIES track.",
        f"- Winning blend policy: `{winner.get('policy', 'none')}`.",
        f"- Winning blend mean_best_tanimoto: `{winner.get('mean_best_tanimoto', 'n/a')}`.",
        f"- Winning blend n_topk_tanimoto_ge_0.7: `{winner.get('n_topk_tanimoto_ge_0.7', 'n/a')}`.",
        f"- Winning blend n_controls_covered: `{winner.get('n_controls_covered', 'n/a')}`.",
        "",
        "## Baseline Comparison",
        (
            "- Retrieval baseline: "
            f"mean_best_tanimoto={comparison_payload['retrieval_baseline'].get('mean_best_tanimoto')}, "
            f"n_topk_tanimoto_ge_0.7={comparison_payload['retrieval_baseline'].get('n_topk_tanimoto_ge_0.7')}, "
            f"n_controls_covered={comparison_payload['retrieval_baseline'].get('n_controls_covered')}."
        ),
        (
            "- SELFIES baseline: "
            f"mean_best_tanimoto={comparison_payload['selfies_baseline'].get('mean_best_tanimoto')}, "
            f"n_topk_tanimoto_ge_0.7={comparison_payload['selfies_baseline'].get('n_topk_tanimoto_ge_0.7')}, "
            f"n_controls_covered={comparison_payload['selfies_baseline'].get('n_controls_covered')}."
        ),
        "",
        "## Hard-Window Diagnostics",
    ]
    for policy_name, payload in sorted(comparison_payload.get("diagnostic_hard_windows", {}).items()):
        actual = payload.get("actual", {})
        ref = payload.get("plan_reference", {})
        lines.append(
            (
                f"- `{policy_name}`: actual coverage={actual.get('n_controls_covered')}, "
                f"actual mean_best_tanimoto={actual.get('mean_best_tanimoto')}; "
                f"plan reference coverage={ref.get('n_controls_covered')}, "
                f"plan reference mean_best_tanimoto={ref.get('mean_best_tanimoto')}."
            )
        )
    lines.extend(
        [
            "",
            "## Decision Note",
            (
                "- This sprint used frozen artifacts only and benchmarked ranking/merge policies without retraining "
                "or regeneration."
            ),
            (
                "- Keep/close gate applied to winning blend policy: "
                f"mean_best_tanimoto >= {PROMOTION_GATES['min_mean_best_tanimoto']}, "
                f"n_topk_tanimoto_ge_0.7 >= {PROMOTION_GATES['min_n_topk_tanimoto_ge_0p7']}, "
                f"n_controls_covered >= {PROMOTION_GATES['min_n_controls_covered']}."
            ),
        ]
    )
    markdown_path.write_text("\n".join(lines) + "\n")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Benchmark fixed SELFIES+retrieval merge policies.")
    parser.add_argument("--selfies-dir", type=str, default="models/x209_selfies_denovo")
    parser.add_argument("--retrieval-dir", type=str, default="models/x209_smiles_predictions_baseline_refresh")
    parser.add_argument("--output-dir", type=str, default="models/x209_selfies_followup")
    parser.add_argument(
        "--artifact-prefix",
        type=str,
        default="x209",
        help="Filename prefix for retrieval/selfies/merged artifacts used within this followup run.",
    )
    parser.add_argument("--selfies-max-per-feature", type=int, default=3)
    parser.add_argument("--known-x", type=str, default="known_X.csv")
    parser.add_argument("--known-x-all", type=str, default="known_X_all_available_structures_simple.csv")
    parser.add_argument("--control-lib-id", type=int, default=209)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument(
        "--coverage-reference-json",
        type=str,
        default="models/x209_reference_baseline.json",
        help="Coverage-preserving reference manifest used when emitting a hit-max reference payload.",
    )
    parser.add_argument(
        "--known-x-hitmax-reference-json",
        type=str,
        default="",
        help="Optional output path for writing a hit-max reference manifest. Disabled by default.",
    )
    return parser


def main() -> None:
    args = get_parser().parse_args()

    repo_root = REPO_ROOT
    selfies_dir = Path(args.selfies_dir).expanduser().resolve()
    retrieval_dir = Path(args.retrieval_dir).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    policies_root = output_dir / "policies"
    policies_root.mkdir(parents=True, exist_ok=True)
    artifact_prefix = str(args.artifact_prefix).strip() or "x209"

    candidate_metrics_csv = resolve_existing_artifact(
        selfies_dir,
        [
            artifact_filename(artifact_prefix, "candidate_metrics_all", ".csv"),
            "selfies_candidate_metrics_all.csv",
            "candidate_metrics_all.csv",
        ],
    )
    selfies_report_csv = resolve_existing_artifact(
        selfies_dir,
        [
            artifact_filename(artifact_prefix, "denovo_report", ".csv"),
            "selfies_denovo_report.csv",
            "denovo_report.csv",
        ],
    )
    selfies_eval_summary_json = selfies_dir / "control_eval_summary.json"
    retrieval_report_csv = resolve_existing_artifact(
        retrieval_dir,
        [
            artifact_filename(artifact_prefix, "retrieval_report", ".csv"),
            "retrieval_report.csv",
            "x209_retrieval_report.csv",
        ],
    )
    retrieval_eval_summary_json = retrieval_dir / "control_eval_summary.json"
    retrieval_merged_csv = resolve_existing_artifact(
        retrieval_dir,
        [
            artifact_filename(artifact_prefix, "merged_predictions", ".csv"),
            "merged_predictions.csv",
            "x209_merged_predictions.csv",
        ],
    )
    retrieval_meta_csv = resolve_existing_artifact(
        retrieval_dir,
        [
            artifact_filename(artifact_prefix, "feature_meta", ".csv"),
            "feature_meta.csv",
            "x209_feature_meta.csv",
        ],
    )

    eval_script = repo_root / "evaluate_denovo_controls.py"
    known_x = Path(args.known_x).expanduser().resolve()
    known_x_all = Path(args.known_x_all).expanduser().resolve()
    selfies_eval_summary_json = ensure_track_metric_summary(
        eval_script=eval_script,
        prediction_dir=selfies_dir,
        summary_json_path=selfies_eval_summary_json,
        known_x=known_x,
        known_x_all=known_x_all,
        control_lib_id=int(args.control_lib_id),
        denovo_csv=selfies_report_csv,
    )
    retrieval_eval_summary_json = ensure_track_metric_summary(
        eval_script=eval_script,
        prediction_dir=retrieval_dir,
        summary_json_path=retrieval_eval_summary_json,
        known_x=known_x,
        known_x_all=known_x_all,
        control_lib_id=int(args.control_lib_id),
        retrieval_csv=retrieval_report_csv,
        merged_csv=retrieval_merged_csv,
    )

    candidate_metrics_df = pd.read_csv(candidate_metrics_csv)
    retrieval_df = pd.read_csv(retrieval_report_csv)
    selfies_report_df = pd.read_csv(selfies_report_csv)
    feature_meta_df = pd.read_csv(retrieval_meta_csv) if retrieval_meta_csv.exists() else pd.DataFrame()

    retrieval_df = retrieval_df.copy()
    retrieval_df["rank"] = pd.to_numeric(retrieval_df["rank"], errors="coerce").fillna(9999).astype(int)
    retrieval_top10_df = retrieval_df[retrieval_df["rank"] <= int(args.top_k)].copy()
    retrieval_top10_df["feature_key"] = retrieval_top10_df["feature_key"].astype(str)

    candidate_metrics_df = candidate_metrics_df.copy()
    candidate_metrics_df["feature_key"] = candidate_metrics_df["feature_key"].astype(str)
    candidate_metrics_df["predicted_smiles"] = candidate_metrics_df["predicted_smiles"].astype(str)
    candidate_metrics_df["abs_mass_error_da"] = pd.to_numeric(
        candidate_metrics_df["abs_mass_error_da"], errors="coerce"
    )

    feature_keys = sorted(set(retrieval_top10_df["feature_key"].astype(str).tolist()))
    rank1_df = retrieval_top10_df[retrieval_top10_df["rank"] == 1].copy()
    rank1_tier_map = {
        str(row.feature_key): str(row.confidence_tier)
        for row in rank1_df[["feature_key", "confidence_tier"]].itertuples(index=False)
    }

    selfies_pool_5da = rank_selfies_pool(candidate_metrics_df, tolerance_da=5.0)
    selfies_pool_1da = rank_selfies_pool(candidate_metrics_df, tolerance_da=1.0)
    selfies_pool_0p5da = rank_selfies_pool(candidate_metrics_df, tolerance_da=0.5)
    pool_alignment = compare_selfies_pool_to_report(selfies_pool_5da, selfies_report_df)

    retrieval_merge_df = map_retrieval_rows_for_merge(retrieval_top10_df)
    selfies_merge_5da = map_selfies_rows_for_merge(selfies_pool_5da, feature_meta_df=feature_meta_df)
    selfies_merge_1da = map_selfies_rows_for_merge(selfies_pool_1da, feature_meta_df=feature_meta_df)
    selfies_merge_0p5da = map_selfies_rows_for_merge(selfies_pool_0p5da, feature_meta_df=feature_meta_df)
    selfies_merge_cap_5da = select_selfies_rows_per_feature(
        selfies_merge_5da, max_per_feature=int(args.selfies_max_per_feature)
    )

    policy_builders: Dict[str, tuple[bool, Callable[[], pd.DataFrame]]] = {
        "retrieval_only": (False, lambda: build_policy_retrieval_only(feature_keys, retrieval_merge_df)),
        "selfies_only_5da": (False, lambda: build_policy_selfies_only(feature_keys, selfies_merge_5da)),
        "pin_high_rank1_joint_rest": (
            False,
            lambda: build_policy_pin_high_rank1_joint_rest(
                feature_keys,
                retrieval_merge_df,
                selfies_merge_cap_5da,
                rank1_tier_map,
            ),
        ),
        "top5_selfies_then_joint_rest": (
            False,
            lambda: build_policy_top5_selfies_then_joint_rest(
                feature_keys,
                retrieval_merge_df,
                selfies_merge_cap_5da,
            ),
        ),
        "high_protected_medium_low_joint": (
            False,
            lambda: build_policy_high_protected_medium_low_joint(
                feature_keys,
                retrieval_merge_df,
                selfies_merge_cap_5da,
                rank1_tier_map,
            ),
        ),
        "retrieval_rank1_plus_best_selfies": (
            False,
            lambda: build_policy_retrieval_rank1_plus_best_selfies(
                feature_keys,
                retrieval_merge_df,
                selfies_merge_5da,
            ),
        ),
        "selfies_only_1da_diag": (True, lambda: build_policy_selfies_only(feature_keys, selfies_merge_1da)),
        "selfies_only_0p5da_diag": (True, lambda: build_policy_selfies_only(feature_keys, selfies_merge_0p5da)),
    }
    ordered_policy_names = CORE_POLICY_ORDER + DIAGNOSTIC_POLICY_ORDER

    policy_results: List[PolicyResult] = []
    for policy_name in ordered_policy_names:
        diagnostic, builder = policy_builders[policy_name]
        print(f"\n=== Policy: {policy_name} ===")
        raw_ranked = builder()
        ranked_with_marks = finalize_policy_rows(raw_ranked, topk=int(args.top_k))

        policy_dir = policies_root / policy_name
        policy_dir.mkdir(parents=True, exist_ok=True)
        ranked_csv = policy_dir / "ranked_policy_rows.csv"
        ranked_with_marks.to_csv(ranked_csv, index=False)

        merged_df = to_merged_predictions_csv_schema(ranked_with_marks, topk=int(args.top_k))
        merged_csv = policy_dir / artifact_filename(
            artifact_prefix, "merged_predictions", ".csv"
        )
        merged_df.to_csv(merged_csv, index=False)
        merged_df.to_csv(policy_dir / "merged_predictions.csv", index=False)
        print(f"  Ranked rows: {len(ranked_with_marks)}")
        print(f"  Merged rows (top-{args.top_k}): {len(merged_df)}")

        summary_json, per_xid_csv = run_control_eval(
            eval_script=eval_script,
            policy_dir=policy_dir,
            merged_csv=merged_csv,
            known_x=known_x,
            known_x_all=known_x_all,
            control_lib_id=int(args.control_lib_id),
        )
        metric = metric_from_policy_eval(summary_json)
        print(
            "  Control metrics:",
            {
                "mean_best_tanimoto": metric.get("mean_best_tanimoto"),
                "n_topk_tanimoto_ge_0.7": metric.get("n_topk_tanimoto_ge_0.7"),
                "n_controls_covered": metric.get("n_controls_covered"),
            },
        )
        policy_results.append(
            PolicyResult(
                name=policy_name,
                diagnostic=diagnostic,
                ranked_csv=ranked_csv,
                merged_csv=merged_csv,
                control_eval_summary_json=summary_json,
                control_eval_per_xid_csv=per_xid_csv,
                metric=metric,
            )
        )

    retrieval_baseline_metric = load_track_metric(retrieval_eval_summary_json, "retrieval")
    selfies_baseline_metric = load_track_metric(selfies_eval_summary_json, "de_novo")
    comparison_payload = build_comparison_payload(
        policy_results=policy_results,
        retrieval_baseline_metric=retrieval_baseline_metric,
        selfies_baseline_metric=selfies_baseline_metric,
        pool_alignment=pool_alignment,
    )

    comparison_json = output_dir / "comparison_summary.json"
    comparison_json.write_text(json.dumps(comparison_payload, indent=2))
    comparison_table_csv = output_dir / "comparison_table.csv"
    pd.DataFrame(comparison_table_rows(comparison_payload)).to_csv(comparison_table_csv, index=False)
    memo_md = output_dir / "followup_memo.md"
    write_followup_memo(memo_md, comparison_payload=comparison_payload)

    hitmax_reference_json = None
    if int(args.selfies_max_per_feature) >= 10 and str(args.known_x_hitmax_reference_json).strip():
        coverage_reference_path = Path(args.coverage_reference_json).expanduser().resolve()
        if coverage_reference_path.exists():
            coverage_reference_payload = json.loads(coverage_reference_path.read_text())
            hitmax_reference_payload = build_known_x_hitmax_reference_payload(
                comparison_payload=comparison_payload,
                coverage_reference_payload=coverage_reference_payload,
                coverage_reference_path=coverage_reference_path,
                output_dir=output_dir,
                selfies_dir=selfies_dir,
                retrieval_dir=retrieval_dir,
                selfies_max_per_feature=int(args.selfies_max_per_feature),
            )
            hitmax_reference_json = Path(args.known_x_hitmax_reference_json).expanduser().resolve()
            hitmax_reference_json.parent.mkdir(parents=True, exist_ok=True)
            hitmax_reference_json.write_text(json.dumps(hitmax_reference_payload, indent=2))

    run_manifest = {
        "selfies_dir": str(selfies_dir),
        "retrieval_dir": str(retrieval_dir),
        "output_dir": str(output_dir),
        "artifact_prefix": artifact_prefix,
        "selfies_max_per_feature": int(args.selfies_max_per_feature),
        "top_k": int(args.top_k),
        "policy_order": ordered_policy_names,
        "policy_results": [
            {
                "policy": r.name,
                "diagnostic": bool(r.diagnostic),
                "ranked_csv": str(r.ranked_csv),
                "merged_csv": str(r.merged_csv),
                "control_eval_summary_json": str(r.control_eval_summary_json),
                "control_eval_per_xid_csv": str(r.control_eval_per_xid_csv),
                "metric": r.metric,
            }
            for r in policy_results
        ],
        "pool_alignment": pool_alignment,
        "comparison_summary_json": str(comparison_json),
        "comparison_table_csv": str(comparison_table_csv),
        "followup_memo_md": str(memo_md),
        "known_x_hitmax_reference_json": str(hitmax_reference_json) if hitmax_reference_json is not None else None,
    }
    (output_dir / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2))

    decision = comparison_payload.get("decision", {})
    print("\n=== Decision ===")
    print(json.dumps(decision, indent=2))
    print(f"Wrote: {comparison_json}")
    print(f"Wrote: {comparison_table_csv}")
    print(f"Wrote: {memo_md}")
    if hitmax_reference_json is not None:
        print(f"Wrote: {hitmax_reference_json}")
    print(f"Wrote: {output_dir / 'run_manifest.json'}")


if __name__ == "__main__":
    main()
