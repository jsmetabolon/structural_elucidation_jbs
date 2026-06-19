#!/usr/bin/env python3
"""Shared X209 retrieval/SELFIES merge policy helpers.

Pipeline context:
    Implements composite SELFIES scoring, per-feature deduplication, and
    multiple merge policies (retrieval-only, tier-pinned joint sort, etc.)
    used by benchmark and production rerank assembly.

Key I/O:
    ``rank_selfies_pool``: mass-window SELFIES ranking per feature.
    ``to_merged_predictions_csv_schema``: policy rows → deliverable CSV columns.
"""

from __future__ import annotations

from typing import Dict, Iterable, List

import numpy as np
import pandas as pd

SELFIES_SCORE_FORMULA_VERSION = (
    "0.55*mass_match_score + 0.25*generation_support_fraction + 0.20*sigmoid(mean_log_prob+1.0)"
)


# ── Scoring primitives ───────────────────────────────────────────────
def sigmoid(x: float | pd.Series) -> float | pd.Series:
    """Element-wise logistic transform."""
    return 1.0 / (1.0 + np.exp(-x))


def compute_selfies_score(
    mass_match_score: float | pd.Series,
    generation_support_fraction: float | pd.Series,
    mean_log_prob: float | pd.Series,
) -> float | pd.Series:
    """Composite score used for policy-level SELFIES ranking.

    See ``SELFIES_SCORE_FORMULA_VERSION`` for weight breakdown.
    """
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
    """Clip numeric scores to ``[0, 1]`` with NaN → 0."""
    return pd.to_numeric(series, errors="coerce").fillna(0.0).clip(lower=0.0, upper=1.0)


def rank_selfies_pool(candidate_metrics_df: pd.DataFrame, tolerance_da: float) -> pd.DataFrame:
    """Rank SELFIES pool for a given mass window and enforce per-feature SMILES uniqueness.

    Args:
        candidate_metrics_df: Per-candidate SELFIES metrics with mass errors.
        tolerance_da: Maximum absolute mass error (Da) to retain.

    Returns:
        Ranked pool with ``selfies_score``, ``normalized_merge_score``, ``selfies_rank``.
    """
    pool = candidate_metrics_df.copy()
    required_columns = {
        "feature_key": pd.Series(dtype=str),
        "predicted_smiles": pd.Series(dtype=str),
        "abs_mass_error_da": pd.Series(dtype=float),
        "mass_match_score": pd.Series(dtype=float),
        "generation_support_fraction": pd.Series(dtype=float),
        "mean_log_prob": pd.Series(dtype=float),
        "predicted_formula": pd.Series(dtype=object),
        "predicted_exact_mass": pd.Series(dtype=float),
        "candidate_latent_cosine": pd.Series(dtype=float),
    }
    for col_name, empty_series in required_columns.items():
        if col_name not in pool.columns:
            pool[col_name] = empty_series
    pool = pool[pd.to_numeric(pool["abs_mass_error_da"], errors="coerce") <= float(tolerance_da)].copy()
    if pool.empty:
        pool["selfies_score"] = pd.Series(dtype=float)
        pool["normalized_merge_score"] = pd.Series(dtype=float)
        pool["selfies_rank"] = pd.Series(dtype=int)
        return pool.reset_index(drop=True)

    pool["selfies_score"] = compute_selfies_score(
        pool["mass_match_score"],
        pool["generation_support_fraction"],
        pool["mean_log_prob"],
    )
    pool["normalized_merge_score"] = normalize_score(pool["selfies_score"])
    pool["abs_mass_error_da"] = pd.to_numeric(pool["abs_mass_error_da"], errors="coerce").fillna(np.inf)
    pool["generation_support_fraction"] = pd.to_numeric(pool["generation_support_fraction"], errors="coerce").fillna(
        0.0
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
    """Compare (feature_key, SMILES) pair sets between pool and report exports."""
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


# ── Track normalization for joint merge ────────────────────────────────
def map_retrieval_rows_for_merge(retrieval_df: pd.DataFrame) -> pd.DataFrame:
    """Add unified merge/sort columns for retrieval track rows."""
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
    out["feature_is_likely_isotopologue"] = out.get("feature_is_likely_isotopologue", False)
    out["feature_isotopologue_reason"] = out.get("feature_isotopologue_reason", None)
    out["isotope_element_hints"] = out.get("isotope_element_hints", np.nan)
    out["isotope_formula_consistent"] = out.get("isotope_formula_consistent", np.nan)
    out["isotope_formula_score"] = out.get("isotope_formula_score", np.nan)
    out["isotope_ratio_score"] = out.get("isotope_ratio_score", np.nan)
    out["isotope_audit_status"] = out.get("isotope_audit_status", None)
    out["isotope_signal_count"] = out.get("isotope_signal_count", np.nan)
    return out


def map_selfies_rows_for_merge(selfies_pool_df: pd.DataFrame, feature_meta_df: pd.DataFrame) -> pd.DataFrame:
    """Add unified merge/sort columns for SELFIES track rows."""
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
    out["feature_is_likely_isotopologue"] = out.get("feature_is_likely_isotopologue", False)
    out["feature_isotopologue_reason"] = out.get("feature_isotopologue_reason", None)
    out["isotope_element_hints"] = out.get("isotope_element_hints", np.nan)
    out["isotope_formula_consistent"] = out.get("isotope_formula_consistent", np.nan)
    out["isotope_formula_score"] = out.get("isotope_formula_score", np.nan)
    out["isotope_ratio_score"] = out.get("isotope_ratio_score", np.nan)
    out["isotope_audit_status"] = out.get("isotope_audit_status", None)
    out["isotope_signal_count"] = out.get("isotope_signal_count", np.nan)

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
    """Sort mixed retrieval/SELFIES rows by score, error, support, and track priority."""
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
    """Assign ``policy_rank`` and ``include_topk`` after policy-specific ordering."""
    if ranked_rows_df.empty:
        return ranked_rows_df
    out = ranked_rows_df.copy()
    out["policy_rank"] = out.groupby("feature_key").cumcount() + 1
    out["include_topk"] = out["policy_rank"] <= int(topk)
    return out


def to_merged_predictions_csv_schema(policy_rows_df: pd.DataFrame, topk: int = 10) -> pd.DataFrame:
    """Project policy rows into the merged predictions deliverable schema."""
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
        "feature_is_likely_isotopologue",
        "feature_isotopologue_reason",
        "isotope_element_hints",
        "isotope_formula_consistent",
        "isotope_formula_score",
        "isotope_ratio_score",
        "isotope_audit_status",
        "isotope_signal_count",
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
    if policy_rows_df.empty:
        return pd.DataFrame(columns=columns)

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

    for col in columns:
        if col not in top_df.columns:
            top_df[col] = np.nan
    return top_df[columns].sort_values(["feature_key", "merged_rank"], kind="mergesort").reset_index(drop=True)


def select_selfies_rows_per_feature(selfies_merge_df: pd.DataFrame, max_per_feature: int) -> pd.DataFrame:
    """Cap SELFIES rows per feature by ``source_rank``."""
    if selfies_merge_df.empty:
        return selfies_merge_df
    return (
        selfies_merge_df.sort_values(["feature_key", "source_rank"], kind="mergesort")
        .groupby("feature_key", as_index=False, group_keys=False)
        .head(int(max_per_feature))
        .reset_index(drop=True)
    )


# ── Merge policy builders ────────────────────────────────────────────
def build_policy_retrieval_only(
    feature_keys: Iterable[str],
    retrieval_merge_df: pd.DataFrame,
) -> pd.DataFrame:
    """Per-feature retrieval-only policy (no SELFIES rows)."""
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
    """Per-feature SELFIES-only policy (no retrieval rows)."""
    chunks: List[pd.DataFrame] = []
    grouped = selfies_merge_df.sort_values(["feature_key", "source_rank", "smiles"], kind="mergesort").groupby(
        "feature_key", sort=False
    )
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
    """Pin rank-1 retrieval when tier is high; joint-sort remainder with SELFIES."""
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
    """Keep retrieval-only ordering for high tier; joint-sort medium/low tiers."""
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
    """Place retrieval rank-1 first, then one eligible SELFIES hit, then remaining retrieval."""
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


def build_policy_top5_selfies_then_joint_rest(
    feature_keys: Iterable[str],
    retrieval_merge_df: pd.DataFrame,
    selfies_cap_df: pd.DataFrame,
) -> pd.DataFrame:
    """Front-load top five SELFIES rows, then joint-sort retrieval with remainder."""
    chunks: List[pd.DataFrame] = []
    retr_group = retrieval_merge_df.groupby("feature_key", sort=False)
    self_group = selfies_cap_df.groupby("feature_key", sort=False)
    empty_retr = retrieval_merge_df.iloc[0:0].copy()
    empty_self = selfies_cap_df.iloc[0:0].copy()

    for feature_key in feature_keys:
        retr = retr_group.get_group(feature_key).copy() if feature_key in retr_group.groups else empty_retr.copy()
        selfs = self_group.get_group(feature_key).copy() if feature_key in self_group.groups else empty_self.copy()
        selfs = selfs.sort_values(["source_rank", "smiles"], kind="mergesort").reset_index(drop=True)
        front_selfies = selfs.head(5).copy()
        remainder_selfies = selfs.iloc[len(front_selfies) :].copy()
        remainder = joint_sort_rows(pd.concat([retr, remainder_selfies], ignore_index=True))
        ranked = pd.concat([front_selfies, remainder], ignore_index=True)
        if not ranked.empty:
            chunks.append(ranked)

    if not chunks:
        return pd.DataFrame()
    return pd.concat(chunks, ignore_index=True)


__all__ = [
    "SELFIES_SCORE_FORMULA_VERSION",
    "compare_selfies_pool_to_report",
    "compute_selfies_score",
    "finalize_policy_rows",
    "joint_sort_rows",
    "map_retrieval_rows_for_merge",
    "map_selfies_rows_for_merge",
    "normalize_score",
    "rank_selfies_pool",
    "select_selfies_rows_per_feature",
    "build_policy_retrieval_only",
    "build_policy_selfies_only",
    "build_policy_pin_high_rank1_joint_rest",
    "build_policy_high_protected_medium_low_joint",
    "build_policy_retrieval_rank1_plus_best_selfies",
    "build_policy_top5_selfies_then_joint_rest",
    "to_merged_predictions_csv_schema",
]
