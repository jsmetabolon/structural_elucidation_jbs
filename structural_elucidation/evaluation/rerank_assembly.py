#!/usr/bin/env python3
"""Frozen-pool reranking for V3 retrieval-anchored de novo assembly candidates."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from structural_elucidation.evaluation.aggregate import best_by_tanimoto, metric_from_rows
from structural_elucidation.evaluation.assess import discover_profile_dirs, parse_lib_dir
from structural_elucidation.evaluation.failure_modes import (
    annotate_pool,
    load_x_controls,
    normalize_prediction_pool,
    pool_for_control,
)
from structural_elucidation.utils.domain import (
    dataframe_to_markdown,
    parse_libs,
    safe_float,
    safe_int,
    write_csv,
    write_json,
)

DEFAULT_V3_ROOT = Path("models/denovo_pipeline_v3/assembly_v3")
DEFAULT_OUTPUT_ROOT = Path("models/denovo_pipeline_v3/assembly_v3_rerank")
GOOD_TANIMOTO = 0.7


@dataclass(frozen=True)
class AssemblyRerankPolicy:
    name: str
    weights: dict[str, float]
    description: str


POLICIES: dict[str, AssemblyRerankPolicy] = {
    "assembly_score_baseline": AssemblyRerankPolicy(
        "assembly_score_baseline",
        {"assembly": 1.0},
        "Reproduce current V3 ordering by assembly_score.",
    ),
    "mass_anchor_balanced": AssemblyRerankPolicy(
        "mass_anchor_balanced",
        {"mass": 0.30, "anchor": 0.25, "assembly": 0.20, "adduct": 0.10, "operator": 0.10, "formula": 0.05},
        "Balance mass/adduct fit, anchor similarity, assembly score, operator support, and formula agreement.",
    ),
    "operator_supported": AssemblyRerankPolicy(
        "operator_supported",
        {"operator": 0.30, "support": 0.15, "mass": 0.20, "anchor": 0.15, "assembly": 0.10, "formula": 0.10},
        "Boost candidates supported by local edit operators and multiple generation modes.",
    ),
    "formula_anchor_guarded": AssemblyRerankPolicy(
        "formula_anchor_guarded",
        {"formula_guard": 0.30, "anchor": 0.25, "mass": 0.20, "operator": 0.10, "assembly": 0.10, "novelty": 0.05},
        "Promote only formula-compatible or strongly anchored candidates.",
    ),
    "novelty_guarded": AssemblyRerankPolicy(
        "novelty_guarded",
        {"novelty_guard": 0.25, "anchor": 0.20, "mass": 0.20, "operator": 0.15, "assembly": 0.10, "formula": 0.10},
        "Avoid near-exact database recovery while rejecting very remote weak-anchor candidates.",
    ),
    "rank_gap_rescue": AssemblyRerankPolicy(
        "rank_gap_rescue",
        {"mass": 0.20, "anchor": 0.20, "operator": 0.20, "formula_guard": 0.20, "support": 0.10, "rank_rescue": 0.10},
        "De-emphasize current rank and promote strong mass/anchor/formula/operator candidates outside top-10.",
    ),
}


def clip01(values: pd.Series) -> pd.Series:
    """Clip a numeric series to the [0, 1] interval."""
    """Weighted feature policy for assembly candidate reranking."""
    return pd.to_numeric(values, errors="coerce").fillna(0.0).clip(lower=0.0, upper=1.0)


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    """Read a numeric column with a default fill value."""
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default).astype(float)


def bool_score(df: pd.DataFrame, col: str) -> pd.Series:
    """Convert a boolean column to 0/1 float scores."""
    if col not in df.columns:
        return pd.Series(0.0, index=df.index, dtype=float)
    return df[col].fillna(False).astype(str).str.lower().isin({"true", "1", "yes", "y"}).astype(float)


def normalize_positive(series: pd.Series) -> pd.Series:
    """Min-max normalize positive values to [0, 1]."""
    values = pd.to_numeric(series, errors="coerce").fillna(0.0)
    max_value = float(values.max()) if len(values) else 0.0
    if max_value <= 0:
        return pd.Series(0.0, index=series.index, dtype=float)
    return (values / max_value).clip(lower=0.0, upper=1.0)


def operator_score(df: pd.DataFrame) -> pd.Series:
    """Map operator-family labels to a rerank prior score."""
    family = df.get("operator_family", pd.Series("", index=df.index)).fillna("").astype(str)
    modes = df.get("assembly_generation_modes", pd.Series("", index=df.index)).fillna("").astype(str)
    scores = pd.Series(0.20, index=df.index, dtype=float)
    scores = scores.mask(family.str.contains("mcs_rgroup_edit", regex=False), 1.00)
    scores = scores.mask(family.str.contains("mass_delta_directed", regex=False), 0.88)
    scores = scores.mask(family.str.contains("rgroup_growth", regex=False), 0.78)
    scores = scores.mask(family.str.contains("mass_constrained_fragment", regex=False), 0.72)
    scores = scores.mask(family.str.contains("formula_preserving", regex=False), 0.66)
    scores = scores.mask(family.str.contains("brics", regex=False), 0.55)
    multi_mode = modes.str.count(r"\|").fillna(0).astype(float).clip(lower=0, upper=3) / 3.0
    return (0.80 * scores + 0.20 * multi_mode).clip(lower=0.0, upper=1.0)


def add_score_features(pool: pd.DataFrame, tolerance_da: float) -> pd.DataFrame:
    """Add derived rerank feature columns to a candidate pool."""
    out = pool.copy()
    abs_err = numeric(out, "abs_mass_error_da", default=float(tolerance_da))
    primary_err = numeric(out, "primary_abs_mass_error_da", default=float(tolerance_da))
    assembly = numeric(out, "assembly_score", default=0.0)
    out["score_assembly"] = normalize_positive(assembly)
    if "mass_match_score" in out.columns:
        out["score_mass"] = clip01(numeric(out, "mass_match_score"))
    else:
        out["score_mass"] = (1.0 - abs_err / max(float(tolerance_da), 1e-9)).clip(lower=0.0, upper=1.0)
    out["score_primary_mass"] = (1.0 - primary_err / max(float(tolerance_da), 1e-9)).clip(lower=0.0, upper=1.0)
    out["score_adduct"] = 0.5 * out["score_mass"] + 0.5 * out["score_primary_mass"]
    out["score_anchor"] = clip01(
        numeric(out, "anchor_tanimoto", default=0.0).combine_first(numeric(out, "retrieval_tanimoto", default=0.0))
    )
    out["score_formula"] = bool_score(out, "formula_anchor_match")
    out["score_support"] = clip01(numeric(out, "generation_support_fraction", default=0.0))
    out["score_operator"] = operator_score(out)
    nearest_db = numeric(out, "nearest_db_tanimoto", default=0.0)
    novelty = numeric(out, "novelty_score", default=np.nan)
    out["score_novelty"] = novelty.where(novelty.notna(), 1.0 - nearest_db).clip(lower=0.0, upper=1.0)
    out["score_formula_guard"] = np.where(
        (out["score_formula"] >= 1.0) | (out["score_anchor"] >= 0.55),
        1.0,
        (0.5 * out["score_anchor"] + 0.5 * out["score_formula"]).clip(lower=0.0, upper=1.0),
    )
    out["score_novelty_guard"] = np.where(
        (nearest_db >= 0.985) | (out["score_anchor"] < 0.20),
        0.25 * out["score_novelty"],
        (0.65 * out["score_novelty"] + 0.35 * out["score_anchor"]).clip(lower=0.0, upper=1.0),
    )
    old_rank = numeric(out, "rank", default=1.0).clip(lower=1.0)
    out["score_rank_rescue"] = np.where(old_rank > 10, 1.0, 0.25)
    return out


def apply_rerank_policy(pool: pd.DataFrame, policy: AssemblyRerankPolicy, tolerance_da: float) -> pd.DataFrame:
    """Apply policy weights and compute a composite rerank score."""
    out = add_score_features(pool, tolerance_da)
    component_map = {
        "assembly": "score_assembly",
        "mass": "score_mass",
        "adduct": "score_adduct",
        "anchor": "score_anchor",
        "operator": "score_operator",
        "support": "score_support",
        "formula": "score_formula",
        "formula_guard": "score_formula_guard",
        "novelty": "score_novelty",
        "novelty_guard": "score_novelty_guard",
        "rank_rescue": "score_rank_rescue",
    }
    score = pd.Series(0.0, index=out.index, dtype=float)
    total = 0.0
    for key, weight in policy.weights.items():
        score = score + float(weight) * numeric(out, component_map[key], default=0.0)
        total += float(weight)
    out["rerank_policy"] = policy.name
    out["rerank_score"] = score / total if total else score
    return out


def rank_topk(scored: pd.DataFrame, top_k: int) -> pd.DataFrame:
    """Sort by rerank score and assign top-k ranks per feature."""
    if scored.empty:
        return scored.copy()
    out = scored.copy()
    out["feature_key"] = out["feature_key"].astype(str)
    out["predicted_smiles"] = out["predicted_smiles"].astype(str)
    out = (
        out.sort_values(
            [
                "feature_key",
                "rerank_score",
                "score_mass",
                "score_anchor",
                "score_operator",
                "score_formula_guard",
                "assembly_score",
                "abs_mass_error_da",
                "predicted_smiles",
            ],
            ascending=[True, False, False, False, False, False, False, True, True],
            kind="mergesort",
        )
        .drop_duplicates(["feature_key", "predicted_smiles"], keep="first")
        .copy()
    )
    out["rank"] = out.groupby("feature_key").cumcount() + 1
    return out[out["rank"].le(int(top_k))].reset_index(drop=True)


def write_denovo_report(df: pd.DataFrame, path: Path) -> None:
    """Write the standard denovo report CSV schema."""
    path.parent.mkdir(parents=True, exist_ok=True)
    out = df.copy()
    if "feature_label" not in out.columns and "feature_key" in out.columns:
        out["feature_label"] = out["feature_key"]
    if "is_valid" not in out.columns:
        out["is_valid"] = True
    preferred = [
        "feature_key",
        "rank",
        "predicted_smiles",
        "predicted_exact_mass",
        "predicted_formula",
        "abs_mass_error_da",
        "primary_abs_mass_error_da",
        "mass_error_hypothesis",
        "mass_hypothesis_adduct",
        "strict_mass_pass",
        "mass_match_score",
        "assembly_score",
        "rerank_policy",
        "rerank_score",
        "operator_family",
        "assembly_generation_modes",
        "anchor_tanimoto",
        "formula_anchor_match",
        "nearest_db_tanimoto",
        "novelty_score",
    ]
    cols = [col for col in preferred if col in out.columns]
    cols.extend([col for col in out.columns if col not in cols])
    out[cols].to_csv(path, index=False)


def first_value(df: pd.DataFrame, col: str, default: Any = "") -> Any:
    """Return the first non-null value in a column or a default."""
    if df.empty or col not in df.columns:
        return default
    value = df.iloc[0].get(col, default)
    return default if pd.isna(value) else value


def eval_reranked(
    *,
    candidate_df: pd.DataFrame,
    ranked_df: pd.DataFrame,
    controls: list[Any],
    lib_id: int,
    ion_mode: str,
    profile: str,
    policy: str,
    top_k: int,
    candidate_path: Path,
    denovo_path: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Evaluate reranked candidates against known controls."""
    candidate_pool = normalize_prediction_pool(candidate_df, candidate_path)
    denovo_pool = normalize_prediction_pool(ranked_df, denovo_path)
    top_pool = denovo_pool[denovo_pool["_rank"].le(int(top_k))].copy()
    rank1_pool = denovo_pool[denovo_pool["_rank"].le(1)].copy()
    scope_rows: dict[str, list[dict[str, Any]]] = {"rank1": [], "top10": [], "full_pool_oracle": []}
    per_rows: list[dict[str, Any]] = []
    for control in controls:
        full = annotate_pool(pool_for_control(candidate_pool, control), control)
        top10 = annotate_pool(pool_for_control(top_pool, control), control)
        rank1 = annotate_pool(pool_for_control(rank1_pool, control), control)
        buckets = {
            "rank1": best_by_tanimoto(rank1),
            "top10": best_by_tanimoto(top10),
            "full_pool_oracle": best_by_tanimoto(full),
        }
        per_row = {
            "policy": str(policy),
            "lib_id": int(lib_id),
            "ion_mode": str(ion_mode),
            "profile": str(profile),
            "holdout_scope": "x_control",
            "control_id": str(control.control_id),
            "n_candidates": int(full.shape[0]),
            "n_top10": int(top10.shape[0]),
            "candidate_metrics_csv": str(candidate_path),
            "denovo_csv": str(denovo_path),
        }
        for scope, best in buckets.items():
            tan = safe_float(first_value(best, "_best_tanimoto", 0.0), 0.0)
            exact = bool(first_value(best, "_exact_match", False))
            scope_rows[scope].append({"best_tanimoto": tan, "exact_match": exact, "n_candidates": int(full.shape[0])})
            per_row[f"{scope}_best_tanimoto"] = round(tan, 4)
            per_row[f"{scope}_smiles"] = first_value(best, "_predicted_smiles", "")
            per_row[f"{scope}_rank"] = first_value(best, "_rank", "")
        per_rows.append(per_row)
    metric_rows = []
    for scope, rows in scope_rows.items():
        metric_rows.append(
            {
                "policy": str(policy),
                "lib_id": int(lib_id),
                "ion_mode": str(ion_mode),
                "profile": str(profile),
                "holdout_scope": "x_control",
                **metric_from_rows(rows, scope),
            }
        )
    return metric_rows, per_rows


def load_controls_for_profile(
    v3_root: Path, known_x: Path, known_x_all: Path, lib_id: int, ion_mode: str, profile: str
) -> list[Any]:
    """Load X controls scoped to one library and ion mode."""
    per_path = v3_root / "aggregate" / "denovo_assembly_v3_per_control.csv"
    wanted: set[str] = set()
    if per_path.exists():
        per = pd.read_csv(per_path)
        if not per.empty:
            mask = per["lib_id"].astype(int).eq(int(lib_id)) & per["profile"].astype(str).eq(str(profile))
            wanted = set(per.loc[mask, "control_id"].astype(str))
    controls = load_x_controls(known_x=known_x, known_x_all=known_x_all, lib_id=int(lib_id), ion_mode=str(ion_mode))
    return [control for control in controls if not wanted or str(control.control_id) in wanted]


def profile_tolerance(profile: str, profile_dir: Path) -> float:
    """Resolve mass tolerance for a profile directory."""
    summary_path = profile_dir / "selfies_run_summary.json"
    if summary_path.exists():
        try:
            payload = json.loads(summary_path.read_text())
            return safe_float(payload.get("mass_tolerance_da"), 1.0 if profile == "strict_1da" else 5.0)
        except Exception:
            pass
    return 1.0 if str(profile).replace("non_x_holdout_", "") == "strict_1da" else 5.0


def gate_policy(
    *,
    policy: str,
    metric: dict[str, Any],
    baseline: dict[str, Any],
    per_control: pd.DataFrame,
) -> bool:
    """Select rerank policy variant based on profile acceptance gates."""
    if str(metric.get("metric_scope")) != "top10":
        return False
    lib_id = safe_int(metric.get("lib_id"), -1)
    mean_delta = safe_float(metric.get("mean_best_tanimoto"), 0.0) - safe_float(baseline.get("mean_best_tanimoto"), 0.0)
    ge07_delta = safe_int(metric.get("n_tanimoto_ge_0p7"), 0) - safe_int(baseline.get("n_tanimoto_ge_0p7"), 0)
    coverage_delta = safe_int(metric.get("n_controls_covered"), 0) - safe_int(baseline.get("n_controls_covered"), 0)
    if coverage_delta < 0:
        return False
    if lib_id == 305:
        if safe_int(metric.get("n_tanimoto_ge_0p7"), 0) < safe_int(baseline.get("n_tanimoto_ge_0p7"), 0):
            return False
        x25003 = per_control[
            per_control["policy"].astype(str).eq(str(policy)) & per_control["control_id"].astype(str).eq("X-25003")
        ]
        if not x25003.empty and not x25003["top10_best_tanimoto"].ge(GOOD_TANIMOTO).any():
            return False
    if lib_id == 209 and safe_int(metric.get("n_tanimoto_ge_0p7"), 0) < safe_int(baseline.get("n_tanimoto_ge_0p7"), 0):
        return False
    return mean_delta >= 0.03 or ge07_delta >= 1


def run_rerank(
    *,
    v3_root: Path,
    output_root: Path,
    known_x: Path,
    known_x_all: Path,
    libs: list[int],
    profiles: list[str],
    top_k: int,
) -> dict[str, Any]:
    """Rerank assembly candidates and write evaluation artifacts."""
    output_root.mkdir(parents=True, exist_ok=True)
    metric_rows: list[dict[str, Any]] = []
    per_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for lib_id, ion_mode, profile, profile_dir in discover_profile_dirs(v3_root):
        if libs and int(lib_id) not in {int(item) for item in libs}:
            continue
        if profiles and str(profile) not in set(profiles):
            continue
        candidate_path = profile_dir / "candidate_metrics_all.csv"
        if not candidate_path.exists():
            continue
        candidate_df = pd.read_csv(candidate_path)
        controls = load_controls_for_profile(v3_root, known_x, known_x_all, int(lib_id), ion_mode, profile)
        tolerance = profile_tolerance(profile, profile_dir)
        baseline_top10: dict[str, Any] | None = None
        profile_out = (
            output_root
            / "runs"
            / "control_targeted_assembly_v3"
            / f"lib_{int(lib_id)}_{ion_mode}"
            / "selfies_profiles"
            / profile
        )
        for policy in POLICIES.values():
            scored = apply_rerank_policy(candidate_df, policy, tolerance)
            ranked = rank_topk(scored, top_k=top_k)
            policy_dir = profile_out / "rerank_profiles" / policy.name
            denovo_path = policy_dir / "denovo_report.csv"
            write_denovo_report(ranked, denovo_path)
            write_denovo_report(ranked, policy_dir / "ranked_topk_candidates.csv")
            metrics, per_control = eval_reranked(
                candidate_df=candidate_df,
                ranked_df=ranked,
                controls=controls,
                lib_id=int(lib_id),
                ion_mode=ion_mode,
                profile=profile,
                policy=policy.name,
                top_k=top_k,
                candidate_path=candidate_path,
                denovo_path=denovo_path,
            )
            metric_rows.extend(metrics)
            per_rows.extend(per_control)
            top10_metric = next((row for row in metrics if row["metric_scope"] == "top10"), {})
            if policy.name == "assembly_score_baseline":
                baseline_top10 = top10_metric
            summary_rows.append({**top10_metric, "description": policy.description})
        if baseline_top10 is None:
            continue
    metrics_df = pd.DataFrame(metric_rows)
    per_df = pd.DataFrame(per_rows)
    summary_df = pd.DataFrame(summary_rows)
    accepted_rows: list[dict[str, Any]] = []
    if not summary_df.empty:
        baselines = {
            (safe_int(row.get("lib_id")), str(row.get("profile"))): row.to_dict()
            for _, row in summary_df[summary_df["policy"].eq("assembly_score_baseline")].iterrows()
        }
        for _, row in summary_df.iterrows():
            if row.get("policy") == "assembly_score_baseline":
                continue
            baseline = baselines.get((safe_int(row.get("lib_id")), str(row.get("profile"))), {})
            passed = gate_policy(
                policy=str(row.get("policy")), metric=row.to_dict(), baseline=baseline, per_control=per_df
            )
            summary_df.loc[row.name, "acceptance_pass"] = bool(passed)
            summary_df.loc[row.name, "top10_mean_delta_vs_baseline"] = round(
                safe_float(row.get("mean_best_tanimoto"), 0.0) - safe_float(baseline.get("mean_best_tanimoto"), 0.0),
                4,
            )
            summary_df.loc[row.name, "top10_ge07_delta_vs_baseline"] = safe_int(
                row.get("n_tanimoto_ge_0p7"), 0
            ) - safe_int(baseline.get("n_tanimoto_ge_0p7"), 0)
            if passed:
                accepted_rows.append(summary_df.loc[row.name].to_dict())
    write_csv(output_root / "rerank_metrics_by_scope.csv", metrics_df)
    write_csv(output_root / "rerank_per_control_metrics.csv", per_df)
    write_csv(output_root / "rerank_policy_summary.csv", summary_df)
    recommendations = {
        "accepted_rows": accepted_rows,
        "next_branch": (
            "Promote the best accepted assembly rerank policy for de novo top-k review."
            if accepted_rows
            else "No assembly rerank policy cleared gates; inspect missed oracle candidate features before more generation."
        ),
    }
    write_json(output_root / "rerank_recommendations.json", recommendations)
    report = "\n".join(
        [
            "# Assembly V3 Frozen-Pool Rerank",
            "",
            f"Recommendation: {recommendations['next_branch']}",
            "",
            "## Top-10 Policy Summary",
            dataframe_to_markdown(summary_df),
            "",
            "## Accepted Rows",
            dataframe_to_markdown(pd.DataFrame(accepted_rows)),
            "",
        ]
    )
    (output_root / "assembly_v3_rerank_report.md").write_text(report)
    payload = {
        "output_root": str(output_root),
        "metrics_csv": str(output_root / "rerank_metrics_by_scope.csv"),
        "per_control_csv": str(output_root / "rerank_per_control_metrics.csv"),
        "summary_csv": str(output_root / "rerank_policy_summary.csv"),
        "recommendations_json": str(output_root / "rerank_recommendations.json"),
        "report_md": str(output_root / "assembly_v3_rerank_report.md"),
        "accepted_rows": accepted_rows,
    }
    write_json(output_root / "assembly_v3_rerank_summary.json", payload)
    return payload


def parse_csv(text: str) -> list[str]:
    """Parse a comma-separated CLI string into a list."""
    return [part.strip() for part in str(text or "").split(",") if part.strip()]


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for assembly reranking."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v3-root", type=Path, default=DEFAULT_V3_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--known-x", type=Path, default=Path("known_X.csv"))
    parser.add_argument("--known-x-all", type=Path, default=Path("known_X_all_available_structures_simple.csv"))
    parser.add_argument("--libs", type=str, default="305,209")
    parser.add_argument("--profiles", type=str, default="strict_1da,balanced")
    parser.add_argument("--top-k", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    """CLI entry point for assembly reranking."""
    args = parse_args()
    payload = run_rerank(
        v3_root=args.v3_root,
        output_root=args.output_root,
        known_x=args.known_x,
        known_x_all=args.known_x_all,
        libs=parse_libs(args.libs),
        profiles=parse_csv(args.profiles),
        top_k=int(args.top_k),
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
