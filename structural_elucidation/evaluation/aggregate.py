#!/usr/bin/env python3
"""Aggregate X-domain hard holdout and retrieval-anchored analog results."""

from __future__ import annotations

from structural_elucidation.paths import REPO_ROOT

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from structural_elucidation.evaluation.failure_modes import (
    GOOD_TANIMOTO,
    annotate_pool,
    load_non_x_holdouts,
    load_x_controls,
    normalize_prediction_pool,
    pool_for_control,
)
from structural_elucidation.evaluation.candidate_pools import find_artifact, load_mass_tolerance, parse_lib_dir
from structural_elucidation.utils.domain import dataframe_to_markdown, read_table, safe_float, safe_int, write_csv, write_json


def best_value(pool: pd.DataFrame, col: str, default: Any = "") -> Any:
    if pool.empty or col not in pool.columns:
        return default
    value = pool.iloc[0].get(col, default)
    return default if pd.isna(value) else value


def best_by_tanimoto(pool: pd.DataFrame) -> pd.DataFrame:
    if pool.empty:
        return pool
    return pool.sort_values(
        ["_best_tanimoto", "_exact_match", "_rank", "_abs_mass_error_da"],
        ascending=[False, False, True, True],
        na_position="last",
    ).head(1)


def metric_from_rows(rows: list[dict[str, Any]], metric_scope: str) -> dict[str, Any]:
    values = [safe_float(row["best_tanimoto"], 0.0) for row in rows]
    covered = [row for row in rows if safe_int(row.get("n_candidates"), 0) > 0]
    return {
        "metric_scope": metric_scope,
        "n_controls": int(len(rows)),
        "n_controls_covered": int(len(covered)),
        "n_exact_match": int(sum(bool(row.get("exact_match")) for row in rows)),
        "n_tanimoto_ge_0p5": int(sum(value >= 0.5 for value in values)),
        "n_tanimoto_ge_0p7": int(sum(value >= GOOD_TANIMOTO for value in values)),
        "n_tanimoto_ge_0p9": int(sum(value >= 0.9 for value in values)),
        "mean_best_tanimoto": round(float(np.mean(values)) if values else 0.0, 4),
        "median_best_tanimoto": round(float(np.median(values)) if values else 0.0, 4),
    }


def evaluate_profile(
    *,
    profile_dir: Path,
    lib_id: int,
    ion_mode: str,
    profile: str,
    holdout_scope: str,
    controls: list[Any],
    top_k: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    candidate_path = find_artifact(profile_dir, int(lib_id), "candidate_metrics_all.csv")
    denovo_path = find_artifact(profile_dir, int(lib_id), "denovo_report.csv")
    if candidate_path is None or not candidate_path.exists():
        return [], []
    candidate_pool = normalize_prediction_pool(pd.read_csv(candidate_path), candidate_path)
    denovo_pool = (
        normalize_prediction_pool(pd.read_csv(denovo_path), denovo_path)
        if denovo_path is not None and denovo_path.exists()
        else pd.DataFrame()
    )
    top_pool = denovo_pool[denovo_pool["_rank"].le(int(top_k))].copy() if not denovo_pool.empty else denovo_pool
    rank1_pool = denovo_pool[denovo_pool["_rank"].le(1)].copy() if not denovo_pool.empty else denovo_pool
    tolerance = load_mass_tolerance(profile, profile_dir, int(lib_id))
    per_rows: list[dict[str, Any]] = []
    scope_rows: dict[str, list[dict[str, Any]]] = {"rank1": [], "top10": [], "full_pool_oracle": []}
    for control in controls:
        full = annotate_pool(pool_for_control(candidate_pool, control), control)
        top10 = annotate_pool(pool_for_control(top_pool, control), control)
        rank1 = annotate_pool(pool_for_control(rank1_pool, control), control)
        buckets = {
            "rank1": best_by_tanimoto(rank1),
            "top10": best_by_tanimoto(top10),
            "full_pool_oracle": best_by_tanimoto(full),
        }
        mass_series = pd.to_numeric(full.get("_abs_mass_error_da"), errors="coerce")
        per_row = {
            "lib_id": int(lib_id),
            "ion_mode": str(ion_mode),
            "profile": str(profile),
            "base_profile": str(profile).replace("non_x_holdout_", ""),
            "holdout_scope": str(holdout_scope),
            "control_id": str(control.control_id),
            "comp_id": "" if control.comp_id is None else int(control.comp_id),
            "n_candidates": int(full.shape[0]),
            "n_top10": int(top10.shape[0]),
            "mass_tolerance_da": float(tolerance),
            "has_candidate_within_1p0_da": bool(mass_series.le(1.0).any()) if not full.empty else False,
            "has_candidate_within_5p0_da": bool(mass_series.le(5.0).any()) if not full.empty else False,
            "candidate_metrics_csv": str(candidate_path),
            "denovo_csv": str(denovo_path) if denovo_path is not None else "",
        }
        for scope, best in buckets.items():
            value = safe_float(best_value(best, "_best_tanimoto", 0.0), 0.0)
            scope_rows[scope].append(
                {
                    "best_tanimoto": value,
                    "exact_match": bool(best_value(best, "_exact_match", False)),
                    "n_candidates": int(full.shape[0]),
                }
            )
            per_row[f"{scope}_best_tanimoto"] = round(value, 4)
            per_row[f"{scope}_smiles"] = best_value(best, "_predicted_smiles", "")
            per_row[f"{scope}_abs_mass_error_da"] = best_value(best, "_abs_mass_error_da", "")
        per_rows.append(per_row)
    metric_rows = []
    for scope, rows in scope_rows.items():
        metric_rows.append(
            {
                "lib_id": int(lib_id),
                "ion_mode": str(ion_mode),
                "profile": str(profile),
                "base_profile": str(profile).replace("non_x_holdout_", ""),
                "holdout_scope": str(holdout_scope),
                **metric_from_rows(rows, scope),
            }
        )
    return metric_rows, per_rows


def discover_profiles(run_root: Path) -> list[tuple[int, str, str, Path]]:
    rows: list[tuple[int, str, str, Path]] = []
    for lib_dir in sorted(run_root.glob("lib_*_*")):
        if not lib_dir.is_dir():
            continue
        lib_id, ion_mode = parse_lib_dir(lib_dir)
        profile_root = lib_dir / "selfies_profiles"
        if not profile_root.exists():
            continue
        for profile_dir in sorted(profile_root.iterdir()):
            if profile_dir.is_dir():
                rows.append((int(lib_id), str(ion_mode), profile_dir.name, profile_dir))
    return rows


def load_baseline_metrics(ablation_root: Path) -> dict[tuple[str, str, str], dict[str, Any]]:
    path = ablation_root / "baseline_candidate_pool_quality" / "candidate_pool_quality_metrics.csv"
    if not path.exists():
        return {}
    rows = read_table(path).to_dict(orient="records")
    out = {}
    for row in rows:
        key = (str(row.get("lib_id")), str(row.get("source_name")), str(row.get("metric_scope")))
        out[key] = row
    return out


def build_deltas(metrics: pd.DataFrame, baseline: dict[tuple[str, str, str], dict[str, Any]]) -> pd.DataFrame:
    rows = []
    if metrics.empty:
        return pd.DataFrame()
    for _, row in metrics.iterrows():
        if str(row.get("holdout_scope")) != "x_control":
            continue
        key = (str(row.get("lib_id")), str(row.get("base_profile")), str(row.get("metric_scope")))
        base = baseline.get(key)
        if base is None:
            continue
        delta_mean = safe_float(row.get("mean_best_tanimoto"), 0.0) - safe_float(base.get("mean_best_tanimoto"), 0.0)
        delta_ge07 = safe_int(row.get("n_tanimoto_ge_0p7"), 0) - safe_int(base.get("n_tanimoto_ge_0p7"), 0)
        delta_covered = safe_int(row.get("n_controls_covered"), 0) - safe_int(base.get("n_controls_covered"), 0)
        rows.append(
            {
                "lib_id": int(row.get("lib_id")),
                "profile": row.get("base_profile"),
                "metric_scope": row.get("metric_scope"),
                "baseline_mean_best_tanimoto": round(safe_float(base.get("mean_best_tanimoto"), 0.0), 4),
                "analog_mean_best_tanimoto": round(safe_float(row.get("mean_best_tanimoto"), 0.0), 4),
                "delta_mean_best_tanimoto": round(delta_mean, 4),
                "baseline_ge07": safe_int(base.get("n_tanimoto_ge_0p7"), 0),
                "analog_ge07": safe_int(row.get("n_tanimoto_ge_0p7"), 0),
                "delta_ge07": int(delta_ge07),
                "baseline_covered": safe_int(base.get("n_controls_covered"), 0),
                "analog_covered": safe_int(row.get("n_controls_covered"), 0),
                "delta_covered": int(delta_covered),
                "acceptance_pass": bool((delta_mean >= 0.03 or delta_ge07 >= 1) and delta_covered >= 0),
            }
        )
    return pd.DataFrame(rows)


def build_hard_x_gap(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return pd.DataFrame()
    rows = []
    for (lib_id, profile, scope), group in metrics.groupby(["lib_id", "base_profile", "metric_scope"], dropna=False):
        x = group[group["holdout_scope"].eq("x_control")]
        hard = group[group["holdout_scope"].eq("non_x_known")]
        if x.empty or hard.empty:
            continue
        x_row = x.iloc[0]
        hard_row = hard.iloc[0]
        rows.append(
            {
                "lib_id": int(lib_id),
                "profile": profile,
                "metric_scope": scope,
                "x_mean_best_tanimoto": safe_float(x_row.get("mean_best_tanimoto"), 0.0),
                "hard_mean_best_tanimoto": safe_float(hard_row.get("mean_best_tanimoto"), 0.0),
                "hard_minus_x_mean": round(safe_float(hard_row.get("mean_best_tanimoto"), 0.0) - safe_float(x_row.get("mean_best_tanimoto"), 0.0), 4),
                "x_ge07": safe_int(x_row.get("n_tanimoto_ge_0p7"), 0),
                "hard_ge07": safe_int(hard_row.get("n_tanimoto_ge_0p7"), 0),
                "hard_minus_x_ge07": safe_int(hard_row.get("n_tanimoto_ge_0p7"), 0) - safe_int(x_row.get("n_tanimoto_ge_0p7"), 0),
            }
        )
    return pd.DataFrame(rows)


def render_report(metrics: pd.DataFrame, deltas: pd.DataFrame, gap: pd.DataFrame, recommendations: dict[str, Any]) -> str:
    display_metrics = metrics[
        [
            col
            for col in [
                "lib_id",
                "profile",
                "holdout_scope",
                "metric_scope",
                "n_controls",
                "n_controls_covered",
                "mean_best_tanimoto",
                "n_tanimoto_ge_0p7",
            ]
            if col in metrics.columns
        ]
    ]
    lines = [
        "# X-Domain Analog V1 Aggregate",
        "",
        f"Recommendation: {recommendations.get('next_branch', '')}",
        "",
        "## Analog Metrics",
        dataframe_to_markdown(display_metrics),
        "",
        "## X Deltas Versus Baseline SELFIES",
        dataframe_to_markdown(deltas),
        "",
        "## Hard Holdout Versus X Gap",
        dataframe_to_markdown(gap),
        "",
    ]
    return "\n".join(lines)


def aggregate(
    *,
    output_root: Path,
    run_root: Path,
    ablation_root: Path,
    hard_holdout_dir: Path,
    known_x: Path,
    known_x_all: Path,
    top_k: int,
) -> dict[str, Any]:
    out_dir = output_root / "aggregate"
    out_dir.mkdir(parents=True, exist_ok=True)
    metric_rows: list[dict[str, Any]] = []
    per_rows: list[dict[str, Any]] = []
    for lib_id, ion_mode, profile, profile_dir in discover_profiles(run_root):
        if str(profile).startswith("non_x_holdout"):
            controls = load_non_x_holdouts(hard_holdout_dir, int(lib_id), str(ion_mode))
            holdout_scope = "non_x_known"
        else:
            controls = load_x_controls(known_x=known_x, known_x_all=known_x_all, lib_id=int(lib_id), ion_mode=str(ion_mode))
            holdout_scope = "x_control"
        metrics, per_control = evaluate_profile(
            profile_dir=profile_dir,
            lib_id=int(lib_id),
            ion_mode=str(ion_mode),
            profile=str(profile),
            holdout_scope=holdout_scope,
            controls=controls,
            top_k=top_k,
        )
        metric_rows.extend(metrics)
        per_rows.extend(per_control)
    metrics_df = pd.DataFrame(metric_rows).sort_values(["lib_id", "profile", "holdout_scope", "metric_scope"]) if metric_rows else pd.DataFrame()
    per_df = pd.DataFrame(per_rows)
    deltas_df = build_deltas(metrics_df, load_baseline_metrics(ablation_root))
    gap_df = build_hard_x_gap(metrics_df)
    accepted = deltas_df[deltas_df["acceptance_pass"].astype(bool)].to_dict(orient="records") if not deltas_df.empty else []
    oracle_accepted = [
        row for row in accepted if str(row.get("metric_scope")) == "full_pool_oracle"
    ]
    topk_accepted = [
        row for row in accepted if str(row.get("metric_scope")) == "top10"
    ]
    recommendations = {
        "accepted_rows": accepted,
        "next_branch": (
            "Analog full-pool oracle improved; evaluate reranking if top-10 did not also improve."
            if oracle_accepted and not topk_accepted
            else "Analog top-10 cleared a gate; inspect retrieval-preserving delivery policies."
            if topk_accepted
            else "Analog did not clear X gates; keep retrieval authoritative and revisit X feature/domain representation."
        ),
    }
    write_csv(out_dir / "x_domain_analog_v1_metrics.csv", metrics_df)
    write_csv(out_dir / "x_domain_analog_v1_per_control.csv", per_df)
    write_csv(out_dir / "x_domain_analog_v1_deltas.csv", deltas_df)
    write_csv(out_dir / "hard_holdout_vs_x_gap.csv", gap_df)
    write_json(out_dir / "next_branch_recommendations.json", recommendations)
    (out_dir / "x_domain_analog_v1_report.md").write_text(render_report(metrics_df, deltas_df, gap_df, recommendations))
    payload = {
        "output_root": str(output_root),
        "run_root": str(run_root),
        "metric_rows": int(metrics_df.shape[0]),
        "per_control_rows": int(per_df.shape[0]),
        "accepted_rows": int(len(accepted)),
        "outputs": {
            "report_md": str(out_dir / "x_domain_analog_v1_report.md"),
            "deltas_csv": str(out_dir / "x_domain_analog_v1_deltas.csv"),
            "hard_holdout_vs_x_gap_csv": str(out_dir / "hard_holdout_vs_x_gap.csv"),
            "recommendations_json": str(out_dir / "next_branch_recommendations.json"),
            "metrics_csv": str(out_dir / "x_domain_analog_v1_metrics.csv"),
            "per_control_csv": str(out_dir / "x_domain_analog_v1_per_control.csv"),
        },
    }
    write_json(out_dir / "x_domain_analog_v1_aggregate_summary.json", payload)
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=Path("models/denovo_x_domain_analog_v1"))
    parser.add_argument("--run-root", type=Path, default=Path("models/denovo_x_domain_analog_v1/runs/retrieval_anchored_analog_v1"))
    parser.add_argument("--ablation-root", type=Path, default=Path("models/decoder_ablation_experiments"))
    parser.add_argument("--hard-holdout-dir", type=Path, default=Path("models/denovo_x_domain_analog_v1/hard_holdouts"))
    parser.add_argument("--known-x", type=Path, default=Path("known_X.csv"))
    parser.add_argument("--known-x-all", type=Path, default=Path("known_X_all_available_structures_simple.csv"))
    parser.add_argument("--top-k", type=int, default=10)
    return parser.parse_args()


def resolve(path: Path, repo_root: Path) -> Path:
    path = path.expanduser()
    return path if path.is_absolute() else repo_root / path


def main() -> None:
    args = parse_args()
    repo_root = REPO_ROOT
    payload = aggregate(
        output_root=resolve(args.output_root, repo_root),
        run_root=resolve(args.run_root, repo_root),
        ablation_root=resolve(args.ablation_root, repo_root),
        hard_holdout_dir=resolve(args.hard_holdout_dir, repo_root),
        known_x=resolve(args.known_x, repo_root),
        known_x_all=resolve(args.known_x_all, repo_root),
        top_k=int(args.top_k),
    )
    print(json.dumps(payload["outputs"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
