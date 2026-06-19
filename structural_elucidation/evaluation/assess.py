#!/usr/bin/env python3
"""Assess V3 assembly rank gaps before frozen-pool reranking."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd

from structural_elucidation.evaluation.aggregate import best_by_tanimoto
from structural_elucidation.evaluation.failure_modes import (
    annotate_pool,
    load_x_controls,
    normalize_prediction_pool,
    pool_for_control,
)
from structural_elucidation.utils.domain import dataframe_to_markdown, safe_float, safe_int, write_csv, write_json

DEFAULT_V3_ROOT = Path("models/denovo_pipeline_v3/assembly_v3")
GOOD_TANIMOTO = 0.7


def first_value(df: pd.DataFrame, col: str, default: Any = "") -> Any:
    """Return the first non-null value in a column or a default."""
    if df.empty or col not in df.columns:
        return default
    value = df.iloc[0].get(col, default)
    return default if pd.isna(value) else value


def parse_lib_dir(lib_dir: Path) -> tuple[int, str]:
    """Parse ``lib_<id>_<ion_mode>`` directory names into components."""
    parts = lib_dir.name.split("_")
    if len(parts) < 3 or parts[0] != "lib":
        raise ValueError(f"Cannot parse lib directory: {lib_dir}")
    return int(parts[1]), "_".join(parts[2:])


def discover_profile_dirs(v3_root: Path) -> list[tuple[int, str, str, Path]]:
    """List V3 assembly profile directories under the run root."""
    run_root = v3_root / "runs" / "control_targeted_assembly_v3"
    rows: list[tuple[int, str, str, Path]] = []
    for lib_dir in sorted(run_root.glob("lib_*_*")):
        if not lib_dir.is_dir():
            continue
        lib_id, ion_mode = parse_lib_dir(lib_dir)
        for profile_dir in sorted((lib_dir / "selfies_profiles").glob("*")):
            if profile_dir.is_dir():
                rows.append((lib_id, ion_mode, profile_dir.name, profile_dir))
    return rows


def load_evaluated_control_ids(v3_root: Path, lib_id: int, profile: str) -> set[str]:
    """Load control ids already present in the V3 aggregate audit."""
    path = v3_root / "aggregate" / "denovo_assembly_v3_per_control.csv"
    if not path.exists():
        return set()
    df = pd.read_csv(path)
    if df.empty:
        return set()
    mask = df["lib_id"].astype(int).eq(int(lib_id)) & df["profile"].astype(str).eq(str(profile))
    return set(df.loc[mask, "control_id"].astype(str))


def top10_miss_reason(oracle_rank: int, oracle_tanimoto: float, top10_tanimoto: float) -> str:
    """Classify why a high-oracle candidate may be missing from top-k."""
    if oracle_tanimoto < GOOD_TANIMOTO:
        return "no_high_oracle_candidate"
    if oracle_rank <= 10:
        return "high_oracle_already_top10"
    if top10_tanimoto >= GOOD_TANIMOTO:
        return "top10_has_different_high_hit"
    return "high_oracle_misranked_outside_top10"


def assess_profile(
    *,
    v3_root: Path,
    lib_id: int,
    ion_mode: str,
    profile: str,
    profile_dir: Path,
    known_x: Path,
    known_x_all: Path,
    top_k: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Compare full-pool oracle vs top-k ranks for one V3 profile."""
    candidate_path = profile_dir / "candidate_metrics_all.csv"
    denovo_path = profile_dir / "denovo_report.csv"
    if not candidate_path.exists():
        return [], []
    candidate_pool = normalize_prediction_pool(pd.read_csv(candidate_path), candidate_path)
    denovo_pool = (
        normalize_prediction_pool(pd.read_csv(denovo_path), denovo_path) if denovo_path.exists() else pd.DataFrame()
    )
    top_pool = denovo_pool[denovo_pool["_rank"].le(int(top_k))].copy() if not denovo_pool.empty else denovo_pool
    wanted_ids = load_evaluated_control_ids(v3_root, lib_id, profile)
    controls = [
        control
        for control in load_x_controls(known_x=known_x, known_x_all=known_x_all, lib_id=lib_id, ion_mode=ion_mode)
        if not wanted_ids or str(control.control_id) in wanted_ids
    ]
    per_control: list[dict[str, Any]] = []
    operator_rows: list[dict[str, Any]] = []
    for control in controls:
        full = annotate_pool(pool_for_control(candidate_pool, control), control)
        top10 = annotate_pool(pool_for_control(top_pool, control), control)
        best_full = best_by_tanimoto(full)
        best_top = best_by_tanimoto(top10)
        oracle_rank = safe_int(first_value(best_full, "_rank", 999999), 999999)
        oracle_tan = safe_float(first_value(best_full, "_best_tanimoto", 0.0), 0.0)
        top10_tan = safe_float(first_value(best_top, "_best_tanimoto", 0.0), 0.0)
        per_control.append(
            {
                "lib_id": int(lib_id),
                "ion_mode": str(ion_mode),
                "profile": str(profile),
                "control_id": str(control.control_id),
                "n_candidates": int(full.shape[0]),
                "n_top10": int(top10.shape[0]),
                "oracle_tanimoto": round(oracle_tan, 4),
                "oracle_current_rank": oracle_rank if oracle_rank < 999999 else "",
                "top10_best_tanimoto": round(top10_tan, 4),
                "rank_gap": max(0, oracle_rank - int(top_k)) if oracle_rank < 999999 else "",
                "top10_miss_reason": top10_miss_reason(oracle_rank, oracle_tan, top10_tan),
                "oracle_operator_family": first_value(best_full, "operator_family", ""),
                "oracle_assembly_score": first_value(best_full, "assembly_score", ""),
                "oracle_anchor_tanimoto": first_value(best_full, "anchor_tanimoto", ""),
                "oracle_abs_mass_error_da": first_value(best_full, "_abs_mass_error_da", ""),
                "oracle_formula_anchor_match": first_value(best_full, "formula_anchor_match", ""),
                "oracle_nearest_db_tanimoto": first_value(best_full, "nearest_db_tanimoto", ""),
                "oracle_smiles": first_value(best_full, "_predicted_smiles", ""),
                "top10_smiles": first_value(best_top, "_predicted_smiles", ""),
            }
        )
        if not full.empty and "operator_family" in full.columns:
            for operator, group in full.groupby(full["operator_family"].fillna("").astype(str), dropna=False):
                best = best_by_tanimoto(group)
                operator_rows.append(
                    {
                        "lib_id": int(lib_id),
                        "profile": str(profile),
                        "control_id": str(control.control_id),
                        "operator_family": str(operator),
                        "n_candidates": int(group.shape[0]),
                        "best_tanimoto": round(safe_float(first_value(best, "_best_tanimoto", 0.0), 0.0), 4),
                        "best_rank": first_value(best, "_rank", ""),
                        "best_abs_mass_error_da": first_value(best, "_abs_mass_error_da", ""),
                    }
                )
    return per_control, operator_rows


def assess_v3_outputs(
    *,
    v3_root: Path,
    known_x: Path,
    known_x_all: Path,
    out_dir: Path,
    top_k: int,
) -> dict[str, Any]:
    """Run rank-gap assessment across all discovered V3 profiles."""
    out_dir.mkdir(parents=True, exist_ok=True)
    per_rows: list[dict[str, Any]] = []
    operator_rows: list[dict[str, Any]] = []
    for lib_id, ion_mode, profile, profile_dir in discover_profile_dirs(v3_root):
        rows, ops = assess_profile(
            v3_root=v3_root,
            lib_id=lib_id,
            ion_mode=ion_mode,
            profile=profile,
            profile_dir=profile_dir,
            known_x=known_x,
            known_x_all=known_x_all,
            top_k=top_k,
        )
        per_rows.extend(rows)
        operator_rows.extend(ops)
    per_df = pd.DataFrame(per_rows)
    op_df = pd.DataFrame(operator_rows)
    skipped_path = v3_root / "skipped_controls.csv"
    skipped_df = pd.read_csv(skipped_path) if skipped_path.exists() else pd.DataFrame()
    summary_rows: list[dict[str, Any]] = []
    if not per_df.empty:
        for (lib_id, profile), group in per_df.groupby(["lib_id", "profile"], dropna=False):
            summary_rows.append(
                {
                    "lib_id": int(lib_id),
                    "profile": str(profile),
                    "n_controls": int(group.shape[0]),
                    "n_misranked_high_oracle": int(
                        group["top10_miss_reason"].eq("high_oracle_misranked_outside_top10").sum()
                    ),
                    "n_oracle_ge_0p7": int(group["oracle_tanimoto"].ge(GOOD_TANIMOTO).sum()),
                    "n_top10_ge_0p7": int(group["top10_best_tanimoto"].ge(GOOD_TANIMOTO).sum()),
                    "mean_oracle_tanimoto": round(float(group["oracle_tanimoto"].mean()), 4),
                    "mean_top10_tanimoto": round(float(group["top10_best_tanimoto"].mean()), 4),
                }
            )
    summary_df = pd.DataFrame(summary_rows)
    write_csv(out_dir / "v3_rank_gap_per_control.csv", per_df)
    write_csv(out_dir / "v3_operator_family_summary.csv", op_df)
    write_csv(out_dir / "v3_rank_gap_summary.csv", summary_df)
    write_csv(out_dir / "v3_skipped_controls_audit.csv", skipped_df)
    report = "\n".join(
        [
            "# V3 Assembly Output Assessment",
            "",
            "## Rank Gap Summary",
            dataframe_to_markdown(summary_df),
            "",
            "## High-Oracle Top-10 Misses",
            dataframe_to_markdown(
                per_df[per_df["top10_miss_reason"].eq("high_oracle_misranked_outside_top10")]
                if not per_df.empty
                else per_df
            ),
            "",
            "## Skipped Controls",
            dataframe_to_markdown(skipped_df),
            "",
        ]
    )
    (out_dir / "v3_output_assessment_report.md").write_text(report)
    payload = {
        "out_dir": str(out_dir),
        "n_per_control_rows": int(per_df.shape[0]),
        "n_operator_rows": int(op_df.shape[0]),
        "outputs": {
            "per_control": str(out_dir / "v3_rank_gap_per_control.csv"),
            "operator_summary": str(out_dir / "v3_operator_family_summary.csv"),
            "summary": str(out_dir / "v3_rank_gap_summary.csv"),
            "skipped": str(out_dir / "v3_skipped_controls_audit.csv"),
            "report": str(out_dir / "v3_output_assessment_report.md"),
        },
    }
    write_json(out_dir / "v3_output_assessment_summary.json", payload)
    return payload


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for V3 output assessment."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v3-root", type=Path, default=DEFAULT_V3_ROOT)
    parser.add_argument("--known-x", type=Path, default=Path("known_X.csv"))
    parser.add_argument("--known-x-all", type=Path, default=Path("known_X_all_available_structures_simple.csv"))
    parser.add_argument("--out-dir", type=Path, default=Path("models/denovo_pipeline_v3/assembly_v3_assessment"))
    parser.add_argument("--top-k", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    """CLI entry point for V3 assembly rank-gap assessment."""
    args = parse_args()
    payload = assess_v3_outputs(
        v3_root=args.v3_root,
        known_x=args.known_x,
        known_x_all=args.known_x_all,
        out_dir=args.out_dir,
        top_k=int(args.top_k),
    )
    print(payload)


if __name__ == "__main__":
    main()
