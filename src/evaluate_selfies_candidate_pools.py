#!/usr/bin/env python3
"""Evaluate SELFIES candidate pools without requiring a full bundle manifest.

This is intended for decoder ablations: train a candidate decoder, generate one
or more profile directories, then compare rank-1, top-k, and full-pool oracle
quality against known controls.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from build_post_run_decision_package import (
    METRIC_COLUMNS,
    LibContext,
    dataframe_to_markdown,
    evaluate_full_pool_oracle,
    evaluate_prediction_path,
    mass_coverage_row,
)
from evaluate_denovo_controls import HAS_RDKIT, load_controls


PROFILE_TOLERANCE_DEFAULTS = {
    "strict_1da": 1.0,
    "balanced": 5.0,
    "coverage_explore": 5.0,
}


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def find_artifact(profile_dir: Path, lib_id: int, suffix: str) -> Path | None:
    direct = first_existing(
        [
            profile_dir / f"merged_{lib_id}_{suffix}",
            profile_dir / f"selfies_{suffix}",
            profile_dir / suffix,
        ]
    )
    if direct is not None:
        return direct
    matches = sorted(profile_dir.glob(f"*{suffix}"))
    return matches[0] if matches else None


def parse_lib_dir(lib_dir: Path) -> tuple[int, str]:
    parts = lib_dir.name.split("_", 2)
    if len(parts) != 3 or parts[0] != "lib":
        raise ValueError(f"Expected lib directory named lib_<id>_<mode>, got {lib_dir}")
    return int(parts[1]), str(parts[2])


def discover_profile_dirs(run_root: Path) -> list[tuple[str, Path, int, str, str, Path]]:
    """Return experiment, lib_dir, lib_id, ion_mode, profile, profile_dir rows."""
    rows: list[tuple[str, Path, int, str, str, Path]] = []
    for lib_dir in sorted(run_root.rglob("lib_*_*")):
        if not lib_dir.is_dir():
            continue
        lib_id, ion_mode = parse_lib_dir(lib_dir)
        profile_root = lib_dir / "selfies_profiles"
        if not profile_root.exists():
            continue
        rel_parts = lib_dir.relative_to(run_root).parts
        experiment = rel_parts[0] if len(rel_parts) > 1 else run_root.name
        for profile_dir in sorted(profile_root.iterdir()):
            if not profile_dir.is_dir():
                continue
            rows.append((experiment, lib_dir, lib_id, ion_mode, profile_dir.name, profile_dir))
    return rows


def load_mass_tolerance(profile: str, profile_dir: Path, lib_id: int) -> float:
    summary_path = find_artifact(profile_dir, lib_id, "selfies_run_summary.json")
    payload = read_json(summary_path) if summary_path is not None else {}
    raw = payload.get("mass_tolerance_da")
    try:
        if raw is not None:
            return float(raw)
    except Exception:
        pass
    return float(PROFILE_TOLERANCE_DEFAULTS.get(str(profile), 5.0))


def evaluate_candidate_pools(
    *,
    run_root: Path,
    known_x: Path,
    known_x_all: Path,
    output_dir: Path,
    top_k: int = 10,
) -> dict[str, Any]:
    if not HAS_RDKIT:
        raise RuntimeError("RDKit is required for structural candidate-pool evaluation.")

    output_dir.mkdir(parents=True, exist_ok=True)
    metric_rows: list[dict[str, Any]] = []
    per_control_rows: list[dict[str, Any]] = []
    oracle_rows: list[dict[str, Any]] = []
    mass_rows: list[dict[str, Any]] = []

    discovered = discover_profile_dirs(run_root)
    for experiment, lib_dir, lib_id, ion_mode, profile, profile_dir in discovered:
        ctx = LibContext(
            lib_dir=lib_dir,
            lib_id=int(lib_id),
            ion_mode=str(ion_mode),
            manifest={"lib_id": int(lib_id), "ion_mode": str(ion_mode)},
        )
        controls, _scope_summary = load_controls(
            known_x,
            known_x_all,
            control_lib_id=int(lib_id),
        )
        denovo_path = find_artifact(profile_dir, int(lib_id), "denovo_report.csv")
        candidate_path = find_artifact(profile_dir, int(lib_id), "candidate_metrics_all.csv")
        mass_tolerance_da = load_mass_tolerance(profile, profile_dir, int(lib_id))

        if denovo_path is not None:
            rows, per_rows = evaluate_prediction_path(
                ctx=ctx,
                controls=controls,
                pred_path=denovo_path,
                source_type="selfies_profile",
                source_name=str(profile),
                top_k=int(top_k),
            )
            for row in rows:
                row["experiment"] = str(experiment)
            for row in per_rows:
                row["experiment"] = str(experiment)
            metric_rows.extend(rows)
            per_control_rows.extend(per_rows)

        if candidate_path is not None:
            candidate_df = pd.read_csv(candidate_path)
            mass_row = mass_coverage_row(ctx, str(profile), profile_dir, candidate_df)
            mass_row["experiment"] = str(experiment)
            mass_rows.append(mass_row)

            summary, rows = evaluate_full_pool_oracle(
                ctx=ctx,
                controls=controls,
                profile=str(profile),
                profile_dir=profile_dir,
                candidate_df=candidate_df,
                denovo_path=denovo_path,
                mass_tolerance_da=float(mass_tolerance_da),
            )
            metric_row = {
                **{
                    col: None
                    for col in METRIC_COLUMNS
                },
                "lib_id": int(lib_id),
                "ion_mode": str(ion_mode),
                "source_type": "selfies_profile",
                "source_name": str(profile),
                "lane": "",
                "policy": "",
                "metric_scope": "full_pool_oracle",
                "n_controls": int(summary.get("n_controls", 0) or 0),
                "n_controls_covered": int(summary.get("n_controls_covered", 0) or 0),
                "n_exact_match": int(summary.get("n_exact_match", 0) or 0),
                "exact_match_rate": float(summary.get("exact_match_rate", 0.0) or 0.0),
                "n_tanimoto_ge_0p5": int(summary.get("n_oracle_tanimoto_ge_0.5", 0) or 0),
                "n_tanimoto_ge_0p7": int(summary.get("n_oracle_tanimoto_ge_0.7", 0) or 0),
                "n_tanimoto_ge_0p9": int(summary.get("n_oracle_tanimoto_ge_0.9", 0) or 0),
                "tanimoto_ge_0p7_rate": float(summary.get("oracle_tanimoto_ge_0.7_rate", 0.0) or 0.0),
                "mean_best_tanimoto": float(summary.get("mean_best_tanimoto", 0.0) or 0.0),
                "median_best_tanimoto": float(summary.get("median_best_tanimoto", 0.0) or 0.0),
                "path": str(candidate_path),
                "experiment": str(experiment),
            }
            metric_rows.append(metric_row)
            for row in rows:
                row["experiment"] = str(experiment)
            oracle_rows.extend(rows)

    metrics_df = pd.DataFrame(metric_rows)
    if not metrics_df.empty:
        ordered = ["experiment", *METRIC_COLUMNS]
        metrics_df = metrics_df.reindex(
            columns=[*ordered, *[col for col in metrics_df.columns if col not in ordered]]
        )
    per_control_df = pd.DataFrame(per_control_rows)
    oracle_df = pd.DataFrame(oracle_rows)
    mass_df = pd.DataFrame(mass_rows)

    metrics_csv = output_dir / "candidate_pool_quality_metrics.csv"
    per_control_csv = output_dir / "candidate_pool_selected_per_control.csv"
    oracle_csv = output_dir / "candidate_pool_oracle_per_control.csv"
    mass_csv = output_dir / "candidate_pool_mass_coverage.csv"
    metrics_df.to_csv(metrics_csv, index=False)
    per_control_df.to_csv(per_control_csv, index=False)
    oracle_df.to_csv(oracle_csv, index=False)
    mass_df.to_csv(mass_csv, index=False)

    report = render_report(
        run_root=run_root,
        metrics_df=metrics_df,
        mass_df=mass_df,
        discovered_count=len(discovered),
        output_paths={
            "metrics": metrics_csv,
            "per_control": per_control_csv,
            "oracle": oracle_csv,
            "mass": mass_csv,
        },
    )
    report_md = output_dir / "candidate_pool_quality_report.md"
    report_md.write_text(report)

    payload = {
        "run_root": str(run_root),
        "known_x": str(known_x),
        "known_x_all": str(known_x_all),
        "top_k": int(top_k),
        "discovered_profiles": int(len(discovered)),
        "outputs": {
            "metrics_csv": str(metrics_csv),
            "per_control_csv": str(per_control_csv),
            "oracle_csv": str(oracle_csv),
            "mass_csv": str(mass_csv),
            "report_md": str(report_md),
        },
    }
    (output_dir / "candidate_pool_quality_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    return payload


def render_report(
    *,
    run_root: Path,
    metrics_df: pd.DataFrame,
    mass_df: pd.DataFrame,
    discovered_count: int,
    output_paths: dict[str, Path],
) -> str:
    lines = [
        "# SELFIES Candidate Pool Quality",
        "",
        f"Run root: `{run_root}`",
        f"Discovered profile directories: {int(discovered_count)}",
        "",
    ]
    if metrics_df.empty:
        lines.append("No metrics were produced.")
    else:
        display = metrics_df[
            metrics_df["metric_scope"].isin(["rank1", "top10", "full_pool_oracle"])
        ].copy()
        keep = [
            "experiment",
            "lib_id",
            "source_name",
            "metric_scope",
            "n_controls_covered",
            "n_exact_match",
            "n_tanimoto_ge_0p7",
            "mean_best_tanimoto",
        ]
        for col in ["mean_best_tanimoto"]:
            display[col] = pd.to_numeric(display[col], errors="coerce").round(4)
        lines.extend(
            [
                "## Rank And Oracle Metrics",
                "",
                dataframe_to_markdown(display.reindex(columns=keep).sort_values(keep[:4])),
                "",
            ]
        )

    if not mass_df.empty:
        mass_display = mass_df[
            [
                "experiment",
                "lib_id",
                "profile",
                "candidate_rows",
                "candidate_features",
                "within_0p5_da_rate",
                "within_1p0_da_rate",
                "within_2p0_da_rate",
                "within_5p0_da_rate",
            ]
        ].copy()
        for col in [c for c in mass_display.columns if c.endswith("_rate")]:
            mass_display[col] = pd.to_numeric(mass_display[col], errors="coerce").round(4)
        lines.extend(
            [
                "## Mass Coverage",
                "",
                dataframe_to_markdown(mass_display.sort_values(["experiment", "lib_id", "profile"])),
                "",
            ]
        )

    lines.extend(
        [
            "## Outputs",
            "",
            *[f"- `{name}`: `{path}`" for name, path in output_paths.items()],
            "",
        ]
    )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate generated SELFIES profile dirs by rank-1, top-k, and full-pool oracle."
    )
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--known-x", type=Path, default=Path("known_X.csv"))
    parser.add_argument(
        "--known-x-all",
        type=Path,
        default=Path("known_X_all_available_structures_simple.csv"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-k", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    evaluate_candidate_pools(
        run_root=args.run_root.expanduser().resolve(),
        known_x=args.known_x.expanduser().resolve(),
        known_x_all=args.known_x_all.expanduser().resolve(),
        output_dir=args.output_dir.expanduser().resolve(),
        top_k=int(args.top_k),
    )


if __name__ == "__main__":
    main()
