#!/usr/bin/env python3
"""Run deconjugation-first biotransformation-guided 209 de novo generation."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from aggregate_x_domain_analog_v1 import best_by_tanimoto, metric_from_rows
from biotransformation_guided_v1_utils import (
    TRACK_ANNOTATION_ASSISTED,
    TRACK_CONTROL_BLIND,
    candidate_mass_score,
    profile_tolerance,
    reattach_transform,
)
from denovo_improvement_utils import canonicalize_smiles, exact_mass_from_smiles, formula_from_smiles, inchikey_from_smiles
from evaluate_denovo_failure_modes import annotate_pool, load_x_controls, normalize_prediction_pool, pool_for_control
from nce_evidence_utils import forward_spectrum_score, load_excluded_structure_keys
from protected_209_fusion_v1_utils import PROTECTED_209_CONTROLS, GOOD_TANIMOTO
from rerank_denovo_assembly_v3 import eval_reranked, write_denovo_report
from run_retrieval_anchored_denovo_assembly_v4 import resolve
from x_domain_analog_utils import dataframe_to_markdown, safe_float, safe_int, write_csv, write_json


DEFAULT_OUTPUT_ROOT = Path("models/biotransformation_guided_denovo_v1")
RUN_NAME = "biotransformation_guided_denovo_v1"


def parse_csv(text: Any) -> list[str]:
    return [part.strip() for part in str(text or "").split(",") if part.strip()]


def inchikey_connectivity(value: Any) -> str:
    return str(value or "").split("-", 1)[0]


def load_exclusion_keys(repo_root: Path) -> tuple[set[str], set[str], set[str]]:
    paths = [
        repo_root / "known_X_all_available_structures_simple.csv",
        repo_root / "models/x_known_curation_v1/curated_x_structures_broad.csv",
        repo_root / "models/x_known_curation_v1/known_x_evidence_catalog.csv",
        repo_root / "models/denovo_generation_improvement_v1/holdouts/non_x_holdout_209.csv",
        repo_root / "models/denovo_x_domain_analog_v1/hard_holdouts/x_like_hard_holdout_all.csv",
    ]
    excluded = load_excluded_structure_keys(paths)
    inchikeys = {str(value) for value in excluded.get("inchikey", pd.Series(dtype=str)).astype(str) if str(value)}
    return (
        {str(value) for value in excluded.get("canonical_smiles", pd.Series(dtype=str)).astype(str) if str(value)},
        inchikeys,
        {inchikey_connectivity(value) for value in inchikeys if inchikey_connectivity(value)},
    )


def observed_bins(row: pd.Series, column: str) -> list[int]:
    try:
        parsed = json.loads(str(row.get(column, "[]") or "[]"))
    except Exception:
        parsed = []
    return [safe_int(value, 0) for value in parsed if safe_int(value, 0)]


def generate_products(core: pd.Series, limit: int) -> list[tuple[str, str]]:
    transform = str(core.get("transform_name", ""))
    core_smiles = str(core.get("core_smiles", "") or "")
    if transform == "peptide":
        canonical = canonicalize_smiles(core_smiles) or ""
        return [(canonical, "peptide_composition")] if canonical else []
    if transform == "none":
        return []
    return reattach_transform(core_smiles, transform, limit=limit)


def candidate_row(core: pd.Series, smiles: str, operators: str, profile: str, track: str, idx: int) -> dict[str, Any]:
    canonical = canonicalize_smiles(smiles) or ""
    mass = exact_mass_from_smiles(canonical)
    target_mass = safe_float(core.get("neutral_mass"), math.nan)
    tolerance = profile_tolerance(profile)
    abs_error = abs(mass - target_mass) if np.isfinite(mass) and np.isfinite(target_mass) else math.nan
    mass_score = candidate_mass_score(mass, target_mass, tolerance)
    formula = formula_from_smiles(canonical)
    target_formula = str(core.get("target_formula", "") or "")
    formula_match = bool(target_formula and formula == target_formula)
    formula_score = 1.0 if formula_match else 0.35 if not target_formula else 0.0
    forward_score = forward_spectrum_score(
        canonical,
        observed_bins(core, "observed_ion_bins_json"),
        observed_bins(core, "observed_loss_bins_json"),
        safe_int(core.get("lib_id"), 209),
    )
    source_score = safe_float(core.get("core_source_score"), 0.0)
    evidence_confidence = safe_float(core.get("evidence_confidence"), 0.0)
    generation_score = float(
        np.clip(
            0.30 * mass_score
            + 0.20 * source_score
            + 0.18 * evidence_confidence
            + 0.16 * formula_score
            + 0.10 * forward_score
            + 0.06 * safe_float(core.get("adduct_priority"), 0.0),
            0.0,
            1.0,
        )
    )
    return {
        "feature_key": str(core.get("feature_key", "")),
        "control_id": str(core.get("control_id", "")),
        "lib_id": 209,
        "ion_mode": "negative",
        "profile": str(profile),
        "track": str(track),
        "predicted_smiles": canonical,
        "predicted_exact_mass": mass,
        "predicted_formula": formula,
        "target_neutral_mass": target_mass,
        "abs_mass_error_da": abs_error,
        "primary_abs_mass_error_da": abs_error,
        "mass_hypothesis_adduct": str(core.get("mass_hypothesis_adduct", "")),
        "mass_error_hypothesis": str(core.get("mass_hypothesis_adduct", "")),
        "strict_mass_pass": bool(np.isfinite(abs_error) and abs_error <= 1.0),
        "mass_match_score": mass_score,
        "target_formula": target_formula,
        "formula_anchor_match": formula_match,
        "formula_score": formula_score,
        "core_smiles": str(core.get("core_smiles", "")),
        "core_formula": str(core.get("source_formula", "") or core.get("core_formula", "")),
        "core_mass": safe_float(core.get("source_exact_mass"), safe_float(core.get("core_mass"), math.nan)),
        "core_mass_error_da": safe_float(core.get("core_mass_error_da"), math.nan),
        "core_source_score": source_score,
        "core_source_type": str(core.get("core_source_type", "")),
        "source_comp_id": core.get("source_comp_id", ""),
        "transform_name": str(core.get("transform_name", "")),
        "transform_delta_mass": safe_float(core.get("transform_delta_mass"), 0.0),
        "peptide_composition": str(core.get("peptide_composition", "")),
        "source_track": str(core.get("evidence_source", "")),
        "evidence_confidence": evidence_confidence,
        "forward_spectrum_score": forward_score,
        "operator_family": f"biotransformation_{core.get('transform_name', '')}",
        "assembly_generation_modes": str(operators),
        "assembly_score": generation_score,
        "generation_score": generation_score,
        "novelty_score": 1.0,
        "_source_order": int(idx),
    }


def rank_candidates(candidates: pd.DataFrame, top_k: int) -> pd.DataFrame:
    if candidates.empty:
        return candidates.copy()
    out = candidates.copy()
    out = (
        out.sort_values(
            [
                "feature_key",
                "generation_score",
                "mass_match_score",
                "formula_anchor_match",
                "core_source_score",
                "forward_spectrum_score",
                "abs_mass_error_da",
                "predicted_smiles",
            ],
            ascending=[True, False, False, False, False, False, True, True],
            kind="mergesort",
        )
        .drop_duplicates(["feature_key", "predicted_smiles"], keep="first")
        .copy()
    )
    out["rank"] = out.groupby("feature_key").cumcount() + 1
    return out[out["rank"].le(int(top_k))].reset_index(drop=True)


def evaluate_outputs(
    *,
    candidate_df: pd.DataFrame,
    ranked_df: pd.DataFrame,
    known_x: Path,
    known_x_all: Path,
    profile: str,
    track: str,
    candidate_path: Path,
    denovo_path: Path,
    top_k: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    ranked_for_eval = ranked_df.copy()
    if not ranked_for_eval.empty and "_rank" not in ranked_for_eval.columns and "rank" in ranked_for_eval.columns:
        ranked_for_eval["_rank"] = ranked_for_eval["rank"]
    controls = load_x_controls(known_x=known_x, known_x_all=known_x_all, lib_id=209, ion_mode="negative")
    try:
        metrics, per_control = eval_reranked(
            candidate_df=candidate_df,
            ranked_df=ranked_for_eval,
            controls=controls,
            lib_id=209,
            ion_mode="negative",
            profile=profile,
            policy=f"{track}_biotransformation_guided",
            top_k=top_k,
            candidate_path=candidate_path,
            denovo_path=denovo_path,
        )
    except KeyError as exc:
        if str(exc).strip("'\"") != "_rank":
            raise
        metrics, per_control = [], []
    for row in metrics:
        row["track"] = str(track)
        row["run_name"] = RUN_NAME
    for row in per_control:
        row["track"] = str(track)
        row["run_name"] = RUN_NAME
    return metrics, per_control


def run_generation(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parent
    output_root = resolve(args.output_root, repo_root)
    core_path = resolve(args.core_candidates_csv, repo_root)
    if not core_path.exists():
        raise FileNotFoundError(f"Missing core candidates CSV: {core_path}")
    cores = pd.read_csv(core_path)
    tracks = set(parse_csv(args.tracks))
    profiles = parse_csv(args.profiles)
    controls = set(parse_csv(args.controls))
    if tracks:
        cores = cores[cores["track"].astype(str).isin(tracks)].copy()
    if controls:
        cores = cores[cores["control_id"].astype(str).isin(controls)].copy()
    if not cores.empty:
        cores = cores.sort_values(["control_id", "track", "evidence_confidence", "core_source_score"], ascending=[True, True, False, False], kind="mergesort")
    excl_canonical, excl_inchikey, excl_connectivity = load_exclusion_keys(repo_root)
    all_metric_rows: list[dict[str, Any]] = []
    all_per_rows: list[dict[str, Any]] = []
    generation_summary_rows: list[dict[str, Any]] = []
    for track in sorted(tracks or {TRACK_CONTROL_BLIND, TRACK_ANNOTATION_ASSISTED}):
        track_cores = cores[cores["track"].astype(str).eq(str(track))].copy()
        for profile in profiles:
            candidate_rows: list[dict[str, Any]] = []
            novelty_rows: list[dict[str, Any]] = []
            seen_per_control: dict[str, int] = {}
            for idx, (_, core) in enumerate(track_cores.iterrows(), start=1):
                control_id = str(core.get("control_id", ""))
                if seen_per_control.get(control_id, 0) >= int(args.max_core_rows_per_control):
                    continue
                products = generate_products(core, limit=int(args.max_products_per_core))
                if not products:
                    continue
                seen_per_control[control_id] = seen_per_control.get(control_id, 0) + 1
                for smiles, operators in products:
                    canonical = canonicalize_smiles(smiles) or ""
                    inchikey = inchikey_from_smiles(canonical) if canonical else ""
                    rejected = (
                        not canonical
                        or canonical in excl_canonical
                        or inchikey in excl_inchikey
                        or inchikey_connectivity(inchikey) in excl_connectivity
                    )
                    novelty_rows.append(
                        {
                            "track": track,
                            "profile": profile,
                            "control_id": control_id,
                            "predicted_smiles": canonical,
                            "inchikey": inchikey,
                            "transform_name": core.get("transform_name", ""),
                            "rejected_exact_excluded_structure": bool(rejected),
                        }
                    )
                    if rejected:
                        continue
                    row = candidate_row(core, canonical, operators, profile, track, idx)
                    if safe_float(row.get("abs_mass_error_da"), math.inf) > profile_tolerance(profile):
                        continue
                    candidate_rows.append(row)
                    if len([r for r in candidate_rows if r["control_id"] == control_id]) >= int(args.max_candidates_per_control):
                        break
            candidate_df = pd.DataFrame(candidate_rows)
            if not candidate_df.empty:
                candidate_df = candidate_df.sort_values(["feature_key", "generation_score", "abs_mass_error_da"], ascending=[True, False, True], kind="mergesort")
                candidate_df["rank"] = candidate_df.groupby("feature_key").cumcount() + 1
            ranked_df = rank_candidates(candidate_df, int(args.top_k))
            profile_dir = (
                output_root
                / "runs"
                / RUN_NAME
                / f"track_{track}"
                / "lib_209_negative"
                / "selfies_profiles"
                / profile
            )
            candidate_path = profile_dir / "candidate_metrics_all.csv"
            denovo_path = profile_dir / "denovo_report.csv"
            novelty_path = profile_dir / "novelty_audit.csv"
            write_csv(candidate_path, candidate_df)
            write_denovo_report(ranked_df, denovo_path)
            write_csv(novelty_path, pd.DataFrame(novelty_rows))
            if ranked_df.empty:
                metrics, per_control = [], []
            else:
                metrics, per_control = evaluate_outputs(
                    candidate_df=candidate_df,
                    ranked_df=ranked_df,
                    known_x=resolve(args.known_x, repo_root),
                    known_x_all=resolve(args.known_x_all, repo_root),
                    profile=profile,
                    track=track,
                    candidate_path=candidate_path,
                    denovo_path=denovo_path,
                    top_k=int(args.top_k),
                )
            all_metric_rows.extend(metrics)
            all_per_rows.extend(per_control)
            generation_summary_rows.append(
                {
                    "track": track,
                    "profile": profile,
                    "n_candidate_rows": int(candidate_df.shape[0]),
                    "n_topk_rows": int(ranked_df.shape[0]),
                    "n_controls_with_candidates": int(candidate_df["control_id"].nunique()) if not candidate_df.empty else 0,
                    "candidate_metrics_csv": str(candidate_path),
                    "denovo_csv": str(denovo_path),
                    "novelty_audit_csv": str(novelty_path),
                }
            )
    metrics_df = pd.DataFrame(all_metric_rows)
    per_df = pd.DataFrame(all_per_rows)
    summary_df = pd.DataFrame(generation_summary_rows)
    write_csv(output_root / "biotransformation_generation_metrics.csv", metrics_df)
    write_csv(output_root / "biotransformation_generation_per_control.csv", per_df)
    write_csv(output_root / "biotransformation_generation_summary.csv", summary_df)
    protected = []
    if not per_df.empty:
        for track in sorted(per_df["track"].dropna().astype(str).unique()):
            scoped = per_df[per_df["track"].astype(str).eq(track)]
            recovered = set(
                scoped[pd.to_numeric(scoped["top10_best_tanimoto"], errors="coerce").fillna(0.0).ge(GOOD_TANIMOTO)]["control_id"].astype(str)
            )
            protected.append(
                {
                    "track": track,
                    "protected_controls_recovered": "|".join(sorted(set(PROTECTED_209_CONTROLS) & recovered)),
                    "n_protected_recovered": int(len(set(PROTECTED_209_CONTROLS) & recovered)),
                    "missing_protected_controls": "|".join(sorted(set(PROTECTED_209_CONTROLS) - recovered)),
                }
            )
    write_csv(output_root / "biotransformation_protected_control_audit.csv", pd.DataFrame(protected))
    report = "\n".join(
        [
            "# Biotransformation-Guided Generation V1",
            "",
            "## Generation Summary",
            dataframe_to_markdown(summary_df),
            "",
            "## Metrics",
            dataframe_to_markdown(metrics_df),
            "",
            "## Protected Control Audit",
            dataframe_to_markdown(pd.DataFrame(protected)),
            "",
        ]
    )
    (output_root / "biotransformation_generation_report.md").write_text(report)
    payload = {
        "output_root": str(output_root),
        "metrics_csv": str(output_root / "biotransformation_generation_metrics.csv"),
        "per_control_csv": str(output_root / "biotransformation_generation_per_control.csv"),
        "summary_csv": str(output_root / "biotransformation_generation_summary.csv"),
        "report_md": str(output_root / "biotransformation_generation_report.md"),
    }
    write_json(output_root / "biotransformation_generation_summary.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--core-candidates-csv", type=Path, default=DEFAULT_OUTPUT_ROOT / "core_search_candidates.csv")
    parser.add_argument("--known-x", type=Path, default=Path("known_X.csv"))
    parser.add_argument("--known-x-all", type=Path, default=Path("known_X_all_available_structures_simple.csv"))
    parser.add_argument("--tracks", default=f"{TRACK_CONTROL_BLIND},{TRACK_ANNOTATION_ASSISTED}")
    parser.add_argument("--profiles", default="strict_1da,balanced")
    parser.add_argument("--controls", default="")
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--max-core-rows-per-control", type=int, default=80)
    parser.add_argument("--max-products-per-core", type=int, default=80)
    parser.add_argument("--max-candidates-per-control", type=int, default=700)
    run_generation(parser.parse_args())


if __name__ == "__main__":
    main()
