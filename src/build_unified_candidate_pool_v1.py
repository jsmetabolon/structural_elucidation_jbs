#!/usr/bin/env python3
# flake8: noqa
"""Build one deduplicated cross-source candidate pool for unified reranking."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from denovo_improvement_utils import canonicalize_smiles
from run_retrieval_anchored_denovo_assembly_v4 import resolve
from x_domain_analog_utils import dataframe_to_markdown, write_csv, write_json


DEFAULT_EXTRA_INPUTS = (
    "models/mlims_anchored_candidates_v1/runs/mlims_anchored_v1/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/mlims_strict_overlay_v1/candidate_metrics_all.csv",
    "models/targeted_combinatorial_v1/candidate_metrics_all.csv",
    "models/mass_aware_combinatorial_v1/candidate_metrics_all.csv",
    "models/biotransformation_gated_v2/runs/biotransformation_gated_v2/track_control_blind/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
    "models/biotransformation_gated_v2/runs/biotransformation_gated_v2/track_annotation_assisted/lib_209_negative/selfies_profiles/strict_1da/candidate_metrics_all.csv",
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


def prefer_verifier_csv(path: Path, repo_root: Path) -> Path:
    resolved = resolve(path, repo_root)
    if resolved.name.endswith("_with_verifier.csv"):
        return resolved
    verifier_path = resolved.with_name(resolved.stem + "_with_verifier.csv")
    return verifier_path if verifier_path.exists() else resolved


def default_extra_inputs_for_scope(lib_id: int, ion_mode: str) -> str:
    scope = f"lib_{int(lib_id)}_{str(ion_mode)}"
    if (int(lib_id), str(ion_mode)) == (209, "negative"):
        return ",".join(DEFAULT_EXTRA_INPUTS)
    return ",".join(
        [
            f"models/ec2_heavy_mlims_anchored_v1_lib{int(lib_id)}/runs/mlims_anchored_v1/track_control_blind/{scope}/selfies_profiles/strict_1da/candidate_metrics_all.csv",
            f"models/ec2_heavy_mlims_strict_overlay_v1_lib{int(lib_id)}/candidate_metrics_all.csv",
            f"models/ec2_heavy_targeted_combinatorial_v1_lib{int(lib_id)}/candidate_metrics_all.csv",
            f"models/ec2_heavy_mass_aware_combinatorial_v1_lib{int(lib_id)}/candidate_metrics_all.csv",
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
        ]
    )


def source_flags(labels: str) -> dict[str, int]:
    text = str(labels or "")
    return {
        "is_peptide": int("peptide_first" in text or "peptide_ranker" in text),
        "is_library_hmdb": int("library_anchored_v1" in text),
        "is_library_mlims": int("mlims_anchored" in text),
        "is_nonpep_gated": int("nonpeptide_gated" in text),
        "is_multi_conjugate": int("multi_conjugate" in text),
        "is_nce_direct": int("nce_direct" in text),
        "is_biotransformation": int("biotransformation" in text),
        "is_truth_blind_biotransformation": int(
            any(
                label in text
                for label in [
                    "sulfate_position_v1",
                    "taurine_variant_v1",
                    "lactoyl_amino_acid_v1",
                    "n_acetyl_v1",
                    "glucuronide_v1",
                    "glycine_v1",
                    "glutathione_v1",
                    "cysteinyl_v1",
                    "hydroxyl_v1",
                    "dehydro_v1",
                    "demethyl_v1",
                    "methylation_v1",
                    "wgcna_seeded_v1",
                ]
            )
        ),
    }


def load_inputs(scored_path: Path, extra_paths: list[Path], repo_root: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    if scored_path.exists():
        frames.append(pd.read_csv(scored_path, low_memory=False))
    for rel in extra_paths:
        path = prefer_verifier_csv(rel, repo_root)
        if not path.exists():
            continue
        try:
            df = pd.read_csv(path, low_memory=False)
        except pd.errors.EmptyDataError:
            continue
        if df.empty:
            continue
        if "source_label" not in df.columns:
            df["source_label"] = "mlims_anchored_v1" if "mlims_anchored" in str(path) else path.parent.name
        frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def truthy_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index, dtype=bool)
    raw = df[col]
    if pd.api.types.is_bool_dtype(raw):
        return raw.fillna(False).astype(bool)
    numeric_values = pd.to_numeric(raw, errors="coerce")
    text_values = raw.fillna("").astype(str).str.strip().str.lower()
    return numeric_values.fillna(0.0).gt(0.0) | text_values.isin({"true", "t", "yes", "y", "1"})


def normalize(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if out.empty:
        return out
    if "predicted_smiles" not in out.columns and "smiles" in out.columns:
        out["predicted_smiles"] = out["smiles"]
    if "control_id" not in out.columns:
        out["control_id"] = out.get("feature_key", pd.Series("", index=out.index)).astype(str).str.extract(r"(X-\d+)", expand=False).fillna("")
    if "source_label" not in out.columns:
        out["source_label"] = out.get("candidate_source", pd.Series("unknown", index=out.index)).fillna("unknown")
    existing = (out["canonical_2d_smiles"] if "canonical_2d_smiles" in out.columns else pd.Series("", index=out.index)).fillna("").astype(str)
    needs = existing.eq("") & out["predicted_smiles"].notna()
    if needs.any():
        existing.loc[needs] = out.loc[needs, "predicted_smiles"].map(lambda value: canonicalize_smiles(value) or "")
    out["canonical_2d_smiles"] = existing.fillna("").astype(str)
    out = out[out["canonical_2d_smiles"].ne("") & out["control_id"].astype(str).str.startswith("X-")].copy()
    out["_rank"] = numeric(out, "rank", default=999999)
    out["_score_for_keep"] = numeric(out, "predicted_spectrum_score_v2", default=0.0) + 0.25 * numeric(out, "rerank_score", default=0.0) + 0.15 * numeric(out, "generation_score", default=0.0)
    return out


def restrict_controls(df: pd.DataFrame, allowed_controls: set[str]) -> pd.DataFrame:
    if df.empty or not allowed_controls:
        return df
    return df[df["control_id"].astype(str).isin(allowed_controls)].copy()


def all_x_controls_for_scope(path: Path, lib_id: int) -> set[str]:
    if not path.exists():
        return set()
    df = pd.read_csv(path, low_memory=False)
    if "lib_id" in df.columns:
        df = df[pd.to_numeric(df["lib_id"], errors="coerce").eq(int(lib_id))].copy()
    series = df.get("x_compound_id", df.get("control_id", pd.Series(dtype=str))).dropna().astype(str)
    return {value for value in series if value.startswith("X-")}


def s0_controls_for_scope(path: Path, lib_id: int) -> set[str]:
    if not path.exists():
        return set()
    df = pd.read_csv(path, low_memory=False)
    if "lib_id" in df.columns:
        df = df[pd.to_numeric(df["lib_id"], errors="coerce").eq(int(lib_id))].copy()
    if "evidence_tier" not in df.columns or "x_compound_id" not in df.columns:
        return set()
    tiers = df["evidence_tier"].astype(str)
    return set(df[tiers.str.contains("S0", na=False)]["x_compound_id"].astype(str))


def dedup(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    work = df.copy()
    score_cols = ["predicted_spectrum_score_v2", "mass_match_score", "formula_prior_score", "evidence_confidence", "ms2_match_score", "nce_pipeline_score", "forward_spectrum_score", "rerank_score", "generation_score"]
    for col in score_cols:
        work[col] = numeric(work, col, default=0.0)
    keys = ["control_id", "canonical_2d_smiles"]
    label_agg = (
        work.groupby(keys, dropna=False)
        .agg(
            source_labels=("source_label", lambda values: "|".join(sorted({str(v) for v in values if str(v)}))),
            n_source_rows=("source_label", "count"),
            min_source_rank=("_rank", "min"),
        )
        .reset_index()
    )
    max_scores = work.groupby(keys, dropna=False)[score_cols].max().reset_index()
    best = work.sort_values(["_score_for_keep", "predicted_spectrum_score_v2", "_rank"], ascending=[False, False, True], kind="mergesort").drop_duplicates(keys, keep="first")
    drop_cols = ["source_labels", "n_source_rows", "min_source_rank"] + score_cols
    best = best.drop(columns=[col for col in drop_cols if col in best.columns], errors="ignore")
    out = best.merge(label_agg, on=keys, how="left").merge(max_scores, on=keys, how="left")
    flags = out["source_labels"].map(source_flags).apply(pd.Series)
    return pd.concat([out.reset_index(drop=True), flags.reset_index(drop=True)], axis=1)


def join_isotope_signals(pool: pd.DataFrame, isotope_signals_csv: Path | None) -> tuple[pd.DataFrame, dict[str, int]]:
    out = pool.copy()
    stats = {"n_isotope_signal_rows": 0, "n_isotope_joined_rows": 0}
    isotope_cols = ["isotope_formula_score", "isotope_ratio_score", "isotope_audit_status", "feature_is_likely_isotopologue", "feature_isotopologue_reason", "isotope_signal_count", "c_predict", "carbon_count_agreement"]
    if isotope_signals_csv is not None and isotope_signals_csv.exists():
        signals = pd.read_csv(isotope_signals_csv, low_memory=False)
        required = {"control_id", "canonical_2d_smiles"}
        if required.issubset(set(signals.columns)):
            keep_cols = ["control_id", "canonical_2d_smiles"] + [col for col in isotope_cols if col in signals.columns]
            signals = signals[keep_cols].copy()
            signals["control_id"] = signals["control_id"].astype(str)
            signals["canonical_2d_smiles"] = signals["canonical_2d_smiles"].map(lambda value: canonicalize_smiles(value) or str(value or ""))
            signals = signals.drop_duplicates(["control_id", "canonical_2d_smiles"], keep="first")
            out = out.drop(columns=[col for col in isotope_cols if col in out.columns], errors="ignore")
            out = out.merge(signals, on=["control_id", "canonical_2d_smiles"], how="left")
            stats["n_isotope_signal_rows"] = int(signals.shape[0])
            stats["n_isotope_joined_rows"] = int(out["isotope_audit_status"].notna().sum()) if "isotope_audit_status" in out.columns else 0
    out["isotope_formula_score"] = pd.to_numeric(out.get("isotope_formula_score", pd.Series(0.5, index=out.index)), errors="coerce").fillna(0.5)
    out["isotope_ratio_score"] = pd.to_numeric(out.get("isotope_ratio_score", pd.Series(0.5, index=out.index)), errors="coerce").fillna(0.5)
    out["c_predict"] = pd.to_numeric(out.get("c_predict", pd.Series(np.nan, index=out.index)), errors="coerce")
    out["carbon_count_agreement"] = pd.to_numeric(out.get("carbon_count_agreement", pd.Series(0.5, index=out.index)), errors="coerce").fillna(0.5)
    out["isotope_audit_status"] = out.get("isotope_audit_status", pd.Series("unavailable", index=out.index)).fillna("unavailable").astype(str)
    out["feature_is_likely_isotopologue"] = truthy_series(out, "feature_is_likely_isotopologue")
    if "feature_isotopologue_reason" not in out.columns:
        out["feature_isotopologue_reason"] = ""
    out["feature_isotopologue_reason"] = out["feature_isotopologue_reason"].fillna("").astype(str)
    out["isotope_signal_count"] = pd.to_numeric(out.get("isotope_signal_count", pd.Series(0, index=out.index)), errors="coerce").fillna(0).astype(int)
    return out, stats


def run(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parent
    output_root = resolve(args.output_root, repo_root)
    output_root.mkdir(parents=True, exist_ok=True)
    current = pd.read_csv(resolve(args.current_integrated_best, repo_root), low_memory=False)
    if "lib_id" in current.columns:
        current = current[pd.to_numeric(current["lib_id"], errors="coerce").eq(int(args.lib_id))].copy()
    if "ion_mode" in current.columns:
        current = current[current["ion_mode"].astype(str).eq(str(args.ion_mode))].copy()
    base_pool_path = resolve(args.base_unified_pool_csv, repo_root) if getattr(args, "base_unified_pool_csv", None) else None
    if base_pool_path is not None:
        raw = pd.read_csv(base_pool_path, low_memory=False)
        norm = raw.copy()
        pool = raw.copy()
        if "control_id" not in pool.columns:
            pool["control_id"] = pool.get("feature_key", pd.Series("", index=pool.index)).astype(str).str.extract(r"(X-\d+)", expand=False).fillna("")
        pool["control_id"] = pool["control_id"].fillna("").astype(str)
        if "canonical_2d_smiles" not in pool.columns:
            if "predicted_smiles" in pool.columns:
                pool["canonical_2d_smiles"] = pool["predicted_smiles"].map(lambda value: canonicalize_smiles(value) or "")
            else:
                pool["canonical_2d_smiles"] = ""
        pool["canonical_2d_smiles"] = pool["canonical_2d_smiles"].fillna("").astype(str)
        pool = pool[pool["control_id"].ne("") & pool["canonical_2d_smiles"].ne("")].copy()
    else:
        extra_text = str(args.extra_inputs or "")
        if extra_text == ",".join(DEFAULT_EXTRA_INPUTS) and (int(args.lib_id), str(args.ion_mode)) != (209, "negative"):
            extra_text = default_extra_inputs_for_scope(int(args.lib_id), str(args.ion_mode))
        raw = load_inputs(resolve(args.scored_candidates, repo_root), [Path(part.strip()) for part in extra_text.split(",") if part.strip()], repo_root)
        allowed_controls = set(current["control_id"].astype(str)) if "control_id" in current.columns else set()
        if bool(args.include_all_x):
            allowed_controls |= all_x_controls_for_scope(resolve(args.all_x_catalog, repo_root), int(args.lib_id))
        norm = restrict_controls(normalize(raw), allowed_controls)
        pool = dedup(norm)
    isotope_signals_path = resolve(args.isotope_signals_csv, repo_root) if getattr(args, "isotope_signals_csv", None) else None
    pool, isotope_stats = join_isotope_signals(pool, isotope_signals_path)
    n_isotope_formula_impossible_dropped = 0
    if bool(getattr(args, "isotope_filter_impossible", True)) and not pool.empty:
        impossible_mask = pool["isotope_audit_status"].astype(str).eq("formula_impossible")
        n_isotope_formula_impossible_dropped = int(impossible_mask.sum())
        if n_isotope_formula_impossible_dropped:
            pool = pool[~impossible_mask].copy()
    if not pool.empty:
        sort_cols = [col for col in ["control_id", "predicted_spectrum_score_v2", "mass_match_score", "min_source_rank"] if col in pool.columns]
        ascending_lookup = {"control_id": True, "predicted_spectrum_score_v2": False, "mass_match_score": False, "min_source_rank": True}
        pool = pool.sort_values(sort_cols, ascending=[ascending_lookup[col] for col in sort_cols], kind="mergesort")
        pool["unified_pool_rank"] = pool.groupby("control_id", dropna=False).cumcount() + 1
    write_csv(output_root / "unified_candidates_per_control.csv", pool)
    s0_controls = s0_controls_for_scope(resolve(args.all_x_catalog, repo_root), int(args.lib_id))
    pool_s0 = pool[pool.get("control_id", pd.Series(dtype=str)).astype(str).isin(s0_controls)].copy() if not pool.empty and s0_controls else pd.DataFrame(columns=pool.columns)
    pool_non_s0 = pool[~pool.get("control_id", pd.Series(dtype=str)).astype(str).isin(s0_controls)].copy() if not pool.empty and s0_controls else pd.DataFrame(columns=pool.columns)
    write_csv(output_root / "unified_candidates_per_control_s0.csv", pool_s0)
    write_csv(output_root / "unified_candidates_per_control_non_s0.csv", pool_non_s0)
    winners = current[pd.to_numeric(current.get("top10_best_tanimoto"), errors="coerce").fillna(0.0).ge(0.7)].copy()
    preserve_rows: list[dict[str, Any]] = []
    pool_keys = set(zip(pool.get("control_id", pd.Series(dtype=str)).astype(str), pool.get("canonical_2d_smiles", pd.Series(dtype=str)).astype(str))) if not pool.empty else set()
    for _, row in winners.iterrows():
        key = (str(row.get("control_id", "")), canonicalize_smiles(row.get("top10_smiles", "")) or "")
        preserve_rows.append({"control_id": key[0], "winner_present": key in pool_keys, "winner_smiles": row.get("top10_smiles", "")})
    preserve = pd.DataFrame(preserve_rows)
    write_csv(output_root / "preserve_25_pool_check.csv", preserve)
    n_expected = int(preserve.shape[0])
    n_present = int(preserve["winner_present"].sum()) if not preserve.empty else 0
    summary = {
        "base_unified_pool_csv": str(base_pool_path) if base_pool_path is not None else "",
        "n_base_rows": int(raw.shape[0]) if base_pool_path is not None else 0,
        "n_raw_rows": int(raw.shape[0]) if not raw.empty else 0,
        "n_normalized_rows": int(norm.shape[0]) if not norm.empty else 0,
        "n_unified_rows": int(pool.shape[0]) if not pool.empty else 0,
        "n_controls": int(pool["control_id"].nunique()) if not pool.empty else 0,
        "n_s0_rows": int(pool_s0.shape[0]),
        "n_non_s0_rows": int(pool_non_s0.shape[0]),
        "n_non_s0_controls": int(pool_non_s0["control_id"].nunique()) if not pool_non_s0.empty else 0,
        "n_isotope_signal_rows": int(isotope_stats.get("n_isotope_signal_rows", 0)),
        "n_isotope_joined_rows": int(isotope_stats.get("n_isotope_joined_rows", 0)),
        "n_isotope_formula_impossible_dropped": int(n_isotope_formula_impossible_dropped),
        "n_current_winners_present": n_present,
        "n_current_winners_expected": n_expected,
        "acceptance_pass": bool(n_expected == 0 or n_present >= min(25, n_expected)),
    }
    write_json(output_root / "unified_candidate_pool_summary.json", summary)
    source_label_counts = pool["source_labels"].value_counts().rename_axis("source_labels").reset_index(name="n_rows") if not pool.empty and "source_labels" in pool.columns else pd.DataFrame()
    report = [
        "# Unified Candidate Pool V1",
        "",
        "## Summary",
        dataframe_to_markdown(pd.DataFrame([summary])),
        "",
        "## Source Label Counts",
        dataframe_to_markdown(source_label_counts),
        "",
    ]
    (output_root / "unified_candidate_pool_report.md").write_text("\n".join(report))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=Path("models/unified_candidate_pool_v1"))
    parser.add_argument("--scored-candidates", type=Path, default=Path("models/verifier_scoring_v1/combined_scored_candidates.csv"))
    parser.add_argument("--extra-inputs", default=",".join(DEFAULT_EXTRA_INPUTS))
    parser.add_argument("--current-integrated-best", type=Path, default=Path("models/integrated_denovo_audit_v1/integrated_per_control_best.csv"))
    parser.add_argument("--include-all-x", action="store_true")
    parser.add_argument("--all-x-catalog", type=Path, default=Path("models/x_known_curation_v1/known_x_evidence_catalog.csv"))
    parser.add_argument("--isotope-signals-csv", type=Path, default=None)
    parser.add_argument("--base-unified-pool-csv", type=Path, default=None)
    parser.add_argument("--isotope-filter-impossible", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--lib-id", type=int, default=209)
    parser.add_argument("--ion-mode", default="negative")
    payload = run(parser.parse_args())
    raise SystemExit(0 if payload["acceptance_pass"] else 2)


if __name__ == "__main__":
    main()
