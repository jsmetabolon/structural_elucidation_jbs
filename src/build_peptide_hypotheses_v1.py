#!/usr/bin/env python3
"""Build upstream peptide hypotheses before candidate generation."""

from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from peptide_first_v1_utils import (
    DEFAULT_OUTPUT_ROOT,
    build_peptide_composition_index,
    feature_rows,
    load_adduct_hypotheses,
    match_peptide_composition_index,
    ordered_sequences_for_composition,
    parse_legacy_annotations,
    ri_plausibility,
    spectrum_score_for_sequence,
)
from peptide_ranker_v1_utils import parse_bins
from peptide_spectrum_supervision_v1_utils import event_type_scores, load_calibration
from run_retrieval_anchored_denovo_assembly_v4 import resolve
from x_domain_analog_utils import dataframe_to_markdown, safe_float, safe_int, write_csv, write_json


def score_hypothesis(row: pd.Series, feature: pd.Series, calibration: dict[str, Any] | None) -> dict[str, Any]:
    ion_bins = parse_bins(feature.get("observed_ion_bins_json"))
    loss_bins = parse_bins(feature.get("observed_loss_bins_json"))
    sequences = ordered_sequences_for_composition(row.get("normalized_peptide_composition"), max_permutations=24)
    sequence_scores = [(sequence, spectrum_score_for_sequence(sequence, ion_bins, loss_bins, calibration)) for sequence in sequences]
    best_sequence, best_score = max(sequence_scores, key=lambda item: item[1], default=("", 0.0))
    event_scores = event_type_scores(best_sequence, ion_bins, loss_bins) if best_sequence else {}
    ri_score = ri_plausibility(feature.get("ri"), safe_int(row.get("peptide_length"), 0), str(row.get("modification", "none")))
    return {
        "best_ordered_sequence": best_sequence,
        "n_ordered_permutations": int(len(sequences)),
        "peptide_spectrum_score": round(float(best_score), 6),
        "ri_plausibility": round(float(ri_score), 6),
        "observed_ion_bins_json": feature.get("observed_ion_bins_json", "[]"),
        "observed_loss_bins_json": feature.get("observed_loss_bins_json", "[]"),
        **event_scores,
    }


def hypothesis_confidence(row: pd.Series) -> float:
    return float(
        np.clip(
            0.34 * safe_float(row.get("mass_score"), 0.0)
            + 0.24 * safe_float(row.get("peptide_spectrum_score"), 0.0)
            + 0.18 * safe_float(row.get("adduct_score"), 0.0)
            + 0.12 * safe_float(row.get("ri_plausibility"), 0.5)
            + 0.12 * safe_float(row.get("annotation_score"), 0.0),
            0.0,
            1.0,
        )
    )


def numeric_column(df: pd.DataFrame, column: str, default: float = 0.0) -> pd.Series:
    if column not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[column], errors="coerce").fillna(default)


def lengths_for_breadth(args: argparse.Namespace) -> tuple[int, ...]:
    lengths = tuple(int(part) for part in str(args.lengths).split(",") if part.strip())
    if str(args.composition_breadth) == "v2":
        return tuple(sorted(set(lengths + (5, 6))))
    return lengths


def annotate_v1_1_eligibility(df: pd.DataFrame, *, mass_tolerance_da: float) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out["peptide_neutral_mass_fit"] = numeric_column(out, "mass_error_da", math.nan).abs().le(float(mass_tolerance_da))
    out["calibrated_by_evidence"] = (
        numeric_column(out, "spectrum_b_ion_count")
        + numeric_column(out, "spectrum_y_ion_count")
        + numeric_column(out, "spectrum_b_loss_count")
        + numeric_column(out, "spectrum_y_loss_count")
    ).gt(0)
    out["immonium_evidence"] = (
        numeric_column(out, "spectrum_immonium_count")
        + numeric_column(out, "spectrum_n_terminal_immonium_count")
        + numeric_column(out, "spectrum_c_terminal_immonium_count")
        + numeric_column(out, "spectrum_pyr_specific_count")
    ).gt(0)
    out["neutral_loss_evidence"] = numeric_column(out, "spectrum_neutral_loss_count").gt(0)
    annotation_text = out.get("annotation_text", pd.Series("", index=out.index)).astype(str).str.lower()
    out["legacy_motif_support"] = numeric_column(out, "annotation_score").gt(0.0) | annotation_text.str.contains(
        "peptide|tripeptide|dipeptide|tetrapeptide|leu/ile|pyr|pyroglu",
        regex=True,
    )
    out["high_confidence_residue_composition_support"] = (
        numeric_column(out, "mass_score").ge(0.95)
        & numeric_column(out, "adduct_score").ge(0.70)
        & numeric_column(out, "ri_plausibility").ge(0.35)
    )
    out["peptide_v1_1_evidence_count"] = (
        out["calibrated_by_evidence"].astype(int)
        + out["immonium_evidence"].astype(int)
        + out["neutral_loss_evidence"].astype(int)
        + out["legacy_motif_support"].astype(int)
        + out["high_confidence_residue_composition_support"].astype(int)
    )
    out["peptide_v1_1_eligible"] = out["peptide_neutral_mass_fit"] & out["peptide_v1_1_evidence_count"].gt(0)
    ms1_buddy = out.get("ms1_buddy_supported", pd.Series(False, index=out.index)).fillna(False).astype(str).str.lower().isin(
        {"true", "1", "yes", "y"}
    )
    out["peptide_v1_2_mass_score_rescue"] = numeric_column(out, "mass_score").ge(0.95)
    out["peptide_v1_2_ms1_buddy_rescue"] = ms1_buddy
    out["peptide_v1_2_eligible"] = out["peptide_neutral_mass_fit"] & (
        out["peptide_v1_1_evidence_count"].ge(1)
        | out["peptide_v1_2_mass_score_rescue"]
        | out["peptide_v1_2_ms1_buddy_rescue"]
    )
    return out


def build_control_blind_hypotheses(
    *,
    features: pd.DataFrame,
    adducts: pd.DataFrame,
    composition_index: pd.DataFrame,
    calibration: dict[str, Any] | None,
    mass_tolerance_da: float,
    max_compositions_per_feature: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    feature_map = {str(row.feature_key): row._asdict() for row in features.itertuples(index=False)}
    for adduct in adducts.itertuples(index=False):
        feature = pd.Series(feature_map.get(str(adduct.feature_key), {}), dtype=object)
        if feature.empty:
            continue
        matches = match_peptide_composition_index(composition_index, safe_float(adduct.neutral_mass, math.nan), mass_tolerance_da)
        if matches.empty:
            continue
        matches = matches.head(int(max_compositions_per_feature) * 3)
        for _, match in matches.iterrows():
            extra = score_hypothesis(match, feature, calibration)
            row = {
                **feature.to_dict(),
                **match.to_dict(),
                **extra,
                "track": "control_blind",
                "adduct_label": getattr(adduct, "adduct_label", "[M-H]-"),
                "neutral_mass": round(safe_float(getattr(adduct, "neutral_mass", math.nan), math.nan), 6),
                "adduct_score": round(safe_float(getattr(adduct, "adduct_score", 0.5), 0.5), 6),
                "adduct_rank": safe_int(getattr(adduct, "adduct_rank", 1), 1),
                "formula_prior_score": round(safe_float(getattr(adduct, "formula_prior_score", 0.0), 0.0), 6),
                "mass_adduct_score": round(safe_float(getattr(adduct, "mass_adduct_score", 0.0), 0.0), 6),
                "ms1_buddy_supported": bool(str(getattr(adduct, "ms1_buddy_supported", False)).lower() in {"true", "1", "yes", "y"}),
                "annotation_score": 0.0,
                "annotation_text": "",
                "hypothesis_source": "mass_adduct_spectrum",
            }
            row["hypothesis_confidence"] = round(hypothesis_confidence(pd.Series(row)), 6)
            rows.append(row)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.sort_values(
        ["feature_key", "track", "hypothesis_confidence", "mass_score", "peptide_spectrum_score", "peptide_length"],
        ascending=[True, True, False, False, False, True],
        kind="mergesort",
    )
    out["hypothesis_rank"] = out.groupby(["feature_key", "track"], dropna=False).cumcount() + 1
    out = out[out["hypothesis_rank"].le(int(max_compositions_per_feature))].reset_index(drop=True)
    return annotate_v1_1_eligibility(out, mass_tolerance_da=mass_tolerance_da)


def build_annotation_hypotheses(
    *,
    features: pd.DataFrame,
    adducts: pd.DataFrame,
    legacy: pd.DataFrame,
    calibration: dict[str, Any] | None,
    mass_tolerance_da: float,
    max_compositions_per_feature: int,
) -> pd.DataFrame:
    if legacy.empty:
        return pd.DataFrame()
    best_adduct = (
        adducts.sort_values(["feature_key", "adduct_rank", "adduct_score"], ascending=[True, True, False], kind="mergesort")
        .drop_duplicates("feature_key", keep="first")
        .copy()
    )
    merged = features.merge(legacy, on=["comp_id", "lib_id"], how="inner").merge(best_adduct, on="feature_key", how="left", suffixes=("", "_adduct"))
    rows: list[dict[str, Any]] = []
    for _, feature in merged.iterrows():
        comp = str(feature.get("legacy_composition", "") or "")
        if not comp:
            continue
        mods = [part for part in str(feature.get("legacy_modification", "none") or "none").split("|") if part] or ["none"]
        for mod in mods:
            from peptide_first_v1_utils import composition_tokens, peptide_formula, peptide_mass

            tokens = composition_tokens(comp)
            mass = peptide_mass(tokens, mod)
            neutral = safe_float(feature.get("neutral_mass"), safe_float(feature.get("neutral_mass_fallback"), math.nan))
            error = mass - neutral if np.isfinite(neutral) else math.nan
            denom = 0.15 if np.isfinite(error) else 1.0
            mass_score = max(0.0, 1.0 - abs(error) / denom) if np.isfinite(error) else 0.0
            base = pd.Series(
                {
                    "normalized_peptide_composition": comp,
                    "peptide_length": len(tokens),
                    "modification": mod,
                    "peptide_formula": peptide_formula(tokens, mod),
                    "theoretical_neutral_mass": round(float(mass), 6),
                    "mass_error_da": round(float(error), 6) if np.isfinite(error) else math.nan,
                    "mass_score": round(float(mass_score), 6),
                    "smiles_supported": mod in {"none", "n_acetyl", "glycoside"},
                }
            )
            extra = score_hypothesis(base, feature, calibration)
            row = {
                **feature.to_dict(),
                **base.to_dict(),
                **extra,
                "track": "annotation_assisted",
                "adduct_label": feature.get("adduct_label", "[M-H]-"),
                "neutral_mass": round(neutral, 6) if np.isfinite(neutral) else math.nan,
                "adduct_score": round(safe_float(feature.get("adduct_score"), 0.5), 6),
                "adduct_rank": safe_int(feature.get("adduct_rank"), 1),
                "formula_prior_score": round(safe_float(feature.get("formula_prior_score"), 0.0), 6),
                "mass_adduct_score": round(safe_float(feature.get("mass_adduct_score"), 0.0), 6),
                "ms1_buddy_supported": bool(str(feature.get("ms1_buddy_supported", False)).lower() in {"true", "1", "yes", "y"}),
                "annotation_score": 1.0,
                "annotation_text": feature.get("details", ""),
                "hypothesis_source": "legacy_annotation",
            }
            row["hypothesis_confidence"] = round(hypothesis_confidence(pd.Series(row)), 6)
            rows.append(row)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.sort_values(["feature_key", "track", "hypothesis_confidence", "mass_score"], ascending=[True, True, False, False], kind="mergesort")
    out["hypothesis_rank"] = out.groupby(["feature_key", "track"], dropna=False).cumcount() + 1
    out = out[out["hypothesis_rank"].le(int(max_compositions_per_feature))].reset_index(drop=True)
    return annotate_v1_1_eligibility(out, mass_tolerance_da=mass_tolerance_da)


def render_report(hypotheses: pd.DataFrame, summary: pd.DataFrame, legacy: pd.DataFrame) -> str:
    return "\n".join(
        [
            "# Peptide-First Hypotheses V1",
            "",
            "Control-blind hypotheses use only neutral-mass/adduct, spectrum, and RI evidence. Legacy NCE/MS2SMILES annotations are isolated to `annotation_assisted`.",
            "",
            "## Summary",
            "",
            dataframe_to_markdown(summary),
            "",
            "## Legacy Annotation Rows",
            "",
            dataframe_to_markdown(legacy.head(40)),
            "",
            "## Top Hypotheses",
            "",
            dataframe_to_markdown(hypotheses.head(60)),
            "",
        ]
    )


def run(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parent
    output_root = resolve(args.output_root, repo_root)
    output_root.mkdir(parents=True, exist_ok=True)
    feature_table = pd.read_csv(resolve(args.feature_table, repo_root), low_memory=False)
    features = feature_rows(feature_table, lib_id=int(args.lib_id))
    adducts = load_adduct_hypotheses(
        resolve(args.adduct_evidence_root, repo_root),
        features,
        lib_id=int(args.lib_id),
        top_n=int(args.top_adducts),
    )
    calibration_path = resolve(args.spectrum_calibration_root, repo_root)
    calibration = load_calibration(calibration_path) if calibration_path.exists() else None
    legacy_path = resolve(args.legacy_unknown_list, repo_root)
    legacy_inputs = output_root / "legacy_inputs"
    legacy_inputs.mkdir(parents=True, exist_ok=True)
    if legacy_path.exists():
        target = legacy_inputs / "unknown_list.csv"
        if legacy_path.resolve() != target.resolve():
            shutil.copy2(legacy_path, target)
    legacy = parse_legacy_annotations(legacy_path, lib_id=int(args.lib_id))
    composition_index = build_peptide_composition_index(lengths=lengths_for_breadth(args))
    control_blind = build_control_blind_hypotheses(
        features=features,
        adducts=adducts,
        composition_index=composition_index,
        calibration=calibration,
        mass_tolerance_da=float(args.mass_tolerance_da),
        max_compositions_per_feature=int(args.max_compositions_per_feature),
    )
    assisted = build_annotation_hypotheses(
        features=features,
        adducts=adducts,
        legacy=legacy,
        calibration=calibration,
        mass_tolerance_da=float(args.mass_tolerance_da),
        max_compositions_per_feature=int(args.max_compositions_per_feature),
    )
    hypotheses = pd.concat([control_blind, assisted], ignore_index=True, sort=False) if not assisted.empty else control_blind
    eligibility_policy = str(args.eligibility_policy).lower()
    if eligibility_policy == "v1_1" and not hypotheses.empty:
        hypotheses = hypotheses[hypotheses["peptide_v1_1_eligible"].astype(bool)].reset_index(drop=True)
    elif eligibility_policy == "v1_2" and not hypotheses.empty:
        hypotheses = hypotheses[hypotheses["peptide_v1_2_eligible"].astype(bool)].reset_index(drop=True)
    summary = (
        hypotheses.groupby("track", dropna=False)
        .agg(
            n_hypotheses=("feature_key", "count"),
            n_features=("feature_key", "nunique"),
            mean_confidence=("hypothesis_confidence", "mean"),
            n_smiles_supported=("smiles_supported", "sum"),
        )
        .reset_index()
        if not hypotheses.empty
        else pd.DataFrame()
    )
    write_csv(output_root / "feature_peptide_rows.csv", features)
    write_csv(output_root / "peptide_neutral_mass_adducts.csv", adducts)
    write_csv(output_root / "legacy_peptide_annotations.csv", legacy)
    write_csv(output_root / "peptide_hypotheses.csv", hypotheses)
    write_csv(output_root / "peptide_hypothesis_summary.csv", summary)
    (output_root / "peptide_hypothesis_report.md").write_text(render_report(hypotheses, summary, legacy))
    payload = {
        "output_root": str(output_root),
        "n_features": int(features.shape[0]),
        "n_adduct_rows": int(adducts.shape[0]),
        "n_control_blind_hypotheses": int(control_blind.shape[0]),
        "n_annotation_assisted_hypotheses": int(assisted.shape[0]),
        "eligibility_policy": str(args.eligibility_policy),
        "n_hypotheses_after_eligibility": int(hypotheses.shape[0]),
        "n_legacy_annotations": int(legacy.shape[0]),
        "composition_breadth": str(args.composition_breadth),
        "n_peptide_v1_2_eligible": int(hypotheses.get("peptide_v1_2_eligible", pd.Series(dtype=bool)).astype(bool).sum()) if not hypotheses.empty else 0,
    }
    write_json(output_root / "peptide_hypothesis_summary.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--lib-id", type=int, default=209)
    parser.add_argument("--feature-table", type=Path, default=Path("merged_encoder/data/ions_unknown_merged.csv"))
    parser.add_argument("--adduct-evidence-root", type=Path, default=Path("models/upstream_adduct_evidence_v1"))
    parser.add_argument("--spectrum-calibration-root", type=Path, default=Path("models/peptide_spectrum_supervision_v1/calibration"))
    parser.add_argument("--legacy-unknown-list", type=Path, default=DEFAULT_OUTPUT_ROOT / "legacy_inputs/unknown_list.csv")
    parser.add_argument("--lengths", default="2,3,4")
    parser.add_argument("--composition-breadth", choices=["v1", "v2"], default="v1")
    parser.add_argument("--top-adducts", type=int, default=3)
    parser.add_argument("--mass-tolerance-da", type=float, default=0.05)
    parser.add_argument("--max-compositions-per-feature", type=int, default=12)
    parser.add_argument("--eligibility-policy", choices=["v1", "v1_1", "v1_2"], default="v1")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
