#!/usr/bin/env python3
"""Build biotransformation hypotheses for 209 deconjugation-first generation."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from biotransformation_guided_v1_utils import (
    TRACK_ANNOTATION_ASSISTED,
    TRACK_CONTROL_BLIND,
    TRANSFORMS,
    background_ion_matches,
    evidence_track_allowed,
    formula_subtract,
    neutral_mass_hypotheses,
    parse_json_int_list,
    peptide_compositions_for_mass,
    summarize_transform_counts,
    transform_by_alias,
    transform_evidence_score,
)
from denovo_improvement_utils import formula_from_smiles
from nce_evidence_utils import formula_mass, lib_ion_mode
from run_retrieval_anchored_denovo_assembly_v4 import resolve
from x_domain_analog_utils import dataframe_to_markdown, safe_float, safe_int, write_csv, write_json


DEFAULT_OUTPUT_ROOT = Path("models/biotransformation_guided_denovo_v1")
DEFAULT_UNKNOWN_LIST = Path("shared/legacy_ms2smiles/unknown_list.csv")
DEFAULT_FEATURE_EVIDENCE = (
    "models/nce_evidence_promotion_v4/feature_evidence_v4.csv,"
    "models/nce_evidence_to_generation_v3/evidence_bundle/feature_evidence_bundle.csv,"
    "models/nce_workflow_leverage_v2/feature_evidence.csv,"
    "models/spectrum_preprocessing_sensitivity_v1/evidence/multiview_spectrum_evidence.csv"
)


def parse_csv(text: Any) -> list[str]:
    return [part.strip() for part in str(text or "").split(",") if part.strip()]


def parse_formula_from_text(text: Any) -> str:
    match = re.search(r"\b(C\d*H\d*(?:N\d*)?(?:O\d*)?(?:P\d*)?(?:S\d*)?(?:Cl\d*)?(?:Br\d*)?(?:F\d*)?(?:I\d*)?)\b", str(text or ""))
    return match.group(1) if match else ""


def parse_core_formula_from_text(text: Any, transform_name: str) -> str:
    lower = str(text or "")
    if transform_name in {"sulfate", "disulfate", "glucuronide", "glycoside", "phosphate"}:
        formula_pattern = r"C\d*H\d*(?:N\d*)?(?:O\d*)?(?:P\d*)?(?:S\d*)?(?:Cl\d*)?(?:Br\d*)?(?:F\d*)?(?:I\d*)?"
        pattern = rf"{transform_name.split('_')[0]}(?:\w*)\s+of\s+({formula_pattern})"
        match = re.search(pattern, lower, flags=re.IGNORECASE)
        if match:
            return match.group(1)
    return ""


def load_target_controls(path: Path) -> set[str]:
    if not path or str(path) == "." or not path.exists():
        return set()
    raw = pd.read_csv(path, low_memory=False)
    for col in ("control_id", "x_compound_id", "Unknown"):
        if col in raw.columns:
            return {str(value).strip() for value in raw[col] if str(value).strip().startswith("X-")}
    return set()


def load_verified_targets(
    known_x: Path,
    known_x_all: Path,
    libs: set[int],
    ion_mode: str,
    *,
    include_unverified_x: bool = False,
    target_controls: set[str] | None = None,
) -> pd.DataFrame:
    known = pd.read_csv(known_x)
    known_all = pd.read_csv(known_x_all)
    s0 = set(known_all["x_compound_id"].astype(str)) if "x_compound_id" in known_all.columns else set()
    smiles_lookup = (
        known_all[["x_compound_id", "smiles"]].rename(columns={"x_compound_id": "control_id", "smiles": "known_smiles"}).copy()
        if {"x_compound_id", "smiles"}.issubset(known_all.columns)
        else pd.DataFrame(columns=["control_id", "known_smiles"])
    )
    work = known[pd.to_numeric(known.get("LIB_ID"), errors="coerce").isin(libs)].copy()
    work["control_id"] = work["Unknown"].astype(str)
    work["lib_id"] = pd.to_numeric(work.get("LIB_ID"), errors="coerce").fillna(0).astype(int)
    work["ion_mode"] = str(ion_mode)
    if not include_unverified_x:
        work = work[work["control_id"].isin(s0)].copy()
    if target_controls:
        work = work[work["control_id"].isin(target_controls)].copy()
    for col in ["m/z", "Formula", "External Information (Class or modification)", "Tentative information", "Cleaned External Info", "KNOWN", "SubLevel"]:
        if col not in work.columns:
            work[col] = ""
    work = work.merge(smiles_lookup, on="control_id", how="left")
    return work.drop_duplicates("control_id").reset_index(drop=True)


def load_unknown_list(path: Path, libs: set[int]) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["control_id"])
    df = pd.read_csv(
        path,
        header=None,
        names=["legacy_comp_id", "lib_id", "legacy_parent_mz", "legacy_adduct_hint", "legacy_note", "legacy_unknown"],
        dtype=str,
        on_bad_lines="skip",
    )
    df["xnum"] = df["legacy_unknown"].astype(str).str.extract(r"Unknown\\s*-\\s*(\\d+)")[0]
    df = df[df["xnum"].notna()].copy()
    df["control_id"] = "X-" + df["xnum"].astype(str)
    df["lib_id"] = pd.to_numeric(df["lib_id"], errors="coerce")
    return df[df["lib_id"].isin(libs)].copy()


def load_feature_evidence(repo_root: Path, paths_text: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for raw in parse_csv(paths_text):
        path = resolve(Path(raw), repo_root)
        if not path.exists():
            continue
        try:
            frame = pd.read_csv(path, low_memory=False)
        except Exception:
            continue
        if frame.empty:
            continue
        if "feature_key" not in frame.columns:
            frame["feature_key"] = ""
        frame["control_id"] = frame["feature_key"].astype(str).str.extract(r"(X-\\d+)")[0]
        if "comp_id" in frame.columns:
            frame["_comp_id"] = pd.to_numeric(frame["comp_id"], errors="coerce")
        frame["_feature_evidence_path"] = str(path)
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["control_id"])
    out = pd.concat(frames, ignore_index=True, sort=False)
    return out.drop_duplicates(["control_id", "feature_key", "_feature_evidence_path"], keep="first")


def best_feature_row(control: pd.Series, evidence: pd.DataFrame) -> pd.Series:
    if evidence.empty:
        return pd.Series(dtype=object)
    by_xid = evidence[evidence.get("control_id", pd.Series(dtype=str)).astype(str).eq(str(control.get("control_id")))].copy()
    if by_xid.empty and "CompID" in control.index and "_comp_id" in evidence.columns:
        by_xid = evidence[evidence["_comp_id"].eq(safe_float(control.get("CompID"), -1))].copy()
    if by_xid.empty:
        return pd.Series(dtype=object)
    return by_xid.iloc[0]


def transform_sources_for_control(control: pd.Series, unknown: pd.Series, feature: pd.Series) -> list[dict[str, Any]]:
    sources: list[dict[str, Any]] = []
    annotation_text = " | ".join(
        str(value or "")
        for value in [
            control.get("External Information (Class or modification)", ""),
            control.get("Tentative information", ""),
            control.get("Cleaned External Info", ""),
            control.get("KNOWN", ""),
            control.get("SubLevel", ""),
            unknown.get("legacy_note", "") if not unknown.empty else "",
        ]
    )
    for name in transform_by_alias(annotation_text):
        sources.append({"transform_name": name, "evidence_source": "annotation_unknown_list_known_x", "raw_text": annotation_text, "base_score": 0.92})
    feature_text = " | ".join(
        str(feature.get(col, "") or "")
        for col in [
            "motif_top_label",
            "motif_labels_json",
            "legacy_top_motif_v2",
            "subsmiles_motifs_json",
            "recursive_top_labels_json",
        ]
        if col in feature.index
    )
    for name in transform_by_alias(feature_text):
        sources.append({"transform_name": name, "evidence_source": "feature_motif_model", "raw_text": feature_text, "base_score": 0.65})
    for name, spec in TRANSFORMS.items():
        score = transform_evidence_score(spec, feature) if not feature.empty else 0.0
        if score >= 0.55:
            sources.append({"transform_name": name, "evidence_source": "diagnostic_ion_loss", "raw_text": "", "base_score": score})
    if not sources:
        sources.append({"transform_name": "none", "evidence_source": "no_transform_evidence", "raw_text": "", "base_score": 0.0})
    return sources


def build_evidence(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parent
    output_root = resolve(args.output_root, repo_root)
    output_root.mkdir(parents=True, exist_ok=True)
    known_x = resolve(args.known_x, repo_root)
    known_x_all = resolve(args.known_x_all, repo_root)
    libs = {int(value) for value in parse_csv(args.libs)}
    target_controls = load_target_controls(resolve(args.target_controls_csv, repo_root)) if str(args.target_controls_csv) else set()
    targets = load_verified_targets(
        known_x,
        known_x_all,
        libs,
        str(args.ion_mode),
        include_unverified_x=bool(args.include_unverified_x),
        target_controls=target_controls,
    )
    unknowns = load_unknown_list(resolve(args.unknown_list, repo_root), libs)
    features = load_feature_evidence(repo_root, args.feature_evidence_paths)
    rows: list[dict[str, Any]] = []
    for _, control in targets.iterrows():
        xid = str(control.get("control_id"))
        unknown = unknowns[unknowns["control_id"].astype(str).eq(xid)].head(1)
        unknown_row = unknown.iloc[0] if not unknown.empty else pd.Series(dtype=object)
        feature_row = best_feature_row(control, features)
        parent_mz = safe_float(control.get("m/z"), math.nan)
        if not np.isfinite(parent_mz):
            parent_mz = safe_float(unknown_row.get("legacy_parent_mz"), math.nan)
        if not np.isfinite(parent_mz):
            parent_mz = safe_float(feature_row.get("parent_mz"), math.nan)
        feature_key = str(feature_row.get("feature_key", "") or "")
        if not feature_key:
            feature_key = f"lib_209_negative_{xid}_mz{parent_mz:.5f}" if np.isfinite(parent_mz) else f"lib_209_negative_{xid}"
        adduct_hint = str(unknown_row.get("legacy_adduct_hint", "") or "")
        lib_id = safe_int(control.get("lib_id"), next(iter(libs)) if len(libs) == 1 else 0)
        ion_mode = str(args.ion_mode) or lib_ion_mode(lib_id)
        formula = str(control.get("Formula", "") or "").strip()
        if not formula or formula.lower() == "nan":
            formula = parse_formula_from_text(unknown_row.get("legacy_note", ""))
        if not formula or formula.lower() == "nan":
            formula = formula_from_smiles(str(control.get("known_smiles", "") or "")) or ""
        mass_hypotheses = neutral_mass_hypotheses(parent_mz, lib_id, adduct_hint=adduct_hint)
        if not mass_hypotheses:
            fallback_adduct = "[M-H]-" if ion_mode.startswith("neg") else "[M+H]+"
            mass_hypotheses = [{"mass_hypothesis_adduct": fallback_adduct, "neutral_mass": math.nan, "adduct_priority": 0.0, "adduct_source": "missing_parent_mz"}]
        if formula and all(not np.isfinite(safe_float(hyp.get("neutral_mass"), math.nan)) for hyp in mass_hypotheses):
            neutral_mass = formula_mass(formula)
            if np.isfinite(neutral_mass):
                fallback_adduct = "[M-H]-" if ion_mode.startswith("neg") else "[M+H]+"
                mass_hypotheses = [{
                    "mass_hypothesis_adduct": fallback_adduct,
                    "neutral_mass": neutral_mass,
                    "adduct_priority": 1.0,
                    "adduct_source": "known_formula_fallback",
                }]
        ion_bins = parse_json_int_list(feature_row.get("observed_ion_bins_json", "[]")) if not feature_row.empty else []
        background = background_ion_matches(ion_bins, 209)
        for source in transform_sources_for_control(control, unknown_row, feature_row):
            transform_name = str(source["transform_name"])
            tracks = [TRACK_ANNOTATION_ASSISTED] if "annotation" in str(source["evidence_source"]) else [TRACK_CONTROL_BLIND, TRACK_ANNOTATION_ASSISTED]
            for track in tracks:
                if not evidence_track_allowed(str(source["evidence_source"]), track):
                    continue
                for mass_hyp in mass_hypotheses[: int(args.max_adduct_hypotheses)]:
                    neutral_mass = safe_float(mass_hyp.get("neutral_mass"), math.nan)
                    if transform_name == "peptide":
                        peptide = peptide_compositions_for_mass(neutral_mass, tolerance_da=float(args.peptide_tolerance_da), max_len=int(args.max_peptide_length))
                        if peptide.empty:
                            rows.append(
                                {
                                    "track": track,
                                    "control_id": xid,
                                    "feature_key": feature_key,
                                    "lib_id": lib_id,
                                    "ion_mode": ion_mode,
                                    "parent_mz": parent_mz,
                                    **mass_hyp,
                                    "transform_name": transform_name,
                                    "transform_delta_mass": 0.0,
                                    "core_mass": neutral_mass,
                                    "target_formula": formula,
                                    "core_formula": "",
                                    "peptide_composition": "",
                                    "evidence_source": source["evidence_source"],
                                    "evidence_confidence": float(source["base_score"]),
                                    "raw_evidence_text": source["raw_text"],
                                    **background,
                                }
                            )
                        else:
                            for _, pep in peptide.iterrows():
                                rows.append(
                                    {
                                        "track": track,
                                        "control_id": xid,
                                        "feature_key": feature_key,
                                        "lib_id": lib_id,
                                        "ion_mode": ion_mode,
                                        "parent_mz": parent_mz,
                                        **mass_hyp,
                                        "transform_name": transform_name,
                                        "transform_delta_mass": 0.0,
                                        "core_mass": neutral_mass,
                                        "target_formula": formula,
                                        "core_formula": "",
                                        "peptide_composition": pep.get("peptide_composition", ""),
                                        "peptide_mass_error_da": pep.get("peptide_mass_error_da", ""),
                                        "leu_ile_ambiguous": pep.get("leu_ile_ambiguous", False),
                                        "evidence_source": source["evidence_source"],
                                        "evidence_confidence": float(source["base_score"]),
                                        "raw_evidence_text": source["raw_text"],
                                        **background,
                                    }
                                )
                        continue
                    spec = TRANSFORMS.get(transform_name)
                    if spec is None:
                        continue
                    core_mass = neutral_mass - spec.delta_mass if np.isfinite(neutral_mass) else math.nan
                    core_formula = parse_core_formula_from_text(source["raw_text"], transform_name)
                    if not core_formula and formula:
                        core_formula = formula_subtract(formula, spec.delta_formula)
                    rows.append(
                        {
                            "track": track,
                            "control_id": xid,
                            "feature_key": feature_key,
                            "lib_id": lib_id,
                            "ion_mode": ion_mode,
                            "parent_mz": parent_mz,
                            **mass_hyp,
                            "transform_name": transform_name,
                            "transform_delta_mass": spec.delta_mass,
                            "transform_delta_formula": spec.delta_formula,
                            "core_mass": core_mass,
                            "target_formula": formula,
                            "core_formula": core_formula,
                            "core_formula_mass_error_da": abs(formula_mass(core_formula) - core_mass) if core_formula and np.isfinite(core_mass) else math.nan,
                            "peptide_composition": "",
                            "evidence_source": source["evidence_source"],
                            "evidence_confidence": float(source["base_score"]),
                            "raw_evidence_text": source["raw_text"],
                            **background,
                        }
                    )
    evidence = pd.DataFrame(rows)
    if not evidence.empty:
        evidence = evidence.sort_values(
            ["track", "control_id", "evidence_confidence", "adduct_priority", "transform_name"],
            ascending=[True, True, False, False, True],
            kind="mergesort",
        ).drop_duplicates(["track", "control_id", "mass_hypothesis_adduct", "transform_name", "core_mass", "peptide_composition"], keep="first")
    summary = summarize_transform_counts(evidence)
    write_csv(output_root / "biotransformation_evidence.csv", evidence)
    write_csv(output_root / "biotransformation_evidence_summary.csv", summary)
    report = "\n".join(
        [
            "# Biotransformation Evidence V1",
            "",
            "Evidence is split into control-blind and annotation-assisted tracks.",
            "",
            "## Transform Counts",
            "",
            dataframe_to_markdown(summary),
            "",
        ]
    )
    (output_root / "biotransformation_evidence_report.md").write_text(report)
    payload = {
        "output_root": str(output_root),
        "libs": sorted(libs),
        "ion_mode": str(args.ion_mode),
        "n_verified_controls": int(targets.shape[0]),
        "n_evidence_rows": int(evidence.shape[0]),
        "n_control_blind_controls": int(evidence[evidence["track"].eq(TRACK_CONTROL_BLIND)]["control_id"].nunique()) if not evidence.empty else 0,
        "n_annotation_assisted_controls": int(evidence[evidence["track"].eq(TRACK_ANNOTATION_ASSISTED)]["control_id"].nunique()) if not evidence.empty else 0,
    }
    write_json(output_root / "biotransformation_evidence_summary.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--known-x", type=Path, default=Path("known_X.csv"))
    parser.add_argument("--known-x-all", type=Path, default=Path("known_X_all_available_structures_simple.csv"))
    parser.add_argument("--unknown-list", type=Path, default=DEFAULT_UNKNOWN_LIST)
    parser.add_argument("--feature-evidence-paths", default=DEFAULT_FEATURE_EVIDENCE)
    parser.add_argument("--libs", default="209")
    parser.add_argument("--ion-mode", default="negative")
    parser.add_argument("--target-controls-csv", type=Path, default=Path(""))
    parser.add_argument("--include-unverified-x", action="store_true")
    parser.add_argument("--max-adduct-hypotheses", type=int, default=4)
    parser.add_argument("--peptide-tolerance-da", type=float, default=0.05)
    parser.add_argument("--max-peptide-length", type=int, default=4)
    build_evidence(parser.parse_args())


if __name__ == "__main__":
    main()
