#!/usr/bin/env python3
"""Build upstream adduct and neutral-mass hypotheses for unknown features."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from nce_evidence_utils import PROTON_MASS, parse_json_list
from run_retrieval_anchored_denovo_assembly_v4 import resolve
from x_domain_analog_utils import dataframe_to_markdown, safe_float, safe_int, write_csv, write_json


DEFAULT_OUTPUT_ROOT = Path("models/upstream_adduct_evidence_v1")


@dataclass(frozen=True)
class AdductSpec:
    label: str
    polarity: str
    charge: int
    multiplier: int
    offset_mass: float
    prevalence: float
    source: str
    neutral_loss_possible: bool = False
    note: str = ""


ADDUCTS: tuple[AdductSpec, ...] = (
    AdductSpec("[M-H]-", "negative", -1, 1, -PROTON_MASS, 0.174013748, "msac_default_acd_common", False),
    AdductSpec("[M+Cl]-", "negative", -1, 1, 34.969402, 0.005324162, "msac_default_acd_common", False, "chloride isotope hint expected"),
    AdductSpec("[M+HCOO]-", "negative", -1, 1, 44.998201, 0.002, "acd_common", False, "formate"),
    AdductSpec("[M+CH3COO]-", "negative", -1, 1, 59.013851, 0.000781929, "msac_default_acd_common", False, "acetate"),
    AdductSpec("[M+Br]-", "negative", -1, 1, 78.918885, 0.0002, "acd_common", False, "bromide isotope hint expected"),
    AdductSpec("[M+CF3COO]-", "negative", -1, 1, 112.985586, 0.0001, "acd_common", False, "trifluoroacetate"),
    AdductSpec("[M-2H]2-", "negative", -2, 1, -2.0 * PROTON_MASS, 0.01235448, "msac_default", False),
    AdductSpec("[2M-H]-", "negative", -1, 2, -PROTON_MASS, 0.012588354, "msac_default", False),
    AdductSpec("[M+H]+", "positive", 1, 1, PROTON_MASS, 0.581491808, "msac_default_acd_common", False),
    AdductSpec("[M+NH4]+", "positive", 1, 1, 18.03382, 0.015676622, "msac_default_acd_common", False),
    AdductSpec("[M+Na]+", "positive", 1, 1, 22.989218, 0.05325994, "msac_default_acd_common", False),
    AdductSpec("[M+CH3OH+H]+", "positive", 1, 1, 33.033489, 0.0005, "acd_common", False),
    AdductSpec("[M+K]+", "positive", 1, 1, 38.9632, 0.009104191, "msac_default_acd_common", False),
    AdductSpec("[M+CH3CN+H]+", "positive", 1, 1, 42.003823, 0.0005, "acd_common", False),
    AdductSpec("[M+2H]2+", "positive", 2, 1, 2.0 * PROTON_MASS, 0.043875382, "msac_default_acd_common", False),
    AdductSpec("[M+H+Na]2+", "positive", 2, 1, PROTON_MASS + 22.989218, 0.002030198, "msac_default_acd_common", False),
    AdductSpec("[2M+H]+", "positive", 1, 2, PROTON_MASS, 0.01935662, "msac_default", False),
)


def parse_csv(text: Any) -> list[str]:
    return [part.strip() for part in str(text or "").split(",") if part.strip()]


def neutral_mass_from_adduct(parent_mz: float, spec: AdductSpec) -> float:
    if not np.isfinite(float(parent_mz)) or int(spec.multiplier) <= 0:
        return float("nan")
    neutral = (float(parent_mz) * abs(int(spec.charge)) - float(spec.offset_mass)) / float(spec.multiplier)
    return float(neutral) if neutral > 0 else float("nan")


def adduct_mz_from_neutral(neutral_mass: float, spec: AdductSpec) -> float:
    if not np.isfinite(float(neutral_mass)) or int(spec.charge) == 0:
        return float("nan")
    return float((float(spec.multiplier) * float(neutral_mass) + float(spec.offset_mass)) / abs(int(spec.charge)))


def feature_key(row: pd.Series) -> str:
    for col in ["feature_key", "feature_id", "feature_label"]:
        value = str(row.get(col, "") or "")
        if value and value.lower() != "nan":
            return value
    comp_id = safe_int(row.get("comp_id"), -1)
    lib_id = safe_int(row.get("lib_id"), 0)
    mz = safe_float(row.get("parent_mz", row.get("mz")), math.nan)
    return f"feature_lib{lib_id}_comp{comp_id}_mz{mz:.4f}"


def ion_mode_for_lib(lib_id: int) -> str:
    return "negative" if int(lib_id) in {209, 305} else "positive"


def observed_bins(row: pd.Series, column: str) -> list[int]:
    return [safe_int(value, 0) for value in parse_json_list(row.get(column, "[]")) if safe_int(value, 0)]


def support_score(row: pd.Series, spec: AdductSpec) -> float:
    score = 0.0
    hint = str(row.get("adduct_hypothesis", row.get("mass_hypothesis_adduct", "")) or "")
    if hint and spec.label.replace("[", "").replace("]", "") in hint.replace("[", "").replace("]", ""):
        score += 0.35
    ion_bins = set(observed_bins(row, "observed_ion_bins_json"))
    loss_bins = set(observed_bins(row, "observed_loss_bins_json"))
    note = spec.note.lower()
    if "formate" in note and (45000 in ion_bins or 45000 in loss_bins):
        score += 0.15
    if "acetate" in note and (59014 in ion_bins or 59014 in loss_bins):
        score += 0.15
    if "chloride" in note:
        score += 0.05
    if "bromide" in note:
        score += 0.05
    return float(min(score, 0.5))


def score_adduct(row: pd.Series, spec: AdductSpec) -> float:
    primary = {"[M-H]-", "[M+H]+"}
    prevalence = min(max(float(spec.prevalence), 0.0), 1.0)
    prevalence_score = math.log10(1.0 + 1000.0 * prevalence) / 3.0
    primary_bonus = 0.25 if spec.label in primary else 0.0
    dimer_penalty = -0.10 if int(spec.multiplier) > 1 else 0.0
    charge_penalty = -0.05 if abs(int(spec.charge)) > 1 else 0.0
    return float(np.clip(0.45 * prevalence_score + primary_bonus + support_score(row, spec) + dimer_penalty + charge_penalty, 0.0, 1.0))


def build_rows(features: pd.DataFrame, libs: set[int], ion_mode: str = "") -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    default_lib_id = next(iter(libs)) if len(libs) == 1 else 0
    ion_mode_override = str(ion_mode or "").strip().lower()
    for _, row in features.iterrows():
        lib_id = safe_int(row.get("lib_id"), default_lib_id)
        if libs and lib_id not in libs:
            continue
        parent_mz = safe_float(row.get("parent_mz", row.get("mz")), math.nan)
        if not np.isfinite(parent_mz):
            continue
        polarity = ion_mode_override or str(row.get("ion_mode", "") or "").lower() or ion_mode_for_lib(lib_id)
        key = feature_key(row)
        source_neutral = safe_float(row.get("source_neutral_mass"), math.nan)
        if np.isfinite(source_neutral) and 50 <= source_neutral <= 2500:
            source_adduct = str(row.get("source_parent_adduct", "") or "chrolib_source_adduct")
            rows.append(
                {
                    "feature_key": key,
                    "comp_id": safe_int(row.get("comp_id"), -1),
                    "lib_id": lib_id,
                    "ion_mode": polarity,
                    "parent_mz": parent_mz,
                    "adduct_label": f"CHRO_LIB:{source_adduct}",
                    "charge": "",
                    "multiplier": "",
                    "offset_mass": "",
                    "neutral_mass": float(source_neutral),
                    "recomputed_mz": parent_mz,
                    "adduct_score": 1.0,
                    "prevalence": 1.0,
                    "prevalence_tier": "source",
                    "source_reference": "chrolib_source_parent_adduct",
                    "neutral_loss_possible": False,
                    "note": "source_parent_adduct_neutral_mass",
                }
            )
        for spec in ADDUCTS:
            if spec.polarity != polarity:
                continue
            neutral = neutral_mass_from_adduct(parent_mz, spec)
            if not np.isfinite(neutral) or neutral < 50 or neutral > 2500:
                continue
            score = score_adduct(row, spec)
            rows.append(
                {
                    "feature_key": key,
                    "comp_id": safe_int(row.get("comp_id"), -1),
                    "lib_id": lib_id,
                    "ion_mode": polarity,
                    "parent_mz": parent_mz,
                    "adduct_label": spec.label,
                    "charge": spec.charge,
                    "multiplier": spec.multiplier,
                    "offset_mass": spec.offset_mass,
                    "neutral_mass": neutral,
                    "recomputed_mz": adduct_mz_from_neutral(neutral, spec),
                    "adduct_score": score,
                    "prevalence": spec.prevalence,
                    "prevalence_tier": "primary" if spec.label in {"[M-H]-", "[M+H]+"} else "common" if spec.prevalence >= 0.005 else "extended",
                    "source_reference": spec.source,
                    "neutral_loss_possible": bool(spec.neutral_loss_possible),
                    "note": spec.note,
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.sort_values(["feature_key", "adduct_score", "prevalence", "adduct_label"], ascending=[True, False, False, True], kind="mergesort")
    out = out.drop_duplicates(["feature_key", "adduct_label", "neutral_mass"], keep="first").copy()
    out["adduct_rank"] = out.groupby("feature_key").cumcount() + 1
    return out


def render_report(hypotheses: pd.DataFrame, conflicts: pd.DataFrame) -> str:
    summary = (
        hypotheses[hypotheses["adduct_rank"].le(1)].groupby(["lib_id", "ion_mode", "adduct_label"]).size().reset_index(name="n_top_features")
        if not hypotheses.empty
        else pd.DataFrame()
    )
    return "\n".join(
        [
            "# Upstream Adduct Evidence V1",
            "",
            "Adduct hypotheses use encoded MSAC-style default/common adducts and ACD common adduct masses. They do not use X truth structures or formulas.",
            "",
            "## Top Adduct Summary",
            "",
            dataframe_to_markdown(summary.sort_values(["lib_id", "n_top_features"], ascending=[True, False]) if not summary.empty else summary),
            "",
            "## Conflict Audit",
            "",
            dataframe_to_markdown(conflicts.head(30)),
            "",
        ]
    )


def run(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parent
    output_root = resolve(args.output_root, repo_root)
    output_root.mkdir(parents=True, exist_ok=True)
    features = pd.read_csv(resolve(args.feature_table, repo_root), low_memory=False)
    libs = {int(value) for value in parse_csv(args.libs)}
    hypotheses = build_rows(features, libs, ion_mode=str(getattr(args, "ion_mode", "")))
    top = hypotheses[hypotheses["adduct_rank"].le(int(args.top_n))].copy() if not hypotheses.empty else pd.DataFrame()
    conflicts = pd.DataFrame()
    if not top.empty:
        span = top.groupby("feature_key")["neutral_mass"].agg(["min", "max", "count"]).reset_index()
        span["neutral_mass_span_da"] = span["max"] - span["min"]
        conflicts = span[span["neutral_mass_span_da"].gt(float(args.conflict_span_da))].copy()
    write_csv(output_root / "feature_adduct_hypotheses.csv", hypotheses)
    write_csv(output_root / "neutral_mass_hypotheses.csv", top)
    write_csv(output_root / "adduct_conflict_audit.csv", conflicts)
    (output_root / "adduct_evidence_report.md").write_text(render_report(hypotheses, conflicts))
    payload = {
        "output_root": str(output_root),
        "n_feature_rows": int(features.shape[0]),
        "n_hypothesis_rows": int(hypotheses.shape[0]),
        "n_top_hypothesis_rows": int(top.shape[0]),
        "adduct_source": str(args.adduct_source),
    }
    write_json(output_root / "adduct_evidence_summary.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--feature-table", type=Path, default=Path("merged_encoder/data/ions_unknown_merged.csv"))
    parser.add_argument("--libs", default="209")
    parser.add_argument("--ion-mode", default="")
    parser.add_argument("--adduct-source", default="msac_acd_common")
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--conflict-span-da", type=float, default=50.0)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
