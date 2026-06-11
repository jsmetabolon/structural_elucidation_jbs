#!/usr/bin/env python3
"""Build v21 neutral-mass hypotheses from observed CHRO adduct ladders."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from build_upstream_adduct_evidence_v1 import build_rows as build_fallback_rows
from run_retrieval_anchored_denovo_assembly_v4 import resolve
from x_domain_analog_utils import dataframe_to_markdown, safe_float, safe_int, write_csv, write_json


DEFAULT_OUTPUT_ROOT = Path("models/upstream_adduct_evidence_v21")
DEFAULT_FEATURE_TABLE = Path("merged_encoder/data/ions_unknown_merged_v21.csv")
DEFAULT_LADDER = Path("merged_encoder/data/chrolib_adduct_ladder_v21.csv")
HYPOTHESIS_COLUMNS = [
    "feature_key",
    "comp_id",
    "lib_id",
    "ion_mode",
    "parent_mz",
    "adduct_label",
    "charge",
    "multiplier",
    "offset_mass",
    "neutral_mass",
    "recomputed_mz",
    "adduct_score",
    "prevalence",
    "prevalence_tier",
    "source_reference",
    "neutral_loss_possible",
    "note",
    "adduct_rank",
]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--feature-table", type=Path, default=DEFAULT_FEATURE_TABLE)
    parser.add_argument("--chrolib-adduct-ladder", type=Path, default=DEFAULT_LADDER)
    parser.add_argument("--libs", default="", help="Optional comma-separated library IDs.")
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--conflict-span-da", type=float, default=50.0)
    parser.add_argument("--force", action="store_true", help="Allow overwriting existing v21 output directory files.")
    return parser.parse_args(argv)


def parse_csv(text: Any) -> list[str]:
    return [part.strip() for part in str(text or "").split(",") if part.strip()]


def normalize_ion_mode(value: Any, lib_id: int) -> str:
    text = str(value or "").strip().lower()
    if text in {"0", "0.0", "negative", "neg"}:
        return "negative"
    if text in {"1", "1.0", "positive", "pos"}:
        return "positive"
    return "negative" if int(lib_id) in {209, 305} else "positive"


def xid_from_feature_key(value: Any) -> str:
    match = re.search(r"(X-\d+)", str(value or ""))
    return match.group(1) if match else ""


def feature_rows_from_peaks(features: pd.DataFrame, libs: set[int]) -> pd.DataFrame:
    if features.empty:
        return pd.DataFrame()
    work = features.copy()
    work["lib_id"] = pd.to_numeric(work.get("lib_id"), errors="coerce")
    work["comp_id"] = pd.to_numeric(work.get("comp_id"), errors="coerce")
    work["parent_mz"] = pd.to_numeric(work.get("parent_mz"), errors="coerce")
    if libs:
        work = work[work["lib_id"].isin(libs)].copy()
    if work.empty:
        return pd.DataFrame()
    if "feature_key" not in work.columns:
        raise KeyError("feature table must contain feature_key")
    if "ion_mode_label" in work.columns:
        work["_ion_mode_source"] = work["ion_mode_label"]
    else:
        work["_ion_mode_source"] = work.get("ion_mode", "")
    agg_spec: dict[str, tuple[str, str]] = {
        "comp_id": ("comp_id", "first"),
        "lib_id": ("lib_id", "first"),
        "parent_mz": ("parent_mz", "first"),
        "ion_mode_source": ("_ion_mode_source", "first"),
    }
    if "source_neutral_mass" in work.columns:
        agg_spec["source_neutral_mass"] = ("source_neutral_mass", "first")
    if "source_parent_adduct" in work.columns:
        agg_spec["source_parent_adduct"] = ("source_parent_adduct", "first")
    out = work.groupby("feature_key", dropna=False).agg(**agg_spec).reset_index()
    out["comp_id"] = pd.to_numeric(out["comp_id"], errors="coerce").astype("Int64")
    out["lib_id"] = pd.to_numeric(out["lib_id"], errors="coerce").astype("Int64")
    out["ion_mode"] = [
        normalize_ion_mode(mode, int(lib_id)) for mode, lib_id in zip(out["ion_mode_source"], out["lib_id"])
    ]
    out["control_id"] = out["feature_key"].map(xid_from_feature_key)
    return out.drop(columns=["ion_mode_source"])


def observed_score(adduct_class: str, adduct_name: str) -> float:
    cls = str(adduct_class or "")
    if cls == "primary":
        return 0.97
    if cls in {"alkali", "ammonium"}:
        return 0.92
    if cls in {"dimer", "isotopologue"}:
        return 0.88
    if cls == "solvent":
        return 0.82
    return 0.78


def observed_rows_for_feature(feature: pd.Series, observed: pd.DataFrame) -> list[dict[str, Any]]:
    valid = observed.copy()
    valid["implied_neutral_mass"] = pd.to_numeric(valid.get("implied_neutral_mass"), errors="coerce")
    valid = valid[
        ~valid.get("is_neutral_loss", pd.Series(False, index=valid.index)).astype(bool)
        & valid["implied_neutral_mass"].between(50, 2500, inclusive="both")
    ].copy()
    if valid.empty:
        return []

    masses = valid["implied_neutral_mass"].astype(float).to_numpy()
    consensus = float(np.median(masses))
    labels = ",".join(sorted({str(v) for v in valid["adduct_name"].astype(str) if str(v)}))
    base = {
        "feature_key": str(feature.get("feature_key")),
        "comp_id": safe_int(feature.get("comp_id"), -1),
        "lib_id": safe_int(feature.get("lib_id"), 0),
        "ion_mode": str(feature.get("ion_mode")),
        "parent_mz": safe_float(feature.get("parent_mz"), math.nan),
        "charge": "",
        "multiplier": "",
        "offset_mass": "",
        "recomputed_mz": safe_float(feature.get("parent_mz"), math.nan),
        "prevalence": 1.0,
        "prevalence_tier": "observed",
        "source_reference": "chrolib_observed_adduct_ladder",
        "neutral_loss_possible": False,
    }
    rows: list[dict[str, Any]] = [
        {
            **base,
            "adduct_label": "CHRO_LIB:observed_consensus",
            "neutral_mass": consensus,
            "adduct_score": 1.0,
            "note": f"consensus_observed_adducts={labels}",
        }
    ]
    valid = valid.sort_values(["adduct_class", "adduct_name"], kind="mergesort")
    for _, row in valid.drop_duplicates(["adduct_name", "implied_neutral_mass"]).iterrows():
        rows.append(
            {
                **base,
                "adduct_label": f"CHRO_LIB:{row.get('adduct_name')}",
                "charge": safe_int(row.get("charge"), 0) or "",
                "multiplier": safe_int(row.get("multiplier"), 0) or "",
                "neutral_mass": safe_float(row.get("implied_neutral_mass"), math.nan),
                "adduct_score": observed_score(str(row.get("adduct_class")), str(row.get("adduct_name"))),
                "note": f"observed_adduct_class={row.get('adduct_class')}",
            }
        )
    return rows


def build_observed_hypotheses(features: pd.DataFrame, ladder: pd.DataFrame) -> tuple[pd.DataFrame, set[str]]:
    if features.empty or ladder.empty:
        return pd.DataFrame(columns=HYPOTHESIS_COLUMNS), set()
    ladder = ladder.copy()
    ladder["lib_id"] = pd.to_numeric(ladder.get("lib_id"), errors="coerce").astype("Int64")
    ladder["comp_id"] = pd.to_numeric(ladder.get("comp_id"), errors="coerce").astype("Int64")
    ladder_map = {
        (int(lib_id), int(comp_id)): group.copy()
        for (lib_id, comp_id), group in ladder.dropna(subset=["lib_id", "comp_id"]).groupby(["lib_id", "comp_id"], sort=False)
    }
    rows: list[dict[str, Any]] = []
    observed_feature_keys: set[str] = set()
    for _, feature in features.iterrows():
        key = (safe_int(feature.get("lib_id"), 0), safe_int(feature.get("comp_id"), -1))
        observed = ladder_map.get(key)
        if observed is None or observed.empty:
            continue
        feature_rows = observed_rows_for_feature(feature, observed)
        if not feature_rows:
            continue
        observed_feature_keys.add(str(feature.get("feature_key")))
        rows.extend(feature_rows)
    out = pd.DataFrame(rows, columns=[col for col in HYPOTHESIS_COLUMNS if col != "adduct_rank"])
    return out, observed_feature_keys


def loss_label(tokens: list[str]) -> str:
    upper = [token.upper() for token in tokens if token]
    labels: list[str] = []
    if upper.count("H2O") > 1:
        labels.append("polyol/sugar")
    elif "H2O" in upper:
        labels.append("hydroxyl/dehydration")
    if "HCOOH" in upper or "CO2" in upper:
        labels.append("carboxyl")
    if "NH3" in upper:
        labels.append("amine/amide")
    return "|".join(labels or ["unknown_loss"])


def build_neutral_loss_hints(features: pd.DataFrame, ladder: pd.DataFrame) -> pd.DataFrame:
    if features.empty or ladder.empty:
        return pd.DataFrame()
    ladder = ladder.copy()
    ladder["lib_id"] = pd.to_numeric(ladder.get("lib_id"), errors="coerce").astype("Int64")
    ladder["comp_id"] = pd.to_numeric(ladder.get("comp_id"), errors="coerce").astype("Int64")
    losses = ladder[ladder.get("is_neutral_loss", pd.Series(False, index=ladder.index)).astype(bool)].copy()
    if losses.empty:
        return pd.DataFrame()
    feature_map = {
        (safe_int(row.get("lib_id"), 0), safe_int(row.get("comp_id"), -1)): row
        for _, row in features.iterrows()
    }
    rows: list[dict[str, Any]] = []
    for _, loss in losses.iterrows():
        key = (safe_int(loss.get("lib_id"), 0), safe_int(loss.get("comp_id"), -1))
        feature = feature_map.get(key)
        if feature is None:
            continue
        tokens = [token for token in str(loss.get("neutral_loss_tokens", "") or "").split("|") if token]
        rows.append(
            {
                "feature_key": feature.get("feature_key"),
                "control_id": feature.get("control_id", ""),
                "comp_id": safe_int(feature.get("comp_id"), -1),
                "lib_id": safe_int(feature.get("lib_id"), 0),
                "ion_mode": feature.get("ion_mode", ""),
                "parent_mz": safe_float(feature.get("parent_mz"), math.nan),
                "adduct_name": loss.get("adduct_name", ""),
                "neutral_loss_tokens": "|".join(tokens),
                "functional_group_hint": loss_label(tokens),
                "hint_source": "chrolib_observed_neutral_loss",
            }
        )
    return pd.DataFrame(rows)


def rank_hypotheses(hypotheses: pd.DataFrame) -> pd.DataFrame:
    if hypotheses.empty:
        return pd.DataFrame(columns=HYPOTHESIS_COLUMNS)
    out = hypotheses.copy()
    for col in HYPOTHESIS_COLUMNS:
        if col not in out.columns:
            out[col] = ""
    out["neutral_mass"] = pd.to_numeric(out["neutral_mass"], errors="coerce")
    out["adduct_score"] = pd.to_numeric(out["adduct_score"], errors="coerce").fillna(0.0)
    out["prevalence"] = pd.to_numeric(out["prevalence"], errors="coerce").fillna(0.0)
    out = out[out["neutral_mass"].between(50, 2500, inclusive="both")].copy()
    out = out.sort_values(
        ["feature_key", "adduct_score", "prevalence", "adduct_label"],
        ascending=[True, False, False, True],
        kind="mergesort",
    )
    out = out.drop_duplicates(["feature_key", "adduct_label", "neutral_mass"], keep="first")
    out["adduct_rank"] = out.groupby("feature_key", dropna=False).cumcount() + 1
    return out[HYPOTHESIS_COLUMNS].reset_index(drop=True)


def render_report(hypotheses: pd.DataFrame, conflicts: pd.DataFrame, hints: pd.DataFrame) -> str:
    summary = (
        hypotheses[hypotheses["adduct_rank"].le(1)].groupby(["lib_id", "ion_mode", "adduct_label"]).size().reset_index(name="n_top_features")
        if not hypotheses.empty
        else pd.DataFrame()
    )
    return "\n".join(
        [
            "# Upstream Adduct Evidence V3",
            "",
            "Observed CHRO adduct ladders provide rank-1 neutral-mass hypotheses when available; the v1 inferred adduct prior is fallback only.",
            "",
            "## Top Adduct Summary",
            "",
            dataframe_to_markdown(summary.sort_values(["lib_id", "n_top_features"], ascending=[True, False]) if not summary.empty else summary),
            "",
            "## Conflict Audit",
            "",
            dataframe_to_markdown(conflicts.head(30)),
            "",
            "## Neutral Loss Hints",
            "",
            dataframe_to_markdown(hints.head(30)),
            "",
        ]
    )


def run(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parent
    output_root = resolve(args.output_root, repo_root)
    if output_root.exists() and not args.force:
        existing = list(output_root.glob("*"))
        if existing:
            raise FileExistsError(f"Refusing to overwrite non-empty output root without --force: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)

    raw_features = pd.read_csv(resolve(args.feature_table, repo_root), low_memory=False)
    libs = {int(value) for value in parse_csv(args.libs)}
    features = feature_rows_from_peaks(raw_features, libs)
    ladder_path = resolve(args.chrolib_adduct_ladder, repo_root)
    ladder = pd.read_csv(ladder_path, low_memory=False) if ladder_path.exists() else pd.DataFrame()
    if libs and not ladder.empty:
        ladder = ladder[pd.to_numeric(ladder.get("lib_id"), errors="coerce").isin(libs)].copy()

    observed, observed_feature_keys = build_observed_hypotheses(features, ladder)
    fallback_features = features[~features["feature_key"].astype(str).isin(observed_feature_keys)].copy()
    fallback = build_fallback_rows(fallback_features, libs=libs, ion_mode="") if not fallback_features.empty else pd.DataFrame()
    hypotheses = rank_hypotheses(pd.concat([observed, fallback], ignore_index=True, sort=False))
    top = hypotheses[hypotheses["adduct_rank"].le(int(args.top_n))].copy() if not hypotheses.empty else pd.DataFrame(columns=HYPOTHESIS_COLUMNS)

    conflicts = pd.DataFrame()
    if not top.empty:
        span = top.groupby("feature_key")["neutral_mass"].agg(["min", "max", "count"]).reset_index()
        span["neutral_mass_span_da"] = span["max"] - span["min"]
        conflicts = span[span["neutral_mass_span_da"].gt(float(args.conflict_span_da))].copy()
    hints = build_neutral_loss_hints(features, ladder)

    write_csv(output_root / "feature_adduct_hypotheses.csv", hypotheses)
    write_csv(output_root / "neutral_mass_hypotheses.csv", top)
    write_csv(output_root / "adduct_conflict_audit.csv", conflicts)
    write_csv(output_root / "neutral_loss_hints_v21.csv", hints)
    (output_root / "adduct_evidence_report.md").write_text(render_report(hypotheses, conflicts, hints))
    payload = {
        "output_root": str(output_root),
        "feature_table": str(resolve(args.feature_table, repo_root)),
        "chrolib_adduct_ladder": str(ladder_path),
        "n_feature_rows": int(features.shape[0]),
        "n_observed_feature_rows": int(len(observed_feature_keys)),
        "n_fallback_feature_rows": int(fallback_features.shape[0]),
        "n_hypothesis_rows": int(hypotheses.shape[0]),
        "n_top_hypothesis_rows": int(top.shape[0]),
        "n_neutral_loss_hints": int(hints.shape[0]),
        "adduct_source": "chrolib_observed_ladder_with_v1_fallback",
    }
    write_json(output_root / "adduct_evidence_summary.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
