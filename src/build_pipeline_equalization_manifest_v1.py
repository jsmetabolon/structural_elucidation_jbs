#!/usr/bin/env python3
"""Build full-pipeline equalization manifests for actionable X controls."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import pandas as pd

from build_observed_ms2_lookup_v1 import load_feature_rows_from_peak_table
from run_retrieval_anchored_denovo_assembly_v4 import resolve
from x209_pkl_refresh_utils import load_validated_ions_pkl
from x_domain_analog_utils import write_csv, write_json


LIB_SCOPES = ((209, "negative"), (305, "negative"), (400, "positive"), (402, "positive"))
DEFAULT_OUTPUT_ROOT = Path("models/pipeline_equalization_v1")
DEFAULT_OBSERVED_209 = Path("models/observed_ms2_lookup_v8/per_control_all_x.csv")
DEFAULT_FEATURE_209 = Path("models/peptide_first_denovo_v1_1/feature_peptide_rows.csv")
DEFAULT_MASS_HYPOTHESES = Path("models/upstream_adduct_evidence_v1/neutral_mass_hypotheses.csv")
DEFAULT_MERGED_IONS = Path("merged_encoder/data/ions_unknown_merged.csv")


def clean_text(value: Any) -> str:
    text = "" if value is None else str(value)
    return "" if text.lower() == "nan" else text.strip()


def safe_read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes", "y"})


def actionable_known(known_x: pd.DataFrame, lib_id: int) -> pd.DataFrame:
    scoped = known_x[pd.to_numeric(known_x.get("LIB_ID"), errors="coerce").eq(int(lib_id))].copy()
    if scoped.empty:
        return scoped
    scoped["control_id"] = scoped.get("Unknown", pd.Series("", index=scoped.index)).map(clean_text)
    scoped["comp_id"] = pd.to_numeric(scoped.get("CompID"), errors="coerce")
    text_cols = [
        "Level",
        "SubLevel",
        "KNOWN",
        "External Information (Class or modification)",
        "Tentative information",
        "Cleaned External Info",
    ]
    text = pd.Series("", index=scoped.index)
    for col in text_cols:
        if col in scoped.columns:
            text = text.str.cat(scoped[col].fillna("").astype(str), sep=" ")
    scoped = scoped[scoped["control_id"].str.startswith("X-")].copy()
    scoped = scoped[~text.str.contains("DNI|Artifact", case=False, na=False)].copy()
    scoped["level"] = scoped.get("Level", pd.Series("", index=scoped.index)).map(clean_text)
    scoped["sublevel"] = scoped.get("SubLevel", pd.Series("", index=scoped.index)).map(clean_text)
    scoped["clean_external_info"] = scoped.get("Cleaned External Info", pd.Series("", index=scoped.index)).map(clean_text)
    keep = ["control_id", "comp_id", "level", "sublevel", "clean_external_info"]
    return scoped[keep].drop_duplicates("control_id", keep="first").reset_index(drop=True)


def observed_209(path: Path, lib_id: int, ion_mode: str) -> pd.DataFrame:
    raw = safe_read_csv(path)
    if raw.empty:
        return raw
    scoped = raw[
        pd.to_numeric(raw.get("lib_id"), errors="coerce").eq(int(lib_id))
        & raw.get("ion_mode", pd.Series("", index=raw.index)).astype(str).eq(str(ion_mode))
    ].copy()
    if "ms2_present" in scoped.columns:
        scoped = scoped[bool_series(scoped["ms2_present"])].copy()
    return scoped.reset_index(drop=True)


def observed_from_pkl(path: Path, lib_id: int, ion_mode: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    raw = load_validated_ions_pkl(path)
    return load_feature_rows_from_peak_table(raw, lib_id=int(lib_id), ion_mode=str(ion_mode))


def observed_from_merged(path: Path, lib_id: int, ion_mode: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    peaks = safe_read_csv(path)
    if peaks.empty:
        return pd.DataFrame()
    scoped = peaks[
        pd.to_numeric(peaks.get("lib_id"), errors="coerce").eq(int(lib_id))
        & peaks.get("ion_mode_label", pd.Series("", index=peaks.index)).astype(str).eq(str(ion_mode))
    ].copy()
    return load_feature_rows_from_peak_table(scoped, lib_id=int(lib_id), ion_mode=str(ion_mode))


def combine_observed(*frames: pd.DataFrame) -> pd.DataFrame:
    usable = [frame for frame in frames if frame is not None and not frame.empty]
    if not usable:
        return pd.DataFrame()
    combined = pd.concat(usable, ignore_index=True, sort=False)
    if combined.empty:
        return combined
    sort_cols = [col for col in ("control_id", "feature_key", "n_ms2_pairs") if col in combined.columns]
    if sort_cols:
        combined = combined.sort_values(sort_cols, ascending=[True] * len(sort_cols), kind="mergesort")
    dedupe_cols = ["control_id"]
    if "feature_key" in combined.columns:
        dedupe_cols.append("feature_key")
    return combined.drop_duplicates(dedupe_cols, keep="last").reset_index(drop=True)


def feature_rows_from_observed(observed: pd.DataFrame, lib_id: int, ion_mode: str) -> pd.DataFrame:
    if observed.empty:
        return pd.DataFrame(
            columns=[
                "feature_key",
                "comp_id",
                "control_id",
                "lib_id",
                "ion_mode",
                "parent_mz",
                "neutral_mass_fallback",
                "n_ms2_peaks",
                "observed_ion_bins_json",
                "observed_loss_bins_json",
            ]
        )
    out = observed.copy()
    if "neutral_mass_fallback" not in out.columns:
        out["neutral_mass_fallback"] = pd.to_numeric(out.get("neutral_mass"), errors="coerce")
    if "n_ms2_peaks" not in out.columns:
        out["n_ms2_peaks"] = pd.to_numeric(out.get("n_ms2_pairs"), errors="coerce").fillna(0).astype(int)
    out["lib_id"] = int(lib_id)
    out["ion_mode"] = str(ion_mode)
    keep = [
        "feature_key",
        "comp_id",
        "control_id",
        "lib_id",
        "ion_mode",
        "parent_mz",
        "neutral_mass_fallback",
        "n_ms2_peaks",
        "observed_ion_bins_json",
        "observed_loss_bins_json",
    ]
    for col in keep:
        if col not in out.columns:
            out[col] = ""
    return out[keep].drop_duplicates(["control_id", "feature_key"]).reset_index(drop=True)


def fallback_mass_hypotheses(feature_rows: pd.DataFrame) -> pd.DataFrame:
    if feature_rows.empty:
        return pd.DataFrame(
            columns=[
                "feature_key",
                "control_id",
                "comp_id",
                "lib_id",
                "ion_mode",
                "parent_mz",
                "adduct_label",
                "neutral_mass",
                "adduct_rank",
            ]
        )
    out = feature_rows.copy()
    out["neutral_mass"] = pd.to_numeric(out.get("neutral_mass_fallback"), errors="coerce")
    out = out[out["neutral_mass"].notna()].copy()
    out["adduct_label"] = "feature_neutral_mass_fallback"
    out["adduct_rank"] = 1
    keep = ["feature_key", "control_id", "comp_id", "lib_id", "ion_mode", "parent_mz", "adduct_label", "neutral_mass", "adduct_rank"]
    return out[keep].drop_duplicates(["control_id", "feature_key", "neutral_mass"]).reset_index(drop=True)


def mass_hypotheses_for_scope(path: Path, feature_rows: pd.DataFrame, lib_id: int, ion_mode: str) -> pd.DataFrame:
    base = safe_read_csv(path)
    frames: list[pd.DataFrame] = []
    if not base.empty:
        scoped = base[
            pd.to_numeric(base.get("lib_id"), errors="coerce").eq(int(lib_id))
            & base.get("ion_mode", pd.Series("", index=base.index)).astype(str).eq(str(ion_mode))
        ].copy()
        if not scoped.empty:
            if "control_id" not in scoped.columns:
                scoped["control_id"] = scoped.get("feature_key", pd.Series("", index=scoped.index)).astype(str).str.extract(r"(X-\d+)", expand=False).fillna("")
            frames.append(scoped)
    frames.append(fallback_mass_hypotheses(feature_rows))
    out = pd.concat([frame for frame in frames if not frame.empty], ignore_index=True, sort=False) if frames else pd.DataFrame()
    if out.empty:
        return fallback_mass_hypotheses(feature_rows)
    out["neutral_mass"] = pd.to_numeric(out.get("neutral_mass"), errors="coerce")
    out = out[out["neutral_mass"].notna()].copy()
    for col in ["feature_key", "control_id", "comp_id", "lib_id", "ion_mode", "parent_mz", "adduct_label", "neutral_mass", "adduct_rank"]:
        if col not in out.columns:
            out[col] = ""
    return out.drop_duplicates(["control_id", "feature_key", "neutral_mass"]).reset_index(drop=True)


def pool_controls(path: Path) -> set[str]:
    if not path.exists():
        return set()
    try:
        raw = pd.read_csv(path, usecols=["control_id"], low_memory=False)
    except Exception:
        return set()
    return {clean_text(value) for value in raw["control_id"] if clean_text(value).startswith("X-")}


def pkl_path(args: argparse.Namespace, lib_id: int) -> Path:
    if lib_id == 305:
        return args.lib305_ions_pkl
    if lib_id == 400:
        return args.lib400_ions_pkl
    if lib_id == 402:
        return args.lib402_ions_pkl
    raise KeyError(lib_id)


def build_scope(args: argparse.Namespace, known: pd.DataFrame, lib_id: int, ion_mode: str) -> tuple[dict[str, Any], pd.DataFrame]:
    repo_root = Path(__file__).resolve().parent
    output_root = resolve(args.output_root, repo_root)
    known_scope = actionable_known(known, int(lib_id))
    if lib_id == 209:
        observed = observed_209(resolve(args.observed_209, repo_root), int(lib_id), str(ion_mode))
    else:
        merged_path = resolve(args.merged_ions_csv, repo_root)
        observed = combine_observed(
            observed_from_pkl(resolve(pkl_path(args, int(lib_id)), repo_root), int(lib_id), str(ion_mode)),
            observed_from_merged(merged_path, int(lib_id), str(ion_mode)),
        )
    feature_rows = feature_rows_from_observed(observed, int(lib_id), str(ion_mode))
    mass_hyp = mass_hypotheses_for_scope(resolve(args.neutral_mass_hypotheses, repo_root), feature_rows, int(lib_id), str(ion_mode))

    feature_path = output_root / f"lib{int(lib_id)}_{ion_mode}_feature_rows.csv"
    mass_path = output_root / f"lib{int(lib_id)}_{ion_mode}_neutral_mass_hypotheses.csv"
    observed_path = output_root / f"lib{int(lib_id)}_{ion_mode}_observed_ms2_lookup.csv"
    write_csv(feature_path, feature_rows)
    write_csv(mass_path, mass_hyp)
    write_csv(observed_path, observed)

    observed_controls = set(feature_rows.get("control_id", pd.Series(dtype=str)).astype(str))
    pool_path = resolve(Path(f"models/unified_pool_v20_lib{int(lib_id)}_{ion_mode}/unified_candidates_per_control.csv"), repo_root)
    processed_controls = pool_controls(pool_path)
    feature_summary = (
        feature_rows.sort_values(["control_id", "feature_key"], kind="mergesort")
        .drop_duplicates("control_id", keep="first")
        .set_index("control_id")
        if not feature_rows.empty
        else pd.DataFrame()
    )
    mass_controls = set(mass_hyp.get("control_id", pd.Series(dtype=str)).astype(str))

    actionable_ms2 = known_scope[known_scope["control_id"].isin(observed_controls)].copy()
    rows: list[dict[str, Any]] = []
    for _, row in actionable_ms2.sort_values("control_id", kind="mergesort").iterrows():
        control_id = clean_text(row.get("control_id"))
        feature = feature_summary.loc[control_id] if control_id in feature_summary.index else pd.Series(dtype=object)
        neutral = pd.to_numeric(pd.Series([feature.get("neutral_mass_fallback", math.nan)]), errors="coerce").iloc[0]
        parent = pd.to_numeric(pd.Series([feature.get("parent_mz", math.nan)]), errors="coerce").iloc[0]
        rows.append(
            {
                "control_id": control_id,
                "library_id": int(lib_id),
                "lib_id": int(lib_id),
                "ion_mode": str(ion_mode),
                "level": row.get("level", ""),
                "sublevel": row.get("sublevel", ""),
                "comp_id": int(row.get("comp_id")) if pd.notna(row.get("comp_id")) else "",
                "feature_key": feature.get("feature_key", ""),
                "parent_mz": float(parent) if pd.notna(parent) else "",
                "neutral_mass_estimate": float(neutral) if pd.notna(neutral) else "",
                "has_feature_rows": bool(control_id in observed_controls),
                "has_mass_hypotheses": bool(control_id in mass_controls or pd.notna(neutral)),
                "already_in_unified_pool": bool(control_id in processed_controls),
                "primary_stage": "already_processed" if control_id in processed_controls else "not_processed",
                "recovery_track": "pipeline_equalization_v1",
                "clean_external_info": row.get("clean_external_info", ""),
            }
        )
    all_manifest = pd.DataFrame(rows)
    gap = all_manifest[~all_manifest.get("already_in_unified_pool", pd.Series(dtype=bool)).astype(bool)].copy()
    all_path = output_root / f"lib{int(lib_id)}_{ion_mode}_all_actionable_ms2_manifest.csv"
    gap_path = output_root / f"lib{int(lib_id)}_{ion_mode}_gap_manifest.csv"
    controls_path = output_root / f"lib{int(lib_id)}_{ion_mode}_gap_controls.txt"
    controls_comma_path = output_root / f"lib{int(lib_id)}_{ion_mode}_gap_controls_comma.txt"
    write_csv(all_path, all_manifest)
    write_csv(gap_path, gap)
    controls = sorted(gap.get("control_id", pd.Series(dtype=str)).astype(str).tolist())
    controls_path.write_text("\n".join(controls) + ("\n" if controls else ""))
    controls_comma_path.write_text(",".join(controls) + ("\n" if controls else ""))
    summary = {
        "lib_id": int(lib_id),
        "ion_mode": str(ion_mode),
        "n_known_actionable": int(known_scope["control_id"].nunique()) if not known_scope.empty else 0,
        "n_observed_ms2_controls": int(len(observed_controls)),
        "n_actionable_with_ms2": int(all_manifest["control_id"].nunique()) if not all_manifest.empty else 0,
        "n_existing_unified_pool_controls": int(len(processed_controls)),
        "n_gap_controls": int(len(controls)),
        "n_feature_rows": int(feature_rows.shape[0]),
        "n_mass_hypothesis_rows": int(mass_hyp.shape[0]),
        "gap_manifest_csv": str(gap_path),
        "gap_controls_txt": str(controls_path),
        "feature_rows_csv": str(feature_path),
        "neutral_mass_hypotheses_csv": str(mass_path),
        "observed_ms2_lookup_csv": str(observed_path),
        "base_unified_pool_csv": str(pool_path),
    }
    return summary, observed


def run(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parent
    output_root = resolve(args.output_root, repo_root)
    output_root.mkdir(parents=True, exist_ok=True)
    known = pd.read_csv(resolve(args.known_x, repo_root), low_memory=False)
    summaries: list[dict[str, Any]] = []
    observed_frames: list[pd.DataFrame] = []
    for lib_id, ion_mode in LIB_SCOPES:
        if str(args.libs).strip() and str(lib_id) not in {part.strip() for part in str(args.libs).split(",") if part.strip()}:
            continue
        summary, observed = build_scope(args, known, int(lib_id), str(ion_mode))
        summaries.append(summary)
        if not observed.empty:
            observed_frames.append(observed)
    combined_observed = pd.concat(observed_frames, ignore_index=True, sort=False) if observed_frames else pd.DataFrame()
    write_csv(output_root / "observed_ms2_lookup_all_scopes.csv", combined_observed)
    summary_df = pd.DataFrame(summaries)
    write_csv(output_root / "pipeline_equalization_manifest_summary.csv", summary_df)
    payload = {
        "output_root": str(output_root),
        "n_total_gap_controls": int(summary_df["n_gap_controls"].sum()) if not summary_df.empty else 0,
        "n_total_actionable_with_ms2": int(summary_df["n_actionable_with_ms2"].sum()) if not summary_df.empty else 0,
        "per_scope": summaries,
        "combined_observed_lookup_csv": str(output_root / "observed_ms2_lookup_all_scopes.csv"),
    }
    write_json(output_root / "pipeline_equalization_manifest_summary.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--known-x", type=Path, default=Path("known_X.csv"))
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--libs", default="", help="Optional comma-separated library IDs to build.")
    parser.add_argument("--observed-209", type=Path, default=DEFAULT_OBSERVED_209)
    parser.add_argument("--neutral-mass-hypotheses", type=Path, default=DEFAULT_MASS_HYPOTHESES)
    parser.add_argument("--lib305-ions-pkl", type=Path, default=Path("ions_unknown_or_x_305.pkl"))
    parser.add_argument("--lib400-ions-pkl", type=Path, default=Path("ions_unknown_or_x_400_new.pkl"))
    parser.add_argument("--lib402-ions-pkl", type=Path, default=Path("ions_unknown_or_x_402.pkl"))
    parser.add_argument("--merged-ions-csv", type=Path, default=DEFAULT_MERGED_IONS)
    return parser.parse_args()


def main() -> int:
    run(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
