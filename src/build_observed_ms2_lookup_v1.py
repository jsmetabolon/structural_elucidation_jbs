#!/usr/bin/env python3
# flake8: noqa
"""Build per-control observed MS2 lookup for real predicted-spectrum scoring."""

from __future__ import annotations

import argparse
import json
import pickle
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd

from evaluate_denovo_failure_modes import load_x_controls
from nce_evidence_utils import bin_base, ion_and_loss_bins, neutral_mass_from_precursor, parse_json_list
from nce_local_runtime_utils import DEFAULT_NCE_ROOT, DEFAULT_REFERENCE_DIR, get_cached_mlims_service
from run_retrieval_anchored_denovo_assembly_v4 import resolve
from x_domain_analog_utils import dataframe_to_markdown, safe_float, write_csv, write_json


def parse_bins(value: Any) -> list[int]:
    out: list[int] = []
    for item in parse_json_list(value):
        try:
            intval = int(float(item))
        except Exception:
            continue
        if intval:
            out.append(intval)
    return sorted(set(out))


def parse_float_list(value: Any) -> list[float]:
    out: list[float] = []
    for item in parse_json_list(value):
        try:
            out.append(float(item))
        except Exception:
            continue
    return out


def approximate_ms2_pairs(ion_bins: list[int], intensities: list[float] | None = None) -> list[list[float]]:
    pairs: list[list[float]] = []
    if not ion_bins:
        return pairs
    if intensities is None or len(intensities) != len(ion_bins):
        intensities = [1.0 / float(len(ion_bins))] * len(ion_bins)
    for bin_id, intensity in zip(ion_bins, intensities):
        if bin_id <= 0:
            continue
        try:
            pairs.append([round(float(bin_base(int(bin_id))), 6), round(float(intensity), 6)])
        except Exception:
            continue
    return pairs


def _aligned_intensity_rows(bins_to_raw_intensity: dict[int, float]) -> tuple[list[int], list[float]]:
    bins = sorted(int(bin_id) for bin_id in bins_to_raw_intensity if int(bin_id) > 0)
    raw = [max(0.0, float(bins_to_raw_intensity.get(bin_id, 0.0))) for bin_id in bins]
    max_raw = max(raw) if raw else 0.0
    if max_raw <= 0:
        return bins, [1.0 / max(len(bins), 1)] * len(bins)
    return bins, [min(1.0, max(0.0, value / max_raw)) for value in raw]


def read_csv_or_empty(path: Path) -> pd.DataFrame:
    if not path.exists() or str(path) == "/dev/null":
        return pd.DataFrame()
    try:
        return pd.read_csv(path, low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def read_pickle_dataframe(path: Path) -> pd.DataFrame:
    try:
        with path.open("rb") as handle:
            obj = pickle.load(handle)
        return obj if isinstance(obj, pd.DataFrame) else pd.DataFrame(obj)
    except Exception:
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = Path(tmp) / "ions.csv"
            script = (
                "import pickle, sys, pandas as pd; "
                "obj=pickle.load(open(sys.argv[1], 'rb')); "
                "df=obj if isinstance(obj, pd.DataFrame) else pd.DataFrame(obj); "
                "df.to_csv(sys.argv[2], index=False)"
            )
            subprocess.run(["python3", "-c", script, str(path), str(csv_path)], check=True)
            return pd.read_csv(csv_path, low_memory=False)


def load_feature_rows_from_ions_pkl(path: Path, lib_id: int, ion_mode: str) -> pd.DataFrame:
    raw = read_pickle_dataframe(path)
    if raw.empty:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    work = raw.copy()
    if "control_id" not in work.columns:
        if "compound_name" in work.columns:
            work["control_id"] = work["compound_name"].astype(str)
        elif "feature_key" in work.columns:
            work["control_id"] = work["feature_key"].astype(str).str.extract(r"(X-\d+)", expand=False).fillna("")
        else:
            work["control_id"] = ""
    group_cols = ["control_id"]
    if "feature_key" in work.columns:
        group_cols.append("feature_key")
    for keys, group in work.groupby(group_cols, dropna=False, sort=False):
        first = group.iloc[0]
        control_id = str(first.get("control_id", ""))
        if not control_id.startswith("X-"):
            continue
        parent_mz = safe_float(first.get("parent_mz", first.get("parent_mass", float("nan"))), float("nan"))
        source_neutral_mass = safe_float(first.get("source_neutral_mass"), float("nan"))
        neutral_mass = source_neutral_mass if pd.notna(source_neutral_mass) else neutral_mass_from_precursor(parent_mz, int(lib_id))
        ion_intensities: dict[int, float] = {}
        loss_intensities: dict[int, float] = {}
        for _, peak in group.iterrows():
            mz = safe_float(peak.get("mass", peak.get("mz", float("nan"))), float("nan"))
            if mz <= 0:
                continue
            raw_intensity = safe_float(peak.get("intensity", peak.get("rel_intensity", 1.0)), 1.0)
            ion_bin, loss_bin = ion_and_loss_bins(mz, parent_mz, int(lib_id))
            if ion_bin:
                ion_intensities[int(ion_bin)] = max(float(raw_intensity), ion_intensities.get(int(ion_bin), 0.0))
            if loss_bin:
                loss_intensities[int(loss_bin)] = max(float(raw_intensity), loss_intensities.get(int(loss_bin), 0.0))
        ions, ion_values = _aligned_intensity_rows(ion_intensities)
        losses, loss_values = _aligned_intensity_rows(loss_intensities)
        feature_key = str(first.get("feature_key", f"ions_pkl_{control_id}"))
        rows.append(
            {
                "control_id": control_id,
                "lib_id": int(lib_id),
                "ion_mode": str(ion_mode),
                "feature_key": feature_key,
                "comp_id": first.get("comp_id", ""),
                "parent_mz": parent_mz,
                "neutral_mass": neutral_mass,
                "observed_ion_bins_json": json.dumps(sorted(ions)),
                "observed_loss_bins_json": json.dumps(sorted(losses)),
                "observed_ion_intensities_json": json.dumps([round(float(value), 6) for value in ion_values]),
                "observed_loss_intensities_json": json.dumps([round(float(value), 6) for value in loss_values]),
                "ms2_pairs_json": json.dumps(approximate_ms2_pairs(ions, ion_values)),
                "n_observed_ion_bins": int(len(ions)),
                "n_observed_loss_bins": int(len(losses)),
                "n_ms2_pairs": int(len(ions)),
                "ms2_present": bool(ions or losses),
                "spectrum_source": "ions_pkl",
            }
        )
    return pd.DataFrame(rows)


def load_feature_rows_from_peak_table(peaks: pd.DataFrame, lib_id: int, ion_mode: str) -> pd.DataFrame:
    if peaks.empty or "mass" not in peaks.columns:
        return pd.DataFrame()
    work = peaks.copy()
    if "control_id" not in work.columns:
        if "compound_name" in work.columns:
            work["control_id"] = work["compound_name"].astype(str)
        else:
            work["control_id"] = work.get("feature_key", pd.Series("", index=work.index)).astype(str).str.extract(r"(X-\d+)", expand=False).fillna("")
    rows: list[dict[str, Any]] = []
    group_cols = ["control_id"]
    if "feature_key" in work.columns:
        group_cols.append("feature_key")
    for _keys, group in work.groupby(group_cols, dropna=False, sort=False):
        first = group.iloc[0]
        control_id = str(first.get("control_id", ""))
        if not control_id.startswith("X-"):
            continue
        parent_mz = safe_float(first.get("parent_mz", float("nan")), float("nan"))
        source_neutral_mass = safe_float(first.get("source_neutral_mass"), float("nan"))
        neutral_mass = source_neutral_mass if pd.notna(source_neutral_mass) else neutral_mass_from_precursor(parent_mz, int(lib_id))
        ion_intensities: dict[int, float] = {}
        loss_intensities: dict[int, float] = {}
        for _, peak in group.iterrows():
            mz = safe_float(peak.get("mass", peak.get("mz", float("nan"))), float("nan"))
            if mz <= 0:
                continue
            raw_intensity = safe_float(peak.get("intensity", peak.get("rel_intensity", 1.0)), 1.0)
            ion_bin, loss_bin = ion_and_loss_bins(mz, parent_mz, int(lib_id))
            if ion_bin:
                ion_intensities[int(ion_bin)] = max(float(raw_intensity), ion_intensities.get(int(ion_bin), 0.0))
            if loss_bin:
                loss_intensities[int(loss_bin)] = max(float(raw_intensity), loss_intensities.get(int(loss_bin), 0.0))
        ions, ion_values = _aligned_intensity_rows(ion_intensities)
        losses, loss_values = _aligned_intensity_rows(loss_intensities)
        rows.append(
            {
                "control_id": control_id,
                "lib_id": int(lib_id),
                "ion_mode": str(ion_mode),
                "feature_key": str(first.get("feature_key", f"html_{control_id}")),
                "comp_id": first.get("comp_id", ""),
                "parent_mz": parent_mz,
                "neutral_mass": neutral_mass,
                "observed_ion_bins_json": json.dumps(sorted(ions)),
                "observed_loss_bins_json": json.dumps(sorted(losses)),
                "observed_ion_intensities_json": json.dumps([round(float(value), 6) for value in ion_values]),
                "observed_loss_intensities_json": json.dumps([round(float(value), 6) for value in loss_values]),
                "ms2_pairs_json": json.dumps(approximate_ms2_pairs(ions, ion_values)),
                "n_observed_ion_bins": int(len(ions)),
                "n_observed_loss_bins": int(len(losses)),
                "n_ms2_pairs": int(len(ions)),
                "ms2_present": bool(ions or losses),
                "spectrum_source": "html_curated",
            }
        )
    return pd.DataFrame(rows)


def mlims_unknown_lookup(lib_id: int, ion_mode: str, nce_root: Path, reference_dir: Path) -> dict[str, dict[str, Any]]:
    """Build spectrum-only fallback rows from cached MLIMS unknown MS2 records."""
    try:
        service = get_cached_mlims_service(nce_root=str(Path(str(nce_root)).expanduser()), reference_dir=str(Path(str(reference_dir)).expanduser()))
        unknown = service.get_unknown_ms2(int(lib_id), mass_ratio_lower_limit=1)
    except Exception:
        return {}
    if unknown.empty or "parent_name" not in unknown.columns:
        return {}
    out: dict[str, dict[str, Any]] = {}
    for name, group in unknown.groupby("parent_name", dropna=False):
        control_id = str(name or "")
        if not control_id.startswith("X-"):
            continue
        first = group.iloc[0]
        parent_mz = safe_float(first.get("parent_mass"), float("nan"))
        ions: set[int] = set()
        losses: set[int] = set()
        for _, peak in group.iterrows():
            mz = safe_float(peak.get("mass"), float("nan"))
            if mz <= 0:
                continue
            ion_bin, loss_bin = ion_and_loss_bins(mz, parent_mz, int(lib_id))
            if ion_bin:
                ions.add(int(ion_bin))
            if loss_bin:
                losses.add(int(loss_bin))
        out[control_id] = {
            "control_id": control_id,
            "lib_id": int(lib_id),
            "ion_mode": str(ion_mode),
            "feature_key": f"mlims_unknown_ms2_{control_id}",
            "comp_id": first.get("comp_id", ""),
            "parent_mz": parent_mz,
            "neutral_mass": neutral_mass_from_precursor(parent_mz, int(lib_id)),
            "observed_ion_bins_json": json.dumps(sorted(ions)),
            "observed_loss_bins_json": json.dumps(sorted(losses)),
            "observed_ion_intensities_json": "[]",
            "observed_loss_intensities_json": "[]",
            "ms2_pairs_json": json.dumps(approximate_ms2_pairs(sorted(ions))),
            "n_observed_ion_bins": int(len(ions)),
            "n_observed_loss_bins": int(len(losses)),
            "n_ms2_pairs": int(len(ions)),
            "ms2_present": bool(ions or losses),
            "spectrum_source": "mlims_tables_get_unknown_ms2",
        }
    return out


def run(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parent
    output_root = resolve(args.output_root, repo_root)
    output_root.mkdir(parents=True, exist_ok=True)
    feature_path = resolve(args.feature_rows, repo_root)
    feature_rows = read_csv_or_empty(feature_path)
    used_ions_pkl = False
    used_peak_table = False
    ions_path = resolve(args.ions_pkl, repo_root) if getattr(args, "ions_pkl", None) is not None else None
    if feature_rows.empty and ions_path is not None and ions_path.exists():
        feature_rows = load_feature_rows_from_ions_pkl(ions_path, lib_id=int(args.lib_id), ion_mode=str(args.ion_mode))
        used_ions_pkl = True
    elif not feature_rows.empty and not {"observed_ion_bins_json", "observed_loss_bins_json"}.issubset(set(feature_rows.columns)) and "mass" in feature_rows.columns:
        feature_rows = load_feature_rows_from_peak_table(feature_rows, lib_id=int(args.lib_id), ion_mode=str(args.ion_mode))
        used_peak_table = True
    mlims_fallback = mlims_unknown_lookup(int(args.lib_id), str(args.ion_mode), args.nce_root, args.reference_dir) if bool(args.enable_mlims_fallback) else {}
    controls = load_x_controls(
        known_x=repo_root / "known_X.csv",
        known_x_all=repo_root / "known_X_all_available_structures_simple.csv",
        lib_id=int(args.lib_id),
        ion_mode=str(args.ion_mode),
    )
    if not feature_rows.empty and "control_id" in feature_rows.columns:
        control_ids = sorted({str(value) for value in feature_rows["control_id"].astype(str) if str(value).startswith("X-")})
    else:
        control_ids = [str(control.control_id) for control in controls]
    rows: list[dict[str, Any]] = []
    for control_id in control_ids:
        scoped = feature_rows[feature_rows.get("control_id", pd.Series(dtype=str)).astype(str).eq(str(control_id))]
        if scoped.empty:
            fallback = mlims_fallback.get(str(control_id))
            if fallback is not None:
                rows.append(fallback)
            else:
                rows.append(
                    {
                        "control_id": control_id,
                        "lib_id": int(args.lib_id),
                        "ion_mode": str(args.ion_mode),
                        "feature_key": "",
                        "comp_id": "",
                        "parent_mz": "",
                        "neutral_mass": "",
                        "observed_ion_bins_json": "[]",
                        "observed_loss_bins_json": "[]",
                        "observed_ion_intensities_json": "[]",
                        "observed_loss_intensities_json": "[]",
                        "ms2_pairs_json": "[]",
                        "n_observed_ion_bins": 0,
                        "n_observed_loss_bins": 0,
                        "n_ms2_pairs": 0,
                        "ms2_present": False,
                        "spectrum_source": "missing",
                    }
                )
            continue
        first = scoped.iloc[0]
        ion_bins = parse_bins(first.get("observed_ion_bins_json", "[]"))
        loss_bins = parse_bins(first.get("observed_loss_bins_json", "[]"))
        ion_intensities = parse_float_list(first.get("observed_ion_intensities_json", "[]"))
        loss_intensities = parse_float_list(first.get("observed_loss_intensities_json", "[]"))
        pairs = approximate_ms2_pairs(ion_bins, ion_intensities)
        rows.append(
            {
                "control_id": control_id,
                "lib_id": int(args.lib_id),
                "ion_mode": str(args.ion_mode),
                "feature_key": first.get("feature_key", ""),
                "comp_id": first.get("comp_id", ""),
                "parent_mz": first.get("parent_mz", ""),
                "neutral_mass": first.get("neutral_mass_fallback", first.get("neutral_mass", "")),
                "observed_ion_bins_json": json.dumps(ion_bins),
                "observed_loss_bins_json": json.dumps(loss_bins),
                "observed_ion_intensities_json": json.dumps(ion_intensities) if len(ion_intensities) == len(ion_bins) else "[]",
                "observed_loss_intensities_json": json.dumps(loss_intensities) if len(loss_intensities) == len(loss_bins) else "[]",
                "ms2_pairs_json": json.dumps(pairs),
                "n_observed_ion_bins": int(len(ion_bins)),
                "n_observed_loss_bins": int(len(loss_bins)),
                "n_ms2_pairs": int(len(pairs)),
                "ms2_present": bool(ion_bins or loss_bins),
                "spectrum_source": "feature_peptide_rows",
            }
        )
    out = pd.DataFrame(rows)
    suffix = str(getattr(args, "output_suffix", "") or "")
    scoped_root = output_root / f"lib_{int(args.lib_id)}_{str(args.ion_mode)}{suffix}"
    legacy_root = output_root if int(args.lib_id) == 209 and str(args.ion_mode) == "negative" and not suffix else scoped_root
    legacy_root.mkdir(parents=True, exist_ok=True)
    scoped_root.mkdir(parents=True, exist_ok=True)
    write_csv(legacy_root / "per_control.csv", out)
    if scoped_root != legacy_root:
        write_csv(scoped_root / "per_control.csv", out)
    summary = pd.DataFrame(
        [
            {
                "n_controls": int(out.shape[0]),
                "n_controls_with_ms2": int(out["ms2_present"].sum()) if not out.empty else 0,
                "min_required_with_ms2": int(args.min_controls_with_ms2),
                "acceptance_pass": bool(int(out["ms2_present"].sum()) >= int(args.min_controls_with_ms2)) if not out.empty else False,
                "feature_rows": str(feature_path),
                "ions_pkl": str(ions_path) if ions_path is not None else "",
                "used_ions_pkl": bool(used_ions_pkl),
                "used_peak_table": bool(used_peak_table),
            }
        ]
    )
    write_csv(legacy_root / "summary.csv", summary)
    if scoped_root != legacy_root:
        write_csv(scoped_root / "summary.csv", summary)
    report = [
        "# Observed MS2 Lookup V1",
        "",
        "## Summary",
        dataframe_to_markdown(summary),
        "",
        "## Missing controls",
        dataframe_to_markdown(out[~out["ms2_present"]][["control_id", "feature_key", "parent_mz"]]),
        "",
    ]
    (legacy_root / "report.md").write_text("\n".join(report))
    if scoped_root != legacy_root:
        (scoped_root / "report.md").write_text("\n".join(report))
    payload = summary.iloc[0].to_dict()
    write_json(legacy_root / "summary.json", payload)
    if scoped_root != legacy_root:
        write_json(scoped_root / "summary.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=Path("models/observed_ms2_lookup_v1"))
    parser.add_argument("--feature-rows", type=Path, default=Path("models/peptide_first_denovo_v1_1/feature_peptide_rows.csv"))
    parser.add_argument("--lib-id", type=int, default=209)
    parser.add_argument("--ion-mode", choices=["negative", "positive"], default="negative")
    parser.add_argument("--ions-pkl", type=Path, default=None, help="Optional raw ions_*.pkl path for libs without feature rows")
    parser.add_argument("--output-suffix", type=str, default="", help="Suffix appended to lib_<id>_<mode>/ output dir, e.g. '_html'")
    parser.add_argument("--min-controls-with-ms2", type=int, default=50)
    parser.add_argument("--enable-mlims-fallback", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--nce-root", type=Path, default=DEFAULT_NCE_ROOT)
    parser.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCE_DIR)
    payload = run(parser.parse_args())
    raise SystemExit(0 if payload["acceptance_pass"] else 2)


if __name__ == "__main__":
    main()
