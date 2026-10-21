#!/usr/bin/env python3
"""Inputs for the scoped SELFIES de novo run on compounds with no library hit at 10 ppm.

``slice``: write per-library encoder feature tables holding only the selected
compounds, each on the library the recovery selection placed it on. Compounds
absent from the encoder table are listed in ``missing_ids.txt``.

``inputs``: collect decoded candidates for those compounds from every
``--candidate-roots`` entry, recompute the mass
error against the adduct-derived neutral masses, and write
``lib*_candidate_metrics_all.csv`` plus a gap manifest per library for
``scripts/run_ppm10_extra_job.sh`` with ``EXTRA_START_AT=verifier``.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Callable

import pandas as pd

ION = {209: "negative", 305: "negative", 400: "positive", 402: "positive"}
XID = re.compile(r"(X-\d+)")
MASS_SCORE_TOLERANCE_DA = 1.0


def selected_by_lib(ids: set[str], selection: dict) -> dict[int, set[str]]:
    by_lib: dict[int, set[str]] = {}
    for row in selection.get("assignments", []):
        cid = str(row["control_id"])
        if cid in ids:
            by_lib.setdefault(int(row["lib"]), set()).add(cid)
    return by_lib


def _xid(series: pd.Series) -> pd.Series:
    return series.astype(str).str.extract(XID, expand=False)


def rescore_to_targets(
    candidates: pd.DataFrame,
    targets: pd.DataFrame,
    score: Callable[[float, float], float],
) -> pd.DataFrame:
    """Set each row's mass error to the closest adduct-derived target for its compound.

    Rows whose compound has no target, or whose exact mass is missing, are dropped.
    """

    if candidates.empty or targets.empty:
        return candidates.iloc[0:0].copy()
    by_cid: dict[str, list[float]] = {}
    for cid, mass in zip(targets["control_id"].astype(str), pd.to_numeric(targets["target_neutral_mass"], errors="coerce")):
        if pd.notna(mass):
            by_cid.setdefault(cid, []).append(float(mass))
    out = candidates.copy()
    masses = pd.to_numeric(out["predicted_exact_mass"], errors="coerce")
    errors: list[float] = []
    chosen: list[float] = []
    for cid, mass in zip(out["control_id"].astype(str), masses):
        options = by_cid.get(cid, [])
        if not options or pd.isna(mass):
            errors.append(math.nan)
            chosen.append(math.nan)
            continue
        target = min(options, key=lambda value: abs(float(mass) - value))
        errors.append(abs(float(mass) - target))
        chosen.append(target)
    out["stored_abs_mass_error_da"] = out.get("abs_mass_error_da")
    out["abs_mass_error_da"] = errors
    out["target_neutral_mass"] = chosen
    out["mass_hypothesis_neutral_mass"] = chosen
    out["mass_match_score"] = [
        score(value, MASS_SCORE_TOLERANCE_DA) if not math.isnan(value) else math.nan for value in errors
    ]
    return out[out["abs_mass_error_da"].notna()].reset_index(drop=True)


def cmd_slice(args: argparse.Namespace) -> int:
    ids = set(args.ids.read_text(encoding="utf-8").split())
    by_lib = selected_by_lib(ids, json.loads(args.selection.read_text(encoding="utf-8")))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    found: set[str] = set()
    written: list[str] = []
    for lib, wanted in sorted(by_lib.items()):
        ion = ION[lib]
        table = args.bundle_dir / f"lib_{lib}_{ion}" / "feature_table_prepared.csv"
        frame = pd.read_csv(table, low_memory=False)
        cid = _xid(frame["compound_name"])
        part = frame[cid.isin(wanted) & pd.to_numeric(frame["lib_id"], errors="coerce").eq(lib)]
        found |= set(_xid(part["compound_name"]).dropna())
        if part.empty:
            continue
        part.to_csv(args.out_dir / f"lib{lib}_{ion}_features.csv", index=False)
        written.append(str(lib))
    missing = sorted(ids - found)
    (args.out_dir / "missing_ids.txt").write_text("\n".join(missing) + ("\n" if missing else ""), encoding="utf-8")
    print(f"[denovo-slice] compounds={len(found)} missing={len(missing)} libs={','.join(written) or 'none'}")
    return 0


def cmd_inputs(args: argparse.Namespace) -> int:
    from structural_elucidation.chemistry.conjugates import linear_mass_match_score
    from structural_elucidation.pipeline.candidates import mass_hypotheses

    ids = set(args.ids.read_text(encoding="utf-8").split())
    by_lib = selected_by_lib(ids, json.loads(args.selection.read_text(encoding="utf-8")))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    used: list[str] = []
    counts: dict[str, dict[str, int]] = {}
    for lib, wanted in sorted(by_lib.items()):
        ion = ION[lib]
        frames = []
        for root in args.candidate_roots:
            path = root / f"lib_{lib}_{ion}" / "selfies_denovo" / "candidate_metrics_all.csv"
            if not path.exists() or path.stat().st_size <= 1:
                continue
            frame = pd.read_csv(path, low_memory=False)
            if frame.empty or "feature_key" not in frame.columns:
                continue
            frame["control_id"] = _xid(frame["feature_key"])
            frames.append(frame[frame["control_id"].isin(wanted)])
        if not frames:
            continue
        decoded = pd.concat(frames, ignore_index=True, sort=False)
        decoded = decoded.drop_duplicates(["control_id", "predicted_smiles"], keep="first")
        features = pd.read_csv(args.feature_dir / f"lib{lib}_{ion}_feature_rows.csv", low_memory=False)
        hypotheses = pd.read_csv(args.feature_dir / f"lib{lib}_{ion}_neutral_mass_hypotheses.csv", low_memory=False)
        targets = mass_hypotheses(features, hypotheses, derive_from_adduct=True)
        rescored = rescore_to_targets(decoded, targets, linear_mass_match_score)
        if rescored.empty:
            continue
        rescored["lib_id"] = lib
        rescored["ion_mode"] = ion
        rescored["candidate_source"] = "merged_structuremasst_selfies_v1"
        rescored = rescored.sort_values(["control_id", "abs_mass_error_da"], kind="mergesort")
        rescored["rank"] = rescored.groupby("control_id").cumcount() + 1
        rescored.to_csv(args.out_dir / f"lib{lib}_{ion}_candidate_metrics_all.csv", index=False)
        manifest = pd.DataFrame({"control_id": sorted(set(rescored["control_id"])), "conjugate_hints": ""})
        manifest.to_csv(args.out_dir / f"lib{lib}_{ion}_gap_manifest.csv", index=False)
        ppm = rescored["abs_mass_error_da"] / pd.to_numeric(rescored["predicted_exact_mass"], errors="coerce") * 1e6
        within = set(rescored.loc[ppm.le(10.0), "control_id"])
        counts[str(lib)] = {
            "rows": int(len(rescored)),
            "compounds": int(rescored["control_id"].nunique()),
            "compounds_within_10ppm": len(within),
        }
        used.append(str(lib))
    payload = {"run_libs": ",".join(used), "counts": counts}
    (args.out_dir / "selection.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print("[denovo-inputs]", json.dumps(payload))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    slice_parser = sub.add_parser("slice")
    slice_parser.add_argument("--bundle-dir", type=Path, required=True)
    slice_parser.add_argument("--ids", type=Path, required=True)
    slice_parser.add_argument("--selection", type=Path, required=True)
    slice_parser.add_argument("--out-dir", type=Path, required=True)
    slice_parser.set_defaults(func=cmd_slice)
    inputs_parser = sub.add_parser("inputs")
    inputs_parser.add_argument(
        "--candidate-roots",
        type=Path,
        nargs="+",
        required=True,
        help="roots holding lib_<lib>_<ion>/selfies_denovo/candidate_metrics_all.csv",
    )
    inputs_parser.add_argument("--ids", type=Path, required=True)
    inputs_parser.add_argument("--selection", type=Path, required=True)
    inputs_parser.add_argument("--feature-dir", type=Path, required=True)
    inputs_parser.add_argument("--out-dir", type=Path, required=True)
    inputs_parser.set_defaults(func=cmd_inputs)
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
