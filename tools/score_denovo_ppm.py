"""Score decoded structures that sit within 10 ppm of the corrected neutral mass."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
from pathlib import Path

import pandas as pd

from structural_elucidation.nce.verifier_v2 import predicted_spectrum_score_v2

XID = re.compile(r"(X-\d+)")


def ppm(mass: float, target: float) -> float:
    if target <= 0:
        return math.inf
    return abs(mass - target) / target * 1e6


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--decode-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    reference_dir = os.environ.get("NCE_REFERENCE_DIR", "").strip() or None
    targets = pd.read_csv(args.targets, dtype=str, keep_default_na=False)
    target_mass = {
        row.x_compound_id: float(row.proton_neutral)
        for row in targets.itertuples(index=False)
        if str(row.proton_neutral).strip()
    }
    meta = targets.set_index("x_compound_id")
    frames = []
    for path in args.decode_root.rglob("candidate_metrics_all.csv"):
        frame = pd.read_csv(path, low_memory=False)
        if frame.empty or "predicted_smiles" not in frame.columns:
            continue
        frames.append(frame)
    rows: list[dict[str, object]] = []
    if frames:
        decoded = pd.concat(frames, ignore_index=True, sort=False)
        decoded["x_compound_id"] = decoded["feature_key"].astype(str).str.extract(XID, expand=False)
        decoded = decoded[decoded.x_compound_id.isin(target_mass)]
        decoded["predicted_exact_mass"] = pd.to_numeric(decoded["predicted_exact_mass"], errors="coerce")
        decoded = decoded.dropna(subset=["predicted_exact_mass"])
        decoded = decoded.drop_duplicates(["x_compound_id", "predicted_smiles"])
        for cid, group in decoded.groupby("x_compound_id"):
            target = target_mass.get(str(cid))
            if target is None:
                parent = pd.to_numeric(group["parent_mz"], errors="coerce").dropna()
                if parent.empty:
                    continue
                ion = str(meta.loc[str(cid)].ion_mode) if str(cid) in meta.index else ""
                shift = 1.007276466583 if ion == "negative" else -1.007276466583
                target = float(parent.iloc[0]) + shift
                target_mass[str(cid)] = target
            group = group.copy()
            group["abs_mass_error_da"] = (group["predicted_exact_mass"] - target).abs()
            group["abs_mass_error_ppm"] = group["abs_mass_error_da"] / target * 1e6
            close = group[group.abs_mass_error_ppm.le(10)].sort_values("abs_mass_error_ppm")
            record = meta.loc[str(cid)]
            best = None
            best_score = 0.0
            for product in close.itertuples(index=False):
                score = float(
                    predicted_spectrum_score_v2(
                        product.predicted_smiles,
                        record.ion_json or "[]",
                        record.loss_json or "[]",
                        lib_id=int(float(record.lib_id or 209)),
                        reference_dir=reference_dir,
                    )
                )
                if score > best_score:
                    best_score = score
                    best = product
            if best is None or best_score <= 0:
                rows.append({"x_compound_id": cid, "passed": False, "n_on_mass": int(len(close)), "best_score": best_score})
                continue
            rows.append(
                {
                    "x_compound_id": cid,
                    "passed": True,
                    "lib_id": record.lib_id,
                    "ion_mode": record.ion_mode,
                    "proton_neutral": target,
                    "candidate_smiles": best.predicted_smiles,
                    "candidate_formula": getattr(best, "predicted_formula", ""),
                    "predicted_exact_mass": float(best.predicted_exact_mass),
                    "abs_mass_error_da": float(best.abs_mass_error_da),
                    "abs_mass_error_ppm": float(best.abs_mass_error_ppm),
                    "predicted_spectrum_score_v2": best_score,
                    "n_on_mass": int(len(close)),
                }
            )
    scored_ids = {row["x_compound_id"] for row in rows}
    for cid in target_mass:
        if cid not in scored_ids:
            rows.append({"x_compound_id": cid, "passed": False, "n_on_mass": 0, "best_score": 0.0})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, index=False)
    passed = sum(1 for row in rows if row["passed"])
    print(json.dumps({"compounds": len(rows), "passed": passed, "out": str(args.out)}))


if __name__ == "__main__":
    main()
