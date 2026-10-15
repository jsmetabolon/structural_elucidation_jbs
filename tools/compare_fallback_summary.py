#!/usr/bin/env python3
"""Compare the stored-mass-first summary with the all-adduct-derived summary.

One row per compound in either summary. A compound only in the earlier summary
gets one reason, checked in this order:

1. ``ineligible``: Exclude or Artifact, or not on libraries 209, 305, 400, 402
2. ``no_inputs``: a fallback ID with no feature row and no hypothesis
3. ``fallback_mass_only``: the fallback search has a candidate within 10 ppm,
   but no such candidate has ``predicted_spectrum_score_v2`` above 0
4. ``fallback_base_pool_not_rerun``: every earlier published row came from a
   base-pool generator, which the fallback job does not rerun
5. ``fallback_over_10ppm``: the fallback search has candidates, none within 10 ppm
6. ``fallback_no_candidate``: the fallback search produced no candidate
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from prepare_ppm10_recovery_inputs import LIBS, catalog_targets
from refill_top10_ppm import MASS_REFERENCE_KEYS, abs_mass_error_ppm, mass_reference

GATED = "unified_rerank_all_top10_ppm10.csv"
SCORED = ("unified_candidates_scored.csv", "non_s0_candidates_scored.csv")
LIBRARY_SOURCES = {"library_anchored_candidates_v1"}
NEEDED = {"control_id", "abs_mass_error_da", "predicted_spectrum_score_v2", *MASS_REFERENCE_KEYS}


def read_ids(path: Path) -> set[str]:
    return set(path.read_text(encoding="utf-8").split()) if path.exists() else set()


def gated_rows(root: Path) -> pd.DataFrame:
    frames = []
    for lib in LIBS:
        path = root / f"rerank_lib{lib}" / GATED
        if path.exists() and path.stat().st_size > 12:
            frame = pd.read_csv(path, low_memory=False)
            frame["lib"] = lib
            frames.append(frame)
    out = pd.concat(frames, ignore_index=True, sort=False)
    out["control_id"] = out["control_id"].astype(str)
    return out


def best_by_compound(root: Path) -> pd.DataFrame:
    """Per compound: best ppm, and whether any row within 10 ppm has spectrum above 0."""

    frames = []
    for lib in LIBS:
        for name in SCORED:
            path = root / f"rerank_lib{lib}" / name
            if not path.exists() or path.stat().st_size <= 12:
                continue
            try:
                frame = pd.read_csv(path, low_memory=False, usecols=lambda c: c in NEEDED)
            except pd.errors.EmptyDataError:
                continue
            if "control_id" in frame.columns and not frame.empty:
                frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["control_id", "best_ppm", "n_le10", "n_le10_spectrum"])
    pool = pd.concat(frames, ignore_index=True, sort=False)
    pool["control_id"] = pool["control_id"].astype(str)
    pool["ppm"] = abs_mass_error_ppm(pool.get("abs_mass_error_da", pd.Series(dtype="float64")), mass_reference(pool))
    spectrum = pd.to_numeric(pool.get("predicted_spectrum_score_v2"), errors="coerce").fillna(0.0)
    pool["le10"] = pool["ppm"].le(10.0)
    pool["le10_spectrum"] = pool["le10"] & spectrum.gt(0.0)
    return (
        pool.groupby("control_id")
        .agg(best_ppm=("ppm", "min"), n_le10=("le10", "sum"), n_le10_spectrum=("le10_spectrum", "sum"))
        .reset_index()
    )


def reason(row: pd.Series) -> str:
    if not row["eligible"]:
        return "ineligible"
    if row["no_inputs"]:
        return "no_inputs"
    if row["fallback_n_le10"] > 0:
        return "fallback_mass_only"
    if row["earlier_base_pool_only"]:
        return "fallback_base_pool_not_rerun"
    if pd.notna(row["fallback_best_ppm"]):
        return "fallback_over_10ppm"
    return "fallback_no_candidate"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--new-root", type=Path, required=True)
    parser.add_argument("--earlier-root", type=Path, required=True)
    parser.add_argument("--stored-root", type=Path, required=True)
    parser.add_argument("--fallback-root", type=Path, required=True)
    parser.add_argument("--fallback-selection", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    args = parser.parse_args(argv)

    new = gated_rows(args.new_root)
    earlier = gated_rows(args.earlier_root)
    settled = read_ids(args.new_root / "settled_ids.txt")
    fallback_ids = read_ids(args.new_root / "fallback_ids.txt")
    searched = {str(item["control_id"]) for item in json.loads(args.fallback_selection.read_text())["assignments"]}
    eligible = set(catalog_targets(pd.read_csv(args.catalog, low_memory=False))[0])

    sources = (
        earlier.assign(candidate_source=earlier.get("candidate_source", pd.Series("", index=earlier.index)).fillna("").astype(str))
        .groupby("control_id")["candidate_source"]
        .agg(lambda s: ";".join(sorted(set(s))))
    )
    earlier_spectrum = (
        pd.to_numeric(earlier.get("predicted_spectrum_score_v2"), errors="coerce")
        .fillna(0.0)
        .groupby(earlier["control_id"])
        .max()
    )
    stored = best_by_compound(args.stored_root).set_index("control_id").add_prefix("stored_")
    fallback = best_by_compound(args.fallback_root).set_index("control_id").add_prefix("fallback_")

    new_ids = set(new["control_id"])
    earlier_ids = set(earlier["control_id"])
    ledger = pd.DataFrame({"x_compound_id": sorted(new_ids | earlier_ids)}).set_index("x_compound_id")
    ledger["in_new"] = ledger.index.isin(new_ids)
    ledger["in_earlier"] = ledger.index.isin(earlier_ids)
    ledger["status"] = "both"
    ledger.loc[ledger["in_earlier"] & ~ledger["in_new"], "status"] = "earlier_only"
    ledger.loc[ledger["in_new"] & ~ledger["in_earlier"], "status"] = "new_only"
    ledger["eligible"] = ledger.index.isin(eligible)
    ledger["settled_stored_mass"] = ledger.index.isin(settled)
    ledger["fallback_id"] = ledger.index.isin(fallback_ids)
    ledger["no_inputs"] = ledger["fallback_id"] & ~ledger.index.isin(searched)
    ledger["new_source"] = "none"
    ledger.loc[ledger["in_new"] & ledger["settled_stored_mass"], "new_source"] = "stored_mass"
    ledger.loc[ledger["in_new"] & ~ledger["settled_stored_mass"], "new_source"] = "adduct_fallback"
    ledger = ledger.join(stored).join(fallback)
    for column in ("stored_n_le10", "stored_n_le10_spectrum", "fallback_n_le10", "fallback_n_le10_spectrum"):
        ledger[column] = ledger[column].fillna(0).astype(int)
    ledger["earlier_candidate_sources"] = sources.reindex(ledger.index).fillna("")
    ledger["earlier_max_spectrum"] = earlier_spectrum.reindex(ledger.index)
    ledger["earlier_spectrum_supported"] = ledger["earlier_max_spectrum"].fillna(0.0).gt(0.0)
    ledger["earlier_base_pool_only"] = ledger["earlier_candidate_sources"].map(
        lambda value: bool(value) and not (set(value.split(";")) & LIBRARY_SOURCES)
    )
    ledger["reason"] = ""
    only = ledger["status"].eq("earlier_only")
    ledger.loc[only, "reason"] = ledger.loc[only].apply(reason, axis=1)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    ledger.reset_index().to_csv(args.out_csv, index=False)
    payload = {
        "new_compounds": len(new_ids),
        "earlier_compounds": len(earlier_ids),
        "status": ledger["status"].value_counts().to_dict(),
        "new_source": ledger.loc[ledger["in_new"], "new_source"].value_counts().to_dict(),
        "earlier_only_reason": ledger.loc[only, "reason"].value_counts().to_dict(),
        "earlier_only_stored_mass_only": int((ledger.loc[only, "stored_n_le10"] > 0).sum()),
        "earlier_spectrum_supported": int(ledger.loc[ledger["in_earlier"], "earlier_spectrum_supported"].sum()),
        "earlier_only_spectrum_supported": int(ledger.loc[only, "earlier_spectrum_supported"].sum()),
        "earlier_only_spectrum_supported_by_reason": ledger.loc[only & ledger["earlier_spectrum_supported"], "reason"]
        .value_counts()
        .to_dict(),
    }
    args.out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print("[compare]", json.dumps(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
