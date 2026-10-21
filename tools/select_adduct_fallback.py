#!/usr/bin/env python3
"""Keep stored-mass predictions that ion intensities support. The rest fall back.

A compound is settled when a scored candidate is within 10 ppm and
``predicted_spectrum_score_v2`` is strictly above 0. First library wins in the
order 209, 305, 400, 402. Eligible catalog compounds that are not settled are
written to ``fallback_ids.txt`` for the adduct-derived search.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from prepare_ppm10_recovery_inputs import LIBS, catalog_targets
from refill_top10_ppm import concat_scored, refill_lib

OUT_NAME = "unified_rerank_all_top10_ppm10.csv"
SCORED_NAMES = ("unified_candidates_scored.csv", "non_s0_candidates_scored.csv")


def pool_dir(lib_dir: Path, staging: Path) -> Path:
    """Use the full scored pool, or the ungated top 10 when that pool is absent."""

    if any((lib_dir / name).exists() and (lib_dir / name).stat().st_size > 1 for name in SCORED_NAMES):
        return lib_dir
    top10 = lib_dir / "unified_rerank_all_top10.csv"
    if not top10.exists() or top10.stat().st_size <= 1:
        return lib_dir
    staging.mkdir(parents=True, exist_ok=True)
    link = staging / "unified_candidates_scored.csv"
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(top10.resolve())
    return staging


def claim_first_library(gated: dict[int, pd.DataFrame]) -> dict[int, pd.DataFrame]:
    """Drop a compound from later libraries once an earlier library has a row."""

    seen: set[str] = set()
    claimed: dict[int, pd.DataFrame] = {}
    for lib in LIBS:
        frame = gated.get(lib, pd.DataFrame())
        if frame.empty or "control_id" not in frame.columns:
            claimed[lib] = frame.iloc[0:0].copy() if not frame.empty else frame
            continue
        keep = ~frame["control_id"].astype(str).isin(seen)
        kept = frame.loc[keep].copy()
        seen |= set(kept["control_id"].astype(str))
        claimed[lib] = kept
    return claimed


def split_ids(settled: set[str], eligible: list[str]) -> tuple[list[str], list[str]]:
    """Return settled IDs that are eligible, then eligible IDs that are not settled."""

    eligible_set = set(eligible)
    kept_settled = sorted(settled & eligible_set)
    fallback = [cid for cid in eligible if cid not in settled]
    return kept_settled, fallback


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scored-root", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--ppm-ceiling", type=float, default=10.0)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--min-spectrum-score", type=float, default=0.0)
    args = parser.parse_args(argv)

    gated: dict[int, pd.DataFrame] = {}
    for lib in LIBS:
        lib_dir = args.scored_root / f"rerank_lib{lib}"
        if not lib_dir.is_dir():
            raise SystemExit(f"missing {lib_dir}")
        source = pool_dir(lib_dir, args.out_root / "_pool" / f"rerank_lib{lib}")
        dest = args.out_root / f"rerank_lib{lib}" / OUT_NAME
        try:
            concat_scored(source)
        except SystemExit:
            dest.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(columns=["control_id"]).to_csv(dest, index=False)
            gated[lib] = pd.DataFrame()
            continue
        refill_lib(
            source,
            ppm_ceiling=args.ppm_ceiling,
            top_k=args.top_k,
            out_name=OUT_NAME,
            min_spectrum_score=args.min_spectrum_score,
            out_path=dest,
        )
        gated[lib] = pd.read_csv(dest, low_memory=False) if dest.stat().st_size > 1 else pd.DataFrame()
    claimed = claim_first_library(gated)
    kept, _excluded = catalog_targets(pd.read_csv(args.catalog, low_memory=False))
    eligible = set(kept)
    settled_ids: set[str] = set()
    for frame in claimed.values():
        if not frame.empty and "control_id" in frame.columns:
            settled_ids |= set(frame["control_id"].astype(str))
    settled, fallback = split_ids(settled_ids, kept)
    for lib, frame in claimed.items():
        dest = args.out_root / f"rerank_lib{lib}" / OUT_NAME
        dest.parent.mkdir(parents=True, exist_ok=True)
        if frame.empty or "control_id" not in frame.columns:
            pd.DataFrame(columns=["control_id"]).to_csv(dest, index=False)
            continue
        kept_rows = frame[frame["control_id"].astype(str).isin(eligible)]
        kept_rows.to_csv(dest, index=False)
    args.out_root.mkdir(parents=True, exist_ok=True)
    (args.out_root / "settled_ids.txt").write_text("\n".join(settled) + ("\n" if settled else ""), encoding="utf-8")
    (args.out_root / "fallback_ids.txt").write_text("\n".join(fallback) + ("\n" if fallback else ""), encoding="utf-8")
    payload = {"settled": len(settled), "fallback": len(fallback), "eligible": len(kept)}
    (args.out_root / "fallback_selection.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print("[fallback]", json.dumps(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
