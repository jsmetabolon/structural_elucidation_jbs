#!/usr/bin/env python3
"""Paired McNemar of rank-1 / rank-3 recovery: baseline vs adduct-derive rerank.

Does not change closest_neutral_mass. Run after EC2
``output/pipeline_equalization_adduct_derive`` exists.

Uses ``truth_tanimoto`` on the rerank CSVs when present (threshold 0.7, same as
``tools/audit_mass_policy.py``). Controls missing from either side are skipped.

Usage:
    python3 tools/compare_adduct_derive_rerank.py \\
      --baseline models/pipeline_equalization \\
      --derive output/pipeline_equalization_adduct_derive
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

csv.field_size_limit(10**7)

ROOT = Path(__file__).resolve().parent.parent
TRUTH_THRESHOLD = 0.70


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact binomial p-value for discordant pairs."""
    n = b + c
    if n == 0:
        return 1.0
    tail = sum(math.comb(n, k) for k in range(min(b, c) + 1))
    return min(1.0, 2.0 * tail / (2**n))

LIBS = (209, 305, 400, 402)
MASS_CUT_DA = 0.005


def _mass_class(err: float | None) -> str:
    if err is None:
        return "unknown"
    return "within_5mda" if abs(err) <= MASS_CUT_DA else "gt_5mda"


def mass_class_shift(base_path: Path, derive_path: Path) -> dict[str, int]:
    """Count rank-1 rows whose abs_mass_error_da class flipped vs baseline."""

    def rank1_err(path: Path) -> dict[str, float | None]:
        out: dict[str, float | None] = {}
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            err_col = None
            for name in ("abs_mass_error_da", "primary_abs_mass_error_da", "mass_error_da"):
                if name in (reader.fieldnames or []):
                    err_col = name
                    break
            for row in reader:
                cid = str(row.get("control_id") or "").strip()
                rank = int(float(row["rank"])) if _num(row.get("rank")) else 10**9
                if not cid or rank != 1:
                    continue
                out[cid] = _num(row.get(err_col)) if err_col else None
        return out

    base = rank1_err(base_path)
    deriv = rank1_err(derive_path)
    flipped = same = missing_err = 0
    for cid in set(base) & set(deriv):
        if base[cid] is None or deriv[cid] is None:
            missing_err += 1
            continue
        if _mass_class(base[cid]) != _mass_class(deriv[cid]):
            flipped += 1
        else:
            same += 1
    return {
        "n_shared_rank1": len(set(base) & set(deriv)),
        "class_flipped": flipped,
        "class_same": same,
        "missing_error_col": missing_err,
    }


def _num(value: object) -> float | None:
    try:
        if value is None or str(value).strip() in {"", "nan", "None"}:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def load_topn(path: Path, n: int) -> dict[str, list[float | None]]:
    """control_id -> list of truth_tanimoto by rank (1-indexed positions)."""
    by: dict[str, list[tuple[int, float | None]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "control_id" not in (reader.fieldnames or []):
            raise SystemExit(f"{path} has no control_id column")
        tan_col = "truth_tanimoto" if "truth_tanimoto" in (reader.fieldnames or []) else None
        rank_col = "rank" if "rank" in (reader.fieldnames or []) else None
        for row in reader:
            cid = str(row.get("control_id") or "").strip()
            if not cid:
                continue
            rank = int(float(row[rank_col])) if rank_col and _num(row.get(rank_col)) else 10**9
            tan = _num(row.get(tan_col)) if tan_col else None
            by.setdefault(cid, []).append((rank, tan))
    out: dict[str, list[float | None]] = {}
    for cid, rows in by.items():
        rows.sort(key=lambda item: item[0])
        out[cid] = [tan for _, tan in rows[:n]]
    return out


def recovered(tans: list[float | None], k: int) -> bool:
    return any((tan or 0.0) >= TRUTH_THRESHOLD for tan in tans[:k])


def pair_counts(
    base: dict[str, list[float | None]],
    deriv: dict[str, list[float | None]],
    k: int,
) -> dict[str, int | float]:
    shared = sorted(set(base) & set(deriv))
    gained = lost = both = neither = 0
    for cid in shared:
        before = recovered(base[cid], k)
        after = recovered(deriv[cid], k)
        if after and not before:
            gained += 1
        elif before and not after:
            lost += 1
        elif before and after:
            both += 1
        else:
            neither += 1
    return {
        "n_shared_controls": len(shared),
        "both_recovered": both,
        "neither": neither,
        "gained": gained,
        "lost": lost,
        "mcnemar_exact_p": round(mcnemar_exact(lost, gained), 6),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=ROOT / "models" / "pipeline_equalization")
    parser.add_argument(
        "--derive",
        type=Path,
        default=ROOT / "output" / "pipeline_equalization_adduct_derive",
    )
    parser.add_argument("--json-out", type=Path, default=None)
    args = parser.parse_args()
    derive_root = args.derive if args.derive.is_absolute() else ROOT / args.derive
    base_root = args.baseline if args.baseline.is_absolute() else ROOT / args.baseline

    missing = []
    payload: dict[str, object] = {
        "truth_threshold": TRUTH_THRESHOLD,
        "libs": {},
        "note": (
            "Do not change closest_neutral_mass until rank-1/rank-3 pairing "
            "is measured on a completed derive run."
        ),
    }
    for lib in LIBS:
        bpath = base_root / f"rerank_lib{lib}" / "unified_rerank_all_top10.csv"
        dpath = derive_root / f"rerank_lib{lib}" / "unified_rerank_all_top10.csv"
        if not bpath.is_file():
            missing.append(str(bpath))
            continue
        if not dpath.is_file():
            missing.append(str(dpath))
            continue
        base = load_topn(bpath, 10)
        deriv = load_topn(dpath, 10)
        payload["libs"][str(lib)] = {
            "rank1": pair_counts(base, deriv, 1),
            "rank3": pair_counts(base, deriv, 3),
        }

    if missing and not payload["libs"]:
        print("Derive (or baseline) rerank CSVs are not present yet. Skip McNemar.")
        print("Missing:")
        for path in missing:
            print(f"  {path}")
        print("After EC2 Job B:")
        print(
            "  python3 tools/compare_adduct_derive_rerank.py "
            "--baseline models/pipeline_equalization "
            "--derive output/pipeline_equalization_adduct_derive"
        )
        print("Export CSVs per docs/ec2_export_handoff.md (cloud VMs cannot SSH to 10.12.5.28).")
        sys.exit(2)

    payload["mass_consistency"] = {}
    for lib in LIBS:
        bpath = base_root / f"rerank_lib{lib}" / "unified_rerank_all_top10.csv"
        dpath = derive_root / f"rerank_lib{lib}" / "unified_rerank_all_top10.csv"
        if bpath.is_file() and dpath.is_file():
            payload["mass_consistency"][str(lib)] = mass_class_shift(bpath, dpath)

    text = json.dumps(payload, indent=2)
    print(text)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
