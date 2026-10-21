#!/usr/bin/env python3
"""Append recovery top-10 rows for compounds the main 10 ppm files do not have.

Existing compounds keep the rows already written. A recovery compound with no
rows in the recovery file stays absent.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

LIBS = ("209", "305", "400", "402")


def replace_ids(main: pd.DataFrame, updates: pd.DataFrame, ids: set[str]) -> pd.DataFrame:
    """Drop ``ids`` from ``main`` and append the update rows that belong to them."""

    if main.empty or "control_id" not in main.columns:
        kept = main
    else:
        kept = main[~main["control_id"].astype(str).isin(ids)].copy()
    if updates.empty or "control_id" not in updates.columns:
        return kept
    add = updates[updates["control_id"].astype(str).isin(ids)].copy()
    if add.empty:
        return kept
    return pd.concat([kept, add], ignore_index=True)


def merge_new_controls(main: pd.DataFrame, extra: pd.DataFrame) -> pd.DataFrame:
    """Return ``main`` plus rows from ``extra`` whose control id is new."""

    if extra.empty:
        return main
    if "control_id" not in extra.columns:
        return main
    have: set[str] = set()
    if not main.empty and "control_id" in main.columns:
        have = set(main["control_id"].astype(str))
    add = extra[~extra["control_id"].astype(str).isin(have)].copy()
    if add.empty:
        return main
    return pd.concat([main, add], ignore_index=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main-root", type=Path, required=True)
    parser.add_argument("--recovery-root", type=Path, required=True)
    parser.add_argument("--libs", default=",".join(LIBS))
    parser.add_argument("--name", default="unified_rerank_all_top10_ppm10.csv")
    args = parser.parse_args(argv)
    added: dict[str, int] = {}
    for lib in [item.strip() for item in args.libs.split(",") if item.strip()]:
        main_path = args.main_root / f"rerank_lib{lib}" / args.name
        extra_path = args.recovery_root / f"rerank_lib{lib}" / args.name
        if not extra_path.exists():
            added[lib] = 0
            continue
        extra = pd.read_csv(extra_path, low_memory=False)
        main = pd.read_csv(main_path, low_memory=False) if main_path.exists() else pd.DataFrame(columns=extra.columns)
        merged = merge_new_controls(main, extra)
        before = set(main["control_id"].astype(str)) if not main.empty and "control_id" in main.columns else set()
        after = set(merged["control_id"].astype(str)) if not merged.empty and "control_id" in merged.columns else set()
        added[lib] = len(after - before)
        if added[lib]:
            backup = main_path.with_name(main_path.name + ".pre_recover")
            if main_path.exists() and not backup.exists():
                main.to_csv(backup, index=False)
            merged.to_csv(main_path, index=False)
        print(f"[merge] lib{lib} new_compounds={added[lib]}")
    print("[merge] TOTAL", json.dumps({"new_compounds": sum(added.values()), "by_lib": added}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
