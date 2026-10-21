#!/usr/bin/env python3
"""Build de novo inputs for compounds published on mass match alone."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

LIBS = (209, 305, 400, 402)
ION = {209: "negative", 305: "negative", 400: "positive", 402: "positive"}
NAME = "unified_rerank_all_top10_ppm10.csv"


def catalog_library(catalog: pd.DataFrame, ids: set[str]) -> dict[str, int]:
    """Assign each ID to the earliest run library that lists it.

    Library order is 209, 305, 400, 402. IDs absent from those libraries are
    left out.
    """

    work = catalog.copy()
    work["lib_id"] = pd.to_numeric(work["lib_id"], errors="coerce")
    work = work[work["lib_id"].isin(LIBS)]
    work["x_compound_id"] = work["x_compound_id"].astype(str)
    work = work[work["x_compound_id"].isin(ids)]
    order = {lib: index for index, lib in enumerate(LIBS)}
    work["ord"] = work["lib_id"].map(order)
    work = work.sort_values(["ord", "x_compound_id"]).drop_duplicates("x_compound_id")
    return {str(row.x_compound_id): int(row.lib_id) for row in work.itertuples()}


def assign_libraries(
    ids: set[str],
    ppm_by_lib: list[tuple[int, set[str]]],
    catalog: pd.DataFrame | None = None,
) -> dict[str, int]:
    """Prefer the library of a published 10 ppm row, then the catalog."""

    assigned: dict[str, int] = {}
    for lib, present in ppm_by_lib:
        for cid in present:
            if cid in ids and cid not in assigned:
                assigned[cid] = lib
    if catalog is not None:
        for cid, lib in catalog_library(catalog, ids - set(assigned)).items():
            assigned.setdefault(cid, lib)
    return assigned


def main() -> int:
    ids = set(Path(sys.argv[1]).read_text(encoding="utf-8").split())
    extras = [Path(item) for item in sys.argv[2].split(",") if item]
    feature_sources = [Path(item) for item in sys.argv[3].split(",")]
    out = Path(sys.argv[4])
    catalog = pd.read_csv(sys.argv[5], low_memory=False) if len(sys.argv) > 5 else None
    out.mkdir(parents=True, exist_ok=True)
    ppm_by_lib: list[tuple[int, set[str]]] = []
    for lib in LIBS:
        present: set[str] = set()
        for root in extras:
            path = root / f"rerank_lib{lib}" / NAME
            if not path.exists() or path.stat().st_size <= 1:
                continue
            frame = pd.read_csv(path, usecols=["control_id"])
            present |= set(frame.control_id.astype(str))
        ppm_by_lib.append((lib, present))
    assigned = assign_libraries(ids, ppm_by_lib, catalog)
    missing_lib = sorted(ids - set(assigned))
    assignments = [{"control_id": cid, "lib": lib} for cid, lib in sorted(assigned.items())]
    selection = {"assignments": assignments, "run_libs": "209,305,400,402"}
    (out / "selection.json").write_text(json.dumps(selection, indent=2), encoding="utf-8")
    (out / "ids.txt").write_text("\n".join(sorted(assigned)) + "\n", encoding="utf-8")
    feature_out = out / "features_union"
    feature_out.mkdir(parents=True, exist_ok=True)
    for lib in LIBS:
        ion = ION[lib]
        for kind in ("feature_rows", "neutral_mass_hypotheses"):
            frames = []
            for source in feature_sources:
                path = source / f"lib{lib}_{ion}_{kind}.csv"
                if path.exists() and path.stat().st_size > 1:
                    frames.append(pd.read_csv(path, low_memory=False))
            if not frames:
                continue
            merged = pd.concat(frames, ignore_index=True, sort=False)
            if "control_id" in merged.columns:
                merged = merged.drop_duplicates("control_id", keep="first")
            merged.to_csv(feature_out / f"lib{lib}_{ion}_{kind}.csv", index=False)
    print(
        json.dumps(
            {
                "assigned": len(assigned),
                "missing_lib": missing_lib,
                "by_lib": {str(lib): sum(v == lib for v in assigned.values()) for lib in LIBS},
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
