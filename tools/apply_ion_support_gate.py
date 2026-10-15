#!/usr/bin/env python3
"""Drop newly added compounds that have no ion-intensity support.

Compounds already in the pre-recovery 10 ppm files are left as they are.
A newly added compound stays only when a scored candidate is within 10 ppm
and predicted_spectrum_score_v2 is strictly above 0.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from refill_top10_ppm import refill_lib
from merge_ppm10_recovery import replace_ids

LIBS = (209, 305, 400, 402)
NAME = "unified_rerank_all_top10_ppm10.csv"
IONS = "unified_rerank_all_top10_ppm10_ions.csv"


def main() -> int:
    base = Path("output")
    main = base / "pipeline_equalization_adduct_derive"
    extras = [
        base / "pipeline_equalization_adduct_derive_recover",
        base / "pipeline_equalization_adduct_derive_library10",
        base / "pipeline_equalization_adduct_derive_denovo",
    ]
    old: set[str] = set()
    for lib in LIBS:
        backup = pd.read_csv(
            main / f"rerank_lib{lib}" / f"{NAME}.pre_recover", usecols=["control_id"]
        )
        old |= set(backup.control_id.astype(str))
    added: set[str] = set()
    for lib in LIBS:
        current = pd.read_csv(main / f"rerank_lib{lib}" / NAME, usecols=["control_id"])
        added |= set(current.control_id.astype(str)) - old
    print(f"added {len(added)}", flush=True)
    for root in extras:
        for lib in LIBS:
            lib_dir = root / f"rerank_lib{lib}"
            scored = lib_dir / "unified_candidates_scored.csv"
            if not scored.exists() or scored.stat().st_size <= 1:
                continue
            summary = refill_lib(
                lib_dir,
                ppm_ceiling=10.0,
                top_k=10,
                out_name=IONS,
                min_spectrum_score=0.0,
            )
            print(
                root.name,
                lib,
                "written",
                summary["n_compounds_written"],
                "empty",
                summary["n_empty"],
                flush=True,
            )
    kept_added: set[str] = set()
    removed: set[str] = set()
    for lib in LIBS:
        main_path = main / f"rerank_lib{lib}" / NAME
        main_df = pd.read_csv(main_path, low_memory=False)
        frames = []
        for root in extras:
            path = root / f"rerank_lib{lib}" / IONS
            if path.exists() and path.stat().st_size > 1:
                frames.append(pd.read_csv(path, low_memory=False))
        updates = (
            pd.concat(frames, ignore_index=True, sort=False)
            if frames
            else pd.DataFrame(columns=main_df.columns)
        )
        if not updates.empty and "control_id" in updates.columns:
            updates = updates.drop_duplicates(["control_id", "predicted_smiles"], keep="first")
        merged = replace_ids(main_df, updates, added)
        before = set(main_df.control_id.astype(str)) & added
        after = set(merged.control_id.astype(str)) & added if not merged.empty else set()
        kept_added |= after
        removed |= before - after
        merged.to_csv(main_path, index=False)
        print(
            "lib",
            lib,
            "rows",
            len(main_df),
            "->",
            len(merged),
            "kept_added",
            len(after),
            "removed",
            len(before - after),
            flush=True,
        )
    (main / "mass_only_ids.txt").write_text("\n".join(sorted(removed)) + "\n", encoding="utf-8")
    (main / "ions_kept_ids.txt").write_text("\n".join(sorted(kept_added)) + "\n", encoding="utf-8")
    published: set[str] = set()
    for lib in LIBS:
        published |= set(
            pd.read_csv(main / f"rerank_lib{lib}" / NAME, usecols=["control_id"]).control_id.astype(str)
        )
    print(
        f"unique {len(published)} removed {len(removed)} kept_added {len(kept_added)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
