#!/usr/bin/env python3
"""Extract + parse the 'Partially Characterized Molecules' (PCM) library set.

These are named reference-library metabolites (metabolon_super = 'Partially
Characterized Molecules') whose structure is unresolved but whose partial
characterization is encoded in the name, e.g. "glucuronide of C12H22O4 (3)*"
(conjugate moiety + core/aglycone formula + isomer index). This tool turns them
into (a) a human-readable task table and (b) a diagnosis CSV that the existing
``structural_elucidation.pipeline.candidates`` stage consumes (it mass-matches
library cores and expands conjugates from ``conjugate_hints``).

Pure standard library. Read-only against the metabolite DB.

Usage:
    python3 tools/extract_pcm_tasks.py [--db PATH] [--out-dir docs/pcm]
"""
from __future__ import annotations

import argparse
import csv
import re
import sqlite3
from pathlib import Path

# Conjugate keywords to detect in PCM names, longest/most-specific first so that
# e.g. "glutathione" is matched before "glutamine"/"glutamate". The hint is the key
# passed to structural_elucidation.chemistry.conjugates.CONJUGATE_REACTIONS.
CONJUGATE_KEYWORDS = [
    "n-acetylglucosamine", "glucuronide", "glutathione", "glutamine",
    "glutamate", "glycine", "carnitine", "taurine", "cysteine", "sulfate",
]
# Reactions already implemented in chemistry/conjugates.py::CONJUGATE_REACTIONS.
KNOWN_REACTIONS = {"glucuronide", "sulfate", "glycine", "carnitine", "taurine"}

# monoisotopic masses for a quick core-mass sanity value (not authoritative)
_MONO = {"C": 12.0, "H": 1.0078250319, "N": 14.0030740052, "O": 15.9949146221,
         "S": 31.97207069, "P": 30.97376151, "Cl": 34.96885271, "F": 18.9984032,
         "Br": 78.9183376, "Na": 22.98976928}
_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")
_CORE_FORMULA = re.compile(r"\b(C\d+H\d+[A-Za-z0-9]*)\b")
_ISOMER_IDX = re.compile(r"\((\d+)\)")


def mono_mass(formula: str) -> float | None:
    if not formula:
        return None
    total = 0.0
    seen = False
    for el, n in _FORMULA_TOKEN.findall(formula):
        if not el:
            continue
        if el not in _MONO:
            return None
        total += _MONO[el] * (int(n) if n else 1)
        seen = True
    return round(total, 5) if seen else None


def parse_name(name: str) -> dict:
    """Return {conjugate, core_formula, isomer_index} parsed from a PCM name.

    Robust to the observed name formats: "glucuronide of C12H22O4 (3)*",
    "glycine conjugate of C10H12O2*", "N-acetylglucosamine conjugate of C24H40O4
    bile acid**", and loose forms like "metabolonic lactone sulfate".
    """
    low = name.lower()
    conjugate = next((k for k in CONJUGATE_KEYWORDS if k in low), "")
    fm = _CORE_FORMULA.search(name)
    idx = _ISOMER_IDX.search(name)
    return {
        "conjugate": conjugate,
        "core_formula": fm.group(1) if fm else "",
        "isomer_index": idx.group(1) if idx else "",
    }


def find_db(explicit: str | None) -> Path:
    if explicit:
        return Path(explicit)
    for c in ("vendor/metabolite_synonyms.sqlite", "_bundle/metabolite_synonyms.sqlite"):
        if Path(c).exists():
            return Path(c)
    raise SystemExit("metabolite DB not found; pass --db")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=None)
    ap.add_argument("--out-dir", type=Path, default=Path("docs/pcm"))
    args = ap.parse_args()

    db = find_db(args.db)
    con = sqlite3.connect(str(db))
    con.row_factory = sqlite3.Row
    rows = [dict(r) for r in con.execute(
        "SELECT chemical_id, canonical_name, smiles, inchikeys, resolution_bucket, resolution_status "
        "FROM metabolites WHERE metabolon_super = 'Partially Characterized Molecules'"
    )]
    args.out_dir.mkdir(parents=True, exist_ok=True)

    tasks = []
    for r in rows:
        name = str(r.get("canonical_name") or "")
        p = parse_name(name)
        conj = p["conjugate"]
        has_smiles = bool(str(r.get("smiles") or "").strip())
        needs_reaction = bool(conj) and conj not in KNOWN_REACTIONS
        tasks.append({
            "chemical_id": r["chemical_id"],
            "canonical_name": name,
            "resolution_bucket": r.get("resolution_bucket") or "",
            "conjugate": conj,
            "conjugate_hint": conj,
            "reaction_available": "" if not conj else ("yes" if not needs_reaction else "NO_ADD_REACTION"),
            "core_formula": p["core_formula"],
            "core_monoisotopic_mass": mono_mass(p["core_formula"]) or "",
            "isomer_index": p["isomer_index"],
            "has_existing_smiles": "yes" if has_smiles else "",
        })

    # (a) human-readable task table
    tasks_path = args.out_dir / "pcm_tasks.csv"
    with tasks_path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(tasks[0]))
        w.writeheader()
        w.writerows(tasks)

    # (b) diagnosis CSV for pipeline.candidates (control_id/primary_stage/conjugate_hints)
    diag_path = args.out_dir / "pcm_diagnosis_seed.csv"
    with diag_path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=[
            "control_id", "primary_stage", "conjugate_hints", "core_formula",
            "core_monoisotopic_mass", "isomer_index", "canonical_name",
        ])
        w.writeheader()
        for t in tasks:
            if t["has_existing_smiles"]:
                continue  # already resolved; skip re-prediction
            w.writerow({
                "control_id": t["chemical_id"],
                "primary_stage": "pcm_unresolved",  # != "recovered" -> targeted by pipeline.candidates
                "conjugate_hints": t["conjugate_hint"],
                "core_formula": t["core_formula"],
                "core_monoisotopic_mass": t["core_monoisotopic_mass"],
                "isomer_index": t["isomer_index"],
                "canonical_name": t["canonical_name"],
            })

    # summary
    from collections import Counter
    conj_counts = Counter(t["conjugate"] or "(none)" for t in tasks)
    need = [t["chemical_id"] for t in tasks if t["reaction_available"] == "NO_ADD_REACTION"]
    with_core = sum(1 for t in tasks if t["core_formula"])
    print(f"[extract-pcm-tasks] db={db}  PCM entries={len(tasks)}")
    print(f"  conjugate breakdown: {dict(conj_counts.most_common())}")
    print(f"  parsed a core formula for {with_core}/{len(tasks)}")
    print(f"  conjugate types missing a reaction (add to conjugates.py): "
          f"{sorted(set(t['conjugate'] for t in tasks if t['reaction_available'] == 'NO_ADD_REACTION'))} "
          f"({len(need)} entries)")
    print(f"  wrote {tasks_path}  and  {diag_path}")


if __name__ == "__main__":
    main()
