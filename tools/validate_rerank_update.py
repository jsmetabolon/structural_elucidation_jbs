#!/usr/bin/env python3
"""Validation harness for a pipeline rerank/scoring update (Track C).

Computes a compact metrics snapshot from a reports tree and compares two snapshots
(baseline vs. post-update), asserting the expected direction of change after the
mass-gate rerank + evidence-aware scoring land:

  - MSI Level-4 count DOWN            (mass-rejected structures recovered)
  - mass-mismatch rate DOWN (-> ~0)   (mass gate working)
  - median headline candidate score UP (evidence-aware score)
  - median evidence_confidence stable/high

Reads each report's ``confidence_inputs.json`` (msi_level, mass_match_score) and
``top_candidates.csv`` rank-1 (candidate_report_score, evidence_confidence). Pure
standard library.

Usage:
    python3 tools/validate_rerank_update.py snapshot <reports_root> --save baseline.json
    python3 tools/validate_rerank_update.py compare  <reports_root> --baseline baseline.json
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

csv.field_size_limit(10 ** 7)


def _num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f != f else f


def _median(xs):
    xs = sorted(x for x in xs if x is not None)
    return xs[len(xs) // 2] if xs else None


def snapshot(root: Path) -> dict:
    msi = Counter()
    mass_scores = []
    report_scores = []
    evidence = []
    n = 0
    for d in sorted(p for p in root.iterdir() if p.is_dir() and p.name.upper().startswith("X-")):
        slug = d.name.lower()
        ci = d / "data" / f"{slug}_confidence_inputs.json"
        if ci.exists():
            try:
                f = json.loads(ci.read_text())
            except (ValueError, OSError):
                f = {}
            n += 1
            lvl = f.get("msi_level")
            if lvl is not None:
                msi[int(lvl)] += 1
            mass_scores.append(_num(f.get("mass_match_score")))
        tc = d / "data" / f"{slug}_top_candidates.csv"
        if tc.exists():
            try:
                with tc.open(newline="") as h:
                    row = next(csv.DictReader(h), {})
            except (OSError, StopIteration):
                row = {}
            report_scores.append(_num(row.get("candidate_report_score")))
            evidence.append(_num(row.get("evidence_confidence")))
    mm = [m for m in mass_scores if m is not None]
    return {
        "n_reports": n,
        "msi_distribution": {str(k): msi[k] for k in sorted(msi)},
        "msi_L4": msi.get(4, 0),
        "mass_mismatch_rate": round(sum(1 for m in mm if m < 0.5) / len(mm), 4) if mm else None,
        "median_candidate_report_score": _median(report_scores),
        "median_evidence_confidence": _median(evidence),
    }


def compare(new: dict, base: dict) -> int:
    def line(label, b, nv, want):
        arrow = "?"
        ok = True
        if b is not None and nv is not None:
            if nv > b:
                arrow = "UP"
            elif nv < b:
                arrow = "DOWN"
            else:
                arrow = "="
            ok = (want == "up" and nv >= b) or (want == "down" and nv <= b) or (want == "any")
        print(f"  {label:42s} {b} -> {nv}  [{arrow}]  {'OK' if ok else 'UNEXPECTED'}")
        return ok
    print("[validate-rerank-update] baseline -> new")
    b_score, n_score = base.get("median_candidate_report_score"), new.get("median_candidate_report_score")
    b_ev, n_ev = base.get("median_evidence_confidence"), new.get("median_evidence_confidence")
    checks = [
        line("MSI Level-4 count (want DOWN)", base.get("msi_L4"), new.get("msi_L4"), "down"),
        line("mass-mismatch rate (want DOWN)", base.get("mass_mismatch_rate"), new.get("mass_mismatch_rate"), "down"),
        line("median headline score (want UP)", b_score, n_score, "up"),
        line("median evidence_confidence (any)", b_ev, n_ev, "any"),
    ]
    print(f"  MSI dist: {base.get('msi_distribution')} -> {new.get('msi_distribution')}")
    ok = all(checks)
    print("  RESULT:", "all expected" if ok else "some changes not in the expected direction")
    return 0 if ok else 1


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["snapshot", "compare"])
    ap.add_argument("reports_root", type=Path)
    ap.add_argument("--save", type=Path, help="snapshot mode: write metrics JSON here")
    ap.add_argument("--baseline", type=Path, help="compare mode: baseline metrics JSON")
    args = ap.parse_args()

    snap = snapshot(args.reports_root)
    if args.mode == "snapshot":
        print(json.dumps(snap, indent=2))
        if args.save:
            args.save.write_text(json.dumps(snap, indent=2))
            print(f"[validate-rerank-update] wrote baseline -> {args.save}")
        return
    if not args.baseline or not args.baseline.exists():
        raise SystemExit("compare mode requires an existing --baseline JSON")
    base = json.loads(args.baseline.read_text())
    raise SystemExit(compare(snap, base))


if __name__ == "__main__":
    main()
