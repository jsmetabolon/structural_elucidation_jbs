#!/usr/bin/env python3
"""Reference prototype for the upstream reranker fixes (handoff items 2 and 5).

The authoritative fix belongs in the inference/rerank pipeline (Track A, their
env). This tool demonstrates the two changes on the *delivered* data (each report's
``top_candidates.csv``) and quantifies their impact, so the pipeline team has a
precise, tested spec:

  Item 2 - hard mass gate: if the rank-1 candidate's mass is rejected
    (mass_match_score < MASS_WEAK) and a lower-ranked candidate fits
    (>= MASS_STRONG), promote the best mass-consistent candidate to rank 1.
  Item 5 - evidence-aware headline score: report ``evidence_confidence`` (library +
    mass + spectral) as the candidate's headline number instead of the raw
    generative ``candidate_report_score``; keep the generative score as an isomer
    tie-breaker.

It is read-only by default (prints projected impact). Pure standard library.

Usage:
    python3 tools/prototype_mass_gate_rerank.py <reports_root>
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
import confidence as C  # noqa: E402

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


def gate_report(rows: list[dict]) -> dict:
    """Apply the mass gate to one candidate list; return before/after summary."""
    if not rows:
        return {}
    r1 = rows[0]
    m1 = _num(r1.get("mass_match_score"))
    old_headline = _num(r1.get("candidate_report_score"))
    gated = r1
    changed = False
    if m1 is not None and m1 < C.MASS_WEAK:
        for row in rows[1:]:
            ms = _num(row.get("mass_match_score"))
            if ms is not None and ms >= C.MASS_STRONG:
                gated = row
                changed = True
                break
    return {
        "changed": changed,
        "mass_rejected": m1 is not None and m1 < C.MASS_WEAK,
        "recovered": changed,
        "old_headline": old_headline,
        "new_headline": _num(gated.get("evidence_confidence")),
        "gated_mass": _num(gated.get("mass_match_score")),
    }


def main() -> None:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("_audit_pull/reports")
    n = 0
    changed = 0
    mass_rejected = 0
    reject_no_alt = 0
    old_scores = []
    new_scores = []
    for d in sorted(p for p in root.iterdir() if p.is_dir() and p.name.upper().startswith("X-")):
        tc = d / "data" / f"{d.name.lower()}_top_candidates.csv"
        if not tc.exists():
            continue
        with tc.open(newline="") as h:
            rows = list(csv.DictReader(h))
        s = gate_report(rows)
        if not s:
            continue
        n += 1
        old_scores.append(s["old_headline"])
        new_scores.append(s["new_headline"])
        if s["mass_rejected"]:
            mass_rejected += 1
            if s["recovered"]:
                changed += 1
            else:
                reject_no_alt += 1
    print(f"[prototype-mass-gate-rerank] reports analyzed: {n}")
    print("\nItem 2 - hard mass gate (projected):")
    print(f"  rank-1 mass-rejected today:            {mass_rejected}")
    print(f"  recoverable by gate (promote alt):     {changed}  "
          f"({100 * changed / max(mass_rejected, 1):.0f}% of mass-rejected)")
    print(f"  mass-rejected with no in-list alt:     {reject_no_alt}  (need conjugate enumeration)")
    print(f"  => projected MSI L4 -> L3+ recovery:   ~{changed} reports "
          f"({100 * changed / max(n, 1):.0f}% of all)")
    print("\nItem 5 - evidence-aware headline score (projected):")
    print(f"  median raw candidate_report_score:     {_median(old_scores)}")
    print(f"  median evidence_confidence (headline): {_median(new_scores)}")
    print("\nAuthoritative fix belongs in the reranker; see docs/pipeline_refinement_handoff.md items 2 and 5.")


if __name__ == "__main__":
    main()
