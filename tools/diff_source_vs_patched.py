#!/usr/bin/env python3
"""Compare from-source rendered reports against the enrich->patch delivered reports.

For every control present in both trees, checks that the confidence-relevant
display elements (MSI badge/line, evidence-confidence widget, relabeled rank score,
Metabolon Level line, mass-consistent-alternative callout, and the internal MS2
panel) are present in both, i.e. that a clean from-source render reaches parity
with the patched bundle. Pure standard library.

Usage:
    python3 tools/diff_source_vs_patched.py <source_root> <patched_root>
"""
from __future__ import annotations

import sys
from pathlib import Path

# marker -> (filename, must-match-presence)
PUBLIC_MARKERS = [
    "MSI Level",
    "MSI identification level:",
    "Evidence confidence",
    "Metabolon identification Level",
    "Mass-consistent alternative",
]
INTERNAL_MARKERS = [
    "Observed fragment peaks",
    "Forward spectral score",
]
# public compliance: these must NOT appear in the public index.html
PUBLIC_FORBIDDEN = ["Observed fragment peaks", "fragment-ion bins", "neutral-loss bins"]


def has(path: Path, marker: str) -> bool:
    if not path.exists():
        return False
    return marker in path.read_text(errors="ignore")


def main() -> None:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/regen")
    pat = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("_audit_pull/reports")
    controls = sorted(
        d.name for d in src.iterdir()
        if d.is_dir() and d.name.upper().startswith("X-") and (pat / d.name).is_dir()
    )
    n = len(controls)
    mism = {}
    leaks = 0
    for cid in controls:
        s_idx, p_idx = src / cid / "index.html", pat / cid / "index.html"
        s_int, p_int = src / cid / "internal.html", pat / cid / "internal.html"
        for m in PUBLIC_MARKERS:
            if has(s_idx, m) != has(p_idx, m):
                mism.setdefault(m, []).append(cid)
        for m in INTERNAL_MARKERS:
            if has(s_int, m) != has(p_int, m):
                mism.setdefault(m, []).append(cid)
        for m in PUBLIC_FORBIDDEN:
            if has(s_idx, m):
                leaks += 1
                break
    print(f"[diff-source-vs-patched] compared {n} controls in both trees")
    print(f"  public MS2-ion leaks in from-source: {leaks}")
    if not mism:
        print("  PARITY OK: all confidence display elements match on every control")
    else:
        for m, ids in sorted(mism.items()):
            print(f"  MISMATCH [{m}]: {len(ids)}/{n}  e.g. {ids[:5]}")
    sys.exit(1 if (mism or leaks) else 0)


if __name__ == "__main__":
    main()
