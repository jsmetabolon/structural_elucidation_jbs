#!/usr/bin/env python3
"""De-crowd the consensus atom-confidence SVG in already-generated reports.

The crowding metric (tools/audit_reports.py) is
``overlap_ratio = 2*radius / min_center_distance``; >1 means atom circles overlap.
Because the metric is scale-invariant, the effective fix is to shrink the atom
glyphs (circle + label) relative to the fixed atom positions so circles no longer
overlap. This edits only the ``<g class="consensus-atom">`` groups inside the
``consensus-box`` (bonds and positions are untouched), in every HTML variant.
Idempotent: once overlap_ratio <= 1 a report is skipped.

Note: for a few molecules with a degenerate 2D layout (atoms collapsed onto each
other, e.g. X-24210) this yields tiny-but-non-overlapping glyphs; a true fix
needs an RDKit re-layout (pipeline/rendering work).

Usage:
    python3 tools/patch_consensus_crowding.py <reports_root> [--target 0.92] [--dry-run]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

CONSENSUS_BOX_RE = re.compile(r'(<div class="consensus-box">)(.*?)(</div>)', re.S)
TRANSLATE_RE = re.compile(r"translate\(([-\d.]+),([-\d.]+)\)")
CIRCLE_R_RE = re.compile(r'<circle r="([\d.]+)"')
ATOM_GROUP_RE = re.compile(r'<g class="consensus-atom".*?</g>', re.S)
ATTR_RE = re.compile(r'(\b(?:r|stroke-width|font-size|y)=")(-?[\d.]+)(")')


def _overlap_ratio(block: str) -> tuple[float, float, float]:
    centers = [(float(x), float(y)) for x, y in TRANSLATE_RE.findall(block)]
    radii = [float(r) for r in CIRCLE_R_RE.findall(block)]
    if len(centers) < 2 or not radii:
        return 0.0, 0.0, 0.0
    radius = radii[0]
    min_d = min(
        ((centers[i][0] - centers[j][0]) ** 2 + (centers[i][1] - centers[j][1]) ** 2) ** 0.5
        for i in range(len(centers)) for j in range(i + 1, len(centers))
    )
    ratio = (2 * radius) / min_d if min_d > 0 else 999.0
    return ratio, radius, min_d


def _scale_group(group: str, f: float) -> str:
    return ATTR_RE.sub(lambda m: f"{m.group(1)}{float(m.group(2)) * f:.2f}{m.group(3)}", group)


def patch_block(block: str, target: float) -> tuple[str, bool]:
    ratio, _r, _d = _overlap_ratio(block)
    if ratio <= 1.0:
        return block, False
    f = target / ratio
    return ATOM_GROUP_RE.sub(lambda m: _scale_group(m.group(0), f), block), True


def patch_report(report_dir: Path, target: float, apply: bool) -> str:
    changed = 0
    for hf in sorted(report_dir.glob("*.html")):
        text = hf.read_text()
        if '<div class="consensus-box">' not in text:
            continue
        new_text, did = "", False

        def repl(m):
            nonlocal did
            body, changed_flag = patch_block(m.group(2), target)
            did = did or changed_flag
            return m.group(1) + body + m.group(3)

        new_text = CONSENSUS_BOX_RE.sub(repl, text)
        if did and new_text != text:
            if apply:
                hf.write_text(new_text)
            changed += 1
    return "patched" if changed else "ok"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("reports_root", type=Path)
    ap.add_argument("--target", type=float, default=0.92)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    counts: dict[str, int] = {}
    for d in sorted(p for p in args.reports_root.iterdir() if p.is_dir()):
        st = patch_report(d, args.target, apply=not args.dry_run)
        counts[st] = counts.get(st, 0) + 1
    print(f"[consensus] {'(dry-run) ' if args.dry_run else ''}{counts}")


if __name__ == "__main__":
    main()
