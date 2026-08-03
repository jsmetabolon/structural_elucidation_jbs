#!/usr/bin/env python3
"""Generate the public capability brief's SVG figures from committed audit JSON.

The brief is customer-facing, so these figures deliberately carry no thresholds,
model parameters, or internal cohort identifiers. They are regenerated from the
same audited artifacts the manuscript uses, which keeps the two documents from
drifting apart.

Canvases are narrower than the manuscript figures in ``tools/paper_figures.py``.
Both are embedded at the same printed width, so a smaller canvas renders larger
type, which matters more for a brief that is read at a glance than for a figure
read alongside a methods section.

The workflow diagram is hand-authored rather than generated; see
``docs/public_brief/figures/fig_workflow_public.svg``.

Usage:
    python3 tools/public_brief_figures.py \
        --structure docs/paper/data/structure/audit_structure_summary.json \
        --presentation docs/paper/data/presentation/audit_presentation_summary.json \
        --out-dir docs/public_brief/figures
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tools.paper_figures import COLORS, GRID, INK, MUTED, WHITE, text  # noqa: E402

# Levels that present a specific structure to the customer, versus those that withhold one.
STRUCTURE_SHOWN_LEVELS = ("1", "2", "3")
STRUCTURE_WITHHELD_LEVELS = ("4", "5")

PRIMARY = COLORS[0]
POSITIVE = COLORS[1]


def recovery_chart(path: Path, structure: dict[str, Any]) -> None:
    """Cumulative structural recovery at rank 1, top 3, and top 10, with Wilson intervals."""

    rank_k = structure["rank_k"]["0.70"]
    bars = (
        ("Top-ranked proposal", rank_k["top1"]),
        ("Within top 3", rank_k["top3"]),
        ("Within top 10", rank_k["top10"]),
    )

    width, height = 660, 340
    left, top = 62, 86
    plot_w, plot_h = width - left - 30, height - top - 62
    slot = plot_w / len(bars)
    bar_w = 92.0

    svg = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
        f"viewBox='0 0 {width} {height}' font-family='Arial, Helvetica, sans-serif'>",
        f"<rect width='{width}' height='{height}' fill='white'/>",
        text(left - 22, 32, "How often the correct structure is found", size=17, anchor="start", bold=True),
        text(
            left - 22,
            54,
            f"Internal controls with independently known structures (n={structure['n_truth_verified']})",
            size=12,
            anchor="start",
            color=MUTED,
        ),
    ]

    for tick in range(0, 101, 25):
        y = top + plot_h - tick / 100.0 * plot_h
        svg.append(f"<line x1='{left}' y1='{y:.1f}' x2='{left + plot_w}' y2='{y:.1f}' stroke='{GRID}'/>")
        svg.append(text(left - 10, y + 5, f"{tick}%", size=11, anchor="end", color=MUTED))

    for position, (label, metric) in enumerate(bars):
        value = float(metric["accuracy"]) * 100.0
        centre = left + slot * (position + 0.5)
        y = top + plot_h - value / 100.0 * plot_h
        svg.append(
            f"<rect x='{centre - bar_w / 2:.1f}' y='{y:.1f}' width='{bar_w:.1f}' "
            f"height='{plot_h - (y - top):.1f}' fill='{PRIMARY}' rx='3'/>"
        )

        low, high = metric["ci95"]
        low_y = top + plot_h - float(low) * 100.0 / 100.0 * plot_h
        high_y = top + plot_h - float(high) * 100.0 / 100.0 * plot_h
        svg += [
            f"<line x1='{centre:.1f}' y1='{high_y:.1f}' x2='{centre:.1f}' y2='{low_y:.1f}' "
            f"stroke='{INK}' stroke-width='1.4'/>",
            f"<line x1='{centre - 7:.1f}' y1='{high_y:.1f}' x2='{centre + 7:.1f}' y2='{high_y:.1f}' stroke='{INK}'/>",
            f"<line x1='{centre - 7:.1f}' y1='{low_y:.1f}' x2='{centre + 7:.1f}' y2='{low_y:.1f}' stroke='{INK}'/>",
        ]
        svg.append(text(centre, high_y - 10, f"{value:.1f}%", size=17, bold=True))
        svg.append(text(centre, y + 24, f"{metric['correct']} of {metric['n']}", size=12, color=WHITE))
        svg.append(text(centre, top + plot_h + 24, label, size=13, bold=True))

    svg.append(
        text(
            left - 22,
            height - 14,
            "Bars show recovery; whiskers show the 95% confidence interval.",
            size=11,
            anchor="start",
            color=MUTED,
        )
    )
    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n")


def coverage_chart(path: Path, presentation: dict[str, Any]) -> None:
    """Split of the delivered set into reports that show a structure and reports that withhold one."""

    distribution = presentation["msi_level_distribution"]
    shown = sum(int(distribution[level]) for level in STRUCTURE_SHOWN_LEVELS)
    withheld = sum(int(distribution[level]) for level in STRUCTURE_WITHHELD_LEVELS)
    total = shown + withheld

    width, height = 660, 250
    left, top = 40, 96
    bar_w, bar_h = width - 2 * left, 78
    shown_w = bar_w * shown / total

    svg = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
        f"viewBox='0 0 {width} {height}' font-family='Arial, Helvetica, sans-serif'>",
        f"<rect width='{width}' height='{height}' fill='white'/>",
        text(left, 34, "Every feature is reported", size=17, anchor="start", bold=True),
        text(
            left,
            56,
            f"All {total:,} delivered features receive a report and biochemical pathway context",
            size=12,
            anchor="start",
            color=MUTED,
        ),
        text(left, 84, "Of those reports:", size=12, anchor="start", color=MUTED),
    ]

    segments = (
        (left, shown_w, POSITIVE, "Structure proposed", shown),
        (left + shown_w, bar_w - shown_w, MUTED, "Structure withheld", withheld),
    )
    for x, segment_w, color, label, count in segments:
        share = count / total * 100.0
        centre = x + segment_w / 2
        svg.append(f"<rect x='{x:.1f}' y='{top}' width='{segment_w:.1f}' height='{bar_h}' fill='{color}'/>")
        svg.append(text(centre, top + 36, f"{share:.1f}%", size=24, bold=True, color=WHITE))
        svg.append(text(centre, top + 60, f"{count:,} features", size=12, color=WHITE))
        svg.append(text(centre, top + bar_h + 24, label, size=13, bold=True))

    svg.append(f"<rect x='{left}' y='{top}' width='{bar_w}' height='{bar_h}' fill='none' stroke='{GRID}'/>")
    svg.append(
        text(
            left,
            height - 14,
            "A structure is shown only where the supporting evidence justifies it.",
            size=12,
            anchor="start",
            color=INK,
        )
    )
    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--structure",
        type=Path,
        default=Path("docs/paper/data/structure/audit_structure_summary.json"),
    )
    parser.add_argument(
        "--presentation",
        type=Path,
        default=Path("docs/paper/data/presentation/audit_presentation_summary.json"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("docs/public_brief/figures"))
    args = parser.parse_args()

    for path in (args.structure, args.presentation):
        if not path.exists():
            raise SystemExit(f"audit artifact not found: {path}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    structure = json.loads(args.structure.read_text())
    presentation = json.loads(args.presentation.read_text())

    recovery_chart(args.out_dir / "fig_recovery_public.svg", structure)
    coverage_chart(args.out_dir / "fig_coverage_public.svg", presentation)
    print(f"[public-brief-figures] wrote 2 figures to {args.out_dir}")


if __name__ == "__main__":
    main()
