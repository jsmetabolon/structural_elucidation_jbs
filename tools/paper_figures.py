#!/usr/bin/env python3
"""Generate publication SVG figures from committed audit JSON artifacts.

Pure standard library: no matplotlib dependency.  Every accuracy panel includes
the evaluation cohort, denominator, and Wilson 95% confidence interval.

Colour, type, and chart furniture come from ``tools/brand.py`` rather than being
chosen per figure, so the manuscript and the public brief cannot be restyled
independently of each other.

Usage:
    python3 tools/paper_figures.py \
        --structure docs/paper/data/structure/audit_structure_summary.json \
        --pathway docs/paper/data/pathway/audit_pathway_summary.json \
        --loo docs/paper/data/pathway/audit_pathway_loo.json \
        --presentation output/_presentation_audit/audit_presentation_summary.json \
        --claims docs/paper/supplement/baseline_claims.json \
        --out-dir docs/paper/figures
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tools.brand import (  # noqa: E402
    AXIS_LABEL,
    AXIS_LABEL_SIZE,
    BAR_PRIMARY,
    BAR_SECONDARY,
    DATA_LABEL,
    DATA_LABEL_SIZE,
    ERROR_BAR,
    FONT_STACK,
    FRAME,
    GRID,
    HAIRLINE,
    INK,
    METHOD_GREY,
    MUTED,
    SERIES,
    WHITE,
    on_fill,
)


def esc(value: Any) -> str:
    return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(
    x: float,
    y: float,
    value: Any,
    *,
    size: float = 11,
    anchor: str = "middle",
    bold: bool = False,
    color: str = INK,
) -> str:
    weight = " font-weight='bold'" if bold else ""
    return (
        f"<text x='{x:.1f}' y='{y:.1f}' font-size='{size:g}' text-anchor='{anchor}' "
        f"fill='{color}'{weight}>{esc(value)}</text>"
    )


def canvas(width: int, height: int) -> list[str]:
    """Open a framed white canvas set in the brand typeface."""

    return [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
        f"viewBox='0 0 {width} {height}' font-family='{FONT_STACK}'>",
        f"<rect width='{width}' height='{height}' fill='{WHITE}'/>",
        f"<rect x='0.5' y='0.5' width='{width - 1}' height='{height - 1}' fill='none' "
        f"stroke='{FRAME}' stroke-width='{HAIRLINE}'/>",
    ]


def axis_text(
    x: float,
    y: float,
    value: Any,
    *,
    anchor: str = "middle",
    size: float = AXIS_LABEL_SIZE,
) -> str:
    """An axis label the way the workbook sets one: small, bold, and in caps.

    The size is adjustable because the brief draws on a narrower canvas that is
    printed to the same width, so holding the pixel size fixed across both
    documents would not hold the printed size fixed. The treatment is what
    carries across; the measurement follows the canvas.
    """

    return text(x, y, str(value).upper(), size=size, anchor=anchor, bold=True, color=AXIS_LABEL)


def value_text(
    x: float,
    y: float,
    value: Any,
    *,
    anchor: str = "middle",
    color: str = DATA_LABEL,
    size: float = DATA_LABEL_SIZE,
) -> str:
    """A data label: same size as the axis labels but regular, so the number reads first."""

    return text(x, y, value, size=size, anchor=anchor, color=color)


def rule(x1: float, y1: float, x2: float, y2: float, *, color: str = GRID, width: float = HAIRLINE) -> str:
    return (
        f"<line x1='{x1:.1f}' y1='{y1:.1f}' x2='{x2:.1f}' y2='{y2:.1f}' "
        f"stroke='{color}' stroke-width='{width:g}'/>"
    )


def error_bar(x: float, low_y: float, high_y: float, *, color: str = ERROR_BAR, cap: float = 5.0) -> list[str]:
    """A capped interval whisker, drawn in grey so it never competes with the bar."""

    return [
        rule(x, high_y, x, low_y, color=color),
        rule(x - cap, high_y, x + cap, high_y, color=color),
        rule(x - cap, low_y, x + cap, low_y, color=color),
    ]


def split_label(label: str, width: int = 15) -> list[str]:
    words = label.split()
    lines, line = [], ""
    for word in words:
        if line and len(line) + 1 + len(word) > width:
            lines.append(line)
            line = word
        else:
            line = f"{line} {word}".strip()
    if line:
        lines.append(line)
    return lines[:3]


def bar_chart(
    path: Path,
    *,
    title: str,
    subtitle: str,
    categories: list[str],
    series: list[dict[str, Any]],
    ymax: float = 100.0,
    ylabel: str = "Accuracy (%)",
) -> None:
    width, height = 900, 520
    left, right, top, bottom = 80, 30, 90, 120
    plot_w, plot_h = width - left - right, height - top - bottom
    group_w = plot_w / len(categories)
    bar_w = min(52.0, group_w * 0.72 / max(1, len(series)))
    svg = canvas(width, height)
    svg += [
        text(left, 30, title, size=18, anchor="start", bold=True),
        text(left, 52, subtitle, size=11, anchor="start", color=MUTED),
    ]
    for tick in range(0, 101, 20):
        y = top + plot_h - (tick / ymax) * plot_h
        svg.append(rule(left, y, left + plot_w, y))
        svg.append(axis_text(left - 10, y + 4, tick, anchor="end"))
    svg.append(rule(left, top, left, top + plot_h))
    svg.append(
        f"<text x='20' y='{top + plot_h / 2:.1f}' font-size='12' fill='{MUTED}' "
        f"text-anchor='middle' transform='rotate(-90 20 {top + plot_h / 2:.1f})'>{esc(ylabel)}</text>"
    )
    for category_index, category in enumerate(categories):
        center = left + group_w * (category_index + 0.5)
        total_w = bar_w * len(series)
        for series_index, item in enumerate(series):
            metric = item["metrics"][category_index]
            value = float(metric["accuracy"]) * 100.0 if metric.get("accuracy") is not None else 0.0
            x = center - total_w / 2 + series_index * bar_w
            y = top + plot_h - value / ymax * plot_h
            svg.append(
                f"<rect x='{x:.1f}' y='{y:.1f}' width='{bar_w - 3:.1f}' "
                f"height='{value / ymax * plot_h:.1f}' fill='{item['color']}'/>"
            )
            center_x = x + (bar_w - 3) / 2
            label_top = y
            ci = metric.get("ci95") or [None, None]
            if ci[0] is not None:
                low_y = top + plot_h - float(ci[0]) * 100.0 / ymax * plot_h
                high_y = top + plot_h - float(ci[1]) * 100.0 / ymax * plot_h
                svg += error_bar(center_x, low_y, high_y)
                label_top = min(y, high_y)
            svg.append(value_text(center_x, label_top - 8, f"{value:.1f}"))
            if metric.get("n") is not None:
                bar_px = top + plot_h - y
                if bar_px >= 26:
                    # Sat at the foot of the bar rather than under its top edge,
                    # which is where the confidence whisker crosses it.
                    svg.append(
                        value_text(center_x, top + plot_h - 9, f"n={metric['n']}", color=on_fill(item["color"]))
                    )
                else:
                    svg.append(value_text(center_x, top + plot_h - 5, f"n={metric['n']}", color=AXIS_LABEL))
        label_y = top + plot_h + 20
        for line_index, line in enumerate(split_label(category)):
            svg.append(axis_text(center, label_y + 13 * line_index, line))
    legend_x = left + plot_w - 180 * len(series)
    for index, item in enumerate(series):
        x = legend_x + 180 * index
        svg.append(f"<rect x='{x:.1f}' y='65' width='12' height='12' fill='{item['color']}'/>")
        svg.append(text(x + 18, 75, item["name"], size=10, anchor="start", color=MUTED))
    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n")


def reliability_chart(path: Path, full: dict[str, Any], ablated: dict[str, Any]) -> None:
    width, height = 680, 560
    left, top, size = 85, 75, 400
    bin_sizes = [int(row["n"]) for row in full["reliability_bins"]] or [0]
    min_bin, max_bin = min(bin_sizes), max(bin_sizes)
    svg = canvas(width, height)
    svg += [
        text(left, 28, "Leave-one-out identification-confidence reliability", size=18, anchor="start", bold=True),
        text(
            left,
            50,
            f"n={full['n']}; full Brier={full['brier_score']:.3f}; no-pathway Brier={ablated['brier_score']:.3f}"
            f"; {min_bin}\u2013{max_bin} controls per bin",
            size=11,
            anchor="start",
            color=MUTED,
        ),
    ]
    for tick in range(0, 6):
        value = tick / 5
        x = left + value * size
        y = top + size - value * size
        svg += [
            rule(x, top, x, top + size),
            rule(left, y, left + size, y),
            axis_text(x, top + size + 20, f"{value:.1f}"),
            axis_text(left - 12, y + 4, f"{value:.1f}", anchor="end"),
        ]
    svg.append(
        f"<line x1='{left}' y1='{top + size}' x2='{left + size}' y2='{top}' "
        f"stroke='{METHOD_GREY}' stroke-dasharray='5 5'/>"
    )
    for series_index, (name, result, color) in enumerate(
        (("Full model", full, SERIES[0]), ("Without pathway axis", ablated, SERIES[1]))
    ):
        points = []
        for row in result["reliability_bins"]:
            x = left + float(row["mean_predicted"]) * size
            y = top + size - float(row["observed"]) * size
            low_y = top + size - float(row["observed_ci95"][0]) * size
            high_y = top + size - float(row["observed_ci95"][1]) * size
            points.append((x, y))
            svg += [
                f"<line x1='{x:.1f}' y1='{high_y:.1f}' x2='{x:.1f}' y2='{low_y:.1f}' "
                f"stroke='{color}' stroke-width='1.5'/>",
                f"<circle cx='{x:.1f}' cy='{y:.1f}' r='3.5' fill='{color}'/>",
            ]
        if len(points) > 1:
            coords = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
            svg.append(f"<polyline points='{coords}' fill='none' stroke='{color}' stroke-width='2'/>")
        legend_y = 105 + series_index * 24
        svg.append(f"<circle cx='525' cy='{legend_y}' r='3.5' fill='{color}'/>")
        svg.append(text(538, legend_y + 4, name, size=10, anchor="start", color=MUTED))
    svg += [
        text(left + size / 2, top + size + 48, "Mean predicted probability", size=12),
        (
            f"<text x='24' y='{top + size / 2:.1f}' font-size='12' fill='{MUTED}' "
            f"text-anchor='middle' transform='rotate(-90 24 {top + size / 2:.1f})'>Observed fraction correct</text>"
        ),
        "</svg>",
    ]
    path.write_text("\n".join(svg) + "\n")


def selective_prediction_chart(path: Path, selective: dict[str, Any]) -> None:
    """Accuracy against the fraction of each cohort retained, most confident first.

    Plotted on one axis because the contrast is the point: two of the three signals
    rank usefully and one does not.
    """

    width, height = 760, 600
    left, top, plot_w, plot_h = 85, 100, 620, 380
    series = (
        ("Pathway super-pathway (by similarity)", selective["pathway_super"], SERIES[0]),
        ("Pathway subpathway (by similarity)", selective["pathway_sub"], SERIES[1]),
        ("Structure rank 1 (by confidence)", selective["structure_rank1"], SERIES[2]),
    )
    grid = [float(value) for value in selective["coverage_grid"]]
    span = max(grid) - min(grid) or 1.0

    def x_of(coverage: float) -> float:
        return left + (coverage - min(grid)) / span * plot_w

    def y_of(fraction: float) -> float:
        return top + plot_h - fraction * plot_h

    svg = canvas(width, height)
    svg += [
        text(left, 30, "Accuracy against the fraction of calls retained", size=18, anchor="start", bold=True),
        text(
            left,
            52,
            "Each cohort sorted most to least confident and scored over cumulative prefixes.",
            size=11,
            anchor="start",
            color=MUTED,
        ),
        text(
            left,
            70,
            "A falling curve means confidence ranks usefully; a flat curve means it does not.",
            size=11,
            anchor="start",
            color=MUTED,
        ),
    ]
    for tick in range(0, 101, 20):
        y = y_of(tick / 100)
        svg += [
            rule(left, y, left + plot_w, y),
            axis_text(left - 12, y + 4, tick, anchor="end"),
        ]
    svg.append(rule(left, top, left, top + plot_h))
    for coverage in grid:
        x = x_of(coverage)
        svg.append(axis_text(x, top + plot_h + 22, f"{coverage:.0%}"))
    for index, (name, block, color) in enumerate(series):
        points = []
        for row in block["points"]:
            if row["accuracy"] is None:
                continue
            x = x_of(float(row["coverage"]))
            y = y_of(float(row["accuracy"]))
            points.append((x, y))
            svg += [
                f"<line x1='{x:.1f}' y1='{y_of(float(row['ci95'][1])):.1f}' x2='{x:.1f}' "
                f"y2='{y_of(float(row['ci95'][0])):.1f}' stroke='{color}' stroke-width='1.4'/>",
                f"<circle cx='{x:.1f}' cy='{y:.1f}' r='3.5' fill='{color}'/>",
            ]
        if len(points) > 1:
            coords = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
            svg.append(f"<polyline points='{coords}' fill='none' stroke='{color}' stroke-width='2'/>")
    for index, (name, block, color) in enumerate(series):
        swatch_y = top + plot_h + 56 + index * 21
        svg += [
            f"<rect x='{left}' y='{swatch_y - 9}' width='11' height='11' fill='{color}'/>",
            text(
                left + 18,
                swatch_y,
                f"{name}, n={block['n']} \u2014 {block['discrimination']:+.0%} "
                "from the most confident fifth to the whole cohort",
                size=11,
                anchor="start",
                color=MUTED,
            ),
        ]
    svg += [
        text(left + plot_w / 2, top + plot_h + 40, "Fraction of the cohort retained", size=12, color=MUTED),
        (
            f"<text x='24' y='{top + plot_h / 2:.1f}' font-size='12' fill='{MUTED}' "
            f"text-anchor='middle' transform='rotate(-90 24 {top + plot_h / 2:.1f})'>"
            "Accuracy or agreement (%)</text>"
        ),
        "</svg>",
    ]
    path.write_text("\n".join(svg) + "\n")


def distribution_chart(path: Path, presentation: dict[str, Any]) -> None:
    width, height = 960, 500
    svg = canvas(width, height)
    svg += [
        text(40, 30, "Identification outputs across all 3,897 reports", size=18, anchor="start", bold=True),
        text(
            40,
            52,
            "Square-root count scale preserves small categories; labels show count and percentage.",
            size=11,
            anchor="start",
            color=MUTED,
        ),
    ]
    panels = [
        (
            "MSI/Schymanski-style reporting level",
            [f"L{i}" for i in range(1, 6)],
            [presentation["msi_level_distribution"].get(str(i), 0) for i in range(1, 6)],
        ),
        (
            "Identification-confidence tier",
            ["Confirmed", "High", "Moderate", "Low", "Tentative"],
            [
                presentation["tier_distribution"].get(key, 0)
                for key in ("Confirmed", "High", "Moderate", "Low", "Tentative")
            ],
        ),
    ]
    total = int(presentation["n_reports"])
    for panel_index, (title, labels, counts) in enumerate(panels):
        x0 = 40 + panel_index * 470
        y0, panel_w, panel_h = 90, 410, 300
        svg.append(text(x0, y0 - 12, title, size=13, anchor="start", bold=True))
        max_root = math.sqrt(max(counts) or 1)
        slot = panel_w / len(labels)
        for index, (label, count) in enumerate(zip(labels, counts)):
            bar_h = math.sqrt(count) / max_root * panel_h * 0.88 if count else 0
            x = x0 + index * slot + 15
            y = y0 + panel_h - bar_h
            svg.append(
                f"<rect x='{x:.1f}' y='{y:.1f}' width='{slot - 25:.1f}' height='{bar_h:.1f}' "
                f"fill='{SERIES[index % len(SERIES)]}'/>"
            )
            svg.append(value_text(x + (slot - 25) / 2, y - 7, f"{count} ({100 * count / total:.1f}%)"))
            svg.append(axis_text(x + (slot - 25) / 2, y0 + panel_h + 18, label))
        svg.append(rule(x0, y0 + panel_h, x0 + panel_w, y0 + panel_h))
    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n")


def mass_mismatch_chart(path: Path, claims: dict[str, dict[str, Any]]) -> None:
    categories = ["All mismatches", "Recoverable", "No in-list alternative", "MSI-L4 deep audit", "L4 recoverable"]
    values = [
        int(claims["rank1_mass_mismatch_all"]["value"]),
        int(claims["rank1_mass_mismatch_recoverable_all"]["value"]),
        int(claims["rank1_mass_mismatch_all"]["value"])
        - int(claims["rank1_mass_mismatch_recoverable_all"]["value"]),
        int(claims["rank1_mass_mismatch_msi_l4_deep_audit"]["value"]),
        int(claims["rank1_mass_mismatch_recoverable_l4"]["value"]),
    ]
    metrics = [{"n": value, "accuracy": value / max(values), "ci95": [None, None]} for value in values]
    bar_chart(
        path,
        title="Rank-1 mass-mismatch cohorts",
        subtitle=(
            "Counts; recoverable means a lower-ranked candidate has strong mass fit. "
            "This is a projection, not proven accuracy."
        ),
        categories=categories,
        series=[{"name": "Reports", "color": BAR_PRIMARY, "metrics": metrics}],
        ymax=100,
        ylabel="Relative count (% of maximum)",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--pathway", type=Path, required=True)
    parser.add_argument("--loo", type=Path, required=True)
    parser.add_argument("--presentation", type=Path, required=True)
    parser.add_argument("--claims", type=Path, required=True)
    parser.add_argument("--selective", type=Path, default=Path("docs/paper/data/selective_prediction.json"))
    parser.add_argument("--out-dir", type=Path, default=Path("docs/paper/figures"))
    args = parser.parse_args()
    if not args.selective.exists():
        raise SystemExit(f"missing {args.selective}; run tools/analyze_selective_prediction.py first")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    structure = json.loads(args.structure.read_text())
    pathway = json.loads(args.pathway.read_text())
    loo = json.loads(args.loo.read_text())
    presentation = json.loads(args.presentation.read_text())
    claim_rows = json.loads(args.claims.read_text())["claims"]
    claims = {row["claim_id"]: row for row in claim_rows}

    rank = structure["rank_k"]["0.70"]
    subgroup = structure["source_subgroups"]
    structure_metrics = [
        rank["top1"],
        rank["top3"],
        rank["top10"],
        subgroup["exact"]["rank1_tanimoto_ge_0p70"],
        subgroup["predicted"]["rank1_tanimoto_ge_0p70"],
    ]
    bar_chart(
        args.out_dir / "fig2_structure_accuracy.svg",
        title="Structural recovery on truth-verified controls",
        subtitle="Correct = truth Tanimoto >= 0.70; Wilson 95% confidence intervals.",
        categories=["Rank 1", "Top 3", "Top 10", "Exact-source rank 1", "Predicted-source rank 1"],
        series=[{"name": "Structural recovery", "color": BAR_PRIMARY, "metrics": structure_metrics}],
    )
    reliability_chart(
        args.out_dir / "fig3_confidence_reliability.svg",
        structure["calibration_full"],
        structure["calibration_without_pathway"],
    )

    known = pathway["accuracy"]
    ref = pathway["reference_db_evaluation"]
    pathway_categories = [
        "Direct lookup consistency",
        "Neighbor vs known_X",
        "All vs known_X",
        "LOO mechanism",
    ]
    pathway_super = [
        ref["by_assignment_source"]["exact_reference_lookup"]["super"],
        known["by_assignment_source"]["neighbor_consensus"]["super"],
        known["overall"]["super"],
        loo["super"],
    ]
    pathway_sub = [
        ref["by_assignment_source"]["exact_reference_lookup"]["sub"],
        known["by_assignment_source"]["neighbor_consensus"]["sub"],
        known["overall"]["sub"],
        loo["sub"],
    ]
    bar_chart(
        args.out_dir / "fig5_pathway_evaluations.svg",
        title="Pathway results by evaluation context",
        subtitle=(
            "Direct lookup is consistency, not prediction accuracy; known_X is "
            "independent but taxonomy-discordant."
        ),
        categories=pathway_categories,
        series=[
            {"name": "Super-pathway", "color": BAR_PRIMARY, "metrics": pathway_super},
            {"name": "Subpathway", "color": BAR_SECONDARY, "metrics": pathway_sub},
        ],
    )
    distribution_chart(args.out_dir / "fig4_identification_distributions.svg", presentation)

    tiers = known["by_pathway_tier"]
    tier_order = ["high", "medium", "low"]
    bar_chart(
        args.out_dir / "figS1_pathway_tier.svg",
        title="Independent-label agreement by pathway-confidence tier",
        subtitle="known_X cohort; Wilson 95% confidence intervals.",
        categories=tier_order,
        series=[
            {"name": "Super-pathway", "color": BAR_PRIMARY, "metrics": [tiers[key]["super"] for key in tier_order]},
            {"name": "Subpathway", "color": BAR_SECONDARY, "metrics": [tiers[key]["sub"] for key in tier_order]},
        ],
    )
    loo_bins = loo["super_accuracy_by_similarity_bin"]
    bin_order = sorted(loo_bins)
    bar_chart(
        args.out_dir / "figS2_loo_similarity.svg",
        title="Neighbor-vote leave-one-out accuracy by similarity",
        subtitle="Intrinsic mechanism validation; does not include production priors/overrides.",
        categories=bin_order,
        series=[{"name": "Super-pathway", "color": BAR_PRIMARY, "metrics": [loo_bins[key] for key in bin_order]}],
    )
    statuses = known["by_assignment_status"]
    status_order = [key for key in ("strong_similarity_vote", "weak_similarity_backfill") if key in statuses]
    bar_chart(
        args.out_dir / "figS3_pathway_status.svg",
        title="Independent-label agreement by assignment status",
        subtitle="known_X cohort; Wilson 95% confidence intervals.",
        categories=["Strong similarity vote", "Weak-similarity backfill"][: len(status_order)],
        series=[
            {
                "name": "Super-pathway",
                "color": BAR_PRIMARY,
                "metrics": [statuses[key]["super"] for key in status_order],
            },
            {"name": "Subpathway", "color": BAR_SECONDARY, "metrics": [statuses[key]["sub"] for key in status_order]},
        ],
    )
    mass_mismatch_chart(args.out_dir / "figS4_mass_mismatch.svg", claims)
    selective_prediction_chart(
        args.out_dir / "figS5_selective_prediction.svg",
        json.loads(args.selective.read_text()),
    )
    print(f"[paper-figures] wrote publication + supplementary SVGs to {args.out_dir}")


if __name__ == "__main__":
    main()
