#!/usr/bin/env python3
"""Generate publication SVG figures from committed audit JSON artifacts.

Pure standard library: no matplotlib dependency.  Every accuracy panel includes
the evaluation cohort, denominator, and Wilson 95% confidence interval.

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
from pathlib import Path
from typing import Any

INK = "#1f2933"
MUTED = "#52606d"
GRID = "#d9e2ec"
COLORS = ("#2563a6", "#2f855a", "#b45309", "#805ad5", "#b83280", "#64748b")


def esc(value: Any) -> str:
    return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(x: float, y: float, value: Any, *, size: int = 11, anchor: str = "middle", bold: bool = False) -> str:
    weight = " font-weight='bold'" if bold else ""
    return (
        f"<text x='{x:.1f}' y='{y:.1f}' font-size='{size}' text-anchor='{anchor}' "
        f"fill='{INK}'{weight}>{esc(value)}</text>"
    )


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
    svg = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
        f"viewBox='0 0 {width} {height}' font-family='Arial, Helvetica, sans-serif'>",
        f"<rect width='{width}' height='{height}' fill='white'/>",
        text(left, 30, title, size=18, anchor="start", bold=True),
        text(left, 52, subtitle, size=11, anchor="start"),
    ]
    for tick in range(0, 101, 20):
        y = top + plot_h - (tick / ymax) * plot_h
        svg.append(f"<line x1='{left}' y1='{y:.1f}' x2='{left + plot_w}' y2='{y:.1f}' stroke='{GRID}'/>")
        svg.append(text(left - 10, y + 4, tick, anchor="end"))
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
                f"height='{value / ymax * plot_h:.1f}' fill='{item['color']}' rx='2'/>"
            )
            ci = metric.get("ci95") or [None, None]
            if ci[0] is not None:
                low_y = top + plot_h - float(ci[0]) * 100.0 / ymax * plot_h
                high_y = top + plot_h - float(ci[1]) * 100.0 / ymax * plot_h
                mid_x = x + (bar_w - 3) / 2
                svg += [
                    f"<line x1='{mid_x:.1f}' y1='{high_y:.1f}' x2='{mid_x:.1f}' "
                    f"y2='{low_y:.1f}' stroke='{INK}' stroke-width='1.2'/>",
                    f"<line x1='{mid_x - 5:.1f}' y1='{high_y:.1f}' x2='{mid_x + 5:.1f}' "
                    f"y2='{high_y:.1f}' stroke='{INK}'/>",
                    f"<line x1='{mid_x - 5:.1f}' y1='{low_y:.1f}' x2='{mid_x + 5:.1f}' "
                    f"y2='{low_y:.1f}' stroke='{INK}'/>",
                ]
            svg.append(text(x + (bar_w - 3) / 2, y - 7, f"{value:.1f}", size=10, bold=True))
            if metric.get("n") is not None:
                svg.append(text(x + (bar_w - 3) / 2, y + 15, f"n={metric['n']}", size=9))
        label_y = top + plot_h + 20
        for line_index, line in enumerate(split_label(category)):
            svg.append(text(center, label_y + 14 * line_index, line, size=10))
    legend_x = left + plot_w - 180 * len(series)
    for index, item in enumerate(series):
        x = legend_x + 180 * index
        svg.append(f"<rect x='{x:.1f}' y='65' width='12' height='12' fill='{item['color']}'/>")
        svg.append(text(x + 18, 75, item["name"], size=10, anchor="start"))
    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n")


def reliability_chart(path: Path, full: dict[str, Any], ablated: dict[str, Any]) -> None:
    width, height = 680, 560
    left, top, size = 85, 75, 400
    svg = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
        f"viewBox='0 0 {width} {height}' font-family='Arial, Helvetica, sans-serif'>",
        f"<rect width='{width}' height='{height}' fill='white'/>",
        text(left, 28, "Leave-one-out identification-confidence reliability", size=18, anchor="start", bold=True),
        text(
            left,
            50,
            f"n={full['n']}; full Brier={full['brier_score']:.3f}; no-pathway Brier={ablated['brier_score']:.3f}",
            size=11,
            anchor="start",
        ),
    ]
    for tick in range(0, 6):
        value = tick / 5
        x = left + value * size
        y = top + size - value * size
        svg += [
            f"<line x1='{x:.1f}' y1='{top}' x2='{x:.1f}' y2='{top + size}' stroke='{GRID}'/>",
            f"<line x1='{left}' y1='{y:.1f}' x2='{left + size}' y2='{y:.1f}' stroke='{GRID}'/>",
            text(x, top + size + 20, f"{value:.1f}", size=10),
            text(left - 12, y + 4, f"{value:.1f}", size=10, anchor="end"),
        ]
    svg.append(
        f"<line x1='{left}' y1='{top + size}' x2='{left + size}' y2='{top}' "
        f"stroke='{MUTED}' stroke-dasharray='5 5'/>"
    )
    for series_index, (name, result, color) in enumerate(
        (("Full model", full, COLORS[0]), ("Without pathway axis", ablated, COLORS[1]))
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
                f"<circle cx='{x:.1f}' cy='{y:.1f}' r='5' fill='{color}'/>",
                text(x + 8, y - 7, f"n={row['n']}", size=9, anchor="start"),
            ]
        if len(points) > 1:
            coords = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
            svg.append(f"<polyline points='{coords}' fill='none' stroke='{color}' stroke-width='2'/>")
        legend_y = 105 + series_index * 24
        svg.append(f"<circle cx='525' cy='{legend_y}' r='5' fill='{color}'/>")
        svg.append(text(538, legend_y + 4, name, size=10, anchor="start"))
    svg += [
        text(left + size / 2, top + size + 48, "Mean predicted probability", size=12),
        (
            f"<text x='24' y='{top + size / 2:.1f}' font-size='12' fill='{MUTED}' "
            f"text-anchor='middle' transform='rotate(-90 24 {top + size / 2:.1f})'>Observed fraction correct</text>"
        ),
        "</svg>",
    ]
    path.write_text("\n".join(svg) + "\n")


def distribution_chart(path: Path, presentation: dict[str, Any]) -> None:
    width, height = 960, 500
    svg = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
        f"viewBox='0 0 {width} {height}' font-family='Arial, Helvetica, sans-serif'>",
        f"<rect width='{width}' height='{height}' fill='white'/>",
        text(40, 30, "Identification outputs across all 3,897 reports", size=18, anchor="start", bold=True),
        text(
            40,
            52,
            "Square-root count scale preserves small categories; labels show count and percentage.",
            size=11,
            anchor="start",
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
            bar_h = math.sqrt(count) / max_root * panel_h if count else 0
            x = x0 + index * slot + 15
            y = y0 + panel_h - bar_h
            svg.append(
                f"<rect x='{x:.1f}' y='{y:.1f}' width='{slot - 25:.1f}' height='{bar_h:.1f}' "
                f"fill='{COLORS[index % len(COLORS)]}' rx='2'/>"
            )
            svg.append(text(x + (slot - 25) / 2, y - 7, f"{count} ({100 * count / total:.1f}%)", size=9, bold=True))
            svg.append(text(x + (slot - 25) / 2, y0 + panel_h + 18, label, size=10))
        svg.append(f"<line x1='{x0}' y1='{y0 + panel_h}' x2='{x0 + panel_w}' y2='{y0 + panel_h}' stroke='{INK}'/>")
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
        series=[{"name": "Reports", "color": COLORS[2], "metrics": metrics}],
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
    parser.add_argument("--out-dir", type=Path, default=Path("docs/paper/figures"))
    args = parser.parse_args()
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
        series=[{"name": "Structural recovery", "color": COLORS[0], "metrics": structure_metrics}],
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
            {"name": "Super-pathway", "color": COLORS[0], "metrics": pathway_super},
            {"name": "Subpathway", "color": COLORS[1], "metrics": pathway_sub},
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
            {"name": "Super-pathway", "color": COLORS[0], "metrics": [tiers[key]["super"] for key in tier_order]},
            {"name": "Subpathway", "color": COLORS[1], "metrics": [tiers[key]["sub"] for key in tier_order]},
        ],
    )
    loo_bins = loo["super_accuracy_by_similarity_bin"]
    bin_order = sorted(loo_bins)
    bar_chart(
        args.out_dir / "figS2_loo_similarity.svg",
        title="Neighbor-vote leave-one-out accuracy by similarity",
        subtitle="Intrinsic mechanism validation; does not include production priors/overrides.",
        categories=bin_order,
        series=[{"name": "Super-pathway", "color": COLORS[0], "metrics": [loo_bins[key] for key in bin_order]}],
    )
    statuses = known["by_assignment_status"]
    status_order = [key for key in ("strong_similarity_vote", "weak_similarity_backfill") if key in statuses]
    bar_chart(
        args.out_dir / "figS3_pathway_status.svg",
        title="Independent-label agreement by assignment status",
        subtitle="known_X cohort; Wilson 95% confidence intervals.",
        categories=["Strong similarity vote", "Weak-similarity backfill"][: len(status_order)],
        series=[
            {"name": "Super-pathway", "color": COLORS[0], "metrics": [statuses[key]["super"] for key in status_order]},
            {"name": "Subpathway", "color": COLORS[1], "metrics": [statuses[key]["sub"] for key in status_order]},
        ],
    )
    mass_mismatch_chart(args.out_dir / "figS4_mass_mismatch.svg", claims)
    print(f"[paper-figures] wrote publication + supplementary SVGs to {args.out_dir}")


if __name__ == "__main__":
    main()
