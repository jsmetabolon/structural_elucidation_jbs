#!/usr/bin/env python3
"""Generate vector (SVG) figures for the manuscript from the audit outputs.

Pure standard library (no matplotlib in this environment). Reads the pathway +
presentation audit JSONs and emits publication-style SVG bar charts plus a
hand-authored workflow schematic into ``docs/paper/figures/``.

Usage:
    python3 tools/paper_figures.py [pathway_summary.json] [presentation_summary.json]
"""
# ruff: noqa: E501  (inline SVG markup is intentionally long)
from __future__ import annotations

import json
import sys
from pathlib import Path

OUT = Path("docs/paper/figures")
INK = "#1b1f24"
MUTED = "#5b6570"
BARS = ["#3b6ea5", "#4a8f6d", "#b06a3b", "#8a6ea5", "#a53b52", "#6a737d"]
GRID = "#d7dce1"


def _esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def bar_chart(path: Path, title: str, cats, vals, *, ymax=None, ylabel="", value_fmt="{:.0f}",
              subtitle="", n_labels=None, ref=None, ref_label=""):
    """Vertical bar chart -> SVG. vals in data units; ymax defaults to max*1.15."""
    W, H = 720, 440
    ml, mr, mt, mb = 70, 24, 74, 96
    pw, ph = W - ml - mr, H - mt - mb
    ymax = ymax if ymax is not None else (max(vals) * 1.15 or 1)
    n = len(cats)
    slot = pw / n
    bw = slot * 0.6
    s = []
    s.append(f"<svg xmlns='http://www.w3.org/2000/svg' width='{W}' height='{H}' viewBox='0 0 {W} {H}' font-family='Helvetica, Arial, sans-serif'>")
    s.append(f"<rect width='{W}' height='{H}' fill='white'/>")
    s.append(f"<text x='{ml}' y='30' font-size='18' font-weight='bold' fill='{INK}'>{_esc(title)}</text>")
    if subtitle:
        s.append(f"<text x='{ml}' y='50' font-size='12' fill='{MUTED}'>{_esc(subtitle)}</text>")
    # y gridlines + ticks (5 divisions)
    for i in range(6):
        yv = ymax * i / 5
        y = mt + ph - (yv / ymax) * ph
        s.append(f"<line x1='{ml}' y1='{y:.1f}' x2='{ml+pw}' y2='{y:.1f}' stroke='{GRID}' stroke-width='1'/>")
        s.append(f"<text x='{ml-8}' y='{y+4:.1f}' font-size='11' fill='{MUTED}' text-anchor='end'>{value_fmt.format(yv)}</text>")
    if ylabel:
        s.append(f"<text x='16' y='{mt+ph/2:.1f}' font-size='12' fill='{MUTED}' text-anchor='middle' transform='rotate(-90 16 {mt+ph/2:.1f})'>{_esc(ylabel)}</text>")
    # reference line
    if ref is not None:
        yr = mt + ph - (ref / ymax) * ph
        s.append(f"<line x1='{ml}' y1='{yr:.1f}' x2='{ml+pw}' y2='{yr:.1f}' stroke='{BARS[4]}' stroke-width='1.5' stroke-dasharray='6 4'/>")
        s.append(f"<text x='{ml+pw}' y='{yr-6:.1f}' font-size='11' fill='{BARS[4]}' text-anchor='end'>{_esc(ref_label)}</text>")
    for i, (c, v) in enumerate(zip(cats, vals)):
        x = ml + slot * i + (slot - bw) / 2
        bh = (v / ymax) * ph
        y = mt + ph - bh
        s.append(f"<rect x='{x:.1f}' y='{y:.1f}' width='{bw:.1f}' height='{bh:.1f}' fill='{BARS[i % len(BARS)]}' rx='2'/>")
        s.append(f"<text x='{x+bw/2:.1f}' y='{y-6:.1f}' font-size='12' font-weight='bold' fill='{INK}' text-anchor='middle'>{value_fmt.format(v)}</text>")
        # category labels (wrap on space)
        words = str(c).split(" ")
        ty = mt + ph + 18
        line = ""
        lines = []
        for w in words:
            if len(line + " " + w) > 14:
                lines.append(line.strip())
                line = w
            else:
                line += " " + w
        lines.append(line.strip())
        for j, ln in enumerate(lines[:3]):
            s.append(f"<text x='{x+bw/2:.1f}' y='{ty+j*13:.1f}' font-size='11' fill='{INK}' text-anchor='middle'>{_esc(ln)}</text>")
        if n_labels:
            s.append(f"<text x='{x+bw/2:.1f}' y='{ty+len(lines)*13:.1f}' font-size='10' fill='{MUTED}' text-anchor='middle'>n={n_labels[i]}</text>")
    s.append(f"<line x1='{ml}' y1='{mt+ph}' x2='{ml+pw}' y2='{mt+ph}' stroke='{INK}' stroke-width='1.2'/>")
    s.append("</svg>")
    path.write_text("\n".join(s))


def workflow_svg(path: Path):
    W, H = 900, 300
    boxes = [
        ("Raw LC-MS/MS", "features"),
        ("Discovery", "MS2 feature\ndetection"),
        ("De novo\nstructure prediction", "MS2->SMILES\ngenerative"),
        ("Candidate rerank", "library + mass +\nspectral evidence"),
        ("Pathway consensus", "neighbor vote of\nMetabolon taxonomy"),
        ("Report", "confidence + MSI\nlevel + index"),
    ]
    n = len(boxes)
    bw, bh, gap = 118, 76, (W - 40 - 118 * n) / (n - 1)
    s = [f"<svg xmlns='http://www.w3.org/2000/svg' width='{W}' height='{H}' viewBox='0 0 {W} {H}' font-family='Helvetica, Arial, sans-serif'>"]
    s.append(f"<rect width='{W}' height='{H}' fill='white'/>")
    s.append(f"<text x='20' y='34' font-size='18' font-weight='bold' fill='{INK}'>Structural-elucidation workflow</text>")
    y = 120
    for i, (t, sub) in enumerate(boxes):
        x = 20 + i * (bw + gap)
        s.append(f"<rect x='{x:.1f}' y='{y}' width='{bw}' height='{bh}' rx='8' fill='#eef2f6' stroke='{BARS[0]}' stroke-width='1.5'/>")
        for j, ln in enumerate(t.split("\n")):
            s.append(f"<text x='{x+bw/2:.1f}' y='{y+22+j*15:.0f}' font-size='12.5' font-weight='bold' fill='{INK}' text-anchor='middle'>{_esc(ln)}</text>")
        base = y + 22 + len(t.split(chr(10))) * 15 + 2
        for j, ln in enumerate(sub.split("\n")):
            s.append(f"<text x='{x+bw/2:.1f}' y='{base+j*13:.0f}' font-size='10.5' fill='{MUTED}' text-anchor='middle'>{_esc(ln)}</text>")
        if i < n - 1:
            ax = x + bw
            s.append(f"<line x1='{ax:.1f}' y1='{y+bh/2}' x2='{ax+gap:.1f}' y2='{y+bh/2}' stroke='{MUTED}' stroke-width='1.5' marker-end='url(#a)'/>")
    s.insert(1, "<defs><marker id='a' markerWidth='8' markerHeight='8' refX='6' refY='3' orient='auto'><path d='M0,0 L6,3 L0,6 Z' fill='#5b6570'/></marker></defs>")
    s.append(f"<text x='20' y='{H-20}' font-size='11' fill='{MUTED}'>3,897 X-compound features across 4 chromatographic libraries (2 positive, 2 negative mode).</text>")
    s.append("</svg>")
    path.write_text("\n".join(s))


def main() -> None:
    pathway = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("_audit_pull/reports/audit_pathway_summary.json")
    present = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("output/_presentation_audit/audit_presentation_summary.json")
    pw = json.loads(pathway.read_text())
    pr = json.loads(present.read_text())
    OUT.mkdir(parents=True, exist_ok=True)

    workflow_svg(OUT / "fig1_workflow.svg")

    sim = pw["accuracy"]["super_accuracy_by_similarity_bin"]
    order = ["[0.0,0.3)", "[0.3,0.5)", "[0.5,0.7)", "[0.7,1.0)"]
    bar_chart(
        OUT / "fig2_pathway_by_similarity.svg",
        "Pathway (super) accuracy vs. structural similarity",
        [k for k in order], [sim[k]["accuracy"] * 100 for k in order],
        ymax=100, ylabel="Super-pathway accuracy (%)",
        subtitle="known_X evaluation; dashed = reference-DB exact-structure accuracy",
        n_labels=[sim[k]["n"] for k in order],
        ref=pw["reference_db_evaluation"]["super_accuracy_vs_db"] * 100,
        ref_label="reference-DB 94.5%",
    )

    tier = pw["accuracy"]["super_accuracy_by_tier"]
    torder = ["high", "medium", "low"]
    bar_chart(
        OUT / "fig3_pathway_by_tier.svg",
        "Pathway (super) accuracy by confidence tier",
        ["high", "medium", "low"], [tier[t]["accuracy"] * 100 for t in torder],
        ymax=100, ylabel="Super-pathway accuracy (%)",
        subtitle="known_X evaluation; accuracy is monotonic in the reported tier",
        n_labels=[tier[t]["n"] for t in torder],
    )

    msi = pr["msi_level_distribution"]
    bar_chart(
        OUT / "fig4_msi_distribution.svg",
        "MSI identification level (all 3,897 reports)",
        ["L1 confirmed", "L2 probable", "L3 tentative", "L4 formula", "L5 mass"],
        [msi[k] for k in ["1", "2", "3", "4", "5"]],
        ylabel="Reports", value_fmt="{:.0f}",
        subtitle="Schymanski/MSI level derived from the evidence axes",
    )

    td = pr["tier_distribution"]
    corder = ["Confirmed", "Moderate", "Low", "Tentative"]
    bar_chart(
        OUT / "fig5_confidence_tiers.svg",
        "Identification-confidence tier (all 3,897 reports)",
        corder, [td.get(c, 0) for c in corder],
        ylabel="Reports", value_fmt="{:.0f}",
        subtitle="calibrated P(top structure correct); no compound reaches High (honest by design)",
    )
    print(f"[paper-figures] wrote 5 SVGs to {OUT}")


if __name__ == "__main__":
    main()
