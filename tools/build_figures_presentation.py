#!/usr/bin/env python3
"""Build a template-styled PPTX where every figure is a native, data-editable object.

Data-driven figures become native PowerPoint charts with an embedded Excel
worksheet (right-click -> Edit Data). Fig 1 is rebuilt as editable shapes and
Fig 3 as a native scatter. Confidence intervals are intentionally omitted;
value and n annotations are retained (n via the category label).

Usage:
    python3 tools/build_figures_presentation.py \
        --template /path/to/Presentation_Template.pptx \
        --data docs/paper/data \
        --claims docs/paper/supplement/baseline_claims.json \
        --out figure_data/figures_presentation.pptx
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from pptx import Presentation
from pptx.chart.data import CategoryChartData, XyChartData
from pptx.dml.color import RGBColor
from pptx.enum.chart import XL_CHART_TYPE, XL_LABEL_POSITION, XL_LEGEND_POSITION
from pptx.enum.chart import XL_MARKER_STYLE
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt, Emu

INK = RGBColor(0x1F, 0x29, 0x33)
MUTED = RGBColor(0x52, 0x60, 0x6D)
BOX_FILL = RGBColor(0xEE, 0xF2, 0xF8)
BOX_LINE = RGBColor(0x2C, 0x5A, 0x8C)


def pct(metric: dict[str, Any]) -> float:
    return round(float(metric["accuracy"]) * 100.0, 1)


def load_inputs(data_dir: Path, claims_path: Path) -> dict[str, Any]:
    s = json.loads((data_dir / "structure/audit_structure_summary.json").read_text())
    p = json.loads((data_dir / "pathway/audit_pathway_summary.json").read_text())
    loo = json.loads((data_dir / "pathway/audit_pathway_loo.json").read_text())
    pres = json.loads((data_dir / "presentation/audit_presentation_summary.json").read_text())
    claims = {c["claim_id"]: c for c in json.loads(claims_path.read_text())["claims"]}
    return {"s": s, "p": p, "loo": loo, "pres": pres, "claims": claims}


def delete_all_slides(prs: Presentation) -> None:
    sld_id_lst = prs.slides._sldIdLst
    for sld_id in list(sld_id_lst):
        prs.part.drop_rel(sld_id.get(qn("r:id")))
        sld_id_lst.remove(sld_id)


def set_placeholder_text(placeholder, text: str) -> None:
    tf = placeholder.text_frame
    tf.clear()
    tf.paragraphs[0].text = text


def slide_geometry(prs: Presentation) -> dict[str, int]:
    margin = Inches(0.6)
    top = Inches(1.5)
    bottom_margin = Inches(1.15)
    return {
        "margin": margin,
        "top": top,
        "left": margin,
        "width": prs.slide_width - 2 * margin,
        "height": prs.slide_height - top - bottom_margin,
        "cap_top": prs.slide_height - bottom_margin + Inches(0.15),
    }


def add_figure_slide(prs, layout, title: str, caption: str):
    slide = prs.slides.add_slide(layout)
    phs = {ph.placeholder_format.idx: ph for ph in slide.placeholders}
    geo = slide_geometry(prs)
    if 21 in phs:
        set_placeholder_text(phs[21], title)
    else:
        tb = slide.shapes.add_textbox(geo["margin"], Inches(0.45), geo["width"], Inches(0.9))
        para = tb.text_frame.paragraphs[0]
        para.text = title
        para.font.size = Pt(24)
        para.font.bold = True
        para.font.color.rgb = INK
    cap = slide.shapes.add_textbox(geo["margin"], geo["cap_top"], geo["width"], Inches(0.8))
    cap.text_frame.word_wrap = True
    para = cap.text_frame.paragraphs[0]
    para.text = caption
    para.font.size = Pt(12)
    para.font.color.rgb = MUTED
    para.alignment = PP_ALIGN.CENTER
    return slide, geo


def style_value_labels(plot, number_format: str) -> None:
    plot.has_data_labels = True
    labels = plot.data_labels
    labels.number_format = number_format
    labels.number_format_is_linked = False
    labels.font.size = Pt(10)
    labels.font.bold = True
    labels.font.color.rgb = INK
    try:
        labels.position = XL_LABEL_POSITION.OUTSIDE_END
    except (ValueError, NotImplementedError):
        pass


def add_column_chart(
    slide,
    geo,
    categories: list[str],
    series: list[tuple[str, list[float]]],
    *,
    number_format: str = '0.0',
    y_max: float | None = None,
    show_legend: bool | None = None,
) -> None:
    chart_data = CategoryChartData()
    chart_data.categories = categories
    for name, values in series:
        chart_data.add_series(name, values, number_format)
    gf = slide.shapes.add_chart(
        XL_CHART_TYPE.COLUMN_CLUSTERED,
        geo["left"], geo["top"], geo["width"], geo["height"], chart_data,
    )
    chart = gf.chart
    chart.has_title = False
    if show_legend is None:
        show_legend = len(series) > 1
    chart.has_legend = show_legend
    if show_legend:
        chart.legend.position = XL_LEGEND_POSITION.BOTTOM
        chart.legend.include_in_layout = False
        chart.legend.font.size = Pt(11)
    value_axis = chart.value_axis
    if y_max is not None:
        value_axis.maximum_scale = y_max
        value_axis.minimum_scale = 0
    value_axis.has_major_gridlines = True
    value_axis.tick_labels.font.size = Pt(10)
    chart.category_axis.tick_labels.font.size = Pt(10)
    style_value_labels(chart.plots[0], number_format)


def add_scatter_chart(slide, geo, series: list[tuple[str, list[tuple[float, float]]]]) -> None:
    chart_data = XyChartData()
    for name, points in series:
        s = chart_data.add_series(name, number_format='0.00')
        for x, y in points:
            s.add_data_point(x, y)
    gf = slide.shapes.add_chart(
        XL_CHART_TYPE.XY_SCATTER,
        geo["left"], geo["top"], geo["width"], geo["height"], chart_data,
    )
    chart = gf.chart
    chart.has_title = False
    chart.has_legend = True
    chart.legend.position = XL_LEGEND_POSITION.RIGHT
    chart.legend.include_in_layout = False
    for axis in (chart.value_axis, chart.category_axis):
        axis.minimum_scale = 0.0
        axis.maximum_scale = 1.0
        axis.has_major_gridlines = True
        axis.tick_labels.font.size = Pt(10)
        axis.tick_labels.number_format = '0.0'
        axis.tick_labels.number_format_is_linked = False
    # The last series is the 1:1 reference line: draw a line, hide markers.
    series_list = list(chart.series)
    ideal = series_list[-1]
    ideal.marker.style = XL_MARKER_STYLE.NONE
    ideal.format.line.color.rgb = MUTED
    ideal.format.line.width = Pt(1.25)
    # Data series: markers only, no connecting line.
    for data_series in series_list[:-1]:
        data_series.marker.style = XL_MARKER_STYLE.CIRCLE
        data_series.marker.size = 7
        data_series.format.line.fill.background()


def build_fig1(slide, prs) -> None:
    steps = [
        ("Raw LC-MS/MS", "features"),
        ("Discovery", "MS2 feature detection"),
        ("De novo structure prediction", "MS2\u2192SMILES generative"),
        ("Candidate rerank", "library + mass + spectral evidence"),
        ("Pathway consensus", "neighbor vote of Metabolon taxonomy"),
        ("Report", "confidence + MSI level + index"),
    ]
    n = len(steps)
    gap = Inches(0.18)
    left0 = Inches(0.5)
    total_w = prs.slide_width - 2 * left0
    box_w = int((total_w - gap * (n - 1)) / n)
    box_h = Inches(1.6)
    top = Inches(2.9)
    for i, (head, body) in enumerate(steps):
        x = int(left0 + i * (box_w + gap))
        shape = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, x, top, box_w, box_h)
        shape.fill.solid()
        shape.fill.fore_color.rgb = BOX_FILL
        shape.line.color.rgb = BOX_LINE
        shape.line.width = Pt(1.25)
        shape.shadow.inherit = False
        tf = shape.text_frame
        tf.word_wrap = True
        tf.vertical_anchor = MSO_ANCHOR.MIDDLE
        p0 = tf.paragraphs[0]
        p0.text = head
        p0.alignment = PP_ALIGN.CENTER
        p0.font.size = Pt(12)
        p0.font.bold = True
        p0.font.color.rgb = INK
        p1 = tf.add_paragraph()
        p1.text = body
        p1.alignment = PP_ALIGN.CENTER
        p1.font.size = Pt(9)
        p1.font.color.rgb = MUTED
        if i < n - 1:
            ax = x + box_w
            arrow = slide.shapes.add_shape(
                MSO_SHAPE.RIGHT_ARROW, ax, int(top + box_h / 2 - Inches(0.09)),
                gap, Inches(0.18),
            )
            arrow.fill.solid()
            arrow.fill.fore_color.rgb = BOX_LINE
            arrow.line.fill.background()
            arrow.shadow.inherit = False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--data", type=Path, default=Path("docs/paper/data"))
    parser.add_argument("--claims", type=Path, default=Path("docs/paper/supplement/baseline_claims.json"))
    parser.add_argument("--out", type=Path, default=Path("figure_data/figures_presentation.pptx"))
    args = parser.parse_args()

    d = load_inputs(args.data, args.claims)
    s, p, loo, pres, claims = d["s"], d["p"], d["loo"], d["pres"], d["claims"]
    a = p["accuracy"]
    ref = p["reference_db_evaluation"]

    prs = Presentation(str(args.template))
    delete_all_slides(prs)
    layouts = {l.name: l for l in prs.slide_layouts}
    title_layout = layouts.get("Title Slide 1", prs.slide_layouts[0])
    fig_layout = layouts.get("Blank Slide with Title", prs.slide_layouts[5])

    # Title slide
    title_slide = prs.slides.add_slide(title_layout)
    phs = {ph.placeholder_format.idx: ph for ph in title_slide.placeholders}
    if 21 in phs:
        set_placeholder_text(phs[21], "Evidence-Aware Reporting of LC-MS/MS Metabolite Structure Hypotheses")
    if 22 in phs:
        set_placeholder_text(phs[22], "Figure compendium (data-editable charts)\nMetabolon, Inc.")

    # Fig 1 - editable shapes
    slide, _ = add_figure_slide(
        prs, fig_layout, "Figure 1. Structural-elucidation workflow",
        "3,897 X-compound features across 4 chromatographic libraries (2 positive, 2 negative mode). Boxes are editable shapes.",
    )
    build_fig1(slide, prs)

    # Fig 2 - structural recovery (single series)
    r = s["rank_k"]["0.70"]
    exact = s["source_subgroups"]["exact"]["rank1_tanimoto_ge_0p70"]
    predicted = s["source_subgroups"]["predicted"]["rank1_tanimoto_ge_0p70"]
    cats = [
        f"Rank 1 (n={r['top1']['n']})",
        f"Top 3 (n={r['top3']['n']})",
        f"Top 10 (n={r['top10']['n']})",
        f"Exact-source rank 1 (n={exact['n']})",
        f"Predicted-source rank 1 (n={predicted['n']})",
    ]
    vals = [pct(r["top1"]), pct(r["top3"]), pct(r["top10"]), pct(exact), pct(predicted)]
    slide, geo = add_figure_slide(
        prs, fig_layout, "Figure 2. Structural recovery",
        "Rank-1, top-3, and top-10 recovery on truth-verified controls; correct = truth Tanimoto \u2265 0.70.",
    )
    add_column_chart(slide, geo, cats, [("Structural recovery (%)", vals)], y_max=100)

    # Fig 3 - reliability scatter
    full_pts = [(b["mean_predicted"], b["observed"]) for b in s["calibration_full"]["reliability_bins"]]
    abl_pts = [(b["mean_predicted"], b["observed"]) for b in s["calibration_without_pathway"]["reliability_bins"]]
    slide, geo = add_figure_slide(
        prs, fig_layout, "Figure 3. Identification-confidence reliability",
        f"Leave-one-out reliability (n={s['calibration_full']['n']}; full Brier={s['calibration_full']['brier_score']:.3f}). "
        "Diagonal is the 1:1 ideal-calibration line.",
    )
    add_scatter_chart(
        slide, geo,
        [
            ("Full model", full_pts),
            ("Without pathway axis", abl_pts),
            ("Ideal (y=x)", [(0.0, 0.0), (1.0, 1.0)]),
        ],
    )

    # Fig 4 - distributions (two charts on one slide)
    slide, geo = add_figure_slide(
        prs, fig_layout, "Figure 4. Reporting-level distributions",
        "Report counts across MSI-style reporting levels and identification-confidence tiers (3,897 reports).",
    )
    half_w = int((geo["width"] - Inches(0.4)) / 2)
    msi = pres["msi_level_distribution"]
    msi_cd = CategoryChartData()
    msi_cd.categories = [f"L{i}" for i in range(1, 6)]
    msi_cd.add_series("Reports", [msi.get(str(i), 0) for i in range(1, 6)], '#,##0')
    gf1 = slide.shapes.add_chart(
        XL_CHART_TYPE.COLUMN_CLUSTERED, geo["left"], geo["top"], half_w, geo["height"], msi_cd
    )
    c1 = gf1.chart
    c1.has_title = True
    c1.chart_title.text_frame.text = "MSI/Schymanski reporting level"
    c1.chart_title.text_frame.paragraphs[0].font.size = Pt(12)
    c1.has_legend = False
    style_value_labels(c1.plots[0], '#,##0')
    tier = pres["tier_distribution"]
    tier_order = ["Confirmed", "High", "Moderate", "Low", "Tentative"]
    tier_cd = CategoryChartData()
    tier_cd.categories = tier_order
    tier_cd.add_series("Reports", [tier.get(k, 0) for k in tier_order], '#,##0')
    gf2 = slide.shapes.add_chart(
        XL_CHART_TYPE.COLUMN_CLUSTERED,
        int(geo["left"] + half_w + Inches(0.4)), geo["top"], half_w, geo["height"], tier_cd,
    )
    c2 = gf2.chart
    c2.has_title = True
    c2.chart_title.text_frame.text = "Identification-confidence tier"
    c2.chart_title.text_frame.paragraphs[0].font.size = Pt(12)
    c2.has_legend = False
    style_value_labels(c2.plots[0], '#,##0')

    # Fig 5 - pathway evaluations (two series)
    fig5_cats = [
        "Direct lookup consistency (n=1,459)",
        "Neighbor vs known_X (n=644)",
        "All vs known_X (n=909)",
        "LOO mechanism (n=5,095)",
    ]
    super_vals = [
        pct(ref["by_assignment_source"]["exact_reference_lookup"]["super"]),
        pct(a["by_assignment_source"]["neighbor_consensus"]["super"]),
        pct(a["overall"]["super"]),
        pct(loo["super"]),
    ]
    sub_vals = [
        pct(ref["by_assignment_source"]["exact_reference_lookup"]["sub"]),
        pct(a["by_assignment_source"]["neighbor_consensus"]["sub"]),
        pct(a["overall"]["sub"]),
        pct(loo["sub"]),
    ]
    slide, geo = add_figure_slide(
        prs, fig_layout, "Figure 5. Pathway results by evaluation context",
        "Direct-lookup consistency is not prediction accuracy; known_X is independent but taxonomy-discordant.",
    )
    add_column_chart(
        slide, geo, fig5_cats,
        [("Super-pathway (%)", super_vals), ("Subpathway (%)", sub_vals)], y_max=100,
    )

    # Fig S1 - by pathway-confidence tier
    tiers = ["high", "medium", "low"]
    s1_cats = [f"{t} (n={a['by_pathway_tier'][t]['super']['n']})" for t in tiers]
    slide, geo = add_figure_slide(
        prs, fig_layout, "Figure S1. Pathway agreement by confidence tier",
        "Independent-label agreement (known_X cohort) stratified by pathway-confidence tier.",
    )
    add_column_chart(
        slide, geo, s1_cats,
        [
            ("Super-pathway (%)", [pct(a["by_pathway_tier"][t]["super"]) for t in tiers]),
            ("Subpathway (%)", [pct(a["by_pathway_tier"][t]["sub"]) for t in tiers]),
        ],
        y_max=100,
    )

    # Fig S2 - LOO by similarity (single series)
    bins = loo["super_accuracy_by_similarity_bin"]
    bin_order = sorted(bins)
    s2_cats = [f"{b} (n={bins[b]['n']})" for b in bin_order]
    slide, geo = add_figure_slide(
        prs, fig_layout, "Figure S2. LOO accuracy by similarity",
        "Neighbor-vote leave-one-out super-pathway accuracy by nearest-neighbor similarity bin.",
    )
    add_column_chart(
        slide, geo, s2_cats,
        [("Super-pathway (%)", [pct(bins[b]) for b in bin_order])], y_max=100,
    )

    # Fig S3 - by assignment status (two series)
    statuses = [k for k in ("strong_similarity_vote", "weak_similarity_backfill") if k in a["by_assignment_status"]]
    labels = {"strong_similarity_vote": "Strong similarity vote", "weak_similarity_backfill": "Weak-similarity backfill"}
    s3_cats = [f"{labels[k]} (n={a['by_assignment_status'][k]['super']['n']})" for k in statuses]
    slide, geo = add_figure_slide(
        prs, fig_layout, "Figure S3. Agreement by assignment status",
        "Independent-label agreement for strong-similarity votes vs weak-similarity backfill.",
    )
    add_column_chart(
        slide, geo, s3_cats,
        [
            ("Super-pathway (%)", [pct(a["by_assignment_status"][k]["super"]) for k in statuses]),
            ("Subpathway (%)", [pct(a["by_assignment_status"][k]["sub"]) for k in statuses]),
        ],
        y_max=100,
    )

    # Fig S4 - mass-mismatch cohorts (counts, single series)
    all_mm = int(claims["rank1_mass_mismatch_all"]["value"])
    recover = int(claims["rank1_mass_mismatch_recoverable_all"]["value"])
    s4_cats = ["All mismatches", "Recoverable", "No in-list alternative", "MSI-L4 deep audit", "L4 recoverable"]
    s4_vals = [
        all_mm,
        recover,
        all_mm - recover,
        int(claims["rank1_mass_mismatch_msi_l4_deep_audit"]["value"]),
        int(claims["rank1_mass_mismatch_recoverable_l4"]["value"]),
    ]
    slide, geo = add_figure_slide(
        prs, fig_layout, "Figure S4. Rank-1 mass-mismatch cohorts",
        "Report counts; recoverable = a lower-ranked candidate has strong mass fit. Projection, not proven accuracy.",
    )
    add_column_chart(slide, geo, s4_cats, [("Reports", s4_vals)], number_format='#,##0')

    args.out.parent.mkdir(parents=True, exist_ok=True)
    prs.save(str(args.out))
    print(f"[figures-pptx] wrote {args.out} with {len(prs.slides._sldIdLst)} slides")


if __name__ == "__main__":
    main()
