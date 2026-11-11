#!/usr/bin/env python3
"""Build one PowerPoint slide per X-compound from a summary sheet (the UPMC slide process).

Vendored from x_compound_info_summary_20261014/build_x_compound_slides.py (2026-10-14),
the deck builder behind UPMC_X-compound_slides_20261014.pptx. Layout, template and styling
are unchanged: the template deck's "2 Columns" layout, one slide per compound, three cards
(rank-1 structure, feature summary, most similar named metabolite) and a takeaway band.
The only changes are that paths are arguments and that blank structure, neighbour, formula
or pathway fields print "Not assigned" rather than an empty drawing.

Called by client_deck.py. Standalone:
  build_slides.py --sheet summary.xlsx --template "20260919 external slide deck.pptx" \\
      --reports-root output/reports_inventory_naming_v2 --out slides.pptx

Needs python-pptx, openpyxl, lxml, Pillow; RDKit is optional (crisp drawings, else dashboard PNGs).
Only fields present in the summary sheet are shown. No internal statistics, paths,
or chemist notes are added.
"""
from __future__ import annotations

import re
import tempfile
from pathlib import Path

from lxml import etree
from openpyxl import load_workbook
from PIL import Image, ImageOps
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Emu, Inches, Pt

NOT_ASSIGNED = "Not assigned"
HTML: Path | None = None  # report root holding X-#####/assets and data; set by build()

# Metabolon theme (ppt/theme/theme1.xml of the reference deck)
FONT = "Montserrat"
NAVY = RGBColor(0x07, 0x1C, 0x35)      # title / takeaway text used on the deck's slides
INK = RGBColor(0x00, 0x1D, 0x37)       # accent1, body text
GREY = RGBColor(0x53, 0x53, 0x53)      # lt2, labels
MID = RGBColor(0x7F, 0x7F, 0x7F)       # footer (bg1 lumMod 50%)
CORAL = RGBColor(0xFF, 0x4F, 0x57)     # accent3
ORANGE = RGBColor(0xFC, 0xA7, 0x50)    # accent5
CARD_LINE = RGBColor(0xD3, 0xDE, 0xE6)  # accent6
BAND = RGBColor(0xF2, 0xF2, 0xF2)      # bg1 lumMod 95%
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

FOOTER_TEXT = "2026 Metabolon, Inc. All rights reserved. "
FOOTER_CONF = "Confidential."

SHADOW_XML = (
    '<a:effectLst xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
    '<a:outerShdw blurRad="241300" sx="101000" sy="101000" algn="ctr" rotWithShape="0">'
    '<a:srgbClr val="000000"><a:alpha val="10000"/></a:srgbClr></a:outerShdw></a:effectLst>'
)


# --------------------------------------------------------------------------- data
def read_sheet(path: Path) -> list[dict]:
    ws = load_workbook(path, read_only=True, data_only=True).active
    rows = [r for r in ws.iter_rows(values_only=True) if any(v is not None for v in r)]
    header = [str(h).strip() for h in rows[0]]
    out = []
    for r in rows[1:]:
        rec = dict(zip(header, r))
        if rec.get("Compound ID"):
            out.append(rec)
    return out


def slug(s: str) -> str:
    s = re.sub(r"\s*\(tied .*$", "", s, flags=re.I)  # drop "(tied with ...)" qualifier
    s = re.sub(r"[^a-z0-9]+", "-", s.lower())
    return s.strip("-")


def neighbour_png(xid: str, name: str) -> Path | None:
    want = slug(name)
    for p in sorted((HTML / xid / "assets").glob("neighbor-*.png")):
        file_slug = re.sub(r"-+", "-", p.stem.split("-", 2)[2]).strip("-")
        if file_slug == want:
            return p
    return None


def structure_png(xid: str) -> Path | None:
    p = HTML / xid / "assets" / f"{xid.lower()}_structure.png"
    return p if p.exists() else None


# --------------------------------------------------------------------------- pptx helpers
def clear_slides(prs: Presentation) -> None:
    sld_ids = prs.slides._sldIdLst
    for sld_id in list(sld_ids):
        prs.part.drop_rel(sld_id.rId)
        sld_ids.remove(sld_id)


def set_run(run, text, size, bold=False, color=INK, italic=False):
    run.text = text
    f = run.font
    f.name = FONT
    f.size = Pt(size)
    f.bold = bold
    f.italic = italic
    f.color.rgb = color
    return run


def textbox(slide, left, top, width, height, *, anchor=MSO_ANCHOR.TOP, wrap=True, margins=0):
    tb = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(height))
    tf = tb.text_frame
    tf.word_wrap = wrap
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = Inches(margins)
    return tb, tf


def card(slide, left, top, width, height):
    shp = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(left), Inches(top), Inches(width), Inches(height))
    shp.adjustments[0] = 0.07178
    shp.fill.solid()
    shp.fill.fore_color.rgb = WHITE
    shp.line.color.rgb = CARD_LINE
    shp.line.width = Pt(0.75)
    shp.shadow.inherit = False  # python-pptx writes an empty effectLst; replace with deck's soft shadow
    spPr = shp._element.spPr
    for old in spPr.findall("{http://schemas.openxmlformats.org/drawingml/2006/main}effectLst"):
        spPr.remove(old)
    spPr.append(etree.fromstring(SHADOW_XML))
    shp.text_frame.text = ""
    return shp


def heading(slide, left, top, width, text):
    tb, tf = textbox(slide, left, top, width, 0.4, anchor=MSO_ANCHOR.BOTTOM)
    p = tf.paragraphs[0]
    set_run(p.add_run(), text, 14, bold=True, color=CORAL)
    bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(left), Inches(top + 0.46), Inches(width), Inches(0.045))
    bar.line.fill.background()
    bar.fill.gradient()
    bar.fill.gradient_angle = 0
    stops = bar.fill.gradient_stops
    stops[0].color.rgb = CORAL
    stops[0].position = 0.0
    stops[1].color.rgb = ORANGE
    stops[1].position = 1.0
    bar.shadow.inherit = False
    return top + 0.6


_TMP = Path(tempfile.mkdtemp(prefix="xslides_"))

try:  # optional: crisp high-resolution drawings (spec2mol/.venv has RDKit); else dashboard PNGs
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import rdMolDraw2D
except ImportError:  # pragma: no cover
    Chem = None


def neighbour_smiles(xid: str, name: str) -> str | None:
    if not name:
        return None
    want = slug(name)
    f = HTML / xid / "data" / "similarity_neighbors.csv"
    if not f.exists():
        return None
    import csv

    for r in csv.DictReader(open(f, encoding="utf-8")):
        if slug(r.get("canonical_name", "")) == want:
            return r.get("smiles") or None
    return None


def rdkit_png(smiles: str, stem: str, w_px=1600, h_px=900) -> Path | None:
    if Chem is None or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    d = rdMolDraw2D.MolDraw2DCairo(w_px, h_px)
    opts = d.drawOptions()
    opts.setBackgroundColour((1, 1, 1, 1))
    opts.bondLineWidth = 3.0
    opts.padding = 0.05
    opts.fixedBondLength = 70
    opts.scaleBondWidth = True
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    out = _TMP / f"{stem}.png"
    out.write_bytes(d.GetDrawingText())
    return out


def prepared_png(path: Path) -> Path:
    """Trim surrounding whitespace so the molecule fills the card. For dashboard PNGs
    (520x360) also drop the baked-in caption strip at the bottom. Written to a temp dir."""
    im = Image.open(path).convert("RGB")
    w, h = im.size
    if (w, h) == (520, 360):
        im = im.crop((0, 0, w, int(h * 0.89)))
    bbox = ImageOps.invert(im).getbbox()
    if bbox:
        pad = max(14, w // 60)
        im = im.crop((max(bbox[0] - pad, 0), max(bbox[1] - pad, 0), min(bbox[2] + pad, w), min(bbox[3] + pad, im.size[1])))
    out = _TMP / ("prep_" + path.name)
    im.save(out)
    return out


def picture_fit(slide, path: Path, left, top, box_w, box_h):
    """Place a drawing centred inside the box, preserving aspect ratio. Low-res dashboard
    PNGs are not upscaled past ~1.6x; RDKit renders are high-res and fill the box."""
    path = prepared_png(path)
    pw, ph = Image.open(path).size
    ar = ph / pw
    w = min(box_w, pw / 96 * 1.6)
    if ar < 0.5:  # very elongated chains: let them use the full width but keep legible height
        w = box_w
    h = w * ar
    if h > box_h:
        h = box_h
        w = h / ar
    x = left + (box_w - w) / 2
    y = top + (box_h - h) / 2
    return slide.shapes.add_picture(str(path), Inches(x), Inches(y), Inches(w), Inches(h))


# A molecular formula, or a carbon-count token such as C9 or C20. Digits that are
# not attached to an element symbol (m/z, Tanimoto, SMILES ring numbers) are left alone.
_FORMULA_TOKEN = re.compile(r"(?<![A-Za-z0-9])((?:[A-Z][a-z]?\d*)+)(?![A-Za-z])")


def _formula_parts(token: str) -> list[str] | None:
    """Atom counts for a molecular formula or a carbon-count token such as C9.

    A single letter plus digits is not enough: gibberellin A12 is a name, not a formula.
    """
    if not re.search(r"\d", token):
        return None
    parts = re.findall(r"[A-Z][a-z]?|\d+", token)
    if not parts or "".join(parts) != token:
        return None
    elements = [part for part in parts if not part.isdigit()]
    if len(elements) >= 2:
        return parts
    if elements == ["C"]:
        return parts
    return None


def add_formula_runs(paragraph, text, size, bold=False, color=INK, italic=False):
    """Write text, with atom counts inside chemical formulas set as subscripts."""
    pos = 0
    for match in _FORMULA_TOKEN.finditer(text):
        parts = _formula_parts(match.group(1))
        if not parts:
            continue
        if match.start() > pos:
            set_run(paragraph.add_run(), text[pos:match.start()], size, bold=bold, color=color, italic=italic)
        for part in parts:
            run = set_run(paragraph.add_run(), part, size, bold=bold, color=color, italic=italic)
            if part.isdigit():
                run.font._element.set("baseline", "-25000")
        pos = match.end()
    if pos < len(text) or pos == 0:
        set_run(paragraph.add_run(), text[pos:], size, bold=bold, color=color, italic=italic)


def label_value(tf, label, value, first=False, value_size=12, formulas=False):
    p = tf.paragraphs[0] if first else tf.add_paragraph()
    p.space_before = Pt(0 if first else 10)
    set_run(p.add_run(), label, 9, bold=True, color=GREY)
    p2 = tf.add_paragraph()
    p2.space_before = Pt(1)
    if formulas:
        add_formula_runs(p2, value, value_size, color=INK)
    else:
        set_run(p2.add_run(), value, value_size, color=INK)


def footer(slide):
    tb, tf = textbox(slide, 9.49, 7.10, 3.21, 0.39, anchor=MSO_ANCHOR.MIDDLE)
    tf.margin_left = Inches(0.1)
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.RIGHT
    set_run(p.add_run(), FOOTER_TEXT, 8, bold=True, color=MID)
    set_run(p.add_run(), FOOTER_CONF, 8, bold=True, color=CORAL)


def fmt_mz(v) -> str:
    try:
        return f"{float(v):.4f}"
    except (TypeError, ValueError):
        return str(v)


def fmt_sim(v) -> str:
    try:
        return f"{float(v):.2f}"
    except (TypeError, ValueError):
        return str(v)


# --------------------------------------------------------------------------- slide builder
def build_slide(prs, layout, rec: dict):
    xid = str(rec["Compound ID"]).strip()
    mz = fmt_mz(rec.get("Observed m/z"))
    formula = str(rec.get("Putative formula") or "")
    superp = str(rec.get("Candidate superpathway") or "")
    subp = str(rec.get("Candidate subpathway") or "")
    smiles = str(rec.get("Rank-1 predicted SMILES") or "")
    info = str(rec.get("Rank-1 additional Information") or "")
    nearest = str(rec.get("Most similar named metabolite") or "")
    sim = fmt_sim(rec.get("Tanimoto similarity")) if rec.get("Tanimoto similarity") is not None else NOT_ASSIGNED
    group = str(rec.get("UPF group") or "").strip()
    pathway_line = f"{superp}  ›  {subp}" if superp else NOT_ASSIGNED

    s = prs.slides.add_slide(layout)

    # placeholders: 23 = title band, 21 = subtitle line; 24/4 = unused two-column bodies
    for ph in list(s.placeholders):
        idx = ph.placeholder_format.idx
        if idx == 23:
            # the deck's content slides sit the title in a light-grey band
            ph.fill.solid()
            ph.fill.fore_color.rgb = BAND
            ph.text_frame.text = ""
            ph.text_frame.margin_left = Inches(0.15)
            set_run(ph.text_frame.paragraphs[0].add_run(), f"{xid}: Preliminary Structural Information", 30, bold=True, color=NAVY)
        elif idx == 21:
            ph.text_frame.text = ""
            p = ph.text_frame.paragraphs[0]
            if group:
                set_run(p.add_run(), f"{group}    ", 13, bold=True, color=NAVY)
            set_run(p.add_run(), "Candidate pathway  ", 13, bold=True, color=GREY)
            set_run(p.add_run(), pathway_line, 13, color=INK)
        else:
            ph._element.getparent().remove(ph._element)

    top, h = 1.45, 4.6
    gap = 0.22
    w_a, w_b = 4.15, 3.85
    x_a = 0.45
    x_b = x_a + w_a + gap
    x_c = x_b + w_b + gap
    w_c = 12.70 - x_c
    pad = 0.22

    # --- card A: rank-1 predicted structure
    card(s, x_a, top, w_a, h)
    y = heading(s, x_a + pad, top + 0.12, w_a - 2 * pad, "Rank-1 Predicted Structure")
    img = (rdkit_png(smiles, f"{xid}_rank1") or structure_png(xid)) if smiles else None
    if img:
        picture_fit(s, img, x_a + pad, y + 0.05, w_a - 2 * pad, 2.35)
    else:
        _, ntf = textbox(s, x_a + pad, y + 0.05, w_a - 2 * pad, 2.35, anchor=MSO_ANCHOR.MIDDLE)
        ntf.paragraphs[0].alignment = PP_ALIGN.CENTER
        set_run(ntf.paragraphs[0].add_run(), "No structure assigned", 14, bold=True, color=GREY)
    tb, tf = textbox(s, x_a + pad, y + 2.5, w_a - 2 * pad, h - (y - top) - 2.6)
    label_value(tf, "Rank-1 predicted SMILES", smiles or NOT_ASSIGNED, first=True, value_size=9)
    p = tf.add_paragraph()
    p.space_before = Pt(6)
    note = ("Computational candidate structure; not a confirmed identification." if smiles
            else "Tentative class annotation only; the molecular structure remains unassigned.")
    set_run(p.add_run(), note, 8, italic=True, color=GREY)

    # --- card B: feature summary from the sheet
    card(s, x_b, top, w_b, h)
    y = heading(s, x_b + pad, top + 0.12, w_b - 2 * pad, "Feature Summary")
    tb, tf = textbox(s, x_b + pad, y + 0.05, w_b - 2 * pad, h - (y - top) - 0.2)
    label_value(tf, "Observed m/z", mz, first=True, value_size=14)
    label_value(tf, "Putative formula", formula or NOT_ASSIGNED, value_size=14, formulas=True)
    label_value(tf, "Candidate superpathway", superp or NOT_ASSIGNED, formulas=True)
    label_value(tf, "Candidate subpathway", subp or NOT_ASSIGNED, formulas=True)

    # --- card C: most similar named metabolite
    card(s, x_c, top, w_c, h)
    y = heading(s, x_c + pad, top + 0.12, w_c - 2 * pad, "Most Similar Named Metabolite")
    nsmi = neighbour_smiles(xid, nearest)
    img = (rdkit_png(nsmi, f"{xid}_neighbour") if nsmi else None) or neighbour_png(xid, nearest)
    if img:
        picture_fit(s, img, x_c + pad, y + 0.05, w_c - 2 * pad, 2.35)
    else:
        _, ntf = textbox(s, x_c + pad, y + 0.05, w_c - 2 * pad, 2.35, anchor=MSO_ANCHOR.MIDDLE)
        ntf.paragraphs[0].alignment = PP_ALIGN.CENTER
        set_run(ntf.paragraphs[0].add_run(), "No library neighbour available", 14, bold=True, color=GREY)
    tb, tf = textbox(s, x_c + pad, y + 2.5, w_c - 2 * pad, h - (y - top) - 2.6)
    label_value(tf, "Named metabolite", nearest or NOT_ASSIGNED, first=True, value_size=12, formulas=True)
    label_value(tf, "Tanimoto similarity", sim, value_size=14)

    # --- takeaway band (same style as the deck's bottom banner: grey band, navy bold)
    band = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0.45), Inches(6.28), Inches(12.25), Inches(0.52))
    band.fill.solid()
    band.fill.fore_color.rgb = BAND
    band.line.fill.background()
    band.shadow.inherit = False
    tf = band.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tf.margin_left = tf.margin_right = Inches(0.15)
    tf.margin_top = tf.margin_bottom = 0
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    set_run(p.add_run(), "Rank-1 candidate: ", 18, bold=True, color=CORAL)
    shown = info[:1].upper() + info[1:] if info else ""
    add_formula_runs(p, shown, 18, bold=True, color=NAVY)

    footer(s)
    return s


def build(sheet: Path, template: Path, out: Path, reports_root: Path) -> int:
    global HTML
    HTML = Path(reports_root)
    recs = read_sheet(Path(sheet))
    prs = Presentation(str(template))
    layout = next(l for s in prs.slides for l in [s.slide_layout] if l.name == "2 Columns")
    clear_slides(prs)
    for rec in recs:
        build_slide(prs, layout, rec)
    prs.core_properties.title = "X-compound preliminary structural information"
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    prs.save(str(out))
    return len(prs.slides)


def main():
    import argparse

    ap = argparse.ArgumentParser(description="Build the X-compound slide deck from a summary sheet.")
    ap.add_argument("--sheet", required=True, type=Path)
    ap.add_argument("--template", required=True, type=Path)
    ap.add_argument("--reports-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    n = build(a.sheet, a.template, a.out, a.reports_root)
    print(f"wrote {a.out.name}: {n} slides from {a.sheet.name}")


if __name__ == "__main__":
    main()
