#!/usr/bin/env python3
"""Build the X-compound prediction overview Word document from Markdown.

The output mirrors the uploaded untargeted-pipeline overview: Letter paper,
brand-face body text, 24-point bold title, blue 16-point Heading 1 sections,
compact bullets/numbering, and shaded Consolas command blocks.

Usage:
    python3 tools/build_x_compound_overview_docx.py \
        --source docs/x_compound_prediction_overview.md \
        --output docs/x_compound_prediction_overview.docx

Pass ``--template path/to/untargeted_pipeline_overview.docx`` to reuse the exact
style definitions from the reference document. Without a template, equivalent
styles are created programmatically.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from docx import Document
from docx.document import Document as DocumentType
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tools.brand import MOLECULAR_BLUE  # noqa: E402
from tools.build_manuscript_docx import set_font  # noqa: E402

INLINE = re.compile(r"(\*\*.+?\*\*|`.+?`)")


def clear_body(document: DocumentType) -> None:
    """Remove template content while retaining styles, section properties, and theme."""

    body = document._element.body
    for child in list(body):
        if child.tag != qn("w:sectPr"):
            body.remove(child)


def configure_document(document: DocumentType) -> None:
    """Apply the reference overview's page and typography settings."""

    section = document.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(1)
    section.right_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)

    normal = document.styles["Normal"]
    set_font(normal)
    normal.font.size = Pt(10.5)
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.1

    heading = document.styles["Heading 1"]
    set_font(heading)
    heading.font.size = Pt(16)
    heading.font.bold = True
    heading.font.color.rgb = RGBColor.from_string(MOLECULAR_BLUE.lstrip("#"))
    heading.paragraph_format.space_before = Pt(16)
    heading.paragraph_format.space_after = Pt(8)
    heading.paragraph_format.keep_with_next = True

    for style_name in ("List Bullet", "List Number"):
        style = document.styles[style_name]
        set_font(style)
        style.font.size = Pt(10.5)
        style.paragraph_format.space_after = Pt(5)
        style.paragraph_format.line_spacing = 1.15


def set_cell_shading(cell, fill: str = "F4F6F8") -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shading = tc_pr.find(qn("w:shd"))
    if shading is None:
        shading = OxmlElement("w:shd")
        tc_pr.append(shading)
    shading.set(qn("w:fill"), fill)


def set_cell_margins(cell, *, top: int = 100, start: int = 120, bottom: int = 100, end: int = 120) -> None:
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for margin, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{margin}"))
        if node is None:
            node = OxmlElement(f"w:{margin}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def add_inline_markup(paragraph, value: str) -> None:
    """Render the small Markdown subset used by the source (bold and inline code)."""

    cursor = 0
    for match in INLINE.finditer(value):
        if match.start() > cursor:
            paragraph.add_run(value[cursor : match.start()])
        token = match.group(0)
        if token.startswith("**"):
            run = paragraph.add_run(token[2:-2])
            run.bold = True
        else:
            run = paragraph.add_run(token[1:-1])
            run.font.name = "Consolas"
            run.font.size = Pt(9)
        cursor = match.end()
    if cursor < len(value):
        paragraph.add_run(value[cursor:])


def add_code_block(document: DocumentType, code: list[str]) -> None:
    table = document.add_table(rows=1, cols=1)
    try:
        table.style = "Normal Table"
    except KeyError:
        pass
    table.autofit = False
    cell = table.cell(0, 0)
    cell.width = Inches(6.5)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.TOP
    set_cell_shading(cell)
    set_cell_margins(cell)
    first = cell.paragraphs[0]
    for index, line in enumerate(code):
        paragraph = first if index == 0 else cell.add_paragraph()
        paragraph.paragraph_format.space_after = Pt(0)
        paragraph.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
        run = paragraph.add_run(line)
        run.font.name = "Consolas"
        run.font.size = Pt(8.5)
        run.font.color.rgb = RGBColor(0, 0, 0)
    spacer = document.add_paragraph()
    spacer.paragraph_format.space_after = Pt(0)


def parse_markdown(document: DocumentType, source: str) -> None:
    lines = source.splitlines()
    in_code = False
    code: list[str] = []
    title_seen = False

    for line in lines:
        if line.startswith("```"):
            if in_code:
                add_code_block(document, code)
                code = []
                in_code = False
            else:
                in_code = True
            continue
        if in_code:
            code.append(line)
            continue
        if not line.strip():
            continue
        if line.startswith("# "):
            paragraph = document.add_paragraph()
            paragraph.paragraph_format.space_before = Pt(8)
            paragraph.paragraph_format.space_after = Pt(4)
            paragraph.paragraph_format.line_spacing = 1.0
            run = paragraph.add_run(line[2:].strip())
            run.bold = True
            set_font(run)
            run.font.size = Pt(24)
            run.font.color.rgb = RGBColor(0, 0, 0)
            document.add_paragraph()
            title_seen = True
            continue
        if line.startswith("## "):
            document.add_paragraph(line[3:].strip(), style="Heading 1")
            continue
        if line.startswith("- "):
            paragraph = document.add_paragraph(style="List Bullet")
            add_inline_markup(paragraph, line[2:].strip())
            continue
        if re.match(r"^\d+\.\s", line):
            paragraph = document.add_paragraph(style="List Number")
            add_inline_markup(paragraph, re.sub(r"^\d+\.\s+", "", line))
            continue
        paragraph = document.add_paragraph()
        add_inline_markup(paragraph, line.strip())

    if in_code:
        raise ValueError("Unclosed Markdown code fence")
    if not title_seen:
        raise ValueError("Markdown source must begin with a '# ' title")


def build(source: Path, output: Path, template: Path | None) -> None:
    document = Document(str(template)) if template else Document()
    clear_body(document)
    configure_document(document)
    document.core_properties.title = "X-Compound Structure Prediction Pipeline README"
    document.core_properties.subject = "Overview of Metabolon's X-compound structure prediction workflow"
    document.core_properties.author = "Metabolon"
    document.core_properties.keywords = "X-compounds, structural elucidation, LC-MS/MS, Spec2Mol"
    parse_markdown(document, source.read_text())
    output.parent.mkdir(parents=True, exist_ok=True)
    document.save(str(output))
    print(f"[x-compound-overview] wrote {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, default=Path("docs/x_compound_prediction_overview.md"))
    parser.add_argument("--output", type=Path, default=Path("docs/x_compound_prediction_overview.docx"))
    parser.add_argument("--template", type=Path, default=None)
    args = parser.parse_args()
    if not args.source.exists():
        raise SystemExit(f"source not found: {args.source}")
    if args.template and not args.template.exists():
        raise SystemExit(f"template not found: {args.template}")
    build(args.source, args.output, args.template)


if __name__ == "__main__":
    main()
