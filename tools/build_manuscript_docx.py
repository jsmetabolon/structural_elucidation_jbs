#!/usr/bin/env python3
"""Build a Nature initial-submission Word manuscript from ``docs/paper/manuscript.md``.

Main SVG figures are rasterized at 2x resolution and embedded; Markdown tables,
headings, lists, captions, emphasis, inline code, superscripts, page numbers, and
line numbering are preserved in Word-native form.

Usage:
    python3 tools/build_manuscript_docx.py \
        --source docs/paper/manuscript.md \
        --output docs/paper/manuscript.docx
"""
from __future__ import annotations

import argparse
import re
import sys
import tempfile
from pathlib import Path

import cairosvg
from docx import Document
from docx.document import Document as DocumentType
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tools.figure_formats import png_dimensions  # noqa: E402

# Column width, and the tallest a figure may be before it is scaled down to fit a
# page. Screenshots are much taller than charts, so without the height clamp a
# full-width capture runs off the bottom of the page.
FIGURE_WIDTH_INCHES = 6.35
FIGURE_MAX_HEIGHT_INCHES = 8.0

INLINE = re.compile(
    r"(\*\*.+?\*\*|(?<!\*)\*[^*]+?\*(?!\*)|`.+?`|<sup>.+?</sup>|\[[^\]]+\]\([^)]+\))"
)
FIGURE = re.compile(r"^!\[(?P<alt>[^\]]*)\]\((?P<path>[^)]+)\)$")
TABLE_SEPARATOR = re.compile(r"^\|(?:\s*:?-+:?\s*\|)+$")
SECTION_ORDER = (
    "Abstract",
    "Introduction",
    "Results",
    "Discussion",
    "Limitations",
    "Conclusion",
    "Future work",
    "References",
    "Methods",
    "Data and code availability",
    "Acknowledgments",
    "Funding",
    "Author contributions",
    "Competing interests",
    "Ethics and human subjects",
    # Last on purpose. Unrecognised headings are appended after this tuple, so the
    # glossary would land here anyway; naming it makes the placement a decision
    # rather than a side effect of the fallback.
    "Glossary",
)


def configure(document: DocumentType) -> None:
    section = document.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(1)
    section.right_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)

    normal = document.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(12)
    normal.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    normal.paragraph_format.space_after = Pt(0)

    title = document.styles["Title"]
    title.font.name = "Times New Roman"
    title.font.size = Pt(16)
    title.font.bold = True
    title.font.color.rgb = RGBColor(0, 0, 0)
    title.paragraph_format.space_after = Pt(12)

    heading1 = document.styles["Heading 1"]
    heading1.font.name = "Times New Roman"
    heading1.font.size = Pt(12)
    heading1.font.bold = True
    heading1.font.color.rgb = RGBColor(0, 0, 0)
    heading1.paragraph_format.space_before = Pt(12)
    heading1.paragraph_format.space_after = Pt(6)
    heading1.paragraph_format.keep_with_next = True

    heading2 = document.styles["Heading 2"]
    heading2.font.name = "Times New Roman"
    heading2.font.size = Pt(12)
    heading2.font.bold = True
    heading2.font.italic = True
    heading2.font.color.rgb = RGBColor(0, 0, 0)
    heading2.paragraph_format.space_before = Pt(10)
    heading2.paragraph_format.space_after = Pt(4)
    heading2.paragraph_format.keep_with_next = True

    caption = document.styles["Caption"]
    caption.font.name = "Times New Roman"
    caption.font.size = Pt(9)
    caption.font.italic = True
    caption.font.color.rgb = RGBColor(0x33, 0x33, 0x33)
    caption.paragraph_format.space_after = Pt(8)
    caption.paragraph_format.keep_with_next = False

    for style_name in ("List Bullet", "List Number"):
        style = document.styles[style_name]
        style.font.name = "Times New Roman"
        style.font.size = Pt(12)
        style.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
        style.paragraph_format.space_after = Pt(3)

    add_line_numbers(section)
    add_page_number(section)


def add_line_numbers(section) -> None:
    sect_pr = section._sectPr
    existing = sect_pr.find(qn("w:lnNumType"))
    if existing is not None:
        sect_pr.remove(existing)
    line_numbers = OxmlElement("w:lnNumType")
    line_numbers.set(qn("w:countBy"), "1")
    line_numbers.set(qn("w:restart"), "newPage")
    line_numbers.set(qn("w:distance"), "360")
    sect_pr.append(line_numbers)


def add_page_number(section) -> None:
    paragraph = section.footer.paragraphs[0]
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run()
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instruction = OxmlElement("w:instrText")
    instruction.set(qn("xml:space"), "preserve")
    instruction.text = " PAGE "
    separate = OxmlElement("w:fldChar")
    separate.set(qn("w:fldCharType"), "separate")
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    run._r.extend((begin, instruction, separate, end))


def add_markup(paragraph, value: str) -> None:
    cursor = 0
    for match in INLINE.finditer(value):
        if match.start() > cursor:
            paragraph.add_run(value[cursor : match.start()])
        token = match.group(0)
        if token.startswith("**"):
            run = paragraph.add_run(token[2:-2])
            run.bold = True
        elif token.startswith("*"):
            run = paragraph.add_run(token[1:-1])
            run.italic = True
        elif token.startswith("`"):
            run = paragraph.add_run(token[1:-1])
            run.font.name = "Consolas"
            run.font.size = Pt(9)
        elif token.startswith("<sup>"):
            run = paragraph.add_run(token[5:-6])
            run.font.superscript = True
        else:
            link = re.match(r"\[([^\]]+)\]\(([^)]+)\)", token)
            label, target = link.groups()
            paragraph.add_run(f"{label} ({target})")
        cursor = match.end()
    if cursor < len(value):
        paragraph.add_run(value[cursor:])


def split_table_row(line: str) -> list[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def shade_cell(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shading = tc_pr.find(qn("w:shd"))
    if shading is None:
        shading = OxmlElement("w:shd")
        tc_pr.append(shading)
    shading.set(qn("w:fill"), fill)


def repeat_table_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    header = OxmlElement("w:tblHeader")
    header.set(qn("w:val"), "true")
    tr_pr.append(header)


def column_widths(rows: list[list[str]], columns: int, total_inches: float) -> list[float]:
    """Split the table width in proportion to how much text each column holds.

    Word's autofit gives a glossary's term column as much room as its definition
    column, which wastes half the width and pushes every definition onto an extra
    line. Sizing by content keeps a short first column short. The floor and
    ceiling stop one long cell from collapsing its neighbours.
    """

    body = rows[1:] or rows
    weights = []
    for index in range(columns):
        lengths = [len(row[index]) for row in body if index < len(row)]
        weights.append(max(sum(lengths) / len(lengths), 1.0) if lengths else 1.0)
    total = sum(weights)
    # The floor is set by the longest label a narrow column has to hold on one
    # line, not by aesthetics: below roughly a fifth of the width, glossary terms
    # start wrapping and cost back the lines the resizing just saved.
    floor, ceiling = 0.20, 0.72
    shares = [min(max(weight / total, floor), ceiling) for weight in weights]
    scale = sum(shares)
    return [total_inches * share / scale for share in shares]


def add_table(document: DocumentType, rows: list[list[str]], total_inches: float = 6.35) -> None:
    if not rows:
        return
    width = max(len(row) for row in rows)
    table = document.add_table(rows=len(rows), cols=width)
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    widths = column_widths(rows, width, total_inches)
    for column_index, inches in enumerate(widths):
        table.columns[column_index].width = Inches(inches)
    for row_index, source_row in enumerate(rows):
        row = table.rows[row_index]
        row._tr.get_or_add_trPr()
        if row_index == 0:
            repeat_table_header(row)
        for column_index in range(width):
            cell = row.cells[column_index]
            # Word honours cell width over column width, so both are set.
            cell.width = Inches(widths[column_index])
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            value = source_row[column_index] if column_index < len(source_row) else ""
            paragraph = cell.paragraphs[0]
            paragraph.paragraph_format.space_after = Pt(0)
            paragraph.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            add_markup(paragraph, value)
            for run in paragraph.runs:
                run.font.name = "Times New Roman"
                run.font.size = Pt(8.5)
                if row_index == 0:
                    run.bold = True
                    run.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
            if row_index == 0:
                shade_cell(cell, "1F4E79")
            elif row_index % 2 == 0:
                shade_cell(cell, "EAF2F8")
    paragraph = document.add_paragraph()
    paragraph.paragraph_format.space_after = Pt(0)


def add_figure(
    document: DocumentType,
    source_root: Path,
    figure_path: str,
    alt: str,
    temp_dir: Path,
    *,
    max_height_inches: float = FIGURE_MAX_HEIGHT_INCHES,
) -> None:
    source = source_root / figure_path
    if not source.exists():
        raise FileNotFoundError(f"figure not found: {source}")
    png = temp_dir / f"{source.stem}.png"
    if source.suffix.lower() == ".svg":
        cairosvg.svg2png(url=str(source), write_to=str(png), scale=2.5)
    else:
        png = source
    width_inches = FIGURE_WIDTH_INCHES
    pixel_width, pixel_height = png_dimensions(png)
    printed_height = FIGURE_WIDTH_INCHES * pixel_height / pixel_width
    if printed_height > max_height_inches:
        width_inches = FIGURE_WIDTH_INCHES * max_height_inches / printed_height
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.paragraph_format.space_after = Pt(3)
    paragraph.paragraph_format.keep_with_next = True
    run = paragraph.add_run()
    run.add_picture(str(png), width=Inches(width_inches))
    shape = document.inline_shapes[-1]
    shape._inline.docPr.set("descr", alt)


def add_front_matter(document: DocumentType, line: str) -> None:
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.paragraph_format.space_after = Pt(3)
    add_markup(paragraph, line)


def clean_heading(value: str) -> str:
    heading = re.sub(r"^\d+(?:\.\d+)*\.?\s+", "", value).strip()
    return "Acknowledgements" if heading == "Acknowledgments" else heading


def nature_order(source_text: str) -> str:
    """Move Methods/declarations to Nature order while retaining inline figures."""

    matches = list(re.finditer(r"(?m)^## (.+)$", source_text))
    if not matches:
        return source_text
    preamble = source_text[: matches[0].start()].rstrip()
    sections: list[tuple[str, str]] = []
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(source_text)
        heading = match.group(1).strip()
        body = source_text[match.end() : end].strip()
        sections.append((heading, body))
    by_clean = {clean_heading(heading): (heading, body) for heading, body in sections}
    ordered = []
    used = set()
    for desired in SECTION_ORDER:
        lookup = "Acknowledgements" if desired == "Acknowledgments" else desired
        if lookup in by_clean:
            original, body = by_clean[lookup]
            ordered.append(f"## {original}\n\n{body}")
            used.add(lookup)
    for heading, body in sections:
        cleaned = clean_heading(heading)
        if cleaned not in used:
            ordered.append(f"## {heading}\n\n{body}")
    return preamble + "\n\n" + "\n\n".join(ordered) + "\n"


def render(source: Path, output: Path) -> None:
    source_text = nature_order(source.read_text())
    title_line = next((line[2:].strip() for line in source_text.splitlines() if line.startswith("# ")), "")
    if len(title_line) > 75:
        raise ValueError(f"Nature title exceeds 75 characters ({len(title_line)}): {title_line}")
    abstract = source_text.split("## Abstract", 1)[1].split("---", 1)[0]
    abstract_words = len(abstract.split())
    if abstract_words > 200:
        raise ValueError(f"Nature summary exceeds ~200 words: {abstract_words}")
    main_figures = len(re.findall(r"!\[[^\]]+\]\(figures/fig\d+[^)]*\)", source_text))
    main_tables = len(re.findall(r"^\*\*Table \d+\.", source_text, flags=re.MULTILINE))
    if main_figures + main_tables > 7:
        raise ValueError(
            f"Nature main display-item limit exceeded: {main_figures} figures + "
            f"{main_tables} tables"
        )

    document = Document()
    configure(document)
    document.core_properties.title = title_line
    document.core_properties.subject = "Structural elucidation methods manuscript"
    document.core_properties.author = "Joshua Simpson et al."
    document.core_properties.keywords = "metabolomics, LC-MS/MS, structural elucidation, MSI"

    lines = source_text.splitlines()
    source_root = source.parent
    section = ""
    front_matter = True
    index = 0
    with tempfile.TemporaryDirectory(prefix="paper_docx_") as temporary:
        temp_dir = Path(temporary)
        while index < len(lines):
            line = lines[index].strip()
            if not line or line == "---":
                index += 1
                continue
            if line.startswith("# "):
                paragraph = document.add_paragraph(style="Title")
                paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
                add_markup(paragraph, line[2:])
                index += 1
                continue
            if line.startswith("## "):
                section = clean_heading(line[3:].strip())
                front_matter = False
                document.add_paragraph(section, style="Heading 1")
                index += 1
                continue
            if line.startswith("### "):
                document.add_paragraph(clean_heading(line[4:].strip()), style="Heading 2")
                index += 1
                continue
            figure = FIGURE.match(line)
            if figure:
                add_figure(document, source_root, figure.group("path"), figure.group("alt"), temp_dir)
                index += 1
                continue
            if line.startswith("|") and index + 1 < len(lines) and TABLE_SEPARATOR.match(lines[index + 1].strip()):
                rows = [split_table_row(line)]
                index += 2
                while index < len(lines) and lines[index].strip().startswith("|"):
                    rows.append(split_table_row(lines[index]))
                    index += 1
                add_table(document, rows)
                continue
            if line.startswith("**Table "):
                paragraph = document.add_paragraph(style="Caption")
                paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
                paragraph.paragraph_format.keep_with_next = True
                add_markup(paragraph, line)
                index += 1
                continue
            if line.startswith("*Figure ") and line.endswith("*"):
                paragraph = document.add_paragraph(style="Caption")
                paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
                add_markup(paragraph, line)
                index += 1
                continue
            if line.startswith("- "):
                paragraph = document.add_paragraph(style="List Bullet")
                add_markup(paragraph, line[2:])
                index += 1
                continue
            numbered = re.match(r"^(\d+)\.\s+(.*)", line)
            if numbered and section != "References":
                paragraph = document.add_paragraph(style="List Number")
                add_markup(paragraph, numbered.group(2))
                index += 1
                continue
            if numbered and section == "References":
                paragraph = document.add_paragraph()
                paragraph.paragraph_format.left_indent = Inches(0.25)
                paragraph.paragraph_format.first_line_indent = Inches(-0.25)
                paragraph.paragraph_format.line_spacing = 1.0
                paragraph.paragraph_format.space_after = Pt(4)
                add_markup(paragraph, line)
                index += 1
                continue
            if front_matter:
                add_front_matter(document, line)
                index += 1
                continue
            paragraph = document.add_paragraph()
            add_markup(paragraph, line)
            index += 1

    output.parent.mkdir(parents=True, exist_ok=True)
    document.save(output)
    print(
        f"[manuscript-docx] wrote {output} "
        f"(title={len(title_line)} chars; summary={abstract_words} words; "
        f"display items={main_figures + main_tables})"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("docs/paper/manuscript.md"))
    parser.add_argument("--output", type=Path, default=Path("docs/paper/manuscript.docx"))
    args = parser.parse_args()
    if not args.source.exists():
        raise SystemExit(f"source not found: {args.source}")
    render(args.source, args.output)


if __name__ == "__main__":
    main()
