#!/usr/bin/env python3
"""Build the public capability brief as a marketing-styled Word document.

Markdown parsing is shared with the Nature manuscript builder; only the page
setup, typography, and budget rules differ. The brief uses the brand face at
reading size with single spacing, coloured headings, and no line numbering, and
it caps its own length so the deliverable stays inside the ten-page limit.

Usage:
    python3 tools/build_public_brief_docx.py \
        --source docs/public_brief/public_brief.md \
        --output docs/public_brief/public_brief.docx
"""

from __future__ import annotations

import argparse
import re
import sys
import tempfile
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docx import Document  # noqa: E402
from docx.document import Document as DocumentType  # noqa: E402
from docx.enum.text import WD_ALIGN_PARAGRAPH  # noqa: E402
from docx.shared import Inches, Pt, RGBColor  # noqa: E402

from tools.brand import DARK_MODE_BLUE, MUTED  # noqa: E402
from tools.build_manuscript_docx import (  # noqa: E402
    FIGURE,
    TABLE_SEPARATOR,
    add_figure,
    add_markup,
    add_table,
    set_font,
    split_table_row,
)

# Budgets sized so the rendered brief stays well under the ten-page requirement.
MAX_WORDS = 2600
MAX_FIGURES = 5
MAX_TABLES = 3

# Page width less the narrower margins this document uses.
TEXT_WIDTH_INCHES = 6.7

# Montserrat sets wider than the Calibri this document was budgeted against, and
# the budget is pages rather than words, so the typeface change had to be paid for
# somewhere. Half a point off the body covers it, and Montserrat's larger x-height
# means the result does not read smaller than the Calibri it replaces. Nothing
# else shrank: the headings and captions were measured at their original sizes and
# did not cost a page.
#
# The value is measured rather than chosen. Rendering through LibreOffice, this
# document turns 10 pages somewhere between 10 and 11 points and 11 pages at 12,
# so 10 leaves a full point of headroom before the count moves at all. Absolute
# counts drift between renderers — the environment behind the README's ten-page
# note measures nine here — so keep the headroom rather than trusting any single
# number, and re-run the gate with --check-pages after changing this.
BODY_SIZE_PT = 10.0

# Much tighter than the manuscript's cap, and the value matters: a report
# screenshot printed any taller than this cannot share a page, so Word gives it
# one of its own and the brief spends a tenth of its budget on white space. The
# figure is an illustration of the deliverable rather than something to read
# word by word, and the section around it states in prose everything it shows.
MAX_FIGURE_HEIGHT_INCHES = 4.0


def _rgb(value: str) -> RGBColor:
    return RGBColor.from_string(value.lstrip("#").upper())


# Dark Mode Blue carries both the brand accent and the body ink here. The brand
# has no second dark tone, and setting headings and body in the same colour at
# different weights is how the guide's own pages are laid out.
BRAND = _rgb(DARK_MODE_BLUE)
INK = _rgb(DARK_MODE_BLUE)
MUTED_INK = _rgb(MUTED)


def configure(document: DocumentType) -> None:
    section = document.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(0.9)
    section.right_margin = Inches(0.9)
    section.bottom_margin = Inches(0.9)
    section.left_margin = Inches(0.9)

    normal = document.styles["Normal"]
    set_font(normal)
    normal.font.size = Pt(BODY_SIZE_PT)
    normal.font.color.rgb = INK
    normal.paragraph_format.space_after = Pt(7)
    normal.paragraph_format.line_spacing = 1.15

    title = document.styles["Title"]
    set_font(title)
    title.font.size = Pt(24)
    title.font.bold = True
    title.font.color.rgb = BRAND
    title.paragraph_format.space_after = Pt(4)

    heading1 = document.styles["Heading 1"]
    set_font(heading1)
    heading1.font.size = Pt(15)
    heading1.font.bold = True
    heading1.font.color.rgb = BRAND
    heading1.paragraph_format.space_before = Pt(14)
    heading1.paragraph_format.space_after = Pt(5)
    heading1.paragraph_format.keep_with_next = True

    caption = document.styles["Caption"]
    set_font(caption)
    caption.font.size = Pt(9)
    caption.font.italic = True
    caption.font.color.rgb = MUTED_INK
    caption.paragraph_format.space_after = Pt(10)

    for style_name in ("List Bullet", "List Number"):
        style = document.styles[style_name]
        set_font(style)
        style.font.size = Pt(BODY_SIZE_PT)
        style.paragraph_format.space_after = Pt(4)
        style.paragraph_format.line_spacing = 1.15


def add_tagline(document: DocumentType, value: str) -> None:
    paragraph = document.add_paragraph()
    paragraph.paragraph_format.space_after = Pt(10)
    add_markup(paragraph, value)
    for run in paragraph.runs:
        run.font.size = Pt(12)
        run.font.color.rgb = MUTED_INK


def body_word_count(source_text: str) -> int:
    """Count prose words, excluding figure links, table rows, and captions."""

    words = 0
    for line in source_text.splitlines():
        stripped = line.strip()
        if not stripped or stripped == "---":
            continue
        if stripped.startswith("|") or stripped.startswith("!["):
            continue
        if stripped.startswith("*Figure ") or stripped.startswith("**Table "):
            continue
        words += len(re.sub(r"[#*`]", "", stripped).split())
    return words


def render(source: Path, output: Path) -> dict[str, int]:
    source_text = source.read_text()

    words = body_word_count(source_text)
    figures = len(re.findall(r"^!\[[^\]]*\]\([^)]+\)$", source_text, flags=re.MULTILINE))
    tables = len(re.findall(r"^\*\*Table \d+\.", source_text, flags=re.MULTILINE))
    if words > MAX_WORDS:
        raise ValueError(f"public brief exceeds the word budget: {words} > {MAX_WORDS}")
    if figures > MAX_FIGURES:
        raise ValueError(f"public brief exceeds the figure budget: {figures} > {MAX_FIGURES}")
    if tables > MAX_TABLES:
        raise ValueError(f"public brief exceeds the table budget: {tables} > {MAX_TABLES}")

    document = Document()
    configure(document)
    title_line = next((line[2:].strip() for line in source_text.splitlines() if line.startswith("# ")), "")
    document.core_properties.title = title_line
    document.core_properties.subject = "Public capability brief"
    document.core_properties.author = "Metabolon, Inc."
    document.core_properties.category = "Public"

    lines = source_text.splitlines()
    source_root = source.parent
    index = 0
    tagline_pending = False
    with tempfile.TemporaryDirectory(prefix="public_brief_") as temporary:
        temp_dir = Path(temporary)
        while index < len(lines):
            line = lines[index].strip()
            if not line or line == "---":
                index += 1
                continue
            if line.startswith("# "):
                paragraph = document.add_paragraph(style="Title")
                add_markup(paragraph, line[2:])
                tagline_pending = True
                index += 1
                continue
            if line.startswith("## "):
                document.add_paragraph(line[3:].strip(), style="Heading 1")
                tagline_pending = False
                index += 1
                continue
            figure = FIGURE.match(line)
            if figure:
                add_figure(
                    document,
                    source_root,
                    figure.group("path"),
                    figure.group("alt"),
                    temp_dir,
                    max_height_inches=MAX_FIGURE_HEIGHT_INCHES,
                )
                index += 1
                continue
            if line.startswith("|") and index + 1 < len(lines) and TABLE_SEPARATOR.match(lines[index + 1].strip()):
                rows = [split_table_row(line)]
                index += 2
                while index < len(lines) and lines[index].strip().startswith("|"):
                    rows.append(split_table_row(lines[index]))
                    index += 1
                add_table(document, rows, total_inches=TEXT_WIDTH_INCHES)
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
            if numbered:
                paragraph = document.add_paragraph(style="List Number")
                add_markup(paragraph, numbered.group(2))
                index += 1
                continue
            if tagline_pending:
                add_tagline(document, line)
                index += 1
                continue
            paragraph = document.add_paragraph()
            add_markup(paragraph, line)
            index += 1

    output.parent.mkdir(parents=True, exist_ok=True)
    document.save(output)
    print(
        f"[public-brief-docx] wrote {output} "
        f"(words={words}/{MAX_WORDS}; figures={figures}/{MAX_FIGURES}; tables={tables}/{MAX_TABLES})"
    )
    return {"words": words, "figures": figures, "tables": tables}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("docs/public_brief/public_brief.md"))
    parser.add_argument("--output", type=Path, default=Path("docs/public_brief/public_brief.docx"))
    args = parser.parse_args()
    if not args.source.exists():
        raise SystemExit(f"source not found: {args.source}")
    render(args.source, args.output)


if __name__ == "__main__":
    main()
