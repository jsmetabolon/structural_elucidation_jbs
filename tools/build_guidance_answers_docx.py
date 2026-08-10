#!/usr/bin/env python3
"""Build the writing-guidance answer worksheets as Word documents.

The guidance arrived as a Word template with prompts to fill in, so the answers
go back in the same form. Markdown parsing is shared with the Nature manuscript
builder; only the typography differs. These are internal working documents, so
they use the readable Calibri setting of the other internal Word outputs rather
than the double-spaced submission format.

Usage:
    python3 tools/build_guidance_answers_docx.py
    python3 tools/build_guidance_answers_docx.py --source docs/writing_guidance_manuscript.md
"""
from __future__ import annotations

import argparse
import re
import sys
import tempfile
from pathlib import Path

from docx import Document
from docx.document import Document as DocumentType
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Pt, RGBColor

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tools.build_manuscript_docx import (  # noqa: E402
    TABLE_SEPARATOR,
    add_markup,
    add_table,
    split_table_row,
)

BRAND = RGBColor(0x1F, 0x4E, 0x79)
INK = RGBColor(0x1F, 0x29, 0x33)

DEFAULT_SOURCES = (
    Path("docs/writing_guidance_public_brief.md"),
    Path("docs/writing_guidance_manuscript.md"),
)


def configure(document: DocumentType) -> None:
    normal = document.styles["Normal"]
    normal.font.name = "Calibri"
    normal.font.size = Pt(10.5)
    normal.font.color.rgb = INK
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.15

    title = document.styles["Title"]
    title.font.name = "Calibri"
    title.font.size = Pt(22)
    title.font.bold = True
    title.font.color.rgb = BRAND
    title.paragraph_format.space_after = Pt(10)

    heading1 = document.styles["Heading 1"]
    heading1.font.name = "Calibri"
    heading1.font.size = Pt(15)
    heading1.font.bold = True
    heading1.font.color.rgb = BRAND
    heading1.paragraph_format.space_before = Pt(16)
    heading1.paragraph_format.space_after = Pt(6)
    heading1.paragraph_format.keep_with_next = True

    heading2 = document.styles["Heading 2"]
    heading2.font.name = "Calibri"
    heading2.font.size = Pt(12)
    heading2.font.bold = True
    heading2.font.color.rgb = INK
    heading2.paragraph_format.space_before = Pt(10)
    heading2.paragraph_format.space_after = Pt(4)
    heading2.paragraph_format.keep_with_next = True

    for style_name in ("List Bullet", "List Number"):
        style = document.styles[style_name]
        style.font.name = "Calibri"
        style.font.size = Pt(10.5)
        style.paragraph_format.space_after = Pt(4)
        style.paragraph_format.line_spacing = 1.15


def render(source: Path, output: Path) -> dict[str, int]:
    lines = source.read_text().splitlines()
    document = Document()
    configure(document)

    tables = 0
    index = 0
    with tempfile.TemporaryDirectory(prefix="guidance_docx_"):
        while index < len(lines):
            line = lines[index].strip()
            if not line or line == "---":
                index += 1
                continue
            if line.startswith("# "):
                paragraph = document.add_paragraph(style="Title")
                paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
                add_markup(paragraph, line[2:])
                index += 1
                continue
            if line.startswith("### "):
                document.add_paragraph(line[4:].strip(), style="Heading 2")
                index += 1
                continue
            if line.startswith("## "):
                document.add_paragraph(line[3:].strip(), style="Heading 1")
                index += 1
                continue
            if line.startswith("|") and index + 1 < len(lines) and TABLE_SEPARATOR.match(lines[index + 1].strip()):
                rows = [split_table_row(line)]
                index += 2
                while index < len(lines) and lines[index].strip().startswith("|"):
                    rows.append(split_table_row(lines[index]))
                    index += 1
                add_table(document, rows)
                tables += 1
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
            paragraph = document.add_paragraph()
            add_markup(paragraph, line)
            index += 1

    document.core_properties.title = next(
        (line[2:].strip() for line in lines if line.startswith("# ")), source.stem
    )
    document.core_properties.subject = "White paper writing guidance: filled-in answers"
    document.core_properties.author = "Metabolon, Inc."
    document.core_properties.category = "Internal"
    output.parent.mkdir(parents=True, exist_ok=True)
    document.save(output)
    counts = {"tables": tables, "paragraphs": len(document.paragraphs)}
    print(f"[guidance-docx] wrote {output} (tables={tables}; paragraphs={counts['paragraphs']})")
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, default=None, help="one worksheet; omit to build both")
    parser.add_argument("--output", type=Path, default=None, help="output path, defaulting to the source with .docx")
    args = parser.parse_args()

    sources = [args.source] if args.source else list(DEFAULT_SOURCES)
    for source in sources:
        if not source.exists():
            raise SystemExit(f"source not found: {source}")
        output = args.output if args.output and args.source else source.with_suffix(".docx")
        render(source, output)


if __name__ == "__main__":
    main()
