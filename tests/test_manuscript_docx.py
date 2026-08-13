"""Smoke test for the Nature-formatted manuscript DOCX."""
from pathlib import Path

from docx import Document
from docx.enum.text import WD_LINE_SPACING
from docx.oxml.ns import qn

from tools.build_manuscript_docx import BODY_FONT, HEADING_COLOR, render


def test_nature_manuscript_docx(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parent.parent
    output = tmp_path / "manuscript.docx"
    render(root / "docs/paper/manuscript.md", output)

    document = Document(output)
    assert len(document.core_properties.title) <= 75
    # Against the builder's own switch rather than a named face, so producing a
    # journal-plain submission stays a one-flag change.
    assert document.styles["Normal"].font.name == BODY_FONT
    assert document.styles["Title"].font.color.rgb == HEADING_COLOR
    assert document.styles["Heading 1"].font.color.rgb == HEADING_COLOR
    assert document.styles["Normal"].font.size.pt == 12
    assert document.styles["Normal"].paragraph_format.line_spacing_rule == WD_LINE_SPACING.DOUBLE
    assert len(document.inline_shapes) == 6
    assert len(document.tables) == 1
    headings = [p.text for p in document.paragraphs if p.style.name == "Heading 1"]
    assert headings.index("Results") < headings.index("References") < headings.index("Methods")
    assert headings.index("Methods") < headings.index("Data and code availability")
    line_numbers = document.sections[0]._sectPr.find(qn("w:lnNumType"))
    assert line_numbers is not None
    assert line_numbers.get(qn("w:countBy")) == "1"

