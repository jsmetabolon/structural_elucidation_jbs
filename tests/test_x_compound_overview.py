"""Smoke test for the generated X-compound prediction overview document."""
from pathlib import Path

from docx import Document

from tools.build_x_compound_overview_docx import build


def test_overview_docx_contains_expected_structure(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parent.parent
    output = tmp_path / "overview.docx"
    build(root / "docs/x_compound_prediction_overview.md", output, template=None)

    document = Document(output)
    headings = [paragraph.text for paragraph in document.paragraphs if paragraph.style.name == "Heading 1"]
    assert document.core_properties.title == "X-Compound Structure Prediction Pipeline README"
    assert document.paragraphs[0].text == "X-Compound Structure Prediction Pipeline"
    assert document.paragraphs[0].runs[0].font.size.pt == 24
    assert "Prediction And Evidence Strategy" in headings
    assert "Current Validation Snapshot" in headings
    assert len(document.tables) == 6
    table_text = "\n".join(
        cell.text for table in document.tables for row in table.rows for cell in row.cells
    )
    assert "run_full_pipeline_ec2.sh" in table_text

