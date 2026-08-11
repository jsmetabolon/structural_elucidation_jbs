"""Tests for the public capability brief and its release gate.

The brief is customer-facing, so the leak scanner is the part that actually
matters here: these tests assert not only that the committed brief is clean,
but that the scanner fires on each category of disclosure it is meant to catch.
A guard that never triggers would pass silently forever.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from docx import Document

from tools.brand import FAMILY
from tools.build_public_brief_docx import MAX_FIGURES, MAX_TABLES, MAX_WORDS, render

ROOT = Path(__file__).resolve().parent.parent
BRIEF_DIR = ROOT / "docs/public_brief"
SOURCE = BRIEF_DIR / "public_brief.md"
VERIFIER = ROOT / "tools/verify_public_brief.py"


def run_verifier(brief_dir: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(VERIFIER), "--brief-dir", str(brief_dir)],
        capture_output=True,
        text=True,
        cwd=ROOT,
    )


@pytest.fixture
def brief_copy(tmp_path: Path) -> Path:
    destination = tmp_path / "public_brief"
    shutil.copytree(BRIEF_DIR, destination)
    return destination


def test_public_brief_docx(tmp_path: Path) -> None:
    output = tmp_path / "public_brief.docx"
    counts = render(SOURCE, output)

    assert counts["words"] <= MAX_WORDS
    assert counts["figures"] <= MAX_FIGURES
    assert counts["tables"] <= MAX_TABLES

    document = Document(output)
    assert document.styles["Normal"].font.name == FAMILY
    assert document.core_properties.category == "Public"
    assert len(document.inline_shapes) == counts["figures"]
    assert len(document.tables) == counts["tables"]


def test_committed_brief_passes_the_release_gate() -> None:
    result = run_verifier(BRIEF_DIR)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "[brief-verify] PASS" in result.stdout


@pytest.mark.parametrize(
    ("leak", "expected"),
    [
        ("It relies on Morgan fingerprints.", "fingerprint"),
        ("Candidates pass a 0.05 Da gate.", "threshold"),
        ("The stack follows the Spec2Mol lineage.", "spec2mol"),
        ("Controls were drawn from library 209.", "library number"),
        ("Inference ran on EC2.", "ec2"),
        ("Pipeline v21 produced these results.", "version identifier"),
    ],
)
def test_scanner_catches_disclosure(brief_copy: Path, leak: str, expected: str) -> None:
    source = brief_copy / "public_brief.md"
    source.write_text(source.read_text() + f"\n{leak}\n")

    result = run_verifier(brief_copy)
    assert result.returncode == 1
    assert expected in result.stdout.lower()


def test_scanner_catches_unsupported_statistics(brief_copy: Path) -> None:
    source = brief_copy / "public_brief.md"
    source.write_text(source.read_text() + "\nRecovery reached 99.4% on customer data.\n")

    result = run_verifier(brief_copy)
    assert result.returncode == 1
    assert "99.4%" in result.stdout
    assert "does not resolve to a row" in result.stdout


def test_scanner_catches_overclaims(brief_copy: Path) -> None:
    source = brief_copy / "public_brief.md"
    source.write_text(source.read_text() + "\nResults are a calibrated probability of correctness.\n")

    result = run_verifier(brief_copy)
    assert result.returncode == 1
    assert "overclaim" in result.stdout


def test_brief_does_not_disturb_the_paper_package() -> None:
    """The brief is an additional artifact; the manuscript package must still verify."""

    result = subprocess.run(
        [sys.executable, "tools/verify_paper_package.py"],
        capture_output=True,
        text=True,
        cwd=ROOT,
    )
    assert result.returncode == 0, result.stdout + result.stderr
