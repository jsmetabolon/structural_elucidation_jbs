#!/usr/bin/env python3
"""Verify that the public capability brief is safe to release.

This is the release gate for a customer-facing document, so it checks three
independent things:

1. Disclosure. The brief must contain none of the recipe-level detail that would
   let a reader reconstruct the workflow, and no internal identifiers. This
   mirrors the public-report compliance scan in ``tools/audit_presentation.py``,
   which asserts that observed-MS2 detail never reaches a public page.
2. Accuracy. Every statistic in the brief must resolve to a row in the audited
   claim ledger ``docs/paper/supplement/claim_evidence.csv``, so marketing copy
   cannot drift away from what was actually measured.
3. Restraint. Banned overclaim phrases, and budgets that keep the deliverable
   inside its page limit.

Usage:
    python3 tools/verify_public_brief.py
    python3 tools/verify_public_brief.py --check-pages
    python3 tools/verify_public_brief.py --strict-human-inputs
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from docx import Document

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tools.figure_formats import figure_error  # noqa: E402

MAX_PAGES = 10
MAX_WORDS = 2600
MAX_FIGURES = 5
MAX_TABLES = 3

# Recipe-level detail: publishing any of this would materially help a competitor
# rebuild the pipeline. Matched case-insensitively on word boundaries.
PROPRIETARY_TERMS = (
    # Model families, architecture, and training procedure
    "spec2mol",
    "msnovelist",
    "sirius",
    "canopus",
    "ms2query",
    "csi:fingerid",
    "fingerid",
    "gradient-boosted",
    "logistic classifier",
    "isotonic",
    "platt",
    "encoder",
    "decoder",
    "embedding bank",
    "de novo",
    "checkpoint",
    "random seed",
    # Descriptors, fingerprints, and structure formats
    "morgan",
    "tanimoto",
    "smiles",
    "inchikey",
    "smarts",
    "selfies",
    "fingerprint",
    # Named reranker inputs
    "formula-prior",
    "forward-spectrum",
    "source-rank",
    "biotransformation",
    "nce",
    # Scoring and rubric mechanics
    "vote share",
    "vote weight",
    "similarity floor",
    "mass gate",
    "neighbor vote",
    "neighbour vote",
    "rank score",
    "brier",
    "wilson",
    # Internal identifiers, controlled inputs, and infrastructure
    "known_x",
    "metabolite_synonyms",
    "p1_reference",
    "rerank csv",
    "pipeline_equalization",
    "ec2",
    "cuda",
    "pytorch",
    "wgcna",
    "first-library-wins",
)

# Regexes for shapes of disclosure that a term list cannot catch.
PROPRIETARY_PATTERNS = (
    (r"\b0\.\d+\b", "a bare decimal, which in this domain reads as a threshold or weight"),
    (r"\bv\d{2}\b", "an internal pipeline version identifier"),
    (r"\blibrar(?:y|ies)\s+\d{3}\b", "an internal chromatographic library number"),
    (r"\b\d+(?:,\d{3})*\.\d+\s*s\b", "an internal wall-clock runtime"),
    (r"\bradius-\d\b", "a fingerprint radius"),
    (r"\b\d[\d,]*-bit\b", "a fingerprint bit width"),
    (r"\bk\s*=\s*\d+\b", "a neighbour count"),
    (r"\b\d+(?:\.\d+)?\s*Da\b", "a mass-gate tolerance"),
)

# Language that would promise more than the validation supports.
OVERCLAIM_PHRASES = (
    "calibrated probability",
    "validated accuracy",
    "confirmed identification",
    "state of the art",
    "state-of-the-art",
    "guaranteed",
    "proven correct",
    "always correct",
)

# Ordered-list markers and figure/table numbering are not statistics.
NUMBER_PATTERN = re.compile(r"\b\d[\d,]*(?:\.\d+)?%?")


def load_ledger(path: Path) -> dict[str, dict[str, str]]:
    with path.open() as handle:
        return {row["claim_id"]: row for row in csv.DictReader(handle)}


def allowed_numbers(ledger: dict[str, dict[str, str]]) -> dict[str, str]:
    """Build the set of numeric strings the brief is permitted to assert.

    Includes each claim's value rendered as a percentage and as a raw count, the
    cohort size, and the numerator where the metric expresses one. Distribution
    claims additionally contribute their component sums, which is how the brief
    can cite the structure-shown and structure-withheld totals.
    """

    allowed: dict[str, str] = {}

    def record(value: str, claim_id: str) -> None:
        allowed.setdefault(value, claim_id)

    for claim_id, row in ledger.items():
        raw = row["value"]
        record(row["n"], claim_id)
        record(f"{int(row['n']):,}", claim_id)

        try:
            numeric = float(raw)
        except ValueError:
            numeric = None

        if numeric is not None:
            if 0.0 <= numeric <= 1.0:
                record(f"{numeric * 100:.1f}%", claim_id)
                record(f"{numeric * 100:.0f}%", claim_id)
            else:
                record(f"{int(numeric)}", claim_id)
                record(f"{int(numeric):,}", claim_id)

        # Numerators such as "44/92" quoted in the metric column.
        for numerator in re.findall(r"\b(\d+)\s*/\s*\d+\b", row["metric"]):
            record(numerator, claim_id)
            record(f"{int(numerator):,}", claim_id)

        # Distribution claims: allow the components and any contiguous subtotal.
        if raw.strip().startswith("{"):
            distribution = json.loads(raw)
            counts = {key: int(value) for key, value in distribution.items()}
            total = sum(counts.values())
            keys = sorted(counts)
            for start in range(len(keys)):
                for stop in range(start + 1, len(keys) + 1):
                    subtotal = sum(counts[key] for key in keys[start:stop])
                    record(f"{subtotal:,}", claim_id)
                    record(str(subtotal), claim_id)
                    if total:
                        record(f"{subtotal / total * 100:.1f}%", claim_id)
    return allowed


def statistics_in(text: str) -> list[str]:
    """Numbers the brief asserts, excluding list markers, figure/table numbers, and years."""

    candidates: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        # Drop ordered-list markers and figure/table numbering before scanning.
        stripped = re.sub(r"^\d+\.\s+", "", stripped)
        stripped = re.sub(r"\*?\*?(?:Figure|Table)\s+\d+\.", "", stripped)
        candidates.extend(NUMBER_PATTERN.findall(stripped))
    return candidates


def page_count(docx_path: Path) -> int | None:
    """Render via LibreOffice and count pages; ``None`` when LibreOffice is unavailable."""

    soffice = shutil.which("soffice") or shutil.which("libreoffice")
    if soffice is None:
        return None
    with tempfile.TemporaryDirectory(prefix="brief_pages_") as temporary:
        subprocess.run(
            [soffice, "--headless", "--convert-to", "pdf", "--outdir", temporary, str(docx_path)],
            check=True,
            capture_output=True,
        )
        pdf = Path(temporary) / f"{docx_path.stem}.pdf"
        if not pdf.exists():
            return None
        return len(re.findall(rb"/Type\s*/Page[^s]", pdf.read_bytes()))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--brief-dir", type=Path, default=Path("docs/public_brief"))
    parser.add_argument("--ledger", type=Path, default=Path("docs/paper/supplement/claim_evidence.csv"))
    parser.add_argument("--check-pages", action="store_true", help="render with LibreOffice and assert the page limit")
    parser.add_argument("--strict-human-inputs", action="store_true")
    args = parser.parse_args()

    brief_dir = args.brief_dir
    source_path = brief_dir / "public_brief.md"
    docx_path = brief_dir / "public_brief.docx"
    errors: list[str] = []
    warnings: list[str] = []

    for path in (source_path, docx_path, args.ledger):
        if not path.exists():
            errors.append(f"missing required artifact: {path}")
    if not source_path.exists():
        raise SystemExit("\n".join(errors))

    brief = source_path.read_text()
    lowered = brief.lower()

    # 1. Disclosure.
    for term in PROPRIETARY_TERMS:
        # Trailing (e)s so a plural does not slip past the singular term.
        if re.search(rf"(?<!\w){re.escape(term)}(?:e?s)?(?!\w)", lowered):
            errors.append(f"proprietary term disclosed: '{term}'")
    for pattern, reason in PROPRIETARY_PATTERNS:
        for hit in set(re.findall(pattern, brief, flags=re.IGNORECASE)):
            errors.append(f"proprietary detail disclosed: '{hit}' is {reason}")

    # 2. Accuracy.
    ledger = load_ledger(args.ledger) if args.ledger.exists() else {}
    permitted = allowed_numbers(ledger)
    for value in statistics_in(brief):
        if value not in permitted:
            errors.append(f"statistic '{value}' does not resolve to a row in {args.ledger}")

    # 3. Restraint.
    for phrase in OVERCLAIM_PHRASES:
        if phrase in lowered:
            errors.append(f"overclaim phrase: '{phrase}'")

    words = len(
        [
            word
            for line in brief.splitlines()
            if line.strip() and not line.strip().startswith(("|", "!["))
            for word in re.sub(r"[#*`]", "", line).split()
        ]
    )
    figures = re.findall(r"!\[[^\]]*\]\((figures/[^)]+)\)", brief)
    tables = re.findall(r"^\*\*Table \d+\.", brief, flags=re.MULTILINE)
    if words > MAX_WORDS:
        errors.append(f"brief exceeds the word budget: {words} > {MAX_WORDS}")
    if len(figures) > MAX_FIGURES:
        errors.append(f"brief exceeds the figure budget: {len(figures)} > {MAX_FIGURES}")
    if len(tables) > MAX_TABLES:
        errors.append(f"brief exceeds the table budget: {len(tables)} > {MAX_TABLES}")

    for link in figures:
        path = brief_dir / link
        if not path.exists():
            errors.append(f"broken figure link: {link}")
            continue
        problem = figure_error(path)
        if problem:
            errors.append(f"{link}: {problem}")

    if docx_path.exists():
        word = Document(docx_path)
        if word.styles["Normal"].font.name != "Calibri":
            errors.append("brief DOCX body is not Calibri")
        if len(word.inline_shapes) != len(figures):
            errors.append(f"DOCX embeds {len(word.inline_shapes)} figures; source references {len(figures)}")
        if word.core_properties.category != "Public":
            errors.append("brief DOCX is not marked as Public")

    pages = None
    if args.check_pages and docx_path.exists():
        pages = page_count(docx_path)
        if pages is None:
            warnings.append("LibreOffice unavailable; page count not verified")
        elif pages > MAX_PAGES:
            errors.append(f"brief exceeds the page limit: {pages} > {MAX_PAGES}")

    placeholders = sorted(set(re.findall(r"\[REQUIRED[^\]]*\]", brief)))
    if placeholders:
        message = f"human inputs still required before release: {len(placeholders)} placeholder(s)"
        (errors if args.strict_human_inputs else warnings).append(message)

    print(f"[brief-verify] words: {words}/{MAX_WORDS}")
    print(f"[brief-verify] figures: {len(figures)}/{MAX_FIGURES}; tables: {len(tables)}/{MAX_TABLES}")
    print(f"[brief-verify] statistics checked against ledger: {len(statistics_in(brief))}")
    if pages is not None:
        print(f"[brief-verify] pages: {pages}/{MAX_PAGES}")
    for warning in warnings:
        print(f"[brief-verify] WARN {warning}")
    if errors:
        for error in errors:
            print(f"[brief-verify] ERROR {error}")
        raise SystemExit(1)
    print("[brief-verify] PASS")


if __name__ == "__main__":
    main()
