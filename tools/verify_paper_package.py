#!/usr/bin/env python3
"""Verify manuscript claims, artifacts, citations, figures, and submission blockers."""
from __future__ import annotations

import argparse
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from docx import Document
from docx.enum.text import WD_LINE_SPACING
from docx.oxml.ns import qn


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper-dir", type=Path, default=Path("docs/paper"))
    parser.add_argument("--strict-human-inputs", action="store_true")
    args = parser.parse_args()
    paper = args.paper_dir
    manuscript_path = paper / "manuscript.md"
    errors: list[str] = []
    warnings: list[str] = []

    required = [
        manuscript_path,
        paper / "manuscript.docx",
        paper / "README.md",
        paper / "REVISION_LOG.md",
        paper / "supplement/supplementary_information.md",
        paper / "supplement/claim_evidence.csv",
        paper / "supplement/baseline_claims.json",
        paper / "data/pathway/audit_pathway_summary.json",
        paper / "data/pathway/audit_pathway_loo.json",
        paper / "data/structure/audit_structure_summary.json",
        paper / "data/presentation/audit_presentation_summary.json",
    ]
    for path in required:
        if not path.exists():
            errors.append(f"missing required artifact: {path}")
    if not manuscript_path.exists():
        raise SystemExit("\n".join(errors))

    manuscript = manuscript_path.read_text()
    abstract = manuscript.split("## Abstract", 1)[1].split("---", 1)[0]
    abstract_words = len(abstract.split())
    if abstract_words > 200:
        errors.append(f"Nature summary exceeds ~200 words: {abstract_words}")

    stale = {
        "94.5%": "old circular reference-DB pathway headline",
        "89.7%": "old circular reference-DB subpathway headline",
        "essentially indistinguishable": "unsupported exact-vs-predicted equivalence language",
        "NMD on the Fly": "internal commercial language",
        "fee-for-service": "internal commercial language",
    }
    for phrase, reason in stale.items():
        if phrase in manuscript:
            errors.append(f"stale phrase '{phrase}': {reason}")

    # Figures carry their own text, and an unreferenced one is invisible to the link
    # checks below, so a retracted headline can survive in the directory indefinitely.
    for figure_path in sorted((paper / "figures").glob("*")):
        if not figure_path.is_file():
            continue
        try:
            figure_text = figure_path.read_text()
        except (OSError, UnicodeDecodeError):
            continue
        for phrase, reason in stale.items():
            if phrase in figure_text:
                errors.append(f"stale phrase '{phrase}' in {figure_path.name}: {reason}")

    expected_claims = (
        "44/92",
        "59/92",
        "68/92",
        "336/644",
        "121/644",
        "503/909",
        "214/909",
        "598/3,897",
        "523/598",
        "519/590",
    )
    for value in expected_claims:
        if value not in manuscript:
            errors.append(f"manuscript missing expected audited count: {value}")

    inline = {
        int(number)
        for group in re.findall(r"\[(\d+(?:\s*,\s*\d+)*)\]", manuscript)
        for number in re.findall(r"\d+", group)
    }
    listed = {int(number) for number in re.findall(r"^(\d+)\.\s", manuscript, flags=re.MULTILINE)}
    if inline - listed:
        errors.append(f"inline citations missing from references: {sorted(inline - listed)}")
    if listed - inline:
        errors.append(f"references never cited: {sorted(listed - inline)}")

    figure_links = re.findall(r"!\[[^\]]*\]\((figures/[^)]+)\)", manuscript)
    figure_numbers = []
    for link in figure_links:
        path = paper / link
        if not path.exists():
            errors.append(f"broken figure link: {link}")
            continue
        try:
            ET.parse(path)
        except ET.ParseError as exc:
            errors.append(f"invalid SVG {link}: {exc}")
        match = re.search(r"fig(\d+)", path.name)
        if match:
            figure_numbers.append(int(match.group(1)))
    if figure_numbers != sorted(figure_numbers):
        errors.append(f"main figures are not in numeric order: {figure_numbers}")

    docx_path = paper / "manuscript.docx"
    if docx_path.exists():
        word = Document(docx_path)
        if len(word.core_properties.title) > 75:
            errors.append(f"Nature title exceeds 75 characters: {word.core_properties.title}")
        if word.styles["Normal"].font.name != "Times New Roman":
            errors.append("Nature DOCX body is not Times New Roman")
        if word.styles["Normal"].font.size.pt != 12:
            errors.append("Nature DOCX body is not 12 point")
        if word.styles["Normal"].paragraph_format.line_spacing_rule != WD_LINE_SPACING.DOUBLE:
            errors.append("Nature DOCX body is not double-spaced")
        if len(word.inline_shapes) != 5 or len(word.tables) != 1:
            errors.append(
                f"Nature display-item packaging changed: {len(word.inline_shapes)} figures, "
                f"{len(word.tables)} tables"
            )
        line_numbers = word.sections[0]._sectPr.find(qn("w:lnNumType"))
        if line_numbers is None:
            errors.append("Nature DOCX lacks line numbering")

    human_placeholders = sorted(set(re.findall(r"\[REQUIRED[^\]]*\]", manuscript)))
    if "[additional authors]" in manuscript:
        human_placeholders.append("[additional authors]")
    if "[email]" in manuscript:
        human_placeholders.append("[email]")
    if human_placeholders:
        message = f"human inputs still required: {human_placeholders}"
        (errors if args.strict_human_inputs else warnings).append(message)

    structure_path = paper / "data/structure/audit_structure_summary.json"
    if structure_path.exists():
        structure = json.loads(structure_path.read_text())
        if structure["p0_p1"]["status"] != "not_evaluated":
            warnings.append(f"P0/P1 status changed: {structure['p0_p1']['status']}")
    pathway_path = paper / "data/pathway/audit_pathway_summary.json"
    if pathway_path.exists():
        pathway = json.loads(pathway_path.read_text())
        direct = pathway["reference_db_evaluation"]["by_assignment_source"]["exact_reference_lookup"]
        if direct["super"]["accuracy"] != 1.0:
            errors.append("unexpected direct-lookup consistency value")

        # Coverage counts are quoted verbatim in Results and have drifted from the
        # audit artifact before; assert each one is still present in the prose.
        coverage = pathway["coverage"]
        quoted = {
            "exact reference match": coverage["status_counts"]["exact_reference_match"],
            "strong-similarity vote": coverage["status_counts"]["strong_similarity_vote"],
            "weak-similarity backfill": coverage["status_counts"]["weak_similarity_backfill"],
            "high pathway-confidence tier": coverage["tier_counts"]["high"],
            "medium pathway-confidence tier": coverage["tier_counts"]["medium"],
            "low pathway-confidence tier": coverage["tier_counts"]["low"],
        }
        for label, count in quoted.items():
            if f"{count:,}" not in manuscript:
                errors.append(f"{label} count drifted from the audit artifact: expected {count:,}")
        median = coverage["best_similarity_pctiles"]["p50"]
        if f"median {median}" not in manuscript:
            errors.append(f"best-reference similarity median drifted from the audit artifact: expected {median}")

    print(f"[paper-verify] abstract words: {abstract_words}")
    print(f"[paper-verify] main figures: {len(figure_links)}")
    for warning in warnings:
        print(f"[paper-verify] WARN {warning}")
    if errors:
        for error in errors:
            print(f"[paper-verify] ERROR {error}")
        raise SystemExit(1)
    print("[paper-verify] PASS")


if __name__ == "__main__":
    main()
