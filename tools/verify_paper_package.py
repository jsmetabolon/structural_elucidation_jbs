#!/usr/bin/env python3
"""Verify manuscript claims, artifacts, citations, figures, and submission blockers."""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

from docx import Document
from docx.enum.text import WD_LINE_SPACING
from docx.oxml.ns import qn

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tools.brand import FAMILY, font_is_resolvable  # noqa: E402
from tools.figure_formats import figure_error  # noqa: E402

# Documents whose job is to record what was withdrawn, and which therefore have to
# quote the withdrawn wording. Matched by filename so the exemption survives being
# pointed at a copied tree.
STALE_EXEMPT_NAMES = frozenset(
    {
        "REVISION_LOG.md",
        "pathway_audit_findings.md",
        "publication_review.md",
    }
)


def audit_statistics(
    structure: dict[str, Any],
    pathway: dict[str, Any],
    loo: dict[str, Any],
    presentation: dict[str, Any],
) -> dict[str, str]:
    """Headline values as documents render them, derived from the audit artifacts."""

    rank = structure["rank_k"]["0.70"]
    neighbor = pathway["accuracy"]["by_assignment_source"]["neighbor_consensus"]
    ceiling = pathway["accuracy"]["by_reference_membership"]["reference_structure"]
    msi = {key: int(value) for key, value in presentation["msi_level_distribution"].items()}
    shown = sum(msi[key] for key in ("1", "2", "3"))
    withheld = sum(msi[key] for key in ("4", "5"))
    return {
        "rank-1 recovery": f"{100 * rank['top1']['accuracy']:.1f}%",
        "top-3 recovery": f"{100 * rank['top3']['accuracy']:.1f}%",
        "top-10 recovery": f"{100 * rank['top10']['accuracy']:.1f}%",
        "reports showing a structure": f"{shown:,}",
        "reports withholding a structure": f"{withheld:,}",
        "neighbor-consensus super-pathway agreement": f"{100 * neighbor['super']['accuracy']:.1f}%",
        "neighbor-consensus subpathway agreement": f"{100 * neighbor['sub']['accuracy']:.1f}%",
        "leave-one-out super-pathway accuracy": f"{100 * loo['super']['accuracy']:.1f}%",
        "leave-one-out subpathway accuracy": f"{100 * loo['sub']['accuracy']:.1f}%",
        "super-pathway annotation ceiling": f"{100 * ceiling['super']['accuracy']:.1f}%",
        "subpathway annotation ceiling": f"{100 * ceiling['sub']['accuracy']:.1f}%",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper-dir", type=Path, default=Path("docs/paper"))
    parser.add_argument("--docs-dir", type=Path, default=Path("docs"))
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
    # Retracted claims have twice survived in files nothing checked: an unreferenced
    # figure, and an internal findings document that still led with the withdrawn
    # headline. Scan the whole docs tree and every figure, not just the manuscript.
    scan_targets = {manuscript_path}
    scan_targets.update(path for path in (paper / "figures").glob("*") if path.is_file())
    if args.docs_dir.exists():
        scan_targets.update(args.docs_dir.rglob("*.md"))
    for path in sorted(scan_targets):
        # These records have to quote the retracted numbers in order to document the
        # retraction, so scanning them would fire on the audit trail itself.
        if path.name in STALE_EXEMPT_NAMES:
            continue
        try:
            content = path.read_text()
        except (OSError, UnicodeDecodeError):
            continue
        for phrase, reason in stale.items():
            if phrase in content:
                errors.append(f"stale phrase '{phrase}' in {path}: {reason}")

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
        problem = figure_error(path)
        if problem:
            errors.append(f"{link}: {problem}")
        match = re.search(r"fig(\d+)", path.name)
        if match:
            figure_numbers.append(int(match.group(1)))
    if figure_numbers != sorted(figure_numbers):
        errors.append(f"main figures are not in numeric order: {figure_numbers}")

    # The figures name the brand face, but naming it is not having it. Word files
    # are built by rasterising those SVGs, so a checkout where the family is
    # missing produces a document whose figures are silently set in a fallback and
    # whose every other check still passes.
    if not font_is_resolvable():
        errors.append(f"{FAMILY} is unavailable to the renderer; figures would rasterise in a fallback face")

    docx_path = paper / "manuscript.docx"
    if docx_path.exists():
        word = Document(docx_path)
        if len(word.core_properties.title) > 75:
            errors.append(f"Nature title exceeds 75 characters: {word.core_properties.title}")
        if word.styles["Normal"].font.name != FAMILY:
            errors.append(f"Nature DOCX body is not {FAMILY}")
        if word.styles["Normal"].font.size.pt != 12:
            errors.append("Nature DOCX body is not 12 point")
        if word.styles["Normal"].paragraph_format.line_spacing_rule != WD_LINE_SPACING.DOUBLE:
            errors.append("Nature DOCX body is not double-spaced")
        if len(word.inline_shapes) != 6 or len(word.tables) != 1:
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

    # Three documents quote the same headline statistics, and all three have to move
    # when the audit moves. Only the manuscript was checked before, so the paper README
    # and the operations overview could diverge from it without anything noticing.
    loo_path = paper / "data/pathway/audit_pathway_loo.json"
    presentation_path = paper / "data/presentation/audit_presentation_summary.json"
    if all(path.exists() for path in (structure_path, pathway_path, loo_path, presentation_path)):
        statistics = audit_statistics(
            json.loads(structure_path.read_text()),
            json.loads(pathway_path.read_text()),
            json.loads(loo_path.read_text()),
            json.loads(presentation_path.read_text()),
        )
        quoting_documents = [
            manuscript_path,
            paper / "README.md",
            args.docs_dir / "x_compound_prediction_overview.md",
        ]
        for document in quoting_documents:
            if not document.exists():
                continue
            content = document.read_text()
            for label, value in statistics.items():
                if value not in content:
                    errors.append(f"{document} does not carry the current {label} ({value})")

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
