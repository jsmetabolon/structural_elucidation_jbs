#!/usr/bin/env python3
"""One-command regeneration of manuscript audits, figures, supplement, and public brief.

Full regeneration requires the controlled report bundle and reference database.
``--verify-only`` validates the committed paper package and public brief from a
clean checkout. ``--dry-run`` prints the commands without executing them.

The public brief is rebuilt from the committed audit JSON rather than from the
controlled bundle, so ``--brief-only`` works in any checkout.

Examples:
    python3 tools/rebuild_paper_artifacts.py --verify-only
    python3 tools/rebuild_paper_artifacts.py --brief-only
    python3 tools/rebuild_paper_artifacts.py \
        --reports-root _audit_pull/reports \
        --db _bundle/metabolite_synonyms.sqlite \
        --known-x data/known_X.csv
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path


def run(command: list[str], *, dry_run: bool) -> None:
    print("+", shlex.join(command))
    if not dry_run:
        subprocess.run(command, check=True)


def build_public_brief(paper_dir: Path, brief_dir: Path, *, python: str, dry_run: bool) -> None:
    """Regenerate the public brief's figures, Word file, and release gate.

    Reads the committed audit JSON under ``paper_dir`` rather than the controlled
    bundle, so this runs in any checkout.
    """

    data = paper_dir / "data"
    run(
        [
            python,
            "tools/public_brief_figures.py",
            "--structure",
            str(data / "structure/audit_structure_summary.json"),
            "--presentation",
            str(data / "presentation/audit_presentation_summary.json"),
            "--out-dir",
            str(brief_dir / "figures"),
        ],
        dry_run=dry_run,
    )
    run(
        [
            python,
            "tools/build_public_brief_docx.py",
            "--source",
            str(brief_dir / "public_brief.md"),
            "--output",
            str(brief_dir / "public_brief.docx"),
        ],
        dry_run=dry_run,
    )
    run(
        [python, "tools/verify_public_brief.py", "--brief-dir", str(brief_dir)],
        dry_run=dry_run,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--reports-root", type=Path, default=Path("_audit_pull/reports"))
    parser.add_argument("--db", type=Path, default=Path("_bundle/metabolite_synonyms.sqlite"))
    parser.add_argument("--known-x", type=Path, default=Path("data/known_X.csv"))
    parser.add_argument("--paper-dir", type=Path, default=Path("docs/paper"))
    parser.add_argument("--brief-dir", type=Path, default=Path("docs/public_brief"))
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--brief-only", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    python = sys.executable

    if args.verify_only:
        run(
            [python, "tools/verify_paper_package.py", "--paper-dir", str(args.paper_dir)],
            dry_run=args.dry_run,
        )
        run(
            [python, "tools/verify_public_brief.py", "--brief-dir", str(args.brief_dir)],
            dry_run=args.dry_run,
        )
        return

    if args.brief_only:
        build_public_brief(args.paper_dir, args.brief_dir, python=python, dry_run=args.dry_run)
        return

    missing = [path for path in (args.reports_root, args.db, args.known_x) if not path.exists()]
    if missing:
        raise SystemExit(
            "Full regeneration requires controlled inputs. Missing: "
            + ", ".join(str(path) for path in missing)
            + ". Use --verify-only to validate committed artifacts."
        )

    data = args.paper_dir / "data"
    supplement = args.paper_dir / "supplement"
    figures = args.paper_dir / "figures"
    presentation = data / "presentation/audit_presentation_summary.json"
    pathway = data / "pathway/audit_pathway_summary.json"
    loo = data / "pathway/audit_pathway_loo.json"
    structure = data / "structure/audit_structure_summary.json"
    claims = supplement / "baseline_claims.json"
    selective = data / "selective_prediction.json"

    run(
        [
            python,
            "tools/audit_presentation.py",
            str(args.reports_root),
            "--out-dir",
            str(data / "presentation"),
            "--summary-only",
        ],
        dry_run=args.dry_run,
    )
    run(
        [
            python,
            "tools/audit_pathway_predictions.py",
            "--reports-root",
            str(args.reports_root),
            "--known-x",
            str(args.known_x),
            "--db",
            str(args.db),
            "--out-dir",
            str(data / "pathway"),
            "--loo",
        ],
        dry_run=args.dry_run,
    )
    run(
        [
            python,
            "tools/audit_structure_predictions.py",
            "--reports-root",
            str(args.reports_root),
            "--out-dir",
            str(data / "structure"),
        ],
        dry_run=args.dry_run,
    )
    run(
        [
            python,
            "tools/build_paper_claim_audit.py",
            "--reports-root",
            str(args.reports_root),
            "--presentation-summary",
            str(presentation),
            "--pathway-summary",
            str(pathway),
            "--loo-summary",
            str(loo),
            "--out-dir",
            str(supplement),
        ],
        dry_run=args.dry_run,
    )
    run(
        [
            python,
            "tools/analyze_selective_prediction.py",
            "--structure-controls",
            str(data / "structure/audit_structure_controls.csv"),
            "--pathway-accuracy",
            str(data / "pathway/audit_pathway_accuracy.csv"),
            "--out",
            str(selective),
        ],
        dry_run=args.dry_run,
    )
    run(
        [
            python,
            "tools/paper_figures.py",
            "--structure",
            str(structure),
            "--pathway",
            str(pathway),
            "--loo",
            str(loo),
            "--presentation",
            str(presentation),
            "--claims",
            str(claims),
            "--selective",
            str(selective),
            "--out-dir",
            str(figures),
        ],
        dry_run=args.dry_run,
    )
    run(
        [
            python,
            "tools/build_paper_supplement.py",
            "--structure",
            str(structure),
            "--pathway",
            str(pathway),
            "--loo",
            str(loo),
            "--presentation",
            str(presentation),
            "--claims",
            str(claims),
            "--out",
            str(supplement / "supplementary_information.md"),
        ],
        dry_run=args.dry_run,
    )
    run(
        [
            python,
            "tools/build_manuscript_docx.py",
            "--source",
            str(args.paper_dir / "manuscript.md"),
            "--output",
            str(args.paper_dir / "manuscript.docx"),
        ],
        dry_run=args.dry_run,
    )
    build_public_brief(args.paper_dir, args.brief_dir, python=python, dry_run=args.dry_run)
    run(
        [python, "tools/verify_paper_package.py", "--paper-dir", str(args.paper_dir)],
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
