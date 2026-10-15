#!/usr/bin/env python3
"""Parity report: vendored neighbor vote vs delivered pathway assignments.

The delivered assignments were produced by a legacy EC2-side builder that also
applied undocumented molecular-weight priors, SMARTS rules and amino-acid
overrides. The vendored vote (``pipeline.enrichment.neighbor_vote``) implements
only the documented mechanism. This tool quantifies the difference so the owner
can decide whether the undocumented extras mattered enough to reconstruct.

Inputs (both fetched assets, not in a plain checkout):
- metabolite synonym SQLite (``./fetch_assets.sh``),
- delivered rerank enrichment outputs
  ``models/pipeline_equalization/rerank_lib{lib}/unified_rerank_all_top10_with_pathways.csv``.

For every unique predicted SMILES it re-votes with the vendored mechanism and
compares super/sub labels, status and tier against the delivered columns.

Usage::

    python tools/pathway_vote_parity.py \
        --equalization-root models/pipeline_equalization \
        --out-json docs/paper/data/pathway/vote_parity_report.json \
        --out-md docs/paper/data/pathway/vote_parity_report.md
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
REFERENCE_DB_CANDIDATES = (
    REPO_ROOT / "Metabolite Synonym Database" / "metabolite_synonyms_merged_with_mtx_giga.sqlite",
    REPO_ROOT / "vendor" / "metabolite_synonyms.sqlite",
)
DEFAULT_LIBS = (209, 305, 400, 402)
DELIVERED_COLUMNS = [
    "predicted_pathway",
    "predicted_subpathway",
    "pathway_confidence_tier",
    "pathway_assignment_status",
    "pathway_best_reference_similarity",
]


def norm(value: Any) -> str:
    text = str("" if value is None else value).strip().lower()
    return "" if text in {"nan", "none", "null", "<na>"} else " ".join(text.split())


def find_reference_db() -> Path | None:
    for candidate in REFERENCE_DB_CANDIDATES:
        if candidate.exists():
            return candidate
    return None


def load_delivered(root: Path, libs: tuple[int, ...]) -> pd.DataFrame:
    """Unique predicted SMILES with the delivered assignment (first seen wins)."""

    frames = []
    for lib in libs:
        path = root / f"rerank_lib{lib}" / "unified_rerank_all_top10_with_pathways.csv"
        if not path.exists():
            print(f"[parity] WARN missing {path}")
            continue
        df = pd.read_csv(path, low_memory=False)
        keep = [c for c in ["predicted_smiles", *DELIVERED_COLUMNS] if c in df.columns]
        frames.append(df.loc[:, keep])
    if not frames:
        raise FileNotFoundError(f"no delivered *_with_pathways.csv under {root}")
    delivered = pd.concat(frames, ignore_index=True)
    delivered["predicted_smiles"] = delivered["predicted_smiles"].astype(str).str.strip()
    delivered = delivered[delivered["predicted_smiles"].ne("")]
    return delivered.drop_duplicates(subset=["predicted_smiles"], keep="first").reset_index(drop=True)


def run(args: argparse.Namespace) -> int:
    db_path = args.reference_db or find_reference_db()
    if db_path is None or not Path(db_path).exists():
        print(
            "[parity] SKIPPED: metabolite synonym SQLite not present "
            f"(looked for {[str(p) for p in REFERENCE_DB_CANDIDATES]}). "
            "Fetch it with ./fetch_assets.sh and re-run."
        )
        return 0
    root = args.equalization_root if args.equalization_root.is_absolute() else REPO_ROOT / args.equalization_root
    if not root.exists():
        print(
            f"[parity] SKIPPED: delivered rerank enrichment outputs not present under {root}. "
            "Fetch models/pipeline_equalization (fetch_get_assets.sh / rsync) and re-run."
        )
        return 0

    from structural_elucidation.pipeline.enrichment.pathway import (
        canonicalize_for_pathway,
        generate_pathway_assignment_map,
    )

    delivered = load_delivered(root, tuple(args.libs))
    delivered["canonical_smiles"] = delivered["predicted_smiles"].map(canonicalize_for_pathway)
    delivered = delivered[delivered["canonical_smiles"].ne("")]
    if args.limit:
        delivered = delivered.head(args.limit)

    vendored = generate_pathway_assignment_map(list(delivered["canonical_smiles"]), Path(db_path))
    merged = delivered.merge(vendored, on="canonical_smiles", how="left", suffixes=("_delivered", "_vendored"))

    def agree(col_delivered: str, col_vendored: str, frame: pd.DataFrame) -> float:
        if not len(frame):
            return 0.0
        matches = sum(
            norm(a) == norm(b) for a, b in zip(frame[col_delivered], frame[col_vendored])
        )
        return round(matches / len(frame), 4)

    summary: dict[str, Any] = {
        "n_unique_smiles": int(len(merged)),
        "super_agreement": agree("predicted_pathway_delivered", "predicted_pathway_vendored", merged),
        "sub_agreement": agree("predicted_subpathway_delivered", "predicted_subpathway_vendored", merged),
        "status_agreement": agree(
            "pathway_assignment_status_delivered", "pathway_assignment_status_vendored", merged
        ),
        "tier_agreement": agree("pathway_confidence_tier_delivered", "pathway_confidence_tier_vendored", merged),
        "by_delivered_status": {},
        "delivered_status_counts": {
            str(k): int(v)
            for k, v in Counter(merged["pathway_assignment_status_delivered"].map(norm)).most_common()
        },
    }
    for status, group in merged.groupby(merged["pathway_assignment_status_delivered"].map(norm)):
        summary["by_delivered_status"][str(status)] = {
            "n": int(len(group)),
            "super_agreement": agree("predicted_pathway_delivered", "predicted_pathway_vendored", group),
            "sub_agreement": agree("predicted_subpathway_delivered", "predicted_subpathway_vendored", group),
        }

    diffs = merged[
        [
            norm(a) != norm(b)
            for a, b in zip(merged["predicted_pathway_delivered"], merged["predicted_pathway_vendored"])
        ]
    ]
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2) + "\n")
    diff_path = args.out_json.with_name(args.out_json.stem + "_diffs.csv")
    diffs.to_csv(diff_path, index=False)

    lines = [
        "# Vendored neighbor vote vs delivered assignments",
        "",
        f"- Unique predicted SMILES compared: **{summary['n_unique_smiles']}**",
        f"- Super agreement: **{summary['super_agreement']:.1%}**; sub: {summary['sub_agreement']:.1%}",
        f"- Status agreement: {summary['status_agreement']:.1%}; tier: {summary['tier_agreement']:.1%}",
        "",
        "| delivered status | n | super agreement | sub agreement |",
        "|---|---|---|---|",
    ]
    for status, block in sorted(summary["by_delivered_status"].items()):
        lines.append(
            f"| {status} | {block['n']} | {block['super_agreement']:.1%} | {block['sub_agreement']:.1%} |"
        )
    lines += ["", f"Disagreements itemized in `{diff_path.name}`.", ""]
    args.out_md.write_text("\n".join(lines))
    print(f"[parity] wrote {args.out_json}, {args.out_md}, {diff_path}")
    print(json.dumps(summary, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--equalization-root", type=Path, default=Path("models/pipeline_equalization"))
    parser.add_argument("--reference-db", type=Path, default=None)
    parser.add_argument("--libs", nargs="+", type=int, default=list(DEFAULT_LIBS))
    parser.add_argument("--limit", type=int, default=0, help="cap unique SMILES for a quick pass")
    parser.add_argument(
        "--out-json",
        type=Path,
        default=REPO_ROOT / "docs" / "paper" / "data" / "pathway" / "vote_parity_report.json",
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=REPO_ROOT / "docs" / "paper" / "data" / "pathway" / "vote_parity_report.md",
    )
    return run(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
