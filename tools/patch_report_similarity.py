#!/usr/bin/env python3
"""Patch already-generated compound reports in place.

Replaces the misleading ``score_rank``-based "Tanimoto" header badge and the
"Rank-1 candidate Tanimoto" metric card with the real *nearest-reference
structural similarity*, and refreshes the ``*_evidence_summary`` sidecar to
carry honest fields. This yields output byte-identical to regenerating with the
fixed ``structural_elucidation/reports/compound.py`` -- but without recomputing
the expensive consensus/MCS maps, structure images, or similarity networks
(none of which the fix changes).

All values are taken from each report's own
``data/<slug>_pathway_prediction.csv`` (``pathway_source`` +
``pathway_best_reference_similarity``) and the existing evidence sidecar, so no
model inputs are needed. Idempotent and safe to re-run.

Usage:
    python -m tools.patch_report_similarity output/reports            # apply
    python -m tools.patch_report_similarity output/reports --dry-run  # report only
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
from pathlib import Path

# The old header badge and metric card, as emitted by the pre-fix renderer.
# ``[^<]*`` safely matches the formatted value ("nan", "0.000", "0.125", ...)
# without ever crossing a tag boundary.
BADGE_RE = re.compile(r'<span class="badge ok">Tanimoto [^<]*</span>')
CARD_RE = re.compile(
    r'<div class="metric-card"><div class="metric-label">Rank-1 candidate Tanimoto</div>'
    r'<div class="metric-value">[^<]*</div>'
    r'<div class="metric-note">[^<]*</div></div>'
)


def compute_display(source: str, sim: float) -> tuple[str, str, str]:
    """Mirror render_html: return (badge, note, value_str) for the similarity."""
    if math.isnan(sim):
        sim = 0.0
    if source == "exact":
        badge = "Exact reference match"
        note = "Top candidate structure exactly matches a reference DB metabolite"
    elif sim > 0:
        badge = f"Nearest-ref Tanimoto {sim:.2f}"
        note = "Max Tanimoto of top candidate to nearest reference metabolite"
    else:
        badge = "No reference match"
        note = "No reference structure above the similarity floor"
    return badge, note, f"{sim:.2f}"


def read_pathway(data_dir: Path) -> tuple[str, float] | None:
    """Read (pathway_source, pathway_best_reference_similarity) from the sidecar."""
    cands = sorted(data_dir.glob("*_pathway_prediction.csv"))
    if not cands:
        return None
    with cands[0].open(newline="") as handle:
        try:
            row = next(csv.DictReader(handle))
        except StopIteration:
            return None
    source = (row.get("pathway_source") or "").strip()
    try:
        sim = float(row.get("pathway_best_reference_similarity"))
    except (TypeError, ValueError):
        sim = 0.0
    if math.isnan(sim):
        sim = 0.0
    return source, sim


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    """Byte-compatible with structural_elucidation.reports.compound.write_csv."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def patch_report(report_dir: Path, apply: bool) -> str:
    """Patch one compound report dir. Returns a status string."""
    data_dir = report_dir / "data"
    if not data_dir.is_dir():
        return "no_data"
    pathway = read_pathway(data_dir)
    if pathway is None:
        return "no_pathway"
    source, sim = pathway
    badge, note, value_str = compute_display(source, sim)
    new_badge = f'<span class="badge ok">{html.escape(badge)}</span>'
    new_card = (
        '<div class="metric-card"><div class="metric-label">Top candidate structural similarity</div>'
        f'<div class="metric-value">{value_str}</div>'
        f'<div class="metric-note">{html.escape(note)}</div></div>'
    )

    html_changed = 0
    for html_file in sorted(report_dir.glob("*.html")):
        text = html_file.read_text()
        if 'badge ok">Tanimoto ' not in text and "Rank-1 candidate Tanimoto" not in text:
            continue
        new_text = BADGE_RE.sub(new_badge, text)
        new_text = CARD_RE.sub(new_card, new_text)
        if new_text != text:
            if apply:
                html_file.write_text(new_text)
            html_changed += 1

    # Refresh evidence sidecars (JSON + CSV) with honest fields.
    for json_file in sorted(data_dir.glob("*_evidence_summary.json")):
        obj = json.loads(json_file.read_text())
        old_ev = obj.get("evidence", {}) if isinstance(obj, dict) else {}
        new_ev = {
            "nearest_reference_tanimoto": sim,
            "nearest_reference_source": source,
            "rerank_pool_score": old_ev.get(
                "full_pool_oracle_best_tanimoto", old_ev.get("rerank_pool_score", "")
            ),
            "parent_mz": old_ev.get("parent_mz", ""),
            "mass_error_da": old_ev.get("mass_error_da", ""),
            "lactoyl_hypothesis_score": old_ev.get("lactoyl_hypothesis_score", ""),
            "diagnostic_bin_score": old_ev.get("diagnostic_bin_score", ""),
            "matched_diagnostic_bins_json": old_ev.get("matched_diagnostic_bins_json", ""),
        }
        if isinstance(obj, dict):
            obj["evidence"] = new_ev
            if apply:
                json_file.write_text(json.dumps(obj, indent=2))
        csv_file = json_file.with_suffix(".csv")
        if csv_file.exists() and apply:
            write_csv(csv_file, [{"field": k, "value": v} for k, v in new_ev.items()])

    if html_changed:
        return "patched"
    return "already_current"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports_root", type=Path, help="Master reports dir (one subdir per compound)")
    parser.add_argument("--dry-run", action="store_true", help="Report what would change without writing")
    args = parser.parse_args()

    root = args.reports_root
    if not root.is_dir():
        raise SystemExit(f"not a directory: {root}")

    counts: dict[str, int] = {}
    for report_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        status = patch_report(report_dir, apply=not args.dry_run)
        counts[status] = counts.get(status, 0) + 1

    verb = "would patch" if args.dry_run else "patched"
    print(f"[patch] root={root}")
    for status in sorted(counts):
        print(f"[patch]   {status}: {counts[status]}")
    print(f"[patch] {verb} {counts.get('patched', 0)} report(s)")


if __name__ == "__main__":
    main()
