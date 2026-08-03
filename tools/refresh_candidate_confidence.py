#!/usr/bin/env python3
"""Refresh per-candidate ``candidate_confidence`` in already-enriched reports.

Recomputes each report's ``confidence_inputs.json`` ``candidate_evidence`` array
to add the rank-consistent ``candidate_confidence`` (+ tie flag), reusing the
report's existing confidence ``fields`` for the rank-1 overall-confidence anchor
and its ``<slug>_candidate_score_summary.csv`` for the candidate rows.

No source bundle (feature rows / observed MS2 / reference DB / WGCNA) is required:
the delivered reports tree alone is sufficient to run this and then
``tools/patch_candidate_scores.py``. Use the full ``tools/enrich_report_evidence.py``
instead only when you also want to re-derive the other evidence fields.

Idempotent (overwrites only ``candidate_evidence``). Pure standard library.

Usage:
    python3 tools/refresh_candidate_confidence.py <reports_root>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import confidence as C  # noqa: E402
import enrich_report_evidence as E  # noqa: E402


def refresh(report_dir: Path) -> str:
    data = report_dir / "data"
    if not data.is_dir():
        return "no_data"
    slug = report_dir.name.lower()
    ci = data / f"{slug}_confidence_inputs.json"
    css_path = data / f"{slug}_candidate_score_summary.csv"
    if not ci.exists() or not css_path.exists():
        return "no_sidecars"
    try:
        fields = json.loads(ci.read_text())
    except (ValueError, OSError):
        return "bad_json"
    css = E._read(css_path)
    if not css:
        return "no_candidates"
    # top_candidates.csv carries evidence_confidence (the CSV summary does not), so
    # use it for the SMILES join to preserve the "Formula/class evidence" context;
    # fall back to css when the file is absent.
    tc = E._read(data / f"{slug}_top_candidates.csv")

    # Sublevel evidence: curated Metabolon SubLevel (authoritative) + same-formula
    # multiplicity (the positional-isomer / 3A signal), both reports-tree-only.
    kx = E._read(data / f"{slug}_known_x_row.csv")
    if kx:
        fields["metabolon_sublevel"] = str(kx[0].get("SubLevel", "") or "").strip()
    fields["n_isomers"] = E.same_formula_multiplicity(tc or css)

    conf = C.compute_confidence(fields)
    fields["msi_level"] = conf.msi_level
    fields["msi_label"] = conf.msi_label
    fields["msi_rationale"] = conf.msi_rationale
    fields["msi_sublevel"] = conf.msi_sublevel
    fields["metabolon_sublevel"] = conf.metabolon_sublevel
    fields["candidate_evidence"] = E.candidate_evidence_rows(css, tc or css, conf.raw_evidence)
    ci.write_text(json.dumps(fields, indent=2))
    return "refreshed"


def main() -> None:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("_audit_pull/reports")
    if not root.is_dir():
        raise SystemExit(f"reports root not found: {root}")
    counts: dict[str, int] = {}
    for d in sorted(p for p in root.iterdir() if p.is_dir() and p.name.upper().startswith("X-")):
        r = refresh(d)
        counts[r] = counts.get(r, 0) + 1
    print(f"[refresh-candidate-confidence] {counts}")


if __name__ == "__main__":
    main()
