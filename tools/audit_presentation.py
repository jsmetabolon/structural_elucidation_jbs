#!/usr/bin/env python3
"""Audit presentation + scientific-evidence completeness across all reports.

For every report this records evidence availability per axis, the hybrid
confidence (score/tier/completeness), and presentation defects - including a
compliance check that public variants never expose observed MS2 ions. Emits a
per-report CSV, an aggregate JSON (coverage %), and a worst-case gap gallery.

Usage:
    python3 tools/audit_presentation.py <reports_root> [--out-dir DIR]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
import confidence as C  # noqa: E402

# Strings that would indicate observed-MS2 ion detail leaking into a public page.
PUBLIC_MS2_MARKERS = ("MS2 fragment fit", "Observed MS2 spectrum", "fragment peaks", "matched MS2")


def _num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f != f else f


def audit_report(d: Path) -> dict | None:
    data = d / "data"
    if not data.is_dir():
        return None
    slug = d.name.lower()
    fields = C.read_fields_from_sidecars(data, slug)
    if not fields:
        return None
    conf = C.compute_confidence(fields)
    by = {a.key: a for a in conf.axes}

    name_resolved = bool(fields.get("name_resolved"))
    synonyms_empty = not bool(fields.get("synonyms"))
    has_correlations = bool(fields.get("correlations"))
    has_level = bool(str(fields.get("metabolon_level", "")).strip())

    # presentation defects from the public page
    public_ms2_leak = 0
    mass_nan_public = 0
    pub = d / "index.html"
    if pub.exists():
        text = pub.read_text()
        public_ms2_leak = int(any(m in text for m in PUBLIC_MS2_MARKERS))
        mass_nan_public = int("Observed parent m/z</div><div class=\"metric-value\">nan<" in text)

    mass = by["mass"].score if by["mass"].available else None
    return {
        "control_id": d.name,
        "source": conf.source,
        "tier": conf.tier,
        "overall": conf.overall,
        "completeness": conf.completeness,
        "structure": by["structure"].score,
        "mass_score": mass,
        "mass_mismatch": int(mass is not None and mass < C.MASS_WEAK),
        "ms2_present": int(by["ms2"].available),
        "ms2_score": by["ms2"].score if by["ms2"].available else "",
        "pathway_available": int(by["pathway"].available),
        "observed_mz_available": int(bool(fields.get("observed_mz_available"))),
        "library_neighbors_ge_050": fields.get("library_neighbors_ge_050", 0),
        "name_resolved": int(name_resolved),
        "synonyms_empty": int(synonyms_empty),
        "has_correlations": int(has_correlations),
        "has_metabolon_level": int(has_level),
        "public_ms2_leak": public_ms2_leak,
        "mass_nan_public": mass_nan_public,
        "n_flags": len(conf.flags),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("reports_root", type=Path)
    ap.add_argument("--out-dir", type=Path, default=Path("output/_presentation_audit"))
    args = ap.parse_args()
    rows = []
    for d in sorted(p for p in args.reports_root.iterdir() if p.is_dir()):
        r = audit_report(d)
        if r:
            rows.append(r)
    n = len(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "audit_presentation_findings.csv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    def frac(key):
        return round(sum(r[key] for r in rows) / n, 4) if n else 0.0

    tiers = Counter(r["tier"] for r in rows)
    overalls = sorted(r["overall"] for r in rows)
    summary = {
        "n_reports": n,
        "tier_distribution": dict(tiers),
        "overall_confidence": {
            "median": overalls[n // 2] if n else 0.0,
            "p10": overalls[int(0.1 * n)] if n else 0.0,
            "p90": overalls[int(0.9 * n)] if n else 0.0,
        },
        "evidence_coverage": {
            "observed_mz_available": frac("observed_mz_available"),
            "ms2_spectrum_present": frac("ms2_present"),
            "mass_mismatch": frac("mass_mismatch"),
            "pathway_available": frac("pathway_available"),
        },
        "presentation_gaps": {
            "name_resolved_from_db": frac("name_resolved"),
            "synonyms_empty": frac("synonyms_empty"),
            "has_wgcna_correlations": frac("has_correlations"),
            "has_metabolon_level": frac("has_metabolon_level"),
            "mass_nan_public": sum(r["mass_nan_public"] for r in rows),
        },
        "compliance": {"public_ms2_leaks": sum(r["public_ms2_leak"] for r in rows)},
        "confidence_completeness_median": sorted(r["completeness"] for r in rows)[n // 2] if n else 0.0,
    }
    (args.out_dir / "audit_presentation_summary.json").write_text(json.dumps(summary, indent=2))

    # gap gallery: worst cases per category
    def worst(key, reverse, limit=15):
        return sorted(rows, key=lambda r: (r[key] if isinstance(r[key], (int, float)) else 0), reverse=reverse)[:limit]

    gallery = ["<html><body><h1>Report presentation gap gallery</h1>"]
    for title, items in [
        ("Mass mismatch (top candidate mass disagrees)", [r for r in rows if r["mass_mismatch"]][:20]),
        ("Lowest confidence", worst("overall", False)),
        ("Missing observed m/z", [r for r in rows if not r["observed_mz_available"]][:20]),
        ("No MS2 spectrum", [r for r in rows if not r["ms2_present"]][:20]),
    ]:
        gallery.append(f"<h2>{title} ({len(items)} shown)</h2><ul>")
        for r in items:
            gallery.append(
                f"<li>{r['control_id']} - {r['tier']} ({r['overall']:.2f}); "
                f"struct={r['structure']}, mass={r['mass_score']}, ms2={r['ms2_score']}</li>"
            )
        gallery.append("</ul>")
    gallery.append("</body></html>")
    (args.out_dir / "gap_gallery.html").write_text("\n".join(gallery))

    print(f"[audit-presentation] n={n}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
