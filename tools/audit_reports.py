#!/usr/bin/env python3
"""Audit generated X-compound HTML report bundles for content + visual quality.

Walks a reports root (``output/reports`` by default), and for every
``<control_id>/`` bundle inspects ``public.html`` plus ``manifest.json`` to flag:

- hardcoded template literals leaked from the X-25003 source template
  (``X-25003``/``N-lactoylalanine``/``lib_209_negative``/fabricated confidence
  cards) that are wrong on other compounds;
- a heading/title that does not match the bundle's own control id;
- broken intra-report links (``href``/``src`` to files missing in the bundle);
- the consensus atom-confidence SVG being too crowded (atom circles overlapping)
  - measured directly from the rendered SVG geometry;
- empty/placeholder sections and a display name that is just the control id.

Outputs ``audit_findings.csv``, ``audit_summary.json`` and ``audit_gallery.html``
under ``--out-dir`` (default ``output/_audit``).

Usage::

    python tools/audit_reports.py --reports-root output/reports
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any

# Unambiguous template-leak markers from the X-25003 source template. These never
# occur as legitimate data (a real neighbor named "X-25003" in the similarity
# network is fine and is intentionally NOT matched here).
STALE_PATTERNS = {
    "h1_x25003": re.compile(r"<h1>X-25003</h1>"),
    "title_x25003": re.compile(r"X-25003 Structure Resolution"),
    "x25003_link": re.compile(r"x25003_[a-z_]*\.csv"),
    "to_target_link": re.compile(r"_to_target\.csv"),
    "fabricated_level2c": re.compile(r"Level 2C"),
    "fabricated_10of10": re.compile(r"10/10"),
    "x25003_feature": re.compile(r"X-25003 feature"),
    "lactoyl_subtitle": re.compile(r"interpretation is <strong>N-lactoylalanine"),
    "lib_209_negative": re.compile(r"lib_209_negative structure-resolution"),
}

TITLE_RE = re.compile(r"<title>([^<]*)</title>")
H1_RE = re.compile(r"<h1>([^<]*)</h1>")
LINK_RE = re.compile(r'(?:href|src)="([^"]+)"')
# The *displayed* consensus image lives in the consensus-box; the candidate maps
# are embedded as escaped JSON in a <script> and must not be parsed as geometry.
CONSENSUS_BOX_RE = re.compile(r'<div class="consensus-box">(.*?)</div>', re.S)
TRANSLATE_RE = re.compile(r"translate\(([-\d.]+),([-\d.]+)\)")
CIRCLE_R_RE = re.compile(r'<circle r="([\d.]+)"')


def consensus_crowding(html_text: str) -> dict[str, Any] | None:
    """Crowding of the displayed consensus image.

    overlap_ratio = (2*radius) / min_center_distance; >1 means atom circles
    overlap (the structure is cramped and hard to read).
    """

    box = CONSENSUS_BOX_RE.search(html_text)
    if not box:
        return None
    block = box.group(1)
    centers = [(float(x), float(y)) for x, y in TRANSLATE_RE.findall(block)]
    radii = [float(r) for r in CIRCLE_R_RE.findall(block)]
    if len(centers) < 2 or not radii:
        return None
    radius = radii[0]
    min_d = min(
        ((centers[i][0] - centers[j][0]) ** 2 + (centers[i][1] - centers[j][1]) ** 2) ** 0.5
        for i in range(len(centers))
        for j in range(i + 1, len(centers))
    )
    ratio = round((2 * radius) / min_d, 2) if min_d > 0 else 999.0
    return {"n_atoms": len(centers), "radius": radius, "min_center_dist": round(min_d, 1), "overlap_ratio": ratio}


def broken_links(html_text: str, bundle_dir: Path) -> list[str]:
    """Relative href/src targets that do not exist inside the bundle."""

    missing: list[str] = []
    seen: set[str] = set()
    for link in LINK_RE.findall(html_text):
        if link in seen:
            continue
        seen.add(link)
        if link.startswith(("http://", "https://", "#", "mailto:", "data:", "//")):
            continue
        target = (bundle_dir / link.split("#")[0].split("?")[0]).resolve()
        if not target.exists():
            missing.append(link)
    return missing


def load_lib_map(reports_root: Path) -> dict[str, str]:
    """control_id -> lib from reports_index.csv if available."""

    index = reports_root / "reports_index.csv"
    if not index.exists():
        return {}
    out: dict[str, str] = {}
    with index.open(newline="") as handle:
        for row in csv.DictReader(handle):
            out[str(row.get("control_id", ""))] = str(row.get("lib", ""))
    return out


def audit_one(bundle: Path, lib_map: dict[str, str]) -> dict[str, Any]:
    """Audit a single ``<control_id>/`` bundle."""

    control_id = bundle.name
    page = bundle / "public.html"
    manifest_path = bundle / "manifest.json"
    rec: dict[str, Any] = {"control_id": control_id, "lib": lib_map.get(control_id, "")}

    name = ""
    if manifest_path.exists():
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            target = payload.get("target") or {}
            name = str(target.get("name") or "")
            rec["db_chemical_id"] = str(payload.get("db_chemical_id") or "")
            rec["n_neighbors"] = int(payload.get("n_neighbors") or 0)
        except Exception:
            pass
    rec["name"] = name
    rec["name_is_control_id"] = name.strip() == control_id

    if not page.exists():
        rec["page_missing"] = True
        return rec
    text = page.read_text(encoding="utf-8", errors="replace")

    for label, pattern in STALE_PATTERNS.items():
        rec[f"n[{label}]"] = len(pattern.findall(text))

    title = TITLE_RE.search(text)
    h1 = H1_RE.search(text)
    rec["title"] = title.group(1) if title else ""
    rec["h1"] = h1.group(1) if h1 else ""
    rec["h1_ok"] = control_id in rec["h1"]
    rec["title_ok"] = control_id in rec["title"]

    missing = broken_links(text, bundle)
    rec["n_broken_links"] = len(missing)
    rec["broken_links_sample"] = "; ".join(sorted(set(missing))[:5])

    crowd = consensus_crowding(text)
    rec["consensus_n_atoms"] = crowd["n_atoms"] if crowd else 0
    rec["consensus_overlap_ratio"] = crowd["overlap_ratio"] if crowd else 0.0

    rec["synonyms_empty"] = "No synonyms returned" in text
    return rec


def summarize(records: list[dict[str, Any]], crowd_threshold: float) -> dict[str, Any]:
    """Aggregate issue counts + distributions across all bundles."""

    n = len(records)

    def count(pred) -> int:
        return sum(1 for r in records if pred(r))

    ratios = sorted(r.get("consensus_overlap_ratio", 0.0) for r in records)

    def pct(p: float) -> float:
        if not ratios:
            return 0.0
        return round(ratios[min(len(ratios) - 1, int(p * len(ratios)))], 2)

    worst_crowd = sorted(records, key=lambda r: r.get("consensus_overlap_ratio", 0.0), reverse=True)[:25]
    worst_links = sorted(records, key=lambda r: r.get("n_broken_links", 0), reverse=True)[:25]
    return {
        "n_reports": n,
        "reports_with_template_leak": count(
            lambda r: any(r.get(f"n[{lbl}]", 0) > 0 for lbl in STALE_PATTERNS)
        ),
        "reports_with_wrong_h1": count(lambda r: not r.get("h1_ok", False)),
        "reports_with_wrong_title": count(lambda r: not r.get("title_ok", False)),
        "reports_with_broken_links": count(lambda r: r.get("n_broken_links", 0) > 0),
        "reports_with_fabricated_confidence": count(lambda r: r.get("n[fabricated_level2c]", 0) > 0),
        "reports_with_hardcoded_name_lactoyl": count(lambda r: r.get("n[lactoyl_subtitle]", 0) > 0),
        "reports_name_is_control_id": count(lambda r: r.get("name_is_control_id", False)),
        "reports_synonyms_empty": count(lambda r: r.get("synonyms_empty", False)),
        "reports_crowded_consensus_over_threshold": count(
            lambda r: r.get("consensus_overlap_ratio", 0.0) > crowd_threshold
        ),
        "crowd_threshold": crowd_threshold,
        "consensus_overlap_ratio_pctiles": {"p50": pct(0.5), "p90": pct(0.9), "p99": pct(0.99), "max": pct(1.0)},
        "max_consensus_atoms": max((r.get("consensus_n_atoms", 0) for r in records), default=0),
        "worst_crowded": [
            {"control_id": r["control_id"], "n_atoms": r.get("consensus_n_atoms", 0),
             "overlap_ratio": r.get("consensus_overlap_ratio", 0.0)}
            for r in worst_crowd
        ],
        "worst_broken_links": [
            {"control_id": r["control_id"], "n_broken_links": r.get("n_broken_links", 0),
             "sample": r.get("broken_links_sample", "")}
            for r in worst_links
        ],
    }


def write_gallery(out_dir: Path, reports_root: Path, summary: dict[str, Any]) -> Path:
    """Small HTML index linking the worst-crowded reports for visual review."""

    rel = Path("..") / reports_root.name
    rows = "\n".join(
        f"<tr><td>{i+1}</td><td><a href='{rel}/{r['control_id']}/public.html'>{r['control_id']}</a></td>"
        f"<td>{r['n_atoms']}</td><td>{r['overlap_ratio']}</td></tr>"
        for i, r in enumerate(summary["worst_crowded"])
    )
    html = (
        "<!doctype html><meta charset='utf-8'><title>Report audit - worst crowded</title>"
        "<style>body{font-family:system-ui,sans-serif;padding:24px}table{border-collapse:collapse}"
        "td,th{border-bottom:1px solid #ccc;padding:6px 12px;text-align:left}</style>"
        f"<h1>Worst-crowded consensus images ({summary['n_reports']} reports audited)</h1>"
        f"<p>overlap_ratio &gt; 1 means atom circles overlap. p90={summary['consensus_overlap_ratio_pctiles']['p90']}, "
        f"max={summary['consensus_overlap_ratio_pctiles']['max']}, max atoms={summary['max_consensus_atoms']}.</p>"
        "<table><tr><th>#</th><th>control_id</th><th>atoms</th><th>overlap_ratio</th></tr>"
        f"{rows}</table>"
    )
    path = out_dir / "audit_gallery.html"
    path.write_text(html, encoding="utf-8")
    return path


def main() -> None:
    """CLI entrypoint."""

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reports-root", default="output/reports")
    parser.add_argument("--out-dir", default="output/_audit")
    parser.add_argument("--crowd-threshold", type=float, default=1.0)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    reports_root = Path(args.reports_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    lib_map = load_lib_map(reports_root)

    bundles = sorted(p for p in reports_root.glob("X-*") if p.is_dir())
    if args.limit:
        bundles = bundles[: args.limit]
    print(f"[audit] auditing {len(bundles)} bundles under {reports_root}")

    records: list[dict[str, Any]] = []
    for i, bundle in enumerate(bundles, 1):
        records.append(audit_one(bundle, lib_map))
        if i % 500 == 0:
            print(f"[audit]   {i}/{len(bundles)}")

    fieldnames: list[str] = []
    for rec in records:
        for key in rec:
            if key not in fieldnames:
                fieldnames.append(key)
    csv_path = out_dir / "audit_findings.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    summary = summarize(records, args.crowd_threshold)
    (out_dir / "audit_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    gallery = write_gallery(out_dir, reports_root, summary)

    print(f"[audit] findings: {csv_path}")
    print(f"[audit] summary:  {out_dir / 'audit_summary.json'}")
    print(f"[audit] gallery:  {gallery}")
    print(json.dumps({k: v for k, v in summary.items() if not isinstance(v, (list, dict))}, indent=2))


if __name__ == "__main__":
    main()
