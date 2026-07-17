#!/usr/bin/env python3
"""Rebuild the master ``reports/index.html`` (and ``reports_index.csv``).

Sources lib/name/formula/etc. from the existing ``reports_index.csv`` (so no full
batch run is needed) and joins the MSI level + enriched display name from each
report's ``confidence_inputs.json``, then renders via the shared stdlib
``index_render`` module (adds an MSI column + numeric column sorting).

Pure standard library (imports ``index_render`` directly to avoid the pandas /
package ``__init__`` chain).

Usage:
    python3 tools/rebuild_master_index.py <reports_root>
"""
from __future__ import annotations

import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
import index_render as IR  # noqa: E402

CSV_COLUMNS = [
    "control_id", "lib", "name", "formula", "exact_mass", "db_chemical_id",
    "n_neighbors", "msi_level", "msi_label", "status", "report",
]


def main() -> None:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("_audit_pull/reports")
    csv_path = root / "reports_index.csv"
    if not csv_path.exists():
        raise SystemExit(f"missing {csv_path}")

    rows: list[IR.IndexRow] = []
    with csv_path.open(newline="") as h:
        for r in csv.DictReader(h):
            cid = str(r.get("control_id", "")).strip()
            if not cid:
                continue
            extra = IR._read_confidence_extra(root / cid)
            try:
                n_neighbors = int(float(r.get("n_neighbors") or 0))
            except (TypeError, ValueError):
                n_neighbors = 0
            rows.append(
                IR.IndexRow(
                    control_id=cid,
                    lib=str(r.get("lib", "") or ""),
                    name=str(extra.get("display_name") or r.get("name", "") or ""),
                    formula=str(r.get("formula", "") or ""),
                    exact_mass=str(r.get("exact_mass", "") or ""),
                    db_chemical_id=str(r.get("db_chemical_id", "") or ""),
                    n_neighbors=n_neighbors,
                    status=str(r.get("status", "existing") or "existing"),
                    report=str(r.get("report") or f"{cid}/index.html"),
                    msi_level=str(extra.get("msi_level") or ""),
                    msi_label=str(extra.get("msi_label") or ""),
                    msi_sublevel=str(extra.get("msi_sublevel") or ""),
                    metabolon_sublevel=str(extra.get("metabolon_sublevel") or ""),
                    match_conf=str(extra.get("match_conf") if extra.get("match_conf") not in (None, "") else ""),
                )
            )
    rows.sort(key=lambda x: (x.lib, x.control_id))
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    (root / "index.html").write_text(IR.render_master_index(rows, generated), encoding="utf-8")

    with csv_path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=CSV_COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({
                "control_id": r.control_id, "lib": r.lib, "name": r.name,
                "formula": r.formula, "exact_mass": r.exact_mass,
                "db_chemical_id": r.db_chemical_id, "n_neighbors": r.n_neighbors,
                "msi_level": r.msi_level, "msi_label": r.msi_label,
                "status": r.status, "report": r.report,
            })
    with_msi = sum(1 for r in rows if r.msi_level)
    print(f"[rebuild-master-index] wrote index.html + reports_index.csv for {len(rows)} rows ({with_msi} with MSI)")


if __name__ == "__main__":
    main()
