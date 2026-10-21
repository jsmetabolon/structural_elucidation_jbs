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
    "n_neighbors", "msi_level", "msi_label", "descriptor_kind",
    "descriptor_confidence", "descriptor_parent_formula",
    "descriptor_conjugate_count", "status", "report",
    "naming_policy_version", "name_source", "identity_candidate_rank",
    "library_collision", "identity_eligible", "identity_exclusion_reason",
    "exception_basis", "review_verdict", "excluded_library_matches",
    "raw_rank1_candidate", "identity_candidate",
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
                    identity_status=str(extra.get("identity_status") or ""),
                    identity_confidence=str(extra.get("identity_confidence") or ""),
                    descriptor_kind=str(extra.get("descriptor_kind") or ""),
                    descriptor_confidence=str(extra.get("descriptor_confidence") or ""),
                    descriptor_parent_formula=str(extra.get("descriptor_parent_formula") or ""),
                    descriptor_conjugate_count=str(extra.get("descriptor_conjugate_count") or ""),
                    naming_policy_version=str(extra.get("naming_policy_version") or ""),
                    name_source=str(extra.get("name_source") or ""),
                    identity_candidate_rank=str(extra.get("identity_candidate_rank") or ""),
                    library_collision=str(extra["library_collision"]) if extra.get("library_collision") is not None else "",
                    identity_eligible=str(extra["identity_eligible"]) if extra.get("identity_eligible") is not None else "",
                    identity_exclusion_reason=str(extra.get("identity_exclusion_reason") or ""),
                    exception_basis=str(extra.get("exception_basis") or ""),
                    review_verdict=str(extra.get("review_verdict") or ""),
                    excluded_library_matches=str(extra.get("excluded_library_matches") or ""),
                    raw_rank1_candidate=str(extra.get("raw_rank1_candidate") or ""),
                    identity_candidate=str(extra.get("identity_candidate") or ""),
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
                "descriptor_kind": r.descriptor_kind,
                "descriptor_confidence": r.descriptor_confidence,
                "descriptor_parent_formula": r.descriptor_parent_formula,
                "descriptor_conjugate_count": r.descriptor_conjugate_count,
                "naming_policy_version": r.naming_policy_version,
                "name_source": r.name_source,
                "identity_candidate_rank": r.identity_candidate_rank,
                "library_collision": r.library_collision,
                "identity_eligible": r.identity_eligible,
                "identity_exclusion_reason": r.identity_exclusion_reason,
                "exception_basis": r.exception_basis,
                "review_verdict": r.review_verdict,
                "excluded_library_matches": r.excluded_library_matches,
                "raw_rank1_candidate": r.raw_rank1_candidate,
                "identity_candidate": r.identity_candidate,
                "status": r.status, "report": r.report,
            })
    with_msi = sum(1 for r in rows if r.msi_level)
    print(f"[rebuild-master-index] wrote index.html + reports_index.csv for {len(rows)} rows ({with_msi} with MSI)")


if __name__ == "__main__":
    main()
