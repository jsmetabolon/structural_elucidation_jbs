#!/usr/bin/env python3
"""Audit a completed Naming v2 report set and build its correction exports."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


LEDGER_COLUMNS = [
    "control_id", "old_name", "new_name", "old_smiles", "new_smiles",
    "old_msi", "new_msi", "old_confidence", "new_confidence",
    "old_pathway", "new_pathway", "collision_group", "exception_basis",
    "name_source", "identity_candidate_rank", "identity_exclusion_reason",
]
GLOBALVIZ_COLUMNS = [
    "compound_id", "title", "formula", "library_id", "mode", "summary",
    "top_candidate_name", "top_candidate_smiles", "inchikey", "confidence_level",
    "confidence_score", "confidence_tier", "msi_level", "msi_label",
    "identity_status", "identity_basis", "nearest_reference_tanimoto", "pathway",
    "subpathway", "pathway_confidence", "pathway_source", "pathway_tier",
    "observed_parent_mz", "exact_mass", "descriptor_kind", "descriptor_confidence",
    "public_html", "confidence_inputs_json", "evidence_summary_json",
    "pathway_prediction_csv", "structure_image", "naming_policy_version",
    "name_source", "identity_candidate_rank", "library_collision", "identity_eligible",
    "identity_exclusion_reason", "excluded_library_matches", "raw_rank1_candidate",
    "identity_candidate", "exception_basis", "review_verdict",
]


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return {}


def read_csv_first(path: Path) -> dict:
    try:
        with path.open(newline="") as handle:
            return next(csv.DictReader(handle), {})
    except OSError:
        return {}


def record(root: Path, cid: str) -> dict:
    directory = root / cid
    slug = cid.lower()
    data = directory / "data"
    fields = read_json(data / f"{slug}_confidence_inputs.json")
    evidence = read_json(data / f"{slug}_evidence_summary.json")
    manifest = read_json(directory / "manifest.json")
    pathway = read_csv_first(data / f"{slug}_pathway_prediction.csv")
    target = manifest.get("target") or evidence.get("target") or {}
    return {"fields": fields, "evidence": evidence, "manifest": manifest,
            "pathway": pathway, "target": target}


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-root", type=Path, required=True)
    parser.add_argument("--new-root", type=Path, required=True)
    parser.add_argument("--expected-collisions", type=int, default=903)
    parser.add_argument("--expected-confirmed", type=int, default=15)
    args = parser.parse_args()
    old_ids = sorted(p.parent.name for p in args.old_root.glob("X-*/manifest.json"))
    new_ids = sorted(p.parent.name for p in args.new_root.glob("X-*/manifest.json"))
    if old_ids != new_ids:
        raise SystemExit(f"report inventory differs: old={len(old_ids)} new={len(new_ids)}")
    ledger: list[dict] = []
    export: list[dict] = []
    failures: list[str] = []
    old_collisions = 0
    classified = 0
    unconfirmed = 0
    confirmed = 0
    new_confirmed = 0
    for cid in old_ids:
        old = record(args.old_root, cid)
        new = record(args.new_root, cid)
        of, nf = old["fields"], new["fields"]
        ot, nt = old["target"], new["target"]
        old_name = str(of.get("display_name") or ot.get("name") or "")
        new_name = str(nf.get("display_name") or nt.get("name") or "")
        if nf.get("naming_policy_version") != "v2":
            failures.append(f"{cid}: missing v2 policy marker")
        if of.get("source") == "exact":
            old_collisions += 1
            if nf.get("exception_basis") == "standard_confirmed":
                new_confirmed += 1
            if nf.get("library_collision"):
                classified += 1
            else:
                failures.append(f"{cid}: old catalog collision not classified")
            if of.get("assurance") == "standard_confirmed":
                confirmed += 1
                if nf.get("exception_basis") != "standard_confirmed":
                    failures.append(f"{cid}: confirmed standard exception not retained")
            else:
                unconfirmed += 1
                if nf.get("exception_basis"):
                    failures.append(f"{cid}: unconfirmed collision has an exception")
                if nf.get("identity_status") == "known_pending":
                    failures.append(f"{cid}: unconfirmed collision promoted to known identity")
                if nt.get("smiles") and nt.get("smiles") == ot.get("smiles"):
                    failures.append(f"{cid}: colliding structure still displayed")
                selected = nf.get("identity_candidate")
                level = nf.get("msi_level")
                if selected is None and level == 3 and nf.get("name_source") != "known_x_tentative":
                    failures.append(f"{cid}: Level 3 without an eligible structure or curated hypothesis")
                if selected is None and level not in (3, 4, 5):
                    failures.append(f"{cid}: unresolved collision has invalid MSI level")
                if selected is None and nt.get("smiles"):
                    failures.append(f"{cid}: unresolved collision has a proposed structure")
                names = {str(match.get("name") or "").casefold() for item in
                         nf.get("excluded_library_matches") or [] for match in item.get("matches") or []}
                asserted = [new_name for new_name in (nf.get("display_name"), nf.get("resolved_name"), nt.get("name"))
                            if new_name]
                asserted.extend(str(value) for value in nf.get("synonyms") or [])
                if any(str(name).casefold().removeprefix("possible ").removeprefix("probable ") in names
                       for name in asserted):
                    failures.append(f"{cid}: excluded library name asserted")
                excluded_keys = {str(match.get("connectivity") or "") for item in
                                 nf.get("excluded_library_matches") or [] for match in item.get("matches") or []}
                if str(nt.get("inchikey") or "").split("-")[0] in excluded_keys:
                    failures.append(f"{cid}: excluded library connectivity asserted")
        ledger.append({
            "control_id": cid, "old_name": old_name, "new_name": new_name,
            "old_smiles": ot.get("smiles") or "", "new_smiles": nt.get("smiles") or "",
            "old_msi": of.get("msi_level") or "", "new_msi": nf.get("msi_level") or "",
            "old_confidence": of.get("overall") or of.get("evidence_score") or "",
            "new_confidence": nf.get("overall") or nf.get("evidence_score") or "",
            "old_pathway": old["pathway"].get("predicted_pathway") or "",
            "new_pathway": new["pathway"].get("predicted_pathway") or "",
            "collision_group": (old["evidence"].get("db_detail") or {}).get("chemical_id") or "",
            "exception_basis": nf.get("exception_basis") or "",
            "name_source": nf.get("name_source") or "",
            "identity_candidate_rank": nf.get("identity_candidate_rank") or "",
            "identity_exclusion_reason": nf.get("identity_exclusion_reason") or "",
        })
        slug = cid.lower()
        selected = nf.get("identity_candidate") or {}
        shown_structure = nf.get("msi_level") in (1, 2, 3)
        export.append({
            "compound_id": cid, "title": f"{cid} — {new_name}",
            "formula": nt.get("formula") or "", "library_id": nf.get("lib_id") or "",
            "mode": "", "summary": f"Naming v2 report for {cid}: {new_name}.",
            "top_candidate_name": new_name, "top_candidate_smiles": (nt.get("smiles") or "") if shown_structure else "",
            "inchikey": (nt.get("inchikey") or "") if shown_structure else "", "confidence_level": nf.get("tier") or "",
            "confidence_score": nf.get("overall") or "", "confidence_tier": nf.get("tier") or "",
            "msi_level": nf.get("msi_level") or "", "msi_label": nf.get("msi_label") or "",
            "identity_status": nf.get("identity_status") or "", "identity_basis": nf.get("identity_basis") or "",
            "nearest_reference_tanimoto": (new["evidence"].get("evidence") or {}).get("nearest_reference_tanimoto") or "",
            "pathway": new["pathway"].get("predicted_pathway") or "",
            "subpathway": new["pathway"].get("predicted_subpathway") or "",
            "pathway_confidence": new["pathway"].get("pathway_confidence_tier") or "",
            "pathway_source": new["pathway"].get("pathway_source") or "",
            "pathway_tier": new["pathway"].get("pathway_confidence_tier") or "",
            "observed_parent_mz": nf.get("observed_mz") or "", "exact_mass": nt.get("exact_mass") or "",
            "descriptor_kind": nf.get("descriptor_kind") or "", "descriptor_confidence": nf.get("descriptor_confidence") or "",
            "public_html": f"{cid}/public.html", "confidence_inputs_json": f"{cid}/data/{slug}_confidence_inputs.json",
            "evidence_summary_json": f"{cid}/data/{slug}_evidence_summary.json",
            "pathway_prediction_csv": f"{cid}/data/{slug}_pathway_prediction.csv",
            "structure_image": (
                f"{cid}/assets/{slug}_structure.svg"
                if shown_structure and nt.get("smiles") and (args.new_root / cid / "assets" / f"{slug}_structure.svg").exists()
                else f"{cid}/assets/{slug}_structure.png"
                if shown_structure and nt.get("smiles") and (args.new_root / cid / "assets" / f"{slug}_structure.png").exists()
                else ""
            ),
            "naming_policy_version": nf.get("naming_policy_version") or "",
            "name_source": nf.get("name_source") or "",
            "identity_candidate_rank": nf.get("identity_candidate_rank") or "",
            "library_collision": nf.get("library_collision"), "identity_eligible": nf.get("identity_eligible"),
            "identity_exclusion_reason": nf.get("identity_exclusion_reason") or "",
            "excluded_library_matches": json.dumps(nf.get("excluded_library_matches") or [], sort_keys=True),
            "raw_rank1_candidate": json.dumps(nf.get("raw_rank1_candidate") or {}, sort_keys=True),
            "identity_candidate": json.dumps(selected, sort_keys=True),
            "exception_basis": nf.get("exception_basis") or "",
            "review_verdict": nf.get("review_verdict") or "",
        })
    if old_collisions != args.expected_collisions or confirmed != args.expected_confirmed:
        failures.append(f"baseline drift: collisions={old_collisions} confirmed={confirmed}; expected {args.expected_collisions}/{args.expected_confirmed}")
    if new_confirmed != args.expected_confirmed:
        failures.append(f"v2 authentic-standard exceptions={new_confirmed}; expected {args.expected_confirmed}")
    if classified != old_collisions:
        failures.append(f"classified {classified}/{old_collisions} catalog collisions")
    out = args.new_root
    write_csv(out / "naming_v2_correction_ledger.csv", LEDGER_COLUMNS, ledger)
    write_csv(out / "globalviz_minimum_v2.csv", GLOBALVIZ_COLUMNS, export)
    audit = {"total_reports": len(old_ids), "old_catalog_collisions": old_collisions,
             "classified_collisions": classified, "standard_confirmed": confirmed,
             "v2_standard_exceptions": new_confirmed,
             "unconfirmed_collisions": unconfirmed, "failures": sorted(failures),
             "passed": not failures}
    (out / "naming_v2_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
