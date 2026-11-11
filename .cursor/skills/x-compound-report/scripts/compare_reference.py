#!/usr/bin/env python3
"""Compare a new client exec summary with a past delivery, compound by compound.

``--new`` is the ``exec_summary_by_compound.csv`` written by ``client_report.py``.
``--reference`` is a past ``exec_summary_by_compound.csv`` or an
``executive_summary_top10.xlsx`` (its ``by_compound`` sheet). Only compounds in
both are compared. Rank-1 SMILES are compared by RDKit canonical form so a
re-serialized but identical structure is not reported as a change.

Writes ``<out>/comparison.csv`` (one row per compound and field that differs)
and ``<out>/comparison.md`` (counts and the first differences), and prints the
counts as JSON.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

FIELDS = ("msi_level", "feature_name", "rank1_candidate_smiles", "rank1_candidate_formula",
          "n_candidates", "naming_policy_version", "identity_candidate_rank", "matched_at_rank_1")


def read_rows(path: Path) -> dict[str, dict[str, str]]:
    if path.suffix.lower() in (".xlsx", ".xlsm"):
        import openpyxl

        sheet = openpyxl.load_workbook(path, read_only=True, data_only=True)["by_compound"]
        rows = sheet.iter_rows(values_only=True)
        header = [str(h) for h in next(rows)]
        records = [dict(zip(header, ("" if v is None else str(v) for v in row))) for row in rows]
    else:
        with path.open(newline="", encoding="utf-8") as handle:
            records = list(csv.DictReader(handle))
    return {r["x_compound_id"]: r for r in records if r.get("x_compound_id")}


def canonical(smiles: str) -> str:
    if not smiles:
        return ""
    try:
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol) if mol is not None else smiles
    except ImportError:
        return smiles


def normalize(field: str, value: str) -> str:
    value = (value or "").strip()
    if field == "rank1_candidate_smiles":
        return canonical(value)
    if value.lower() in ("nan", "none"):
        return ""
    try:
        number = float(value)
        return str(int(number)) if number.is_integer() else value
    except ValueError:
        return value


def compare(new: dict[str, dict[str, str]], reference: dict[str, dict[str, str]]) -> dict:
    shared = sorted(set(new) & set(reference))
    fields = [f for f in FIELDS if any(f in r for r in new.values()) and any(f in r for r in reference.values())]
    diffs = []
    identical = 0
    for cid in shared:
        changed = False
        for field in fields:
            a = normalize(field, reference[cid].get(field, ""))
            b = normalize(field, new[cid].get(field, ""))
            if a != b:
                changed = True
                diffs.append({"x_compound_id": cid, "field": field, "reference": a, "new": b})
        identical += not changed
    per_field = {f: sum(1 for d in diffs if d["field"] == f) for f in fields}
    return {"n_new": len(new), "n_reference": len(reference), "n_shared": len(shared),
            "only_in_new": sorted(set(new) - set(reference)), "n_only_in_reference": len(set(reference) - set(new)),
            "fields_compared": fields, "identical_compounds": identical, "differences_by_field": per_field,
            "diffs": diffs}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--new", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    result = compare(read_rows(args.new), read_rows(args.reference))
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["x_compound_id", "field", "reference", "new"])
        writer.writeheader()
        writer.writerows(result["diffs"])
    lines = [f"# Comparison: {args.new} vs {args.reference}", "",
             f"- shared compounds: {result['n_shared']} (new {result['n_new']}, reference {result['n_reference']})",
             f"- identical on every compared field: {result['identical_compounds']}",
             f"- only in new: {', '.join(result['only_in_new']) or 'none'}", "",
             "| field | compounds differing |", "|---|---:|"]
    lines += [f"| {f} | {n} |" for f, n in result["differences_by_field"].items()]
    lines += ["", "| compound | field | reference | new |", "|---|---|---|---|"]
    lines += [f"| {d['x_compound_id']} | {d['field']} | {d['reference'][:80]} | {d['new'][:80]} |"
              for d in result["diffs"][:60]]
    (args.out / "comparison.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "diffs"}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
