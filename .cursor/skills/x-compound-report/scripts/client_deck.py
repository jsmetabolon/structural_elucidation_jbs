#!/usr/bin/env python3
"""Build a client X-compound deck folder (summary CSV + one PowerPoint) from a list of X ids.

This is the UPMC slide process (x_compound_info_summary_20261014/build_x_compound_slides.py,
2026-10-14): one deck, one slide per compound, "2 Columns" layout of the template deck
"20260919 external slide deck.pptx", with the summary sheet beside it. Unlike UPMC, the sheet
is generated from the latest results instead of typed by hand.

Stages:
1. screen   normalize ids, drop duplicates, block compounds ``data/known_X.csv`` marks
            ``Exclude`` / ``Artifact`` / ``DO NOT INCLUDE`` (other clients' proprietary or spiked).
2. sources  read the naming-v2 10 ppm report inventory (``output/reports_inventory_naming_v2``,
            whose X-##### folders link to ``reports_ppm10_full_naming_v2``). June bundle and
            ``reports 4`` are refused. A compound missing from it is listed, never back-filled.
3. sheet    one row per compound: m/z, formula, candidate pathway, rank-1 SMILES, report name,
            nearest named library metabolite + Tanimoto, MSI level (+ label, pathway basis,
            Metabolon level). Written as CSV (client) and xlsx (slide input).
4. slides   ``build_slides.py`` on the template.

Run from anywhere in the elucidation venv (python-pptx, openpyxl, Pillow; RDKit optional):
  client_deck.py --compounds ids.txt --client NIH --out NIH_deliverable \
      [--reports-root output/reports_inventory_naming_v2] [--template "...deck.pptx"]
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from datetime import date, datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[3]
sys.path.insert(0, str(SCRIPT_DIR))
csv.field_size_limit(sys.maxsize)

import client_report as CR  # noqa: E402  (id parsing and the proprietary-compound screen)

DEFAULT_REPORTS_ROOT = REPO_ROOT / "output" / "reports_inventory_naming_v2"
TEMPLATE_NAME = "20260919 external slide deck.pptx"
TEMPLATE_DIR = "x_compound_info_summary_20261014"
STALE_MARKERS = ("reports 4", "reports_4", "se_report_inputs", "handoff")

HEADER = [
    "Compound ID", "Observed m/z", "Putative formula", "Candidate superpathway",
    "Candidate subpathway", "Rank-1 predicted SMILES", "Rank-1 additional Information",
    "Most similar named metabolite", "Tanimoto similarity",
    "MSI level", "MSI label", "Pathway basis", "Metabolon level",
]
SHEET_HEADER = ["UPF group"] + HEADER  # the UPMC sheet's column layout; the group column stays blank


def read_rows(path: Path) -> list[dict]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def find_template(explicit: Path | None) -> Path:
    candidates = [explicit] if explicit else [
        REPO_ROOT.parent / TEMPLATE_DIR / TEMPLATE_NAME,
        REPO_ROOT / TEMPLATE_DIR / TEMPLATE_NAME,
    ]
    for path in candidates:
        if path and path.is_file():
            return path
    raise SystemExit(f"template deck not found; pass --template (looked for {TEMPLATE_NAME} in "
                     f"{[str(c) for c in candidates if c]}). It is not in git (22 MB).")


def check_reports_root(root: Path) -> Path:
    root = root.expanduser().resolve()
    if any(marker in str(root).lower() for marker in STALE_MARKERS):
        raise SystemExit(f"{root} is a stale or hand-off tree; use the naming-v2 10 ppm inventory")
    summary = root / "exec_summary" / "exec_summary_by_compound.csv"
    if not summary.is_file():
        raise SystemExit(f"{summary} not found; pass --reports-root output/reports_inventory_naming_v2")
    return root


def observed_mz(data_dir: Path, xid: str) -> float | None:
    path = data_dir / f"{xid.lower()}_raw_top_candidates.csv"
    if not path.is_file():
        return None
    match = re.search(r"_mz(\d+\.\d+)_", path.read_text(encoding="utf-8"))
    return float(match.group(1)) if match else None


def pathway(data_dir: Path, xid: str) -> tuple[str, str, str]:
    rows = read_rows(data_dir / f"{xid.lower()}_pathway_prediction.csv")
    if not rows:
        return "", "", ""
    row = rows[0]
    if row.get("predicted_pathway"):
        basis = "Curated X annotation" if row.get("pathway_source") == "known_x_annotation" else "Reference match"
        return row["predicted_pathway"], row.get("predicted_subpathway", ""), basis
    if row.get("pathway_vote_super"):
        share = float(row.get("pathway_vote_share") or 0)
        return row["pathway_vote_super"], row.get("pathway_vote_sub", ""), f"Top-candidate vote ({share:.0%})"
    return "", "", ""


def nearest(data_dir: Path) -> tuple[str, float | None]:
    rows = read_rows(data_dir / "similarity_neighbors.csv")
    if not rows:
        return "", None
    top = float(rows[0]["similarity_to_target"])
    tied = [r["canonical_name"] for r in rows[1:] if abs(float(r["similarity_to_target"]) - top) < 1e-9]
    name = rows[0]["canonical_name"]
    if len(tied) == 1:
        name += f" (tied with {tied[0]})"
    elif tied:
        name += f" (tied with {len(tied)} others)"
    return name, round(top, 4)


def build_row(root: Path, xid: str, summary_row: dict) -> list:
    data_dir = root / xid / "data"
    known = (read_rows(data_dir / f"{xid.lower()}_known_x_row.csv") or [{}])[0]
    has_structure = bool(read_rows(data_dir / f"{xid.lower()}_top_candidates.csv"))
    smiles = summary_row.get("rank1_candidate_smiles", "") if has_structure else ""
    formula = summary_row.get("rank1_candidate_formula") or known.get("Formula") or ""
    superpathway, subpathway, basis = pathway(data_dir, xid)
    name, similarity = nearest(data_dir)
    msi = summary_row.get("msi_level", "")
    msi = msi[:-2] if msi.endswith(".0") else msi
    level = known.get("Level", "")
    if known.get("SubLevel"):
        level += f"/{known['SubLevel']}"
    return [xid, observed_mz(data_dir, xid), formula, superpathway, subpathway, smiles,
            summary_row.get("report_display_name", ""), name, similarity,
            msi, summary_row.get("msi_label", ""), basis, level]


def write_sheet(rows: list[list], csv_path: Path, xlsx_path: Path) -> None:
    from openpyxl import Workbook
    from openpyxl.styles import Font

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        for row in rows:
            writer.writerow(["" if v is None else v for v in row])
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "summary"
    sheet.append(SHEET_HEADER)
    for row in rows:
        sheet.append([""] + row)
    font = Font(name="Aptos Narrow", size=12)
    for line in sheet.iter_rows():
        for cell in line:
            cell.font = font
    workbook.save(xlsx_path)


def mtime(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--compounds", nargs="+", required=True, help="ids, or .xlsx/.csv/.tsv/.txt files")
    ap.add_argument("--client", required=True, help="client name used in file names")
    ap.add_argument("--out", required=True, type=Path, help="client folder to create")
    ap.add_argument("--reports-root", type=Path, default=DEFAULT_REPORTS_ROOT)
    ap.add_argument("--template", type=Path, default=None)
    ap.add_argument("--known-x", type=Path, default=CR.KNOWN_X)
    ap.add_argument("--allow-excluded", action="store_true", help="internal QA only")
    ap.add_argument("--dry-run", action="store_true", help="screen and resolve sources; build nothing")
    args = ap.parse_args(argv)

    root = check_reports_root(args.reports_root)
    pairs, unparsed = CR.read_compound_inputs(args.compounds)
    screened = CR.screen(pairs, CR.load_known_x(args.known_x), allow_excluded=args.allow_excluded)
    summary_rows = {r["x_compound_id"]: r for r in read_rows(root / "exec_summary" / "exec_summary_by_compound.csv")}

    ids, excluded, missing = [], [], []
    for item in screened:
        cid = item["control_id"]
        if item["status"] == "excluded":
            excluded.append(cid)
        elif cid not in summary_rows or not (root / cid / "data").is_dir():
            missing.append(cid)
            item["status"], item["source_status"] = "missing", "not in naming-v2 inventory"
        else:
            ids.append(cid)
            item["source_status"] = "naming-v2"
    print(f"build {len(ids)}; excluded {excluded}; missing from naming-v2 {missing}; unparsed {unparsed}")
    if args.dry_run or not ids:
        return 0 if ids or args.dry_run else 1

    template = find_template(args.template)
    out = args.out
    work = out / "_build"
    work.mkdir(parents=True, exist_ok=True)
    safe = re.sub(r"[^A-Za-z0-9]+", "_", args.client).strip("_")
    csv_path = out / f"{safe}_X-compound_summary.csv"
    deck_path = out / f"{safe}_X-compound_slides_{date.today():%Y%m%d}.pptx"
    xlsx_path = work / f"{safe}_X-compound_summary.xlsx"

    write_sheet([build_row(root, cid, summary_rows[cid]) for cid in ids], csv_path, xlsx_path)

    import build_slides

    n_slides = build_slides.build(xlsx_path, template, deck_path, root)

    with (work / "compound_screen.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(screened[0].keys()))
        writer.writeheader()
        writer.writerows(screened)
    exec_csv = root / "exec_summary" / "exec_summary_by_compound.csv"
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "process": "UPMC slide process (x_compound_info_summary_20261014/build_x_compound_slides.py)",
        "client": args.client, "built": ids, "excluded": excluded, "missing": missing, "unparsed": unparsed,
        "slides": n_slides, "csv": csv_path.name, "deck": deck_path.name,
        "sources": {
            "reports_root": {"path": str(root), "mtime": mtime(root)},
            "exec_summary_by_compound": {"path": str(exec_csv), "mtime": mtime(exec_csv)},
            "template": {"path": str(template), "mtime": mtime(template)},
            "scores_json": "none: no scores/<id>.json or run_novel output is used",
        },
    }
    (work / "deck_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (out / "README.md").write_text(
        f"# {args.client} X-compound deliverable\n\n"
        f"| File | What it is |\n| --- | --- |\n"
        f"| `{csv_path.name}` | One row per compound: m/z, formula, candidate pathway, rank-1 SMILES, "
        f"nearest named library metabolite, MSI level. |\n"
        f"| `{deck_path.name}` | One slide per compound. |\n"
        f"| `_build/` | Provenance and build inputs. Internal; do not send. |\n\n"
        f"Structures are computational candidates, not confirmed identifications. "
        f"\"Not assigned\" means the current results do not support a value.\n",
        encoding="utf-8")
    print(f"wrote {out}: {csv_path.name}, {deck_path.name} ({n_slides} slides)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
