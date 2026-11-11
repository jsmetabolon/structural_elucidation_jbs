#!/usr/bin/env python3
"""Build a client X-compound report delivery from a list of X compounds.

Stages (each writes under ``--out``):

1. ``screen``: normalize the supplied ids, drop duplicates, block compounds that
   ``data/known_X.csv`` marks ``Exclude``/``Artifact``/``DO NOT INCLUDE`` (other
   clients' proprietary and spiked compounds), flag ids not in ``known_X``.
2. ``sources``: pick the latest predictions / scoring JSONs / rendered reports
   with ``resolve_latest.py`` unless they are passed explicitly.
3. ``reports``: either copy bundles from an existing rendered delivery
   (``--reports-root``, exact reproduction of what was delivered) or regenerate
   them with ``structural_elucidation.reports.batch`` from rerank CSVs
   (``--predictions-root``, optionally ``--scores-dir``).
4. ``summary``: ``tools/build_exec_summary.py`` + ``tools/build_exec_summary_workbook.py``.
5. ``audit``: ``tools/audit_reports.py`` plus a public-package disclosure check.
6. ``package``: a public-only client folder (master index, each compound's
   ``public.html`` as ``index.html`` with the files it links, the workbook).

Run from the repository root in the elucidation venv (Python 3.10, pandas<2, RDKit).
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[3]
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(REPO_ROOT))

import resolve_latest  # noqa: E402

KNOWN_X = REPO_ROOT / "data" / "known_X.csv"
ID_TOKEN = re.compile(r"\bX\s*[-_ ]?\s*(\d{3,9})\b", re.IGNORECASE)
ID_COLUMNS = ("control_id", "x_compound_id", "unknown", "compound", "compound_id", "x_id", "x compound", "id")
BLOCKED_LEVELS = {"exclude", "artifact"}
SECONDS_PER_COMPOUND = 7  # measured on the June bundle with the synonym DB, one job
PUBLIC_FORBIDDEN = ("internal.html", "_scoring.json", "_build.log", "candidate_score_summary", "observed_ms2")


def canonical_id(number: str) -> str:
    return f"X-{int(number):05d}"


def ids_from_text(text: str) -> list[tuple[str, str]]:
    return [(match.group(0), canonical_id(match.group(1))) for match in ID_TOKEN.finditer(text)]


def ids_from_table(rows: list[list[str]]) -> list[tuple[str, str]]:
    """Use a named id column when present, else the column where most cells are X ids."""
    if not rows:
        return []
    header = [str(c or "").strip().lower() for c in rows[0]]
    column = next((i for i, name in enumerate(header) if name in ID_COLUMNS), None)
    body = rows[1:] if column is not None else rows
    if column is None:
        width = max(len(r) for r in rows)
        hits = [sum(1 for r in rows if i < len(r) and ID_TOKEN.fullmatch(str(r[i] or "").strip()))
                for i in range(width)]
        column = max(range(width), key=lambda i: hits[i])
        if not hits[column]:
            return []
    out = []
    for row in body:
        if column < len(row):
            out.extend(ids_from_text(str(row[column] or "")))
    return out


def read_compound_inputs(values: Iterable[str]) -> tuple[list[tuple[str, str]], list[str]]:
    """Return (raw token, canonical id) pairs and tokens that held no id."""
    pairs: list[tuple[str, str]] = []
    unparsed: list[str] = []
    for value in values:
        path = Path(value).expanduser()
        if path.is_file():
            suffix = path.suffix.lower()
            if suffix in (".xlsx", ".xlsm", ".xls"):
                import openpyxl

                workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
                for sheet in workbook.worksheets:
                    pairs.extend(ids_from_table([list(r) for r in sheet.iter_rows(values_only=True)]))
            elif suffix in (".csv", ".tsv"):
                with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
                    pairs.extend(ids_from_table(list(csv.reader(handle, delimiter="\t" if suffix == ".tsv" else ","))))
            else:
                for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
                    found = ids_from_text(line)
                    if found:
                        pairs.extend(found)
                    elif line.strip() and not line.lstrip().startswith("#"):
                        unparsed.append(line.strip())
        else:
            found = ids_from_text(value)
            if found:
                pairs.extend(found)
            elif value.strip():
                unparsed.append(value.strip())
    return pairs, unparsed


def is_blocked(row: dict[str, str]) -> bool:
    level = (row.get("Level") or "").strip().lower()
    return level in BLOCKED_LEVELS or "do not include" in (row.get("Cleaned External Info") or "").lower()


def load_known_x(path: Path = KNOWN_X) -> dict[str, list[dict[str, str]]]:
    by_id: dict[str, list[dict[str, str]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("Unknown"):
                by_id.setdefault(row["Unknown"].strip(), []).append(row)
    return by_id


def screen(pairs: list[tuple[str, str]], known: dict[str, list[dict[str, str]]], *,
           allow_excluded: bool) -> list[dict[str, Any]]:
    seen: dict[str, dict[str, Any]] = {}
    for raw, cid in pairs:
        if cid in seen:
            seen[cid]["duplicates"] += 1
            continue
        rows = known.get(cid, [])
        levels = sorted({(r.get("Level") or "").strip() for r in rows} - {""})
        blocked_rows = [r for r in rows if is_blocked(r)]
        reason = ""
        if blocked_rows:
            note = (blocked_rows[0].get("Tentative information") or "").strip()
            reason = f"known_X {blocked_rows[0].get('Level') or 'DO NOT INCLUDE'}: {note}".strip()
        status = "excluded" if blocked_rows and not allow_excluded else "included"
        seen[cid] = {
            "control_id": cid, "input_token": raw, "duplicates": 0, "status": status,
            "in_known_x": bool(rows), "known_x_level": ";".join(levels),
            "known_x_sublevel": ";".join(sorted({(r.get("SubLevel") or "").strip() for r in rows} - {""})),
            "known_x_libs": ";".join(sorted({(r.get("LIB_ID") or "").split(".")[0] for r in rows} - {""})),
            "exclusion_reason": reason, "source_status": "",
        }
    return list(seen.values())


def run(cmd: list[str], log: Path, *, check: bool = True) -> int:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a", encoding="utf-8") as handle:
        handle.write(f"\n$ {' '.join(cmd)}\n")
        handle.flush()
        proc = subprocess.run(cmd, cwd=str(REPO_ROOT), stdout=handle, stderr=subprocess.STDOUT)
    if check and proc.returncode != 0:
        raise SystemExit(f"command failed ({proc.returncode}): {' '.join(cmd)}; see {log}")
    return proc.returncode


def copy_bundles(source_root: Path, ids: list[str], dest: Path) -> list[str]:
    copied = []
    for cid in ids:
        src = source_root / cid
        if (src / "manifest.json").is_file():
            target = dest / cid
            if target.exists():
                shutil.rmtree(target)
            shutil.copytree(src, target)
            copied.append(cid)
    return copied


def write_index(reports_dir: Path, predictions: list[str]) -> None:
    from structural_elucidation.reports import batch

    rows = batch.collect_index_rows(reports_dir, [])
    source_index = reports_dir / "_source_reports_index.csv"
    if source_index.is_file():
        with source_index.open(newline="", encoding="utf-8") as handle:
            libs = {r["control_id"]: r.get("lib", "") for r in csv.DictReader(handle)}
        for row in rows:
            row.lib = libs.get(row.control_id, row.lib)
    summary = batch.BatchSummary(generated_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                                 out_dir=str(reports_dir), jobs=0, predictions=predictions,
                                 total_control_ids=len(rows))
    batch.write_master_outputs(reports_dir, rows, summary)


def linked_files(html_path: Path) -> set[str]:
    text = html_path.read_text(encoding="utf-8", errors="replace")
    links = set()
    for target in re.findall(r"""(?:href|src)=["']([^"'#?]+)""", text):
        if re.match(r"^[a-z]+:", target) or target.startswith(("/", "..")):
            continue
        links.add(target)
    return links


def build_package(reports_dir: Path, package: Path, ids: list[str], workbook: Path | None) -> dict[str, Any]:
    if package.exists():
        shutil.rmtree(package)
    package.mkdir(parents=True)
    problems: list[str] = []
    for cid in ids:
        bundle = reports_dir / cid
        public = bundle / "public.html"
        if not public.is_file():
            continue
        (package / cid).mkdir()
        shutil.copy2(public, package / cid / "index.html")
        for link in sorted(linked_files(public)):
            src = bundle / link
            if src.is_file():
                (package / cid / link).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, package / cid / link)
            else:
                problems.append(f"{cid}: public.html links missing file {link}")
    for name in ("index.html", "reports_index.csv"):
        if (reports_dir / name).is_file():
            shutil.copy2(reports_dir / name, package / name)
    if workbook and workbook.is_file():
        shutil.copy2(workbook, package / workbook.name)
    for path in package.rglob("*"):
        if any(marker in path.name for marker in PUBLIC_FORBIDDEN):
            problems.append(f"internal-only file in client package: {path.relative_to(package)}")
    index = package / "index.html"
    if index.is_file():
        for link in linked_files(index):
            if not (package / link).is_file():
                problems.append(f"master index links missing file {link}")
    return {"path": str(package), "n_reports": sum(1 for p in package.glob("X-*/index.html")), "problems": problems}


def msi_counts(exec_dir: Path) -> dict[str, int]:
    path = exec_dir / "exec_summary_by_compound.csv"
    counts: dict[str, int] = {}
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                key = row.get("msi_level") or "blank"
                counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--compounds", nargs="+", required=True,
                        help="Files (.txt/.csv/.tsv/.xlsx) and/or literal ids, e.g. X-12345 'X 2249'")
    parser.add_argument("--client", required=True, help="Client label used for the package folder, e.g. UPMC")
    parser.add_argument("--out", type=Path, required=True, help="Delivery working directory")
    source = parser.add_argument_group("result sources (default: resolve_latest.py over --search)")
    source.add_argument("--search", nargs="+", type=Path, default=[REPO_ROOT / "output", REPO_ROOT / "models"])
    source.add_argument("--reports-root", type=Path, default=None, help="Reuse bundles from a rendered delivery")
    source.add_argument("--predictions-root", type=Path, default=None, help="Root with rerank_lib*/<name>.csv")
    source.add_argument("--predictions-name", default=None, help="Per-lib CSV name (default: chosen by resolver)")
    source.add_argument("--scores-dir", type=Path, default=None, help="Canonical <control_id>.json scores")
    source.add_argument("--mode", choices=("auto", "reuse", "regenerate"), default="auto",
                        help="auto = reuse a rendered v2 bundle covering every id, else regenerate")
    build = parser.add_argument_group("build")
    build.add_argument("--jobs", type=int, default=4)
    build.add_argument("--candidate-tab-n", default="10", help="Candidates per report tab (deliveries used 10)")
    build.add_argument("--observed-lookup-dir", type=Path, default=None)
    build.add_argument("--db", type=Path, default=None, help="metabolite_synonyms.sqlite (default vendor/)")
    build.add_argument("--pure-nce-scores", type=Path, default=None)
    build.add_argument("--allow-excluded", action="store_true",
                       help="Keep known_X Exclude/Artifact compounds (never for an external client)")
    build.add_argument("--max-compounds", type=int, default=60,
                       help="Refuse to regenerate more than this many; hand the printed command to a human")
    build.add_argument("--zip", action="store_true", help="Also write <client>_x_compound_reports.zip")
    build.add_argument("--dry-run", action="store_true", help="Screen and resolve sources only")
    args = parser.parse_args(argv)

    out = args.out.expanduser().resolve()
    work = out / "work"
    log = out / "build.log"
    work.mkdir(parents=True, exist_ok=True)
    started = time.time()

    pairs, unparsed = read_compound_inputs(args.compounds)
    if not pairs:
        raise SystemExit("no X-compound ids found in --compounds")
    rows = screen(pairs, load_known_x(), allow_excluded=args.allow_excluded)
    included = [r["control_id"] for r in rows if r["status"] == "included"]
    (out / "compounds_normalized.txt").write_text("\n".join(included) + "\n", encoding="utf-8")

    latest = resolve_latest.resolve([p.expanduser().resolve() for p in args.search], set(included))
    (out / "latest_results.json").write_text(json.dumps(latest, indent=2) + "\n", encoding="utf-8")
    chosen = latest["chosen"]
    reports_root = args.reports_root or (Path(chosen["reports"]["root"]) if chosen["reports"] else None)
    predictions_root = args.predictions_root or (Path(chosen["predictions"]["root"]) if chosen["predictions"] else None)
    chosen_name = (chosen["predictions"] or {}).get("predictions_name")
    predictions_name = args.predictions_name or chosen_name or "unified_rerank_all_top10.csv"
    fresh_scores = chosen["scores"] if chosen["scores"] and chosen["scores"]["calibration_matches_repo"] else None
    scores_dir = args.scores_dir or (Path(fresh_scores["root"]) if fresh_scores else None)

    mode = args.mode
    if mode == "auto":
        cover = set()
        if reports_root and reports_root.is_dir():
            cover = {p.parent.name for p in reports_root.glob("X-*/manifest.json")}
        mode = "reuse" if args.reports_root or (cover and set(included) <= cover) else "regenerate"
    if mode == "reuse" and not reports_root:
        raise SystemExit("--mode reuse needs --reports-root (no rendered bundle was found)")
    if mode == "regenerate" and not predictions_root:
        raise SystemExit("no predictions found; pass --predictions-root or extend --search")

    reports_dir = work / "reports"
    batch_cmd = [sys.executable, "-m", "structural_elucidation.reports.batch",
                 "--equalization-root", str(predictions_root), "--predictions-name", predictions_name,
                 "--out-dir", str(reports_dir), "--control-ids", f"@{out / 'compounds_normalized.txt'}",
                 "--jobs", str(args.jobs), "--overwrite", "--candidate-tab-n", str(args.candidate_tab_n),
                 "--naming-policy-version", "v2"] if predictions_root else []
    if batch_cmd and scores_dir:
        batch_cmd += ["--scoring-dir", str(scores_dir)]
    if batch_cmd and args.observed_lookup_dir:
        batch_cmd += ["--observed-lookup-dir", str(args.observed_lookup_dir)]
    if batch_cmd and args.db:
        batch_cmd += ["--db", str(args.db)]

    manifest: dict[str, Any] = {
        "client": args.client, "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "repo_head": latest["repo_head"], "repo_calibration_sha256": latest["repo_calibration_sha256"],
        "mode": mode, "reports_root": str(reports_root) if reports_root else None,
        "predictions_root": str(predictions_root) if predictions_root else None,
        "predictions_name": predictions_name, "scores_dir": str(scores_dir) if scores_dir else None,
        "n_input_tokens": len(pairs), "unparsed_tokens": unparsed,
        "n_unique": len(rows), "n_included": len(included),
        "excluded": [{k: r[k] for k in ("control_id", "exclusion_reason")} for r in rows if r["status"] == "excluded"],
        "not_in_known_x": [r["control_id"] for r in rows if not r["in_known_x"]],
        "warnings": list(latest["warnings"]),
        "estimated_regenerate_seconds": len(included) * SECONDS_PER_COMPOUND // max(1, args.jobs),
        "regenerate_command": " ".join(batch_cmd) if batch_cmd else None,
    }
    curated = [r["control_id"] for r in rows
               if r["status"] == "included" and set(r["known_x_level"].split(";")) & {"1", "2"}]
    if scores_dir and curated:
        manifest["warnings"].append(
            f"Scoring JSONs are applied to {len(curated)} compounds with a curated known_X Level 1/2 "
            f"({', '.join(curated[:5])}{'...' if len(curated) > 5 else ''}); novel scoring ignores known_X levels, "
            "so confirm the JSONs were produced for registered X compounds before delivery")
    if scores_dir is None and chosen["scores"]:
        manifest["warnings"].append(f"Scoring JSONs at {chosen['scores']['root']} were not used: calibration does "
                                    "not match this checkout. Pass --scores-dir to force them")

    def finish(status: str) -> int:
        manifest["status"] = status
        manifest["elapsed_seconds"] = round(time.time() - started, 1)
        with (out / "compound_screen.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        (out / "delivery_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        brief = {k: v for k, v in manifest.items() if k not in ("audit", "regenerate_command")}
        if "audit" in manifest:
            brief["audit"] = {k: v for k, v in manifest["audit"].items() if isinstance(v, (int, float, str))}
        print(json.dumps(brief, indent=2))
        print(f"[client_report] manifest: {out / 'delivery_manifest.json'}", file=sys.stderr)
        return 0 if status in ("ok", "dry_run") else 1

    if args.dry_run:
        return finish("dry_run")
    if mode == "regenerate" and len(included) > args.max_compounds:
        manifest["warnings"].append(f"{len(included)} compounds exceeds --max-compounds {args.max_compounds}; "
                                    "run regenerate_command manually, then rerun with --mode reuse "
                                    f"--reports-root {reports_dir}")
        return finish("needs_manual_run")

    in_place = mode == "reuse" and reports_root.resolve() == reports_dir
    if not in_place:
        if reports_dir.exists():
            shutil.rmtree(reports_dir)
        reports_dir.mkdir(parents=True)
    if in_place:
        copied = sorted(cid for cid in included if (reports_dir / cid / "manifest.json").is_file())
        write_index(reports_dir, [str(reports_root)])
    elif mode == "reuse":
        copied = copy_bundles(reports_root, included, reports_dir)
        if (reports_root / "reports_index.csv").is_file():
            shutil.copy2(reports_root / "reports_index.csv", reports_dir / "_source_reports_index.csv")
        write_index(reports_dir, [str(reports_root)])
        (reports_dir / "_source_reports_index.csv").unlink(missing_ok=True)
    else:
        run(batch_cmd, log)
        copied = sorted(p.parent.name for p in reports_dir.glob("X-*/manifest.json"))
    missing = sorted(set(included) - set(copied))
    for row in rows:
        if row["status"] == "included":
            row["source_status"] = "report" if row["control_id"] in copied else "no_report"
    manifest["n_reports"] = len(copied)
    manifest["no_report"] = missing
    if missing:
        manifest["warnings"].append(f"{len(missing)} requested compounds have no report in the chosen source")

    exec_dir = work / "exec_summary"
    summary_cmd = [sys.executable, "tools/build_exec_summary.py", "--reports-root", str(reports_dir),
                   "--out-dir", str(exec_dir)]
    if args.pure_nce_scores:
        summary_cmd += ["--pure-nce-scores", str(args.pure_nce_scores)]
    workbook = None
    if copied:
        run(summary_cmd, log)
        run([sys.executable, "tools/build_exec_summary_workbook.py", "--exec-dir", str(exec_dir)], log)
        workbook = exec_dir / "executive_summary_top10.xlsx"
        run([sys.executable, "tools/audit_reports.py", "--reports-root", str(reports_dir),
             "--out-dir", str(work / "audit")], log, check=False)
    audit_summary = work / "audit" / "audit_summary.json"
    if audit_summary.is_file():
        manifest["audit"] = json.loads(audit_summary.read_text(encoding="utf-8"))
    manifest["msi_level_counts"] = msi_counts(exec_dir)

    safe_client = re.sub(r"[^A-Za-z0-9]+", "_", args.client).strip("_") or "client"
    package = build_package(reports_dir, out / f"{safe_client}_x_compound_reports", copied, workbook)
    manifest["package"] = package
    if package["problems"]:
        manifest["warnings"].extend(package["problems"])
    if args.zip:
        manifest["zip"] = shutil.make_archive(str(out / f"{safe_client}_x_compound_reports"), "zip",
                                              root_dir=out, base_dir=Path(package["path"]).name)
    return finish("ok" if not package["problems"] else "package_problems")


if __name__ == "__main__":
    sys.exit(main())
