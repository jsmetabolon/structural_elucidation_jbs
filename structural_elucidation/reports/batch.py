"""Batch HTML report generation for every X compound in the rerank outputs.

Wraps the single-compound builder (``structural_elucidation.reports.compound``)
and runs it once per ``control_id`` discovered in one or more rerank/prediction
CSVs (typically ``<root>/rerank_lib{209,305,400,402}/unified_rerank_all_top10.csv``).
Each compound is built in an isolated subprocess so one failure never aborts the
batch, results can be produced in parallel (``--jobs``), and partial runs resume
by skipping compounds that already have an ``index.html``.

After per-compound generation it rolls everything up into a master ``index.html``
(searchable/sortable table), ``reports_index.csv``, and ``batch_summary.json``.

CLI examples::

    # Discover the 4 lib CSVs under an equalization root and build all reports
    python -m structural_elucidation.reports.batch \\
        --equalization-root output/pipeline_equalization \\
        --out-dir output/reports --jobs 8

    # Explicit prediction files
    python -m structural_elucidation.reports.batch \\
        --predictions a/unified_rerank_all_top10.csv b/unified_rerank_all_top10.csv \\
        --out-dir output/reports

    # Just enumerate what would be built (fast, no generation)
    python -m structural_elucidation.reports.batch \\
        --equalization-root output/pipeline_equalization --dry-run
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from threading import Lock
from typing import Any

import pandas as pd

from structural_elucidation.paths import REPO_ROOT, resolve
from structural_elucidation.reports.index_render import (
    IndexRow,
    _read_confidence_extra,
    render_master_index,
)

ID_COLUMNS = ("control_id", "x_compound_id")
DEFAULT_PREDICTIONS_NAME = "unified_rerank_all_top10.csv"
DEFAULT_LIBS = ("209", "305", "400", "402")
LIB_ION_MODE = {"209": "negative", "305": "negative", "400": "positive", "402": "positive"}


# ── Discovery + enumeration ───────────────────────────────────────────────────

def lib_from_path(path: Path) -> str:
    """Extract the library id from a ``.../rerank_lib<NNN>/...`` path, else ``""``."""

    match = re.search(r"rerank_lib(\d+)", str(path))
    return match.group(1) if match else ""


def discover_predictions(args: argparse.Namespace) -> list[Path]:
    """Resolve the ordered list of prediction CSVs from CLI inputs.

    Explicit ``--predictions`` paths win; otherwise build the per-lib paths under
    ``--equalization-root`` for each requested lib that exists on disk.
    """

    paths: list[Path] = []
    if args.predictions:
        for raw in args.predictions:
            paths.append(resolve(Path(raw), REPO_ROOT))
    elif args.equalization_root:
        root = resolve(Path(args.equalization_root), REPO_ROOT)
        libs = [lib.strip() for lib in str(args.libs).split(",") if lib.strip()]
        for lib in libs:
            candidate = root / f"rerank_lib{lib}" / args.predictions_name
            if candidate.exists():
                paths.append(candidate)
            else:
                print(f"[batch] WARN: missing predictions for lib {lib}: {candidate}", file=sys.stderr)
    else:
        raise SystemExit("Provide either --predictions or --equalization-root")

    existing: list[Path] = []
    for path in paths:
        if path.exists():
            existing.append(path)
        else:
            print(f"[batch] WARN: predictions file not found: {path}", file=sys.stderr)
    if not existing:
        raise SystemExit("No usable prediction CSVs were found.")
    return existing


def detect_id_column(path: Path) -> str:
    """Return the id column name present in ``path`` header (control_id|x_compound_id)."""

    header = pd.read_csv(path, nrows=0)
    for col in ID_COLUMNS:
        if col in header.columns:
            return col
    raise SystemExit(f"{path} has none of the id columns {ID_COLUMNS}")


def enumerate_control_ids(predictions: list[Path]) -> list[tuple[str, str, Path]]:
    """Return ordered unique ``(control_id, lib, source_csv)`` across all CSVs.

    First occurrence of a control id wins when it appears in multiple files.
    """

    seen: set[str] = set()
    rows: list[tuple[str, str, Path]] = []
    for path in predictions:
        id_col = detect_id_column(path)
        series = pd.read_csv(path, usecols=[id_col])[id_col]
        lib = lib_from_path(path)
        for raw in series:
            if pd.isna(raw):
                continue
            cid = str(raw).strip()
            if not cid or cid.lower() in {"nan", "none"}:
                continue
            if cid in seen:
                continue
            seen.add(cid)
            rows.append((cid, lib, path))
    return rows


def filter_control_ids(
    rows: list[tuple[str, str, Path]],
    only: set[str] | None,
    limit: int | None,
) -> list[tuple[str, str, Path]]:
    """Apply optional ``--control-ids`` allowlist and ``--limit`` cap, preserving order."""

    out = rows
    if only:
        out = [r for r in out if r[0] in only]
    if limit is not None and limit >= 0:
        out = out[:limit]
    return out


# ── Per-compound build ────────────────────────────────────────────────────────

@dataclass
class BuildResult:
    """Outcome of one per-compound build subprocess."""

    control_id: str
    lib: str
    out_dir: Path
    status: str  # generated | skipped | failed
    returncode: int = 0
    error_tail: str = ""


def passthrough_args(args: argparse.Namespace, lib: str) -> list[str]:
    """Build the optional flag list forwarded to the single-compound CLI.

    Only flags the user explicitly set are forwarded so the compound module stays
    the single source of default values.
    """

    extra: list[str] = []
    simple = {
        "--db": args.db,
        "--known-x": args.known_x,
        "--known-x-all": args.known_x_all,
        "--spectrum-tolerance-da": args.spectrum_tolerance_da,
        "--min-similarity": args.min_similarity,
        "--max-neighbors": args.max_neighbors,
        "--candidate-tab-n": args.candidate_tab_n,
        "--consensus-top-n": args.consensus_top_n,
        "--extended-known-render-limit": args.extended_known_render_limit,
        "--extended-known-min-similarity": args.extended_known_min_similarity,
        "--extended-x-min-similarity": args.extended_x_min_similarity,
        "--extended-target-edge-min": args.extended_target_edge_min,
        "--extended-pair-edge-min": args.extended_pair_edge_min,
        "--extended-max-edges-per-node": args.extended_max_edges_per_node,
    }
    for flag, value in simple.items():
        if value is not None:
            extra.extend([flag, str(value)])

    if args.observed_lookup_dir:
        ion = LIB_ION_MODE.get(lib, "")
        lookup = resolve(Path(args.observed_lookup_dir), REPO_ROOT) / f"lib{lib}_{ion}_observed_ms2_lookup.csv"
        if lookup.exists():
            extra.extend(["--peak-csv", str(lookup)])
    elif args.peak_csv:
        extra.extend(["--peak-csv", str(resolve(Path(args.peak_csv), REPO_ROOT))])
    return extra


def build_one(
    control_id: str,
    lib: str,
    predictions: Path,
    out_root: Path,
    args: argparse.Namespace,
) -> BuildResult:
    """Generate one compound bundle via an isolated subprocess; never raises."""

    import subprocess

    out_dir = out_root / control_id
    index_html = out_dir / "index.html"
    if index_html.exists() and not args.overwrite:
        return BuildResult(control_id, lib, out_dir, "skipped")

    cmd = [
        sys.executable,
        "-m",
        "structural_elucidation.reports.compound",
        "--control-id",
        control_id,
        "--predictions",
        str(predictions),
        "--out-dir",
        str(out_dir),
        *passthrough_args(args, lib),
    ]
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            timeout=args.per_compound_timeout if args.per_compound_timeout and args.per_compound_timeout > 0 else None,
        )
    except subprocess.TimeoutExpired:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "_build.log").write_text("TIMEOUT\n", encoding="utf-8")
        return BuildResult(control_id, lib, out_dir, "failed", -1, "timeout")

    if proc.returncode != 0 or not index_html.exists():
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "_build.log").write_text(
            f"$ {' '.join(cmd)}\n\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}\n",
            encoding="utf-8",
        )
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()
        return BuildResult(control_id, lib, out_dir, "failed", proc.returncode, "\n".join(tail[-4:]))
    stale_log = out_dir / "_build.log"
    if stale_log.exists():
        stale_log.unlink()
    return BuildResult(control_id, lib, out_dir, "generated", 0)


# ── Master index + summary ────────────────────────────────────────────────────

def collect_index_rows(out_root: Path, results: list[BuildResult]) -> list[IndexRow]:
    """Read each compound ``manifest.json`` and merge in batch build status."""

    lib_by_id = {r.control_id: r.lib for r in results}
    status_by_id = {r.control_id: r.status for r in results}
    rows: list[IndexRow] = []
    for manifest_path in sorted(out_root.glob("*/manifest.json")):
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        cid = str(payload.get("control_id") or manifest_path.parent.name)
        target = payload.get("target") or {}
        extra = _read_confidence_extra(manifest_path.parent)
        rows.append(
            IndexRow(
                control_id=cid,
                lib=lib_by_id.get(cid, ""),
                name=str(extra.get("display_name") or target.get("name") or ""),
                formula=str(target.get("formula") or ""),
                exact_mass=("" if target.get("exact_mass") in (None, "") else f"{float(target['exact_mass']):.4f}"),
                db_chemical_id=str(payload.get("db_chemical_id") or ""),
                n_neighbors=int(payload.get("n_neighbors") or 0),
                status=status_by_id.get(cid, "existing"),
                report=f"{manifest_path.parent.name}/index.html",
                msi_level=str(extra.get("msi_level") or ""),
                msi_label=str(extra.get("msi_label") or ""),
                identity_status=str(extra.get("identity_status") or ""),
                identity_confidence=str(extra.get("identity_confidence") or ""),
                descriptor_kind=str(extra.get("descriptor_kind") or ""),
                descriptor_confidence=str(extra.get("descriptor_confidence") or ""),
                descriptor_parent_formula=str(extra.get("descriptor_parent_formula") or ""),
                descriptor_conjugate_count=str(extra.get("descriptor_conjugate_count") or ""),
            )
        )
    rows.sort(key=lambda r: (r.lib, r.control_id))
    return rows


def write_reports_index_csv(out_root: Path, rows: list[IndexRow]) -> Path:
    """Write the flat ``reports_index.csv`` table of all generated compounds."""

    df = pd.DataFrame(
        [
            {
                "control_id": r.control_id,
                "lib": r.lib,
                "name": r.name,
                "formula": r.formula,
                "exact_mass": r.exact_mass,
                "db_chemical_id": r.db_chemical_id,
                "n_neighbors": r.n_neighbors,
                "msi_level": r.msi_level,
                "msi_label": r.msi_label,
                "descriptor_kind": r.descriptor_kind,
                "descriptor_confidence": r.descriptor_confidence,
                "descriptor_parent_formula": r.descriptor_parent_formula,
                "descriptor_conjugate_count": r.descriptor_conjugate_count,
                "status": r.status,
                "report": r.report,
            }
            for r in rows
        ]
    )
    path = out_root / "reports_index.csv"
    df.to_csv(path, index=False)
    return path


@dataclass
class BatchSummary:
    """Aggregate counts + failures for ``batch_summary.json``."""

    generated_at: str
    out_dir: str
    jobs: int
    predictions: list[str]
    total_control_ids: int
    generated: int = 0
    skipped_existing: int = 0
    failed: list[dict[str, Any]] = field(default_factory=list)
    elapsed_seconds: float = 0.0


def write_master_outputs(out_root: Path, rows: list[IndexRow], summary: BatchSummary) -> dict[str, str]:
    """Write master index.html, reports_index.csv, and batch_summary.json."""

    out_root.mkdir(parents=True, exist_ok=True)
    index_path = out_root / "index.html"
    index_path.write_text(render_master_index(rows, summary.generated_at), encoding="utf-8")
    csv_path = write_reports_index_csv(out_root, rows)
    summary_path = out_root / "batch_summary.json"
    summary_path.write_text(json.dumps(summary.__dict__, indent=2), encoding="utf-8")
    return {"index": str(index_path), "csv": str(csv_path), "summary": str(summary_path)}


# ── Orchestration ─────────────────────────────────────────────────────────────

def run_batch(args: argparse.Namespace) -> dict[str, Any]:
    """Discover compounds, build each report (parallel), and write master outputs."""

    started = time.time()
    generated_at = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    predictions = discover_predictions(args)
    print(f"[batch] {len(predictions)} prediction CSV(s):")
    for path in predictions:
        print(f"[batch]   - {path}")

    all_rows = enumerate_control_ids(predictions)
    only = None
    if args.control_ids:
        only = _read_control_ids_arg(args.control_ids)
    work = filter_control_ids(all_rows, only, args.limit)

    per_lib: dict[str, int] = {}
    for _, lib, _ in work:
        per_lib[lib] = per_lib.get(lib, 0) + 1
    print(f"[batch] {len(all_rows)} unique control ids discovered; {len(work)} selected for this run")
    for lib in sorted(per_lib):
        print(f"[batch]   lib {lib or '?'}: {per_lib[lib]}")

    out_root = resolve(Path(args.out_dir), REPO_ROOT)

    if args.dry_run:
        print("[batch] dry-run: no reports generated. Planned output dirs:")
        for cid, lib, src in work[: args.dry_run_preview]:
            print(f"[batch]   {out_root / cid}  (lib {lib or '?'}, from {src.name})")
        if len(work) > args.dry_run_preview:
            print(f"[batch]   ... and {len(work) - args.dry_run_preview} more")
        return {
            "dry_run": True,
            "total_control_ids": len(all_rows),
            "selected": len(work),
            "per_lib": per_lib,
            "out_dir": str(out_root),
        }

    if args.index_only:
        rows = collect_index_rows(out_root, [])
        summary = BatchSummary(
            generated_at=generated_at,
            out_dir=str(out_root),
            jobs=args.jobs,
            predictions=[str(p) for p in predictions],
            total_control_ids=len(rows),
        )
        paths = write_master_outputs(out_root, rows, summary)
        print(f"[batch] index-only: wrote master index for {len(rows)} compounds -> {paths['index']}")
        return {"index_only": True, **paths, "n": len(rows)}

    out_root.mkdir(parents=True, exist_ok=True)
    results: list[BuildResult] = []
    lock = Lock()
    done = 0
    total = len(work)

    def task(item: tuple[str, str, Path]) -> BuildResult:
        cid, lib, src = item
        return build_one(cid, lib, src, out_root, args)

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {pool.submit(task, item): item for item in work}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            with lock:
                done += 1
                marker = {"generated": "ok", "skipped": "--", "failed": "XX"}.get(result.status, "??")
                print(f"[batch] [{done}/{total}] {marker} {result.control_id} (lib {result.lib or '?'})")
                if result.status == "failed" and result.error_tail:
                    print(f"[batch]        {result.error_tail.splitlines()[-1] if result.error_tail else ''}")

    summary = BatchSummary(
        generated_at=generated_at,
        out_dir=str(out_root),
        jobs=args.jobs,
        predictions=[str(p) for p in predictions],
        total_control_ids=len(all_rows),
        generated=sum(1 for r in results if r.status == "generated"),
        skipped_existing=sum(1 for r in results if r.status == "skipped"),
        failed=[
            {"control_id": r.control_id, "lib": r.lib, "returncode": r.returncode, "error_tail": r.error_tail}
            for r in results
            if r.status == "failed"
        ],
        elapsed_seconds=round(time.time() - started, 2),
    )
    rows = collect_index_rows(out_root, results)
    paths = write_master_outputs(out_root, rows, summary)

    print(
        f"[batch] done: {summary.generated} generated, {summary.skipped_existing} skipped, "
        f"{len(summary.failed)} failed in {summary.elapsed_seconds}s"
    )
    print(f"[batch] master index: {paths['index']}")
    print(f"[batch] summary: {paths['summary']}")
    if summary.failed:
        print(f"[batch] {len(summary.failed)} failures recorded in batch_summary.json")
    return {"dry_run": False, **paths, "summary_obj": summary.__dict__}


def _read_control_ids_arg(value: str) -> set[str]:
    """Parse ``--control-ids`` as ``@file`` (one per line) or comma-separated list."""

    if value.startswith("@"):
        path = resolve(Path(value[1:]), REPO_ROOT)
        return {line.strip() for line in path.read_text().splitlines() if line.strip()}
    return {item.strip() for item in value.split(",") if item.strip()}


def build_parser() -> argparse.ArgumentParser:
    """Construct the batch CLI argument parser."""

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = parser.add_argument_group("inputs")
    src.add_argument("--predictions", nargs="+", help="One or more rerank/prediction CSV paths")
    src.add_argument("--equalization-root", help="Root holding rerank_lib*/<predictions-name>")
    src.add_argument("--libs", default=",".join(DEFAULT_LIBS), help="Comma libs for --equalization-root")
    src.add_argument("--predictions-name", default=DEFAULT_PREDICTIONS_NAME, help="Per-lib CSV filename")

    out = parser.add_argument_group("outputs / control")
    out.add_argument("--out-dir", default="output/reports", help="Master output dir (one subdir per compound)")
    out.add_argument("--jobs", type=int, default=1, help="Parallel compound builds")
    out.add_argument("--overwrite", action="store_true", help="Rebuild compounds that already have index.html")
    out.add_argument("--limit", type=int, default=None, help="Only build the first N compounds (quick test)")
    out.add_argument("--control-ids", default=None, help="Restrict to comma list or @file of control ids")
    out.add_argument("--dry-run", action="store_true", help="Enumerate compounds and exit without building")
    out.add_argument("--dry-run-preview", type=int, default=20, help="How many planned dirs to print in --dry-run")
    out.add_argument("--index-only", action="store_true", help="Rebuild master index from existing manifests only")
    out.add_argument("--per-compound-timeout", type=int, default=0, help="Seconds per compound (0 = no timeout)")

    # Optional passthrough to the single-compound builder (forwarded only when set).
    fwd = parser.add_argument_group("forwarded to reports.compound (optional)")
    fwd.add_argument("--db", default=None)
    fwd.add_argument("--known-x", default=None)
    fwd.add_argument("--known-x-all", default=None)
    fwd.add_argument("--peak-csv", default=None, help="Single peak CSV applied to every compound")
    fwd.add_argument("--observed-lookup-dir", default=None, help="Dir of lib{N}_{ion}_observed_ms2_lookup.csv")
    fwd.add_argument("--spectrum-tolerance-da", default=None)
    fwd.add_argument("--min-similarity", default=None)
    fwd.add_argument("--max-neighbors", default=None)
    fwd.add_argument("--candidate-tab-n", default=None)
    fwd.add_argument("--consensus-top-n", default=None)
    fwd.add_argument("--extended-known-render-limit", default=None)
    fwd.add_argument("--extended-known-min-similarity", default=None)
    fwd.add_argument("--extended-x-min-similarity", default=None)
    fwd.add_argument("--extended-target-edge-min", default=None)
    fwd.add_argument("--extended-pair-edge-min", default=None)
    fwd.add_argument("--extended-max-edges-per-node", default=None)
    return parser


def main() -> None:
    """CLI entrypoint: parse args and run the batch."""

    args = build_parser().parse_args()
    run_batch(args)


if __name__ == "__main__":
    main()
