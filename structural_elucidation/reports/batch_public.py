"""Batch-generate public compound report bundles from prediction CSVs."""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import pandas as pd

from structural_elucidation.paths import REPO_ROOT, resolve
from structural_elucidation.reports import compound


DEFAULT_EQUALIZATION_ROOT = Path("output/pipeline_equalization")
DEFAULT_INPUT_NAME = "unified_rerank_all_top10_with_pathways_lipid_class_peptide_composition.csv"
DEFAULT_OUT_ROOT = Path("output/public_compound_reports")
DEFAULT_LIBS = "209,305,400,402"
CONTROL_COLUMNS = ("control_id", "x_compound_id", "compound_id")
RANK_COLUMNS = ("rank", "unified_pool_rank", "candidate_rank")
MANIFEST_FIELDS = [
    "slug",
    "control_id",
    "lib_id",
    "ion_mode",
    "status",
    "report_dir",
    "index_html",
    "public_html",
    "predictions_csv",
    "n_prediction_rows",
    "rank1_name",
    "rank1_smiles",
    "rank1_formula",
    "predicted_pathway",
    "predicted_subpathway",
    "lipid_class",
    "lipid_shorthand",
    "peptide_composition_label",
]
FAILURE_FIELDS = ["slug", "control_id", "lib_id", "predictions_csv", "reason", "traceback"]


@dataclass(frozen=True)
class PredictionCsv:
    path: Path
    lib_id: str = ""


@dataclass(frozen=True)
class ReportJob:
    control_id: str
    lib_id: str
    predictions_csv: Path
    rows: pd.DataFrame
    slug: str
    out_dir: Path


def parse_libs(values: list[str] | tuple[str, ...] | str | None) -> list[int]:
    if values is None:
        return []
    raw_values = [values] if isinstance(values, str) else list(values)
    libs: list[int] = []
    for value in raw_values:
        for token in str(value).split(","):
            token = token.strip()
            if token:
                libs.append(int(token))
    return libs


def clean_text(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "null", "<na>"} else text


def infer_lib_id(path: Path, frame: pd.DataFrame | None = None) -> str:
    if frame is not None and "lib_id" in frame.columns:
        libs = sorted({clean_text(value) for value in frame["lib_id"].dropna().unique() if clean_text(value)})
        if len(libs) == 1:
            return libs[0]
    match = re.search(r"(?:rerank_)?lib[_-]?(\d+)|lib_(\d+)", str(path))
    if match:
        return next(group for group in match.groups() if group)
    return ""


def slugify(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", clean_text(value))
    return slug.strip("._") or "report"


def report_slug(lib_id: str, control_id: str) -> str:
    control_slug = slugify(control_id)
    return f"lib{slugify(lib_id)}_{control_slug}" if clean_text(lib_id) else control_slug


def discover_prediction_csvs(equalization_root: Path, libs: list[int], input_name: str) -> list[PredictionCsv]:
    root = resolve(equalization_root, REPO_ROOT)
    csvs: list[PredictionCsv] = []
    for lib_id in libs:
        path = root / f"rerank_lib{int(lib_id)}" / input_name
        if path.exists():
            csvs.append(PredictionCsv(path=path, lib_id=str(int(lib_id))))
    return csvs


def load_prediction_index(prediction_csvs: list[PredictionCsv]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for item in prediction_csvs:
        path = resolve(item.path, REPO_ROOT)
        frame = pd.read_csv(path, low_memory=False)
        control_col = next((column for column in CONTROL_COLUMNS if column in frame.columns), "")
        if not control_col:
            raise ValueError(f"{path} missing one of: {', '.join(CONTROL_COLUMNS)}")
        frame = frame.copy()
        frame["_control_id"] = frame[control_col].map(clean_text)
        frame["_source_csv"] = str(path)
        frame["_source_lib_id"] = item.lib_id or infer_lib_id(path, frame)
        if "lib_id" in frame.columns:
            frame["_job_lib_id"] = frame["lib_id"].map(clean_text)
            frame.loc[frame["_job_lib_id"].eq(""), "_job_lib_id"] = frame["_source_lib_id"]
        else:
            frame["_job_lib_id"] = frame["_source_lib_id"]
        frames.append(frame.loc[frame["_control_id"].str.startswith("X-", na=False)].copy())
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def sort_rows_for_rank(frame: pd.DataFrame) -> pd.DataFrame:
    rank_col = next((column for column in RANK_COLUMNS if column in frame.columns), "")
    if not rank_col:
        return frame.copy()
    out = frame.copy()
    out["_rank_sort"] = pd.to_numeric(out[rank_col], errors="coerce").fillna(10**9)
    return out.sort_values("_rank_sort", kind="mergesort").drop(columns=["_rank_sort"], errors="ignore")


def iter_report_jobs(
    index: pd.DataFrame,
    *,
    out_root: Path,
    control_ids: set[str] | None = None,
    limit: int | None = None,
) -> list[ReportJob]:
    if index.empty:
        return []
    root = resolve(out_root, REPO_ROOT)
    filters = {clean_text(value) for value in control_ids or set() if clean_text(value)}
    jobs: list[ReportJob] = []
    group_cols = ["_source_csv", "_job_lib_id", "_control_id"]
    for (source_csv, lib_id, control_id), group in index.groupby(group_cols, sort=True, dropna=False):
        control_id = clean_text(control_id)
        if filters and control_id not in filters:
            continue
        lib_id = clean_text(lib_id)
        slug = report_slug(lib_id, control_id)
        rows = group.drop(columns=["_source_csv", "_source_lib_id", "_job_lib_id"], errors="ignore").copy()
        if "control_id" not in rows.columns:
            rows["control_id"] = control_id
        if lib_id and "lib_id" not in rows.columns:
            rows["lib_id"] = lib_id
        jobs.append(
            ReportJob(
                control_id=control_id,
                lib_id=lib_id,
                predictions_csv=Path(str(source_csv)),
                rows=sort_rows_for_rank(rows),
                slug=slug,
                out_dir=root / "compounds" / slug,
            )
        )
        if limit is not None and len(jobs) >= int(limit):
            break
    return jobs


def first_present(row: pd.Series, columns: tuple[str, ...] | list[str]) -> str:
    for column in columns:
        if column in row.index:
            value = clean_text(row.get(column))
            if value:
                return value
    return ""


def manifest_row_for_job(job: ReportJob, *, status: str) -> dict[str, Any]:
    rows = sort_rows_for_rank(job.rows)
    rank1 = rows.iloc[0] if not rows.empty else pd.Series(dtype=object)
    return {
        "slug": job.slug,
        "control_id": job.control_id,
        "lib_id": job.lib_id,
        "ion_mode": first_present(rank1, ("ion_mode", "ion_mode_label", "polarity")),
        "status": status,
        "report_dir": str(job.out_dir),
        "index_html": f"compounds/{job.slug}/index.html",
        "public_html": f"compounds/{job.slug}/public.html",
        "predictions_csv": str(job.predictions_csv),
        "n_prediction_rows": int(len(job.rows)),
        "rank1_name": first_present(rank1, ("candidate_name", "base_name", "name", "compound_name", "library_name")),
        "rank1_smiles": first_present(rank1, ("predicted_smiles", "smiles", "canonical_smiles", "candidate_smiles")),
        "rank1_formula": first_present(rank1, ("predicted_formula", "formula", "candidate_formula")),
        "predicted_pathway": first_present(rank1, ("predicted_pathway", "pathway")),
        "predicted_subpathway": first_present(rank1, ("predicted_subpathway", "subpathway")),
        "lipid_class": first_present(rank1, ("lipid_class",)),
        "lipid_shorthand": first_present(rank1, ("lipid_shorthand",)),
        "peptide_composition_label": first_present(rank1, ("peptide_composition_label",)),
    }


def resolve_peak_csv(job: ReportJob, args: argparse.Namespace) -> Path | None:
    if getattr(args, "peak_csv", None):
        peak_csv = resolve(args.peak_csv, REPO_ROOT)
        return peak_csv if peak_csv.exists() else None
    if getattr(args, "peak_root", None):
        peak_root = resolve(args.peak_root, REPO_ROOT)
        for name in (f"{job.slug}.csv", f"{job.control_id}.csv", f"lib{job.lib_id}_{job.control_id}.csv"):
            path = peak_root / name
            if path.exists():
                return path
    return None


def build_renderer_args(job: ReportJob, args: argparse.Namespace, predictions_csv: Path) -> argparse.Namespace:
    return argparse.Namespace(
        control_id=job.control_id,
        predictions=predictions_csv,
        out_dir=job.out_dir,
        db=getattr(args, "db", compound.DEFAULT_DB),
        known_x=getattr(args, "known_x", compound.DEFAULT_KNOWN_X),
        known_x_all=getattr(args, "known_x_all", compound.DEFAULT_KNOWN_X_ALL),
        candidate_pool=getattr(args, "candidate_pool", None),
        hypotheses=getattr(args, "hypotheses", None),
        peak_csv=resolve_peak_csv(job, args),
        target_smiles="",
        target_name="",
        chemical_id="",
        spectrum_tolerance_da=float(getattr(args, "spectrum_tolerance_da", 1.0)),
        min_similarity=float(getattr(args, "min_similarity", 0.28)),
        max_neighbors=int(getattr(args, "max_neighbors", 10)),
        candidate_tab_n=int(getattr(args, "candidate_tab_n", 20)),
        consensus_top_n=int(getattr(args, "consensus_top_n", 20)),
        extended_known_render_limit=int(getattr(args, "extended_known_render_limit", 120)),
        extended_known_min_similarity=float(getattr(args, "extended_known_min_similarity", 0.30)),
        extended_x_min_similarity=float(getattr(args, "extended_x_min_similarity", 0.0)),
        extended_target_edge_min=float(getattr(args, "extended_target_edge_min", 0.32)),
        extended_pair_edge_min=float(getattr(args, "extended_pair_edge_min", 0.42)),
        extended_max_edges_per_node=int(getattr(args, "extended_max_edges_per_node", 6)),
    )


def build_one_report(
    job: ReportJob,
    args: argparse.Namespace,
    *,
    renderer: Callable[[argparse.Namespace], dict[str, Any]] = compound.build_report,
) -> dict[str, Any]:
    if job.out_dir.exists() and (job.out_dir / "manifest.json").exists() and getattr(args, "skip_existing", True) and not getattr(args, "force", False):
        return manifest_row_for_job(job, status="skipped_existing")
    if job.out_dir.exists() and getattr(args, "force", False):
        for path in sorted(job.out_dir.rglob("*"), reverse=True):
            if path.is_file() or path.is_symlink():
                path.unlink()
            elif path.is_dir():
                path.rmdir()
    data_dir = job.out_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    prediction_slice = data_dir / "source_prediction_rows.csv"
    job.rows.to_csv(prediction_slice, index=False)
    renderer(build_renderer_args(job, args, prediction_slice))
    return manifest_row_for_job(job, status="generated")


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def render_landing_page(manifest_rows: list[dict[str, Any]]) -> str:
    generated = [row for row in manifest_rows if row.get("status") in {"generated", "skipped_existing"}]
    body_rows = []
    for row in generated:
        link = html.escape(str(row.get("index_html", "")))
        body_rows.append(
            "<tr>"
            f"<td><a href='{link}'>{html.escape(str(row.get('control_id', '')))}</a></td>"
            f"<td>{html.escape(str(row.get('lib_id', '')))}</td>"
            f"<td>{html.escape(str(row.get('rank1_name', '')))}</td>"
            f"<td class='mono'>{html.escape(str(row.get('rank1_formula', '')))}</td>"
            f"<td>{html.escape(str(row.get('predicted_pathway', '')))}</td>"
            f"<td>{html.escape(str(row.get('predicted_subpathway', '')))}</td>"
            f"<td>{html.escape(str(row.get('lipid_shorthand', '')))}</td>"
            f"<td>{html.escape(str(row.get('peptide_composition_label', '')))}</td>"
            "</tr>"
        )
    table_body = "\n".join(body_rows) or "<tr><td colspan='8'>No reports generated.</td></tr>"
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Public Compound Reports</title>
  <style>
    body {{ margin:0; font-family:Arial, sans-serif; color:#162235; background:#f6f8fb; }}
    main {{ max-width:1200px; margin:0 auto; padding:32px 24px; }}
    h1 {{ margin:0 0 8px; font-size:28px; }}
    p {{ color:#536172; margin:0 0 24px; }}
    table {{ width:100%; border-collapse:collapse; background:white; border:1px solid #dce3ec; }}
    th, td {{ padding:10px 12px; border-bottom:1px solid #e6ebf2; text-align:left; vertical-align:top; }}
    th {{ background:#edf2f7; font-size:12px; text-transform:uppercase; letter-spacing:.04em; color:#415166; }}
    a {{ color:#174ea6; text-decoration:none; }}
    a:hover {{ text-decoration:underline; }}
    .mono {{ font-family:ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size:12px; }}
  </style>
</head>
<body>
  <main>
    <h1>Public Compound Reports</h1>
    <p>{len(generated)} report bundles available.</p>
    <table>
      <thead><tr><th>Control</th><th>Lib</th><th>Rank-1 candidate</th><th>Formula</th><th>Pathway</th><th>Subpathway</th><th>Lipid</th><th>Peptide</th></tr></thead>
      <tbody>{table_body}</tbody>
    </table>
  </main>
</body>
</html>"""


def write_manifest(
    out_root: Path,
    *,
    manifest_rows: list[dict[str, Any]],
    failures: list[dict[str, Any]],
    input_files: list[PredictionCsv],
) -> dict[str, Any]:
    root = resolve(out_root, REPO_ROOT)
    root.mkdir(parents=True, exist_ok=True)
    write_csv(root / "manifest.csv", manifest_rows, MANIFEST_FIELDS)
    write_csv(root / "failures.csv", failures, FAILURE_FIELDS)
    summary = {
        "out_root": str(root),
        "input_files": [{"path": str(item.path), "lib_id": item.lib_id} for item in input_files],
        "n_reports": int(len(manifest_rows)),
        "n_generated": int(sum(1 for row in manifest_rows if row.get("status") == "generated")),
        "n_skipped_existing": int(sum(1 for row in manifest_rows if row.get("status") == "skipped_existing")),
        "n_failures": int(len(failures)),
        "reports": manifest_rows,
        "failures": failures,
    }
    (root / "manifest.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (root / "index.html").write_text(render_landing_page(manifest_rows), encoding="utf-8")
    return summary


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--equalization-root", type=Path, default=DEFAULT_EQUALIZATION_ROOT)
    parser.add_argument("--libs", nargs="+", default=[DEFAULT_LIBS])
    parser.add_argument("--input-name", default=DEFAULT_INPUT_NAME)
    parser.add_argument("--predictions", type=Path, action="append", default=[])
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument("--known-x", type=Path, default=compound.DEFAULT_KNOWN_X)
    parser.add_argument("--known-x-all", type=Path, default=compound.DEFAULT_KNOWN_X_ALL)
    parser.add_argument("--db", type=Path, default=compound.DEFAULT_DB)
    parser.add_argument("--candidate-pool", type=Path, default=None)
    parser.add_argument("--hypotheses", type=Path, default=None)
    parser.add_argument("--peak-csv", type=Path, default=None)
    parser.add_argument("--peak-root", type=Path, default=None)
    parser.add_argument("--control-id", action="append", default=[])
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--fail-fast", action="store_true")
    parser.add_argument("--progress-every", type=int, default=25)
    parser.add_argument("--spectrum-tolerance-da", type=float, default=1.0)
    parser.add_argument("--min-similarity", type=float, default=0.28)
    parser.add_argument("--max-neighbors", type=int, default=10)
    parser.add_argument("--candidate-tab-n", type=int, default=20)
    parser.add_argument("--consensus-top-n", type=int, default=20)
    parser.add_argument("--extended-known-render-limit", type=int, default=120)
    parser.add_argument("--extended-known-min-similarity", type=float, default=0.30)
    parser.add_argument("--extended-x-min-similarity", type=float, default=0.0)
    parser.add_argument("--extended-target-edge-min", type=float, default=0.32)
    parser.add_argument("--extended-pair-edge-min", type=float, default=0.42)
    parser.add_argument("--extended-max-edges-per-node", type=int, default=6)
    return parser.parse_args(argv)


def run(args: argparse.Namespace, *, renderer: Callable[[argparse.Namespace], dict[str, Any]] = compound.build_report) -> dict[str, Any]:
    if args.predictions:
        prediction_csvs = [PredictionCsv(path=resolve(path, REPO_ROOT)) for path in args.predictions]
        prediction_csvs = [PredictionCsv(path=item.path, lib_id=infer_lib_id(item.path)) for item in prediction_csvs]
    else:
        prediction_csvs = discover_prediction_csvs(args.equalization_root, parse_libs(args.libs), args.input_name)
    if not prediction_csvs:
        raise FileNotFoundError("No prediction CSVs found.")

    index = load_prediction_index(prediction_csvs)
    jobs = iter_report_jobs(
        index,
        out_root=args.out_root,
        control_ids=set(args.control_id or []),
        limit=args.limit,
    )
    manifest_rows: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    for idx, job in enumerate(jobs, start=1):
        try:
            manifest_rows.append(build_one_report(job, args, renderer=renderer))
        except Exception as exc:
            failures.append(
                {
                    "slug": job.slug,
                    "control_id": job.control_id,
                    "lib_id": job.lib_id,
                    "predictions_csv": str(job.predictions_csv),
                    "reason": str(exc),
                    "traceback": traceback.format_exc(),
                }
            )
            if args.fail_fast:
                raise
        if args.progress_every and idx % int(args.progress_every) == 0:
            print(f"[batch-public] processed {idx}/{len(jobs)} reports")

    summary = write_manifest(args.out_root, manifest_rows=manifest_rows, failures=failures, input_files=prediction_csvs)
    print(json.dumps({k: summary[k] for k in ("n_reports", "n_generated", "n_skipped_existing", "n_failures")}, sort_keys=True))
    return summary


def main() -> int:
    run(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
