#!/usr/bin/env python3
"""A/B the evidence_confidence backfill through the reranker on one pool.

Trains ``pipeline.rerank`` twice on the same candidate pool -- once with the
backfill enabled and once with ``--no-evidence-backfill`` -- and reports what
moved. Both arms read the *unpopulated* pool so the only difference is the
feature values the model sees; do not feed this an already-backfilled pool.

``rerank.run`` reads the two truth CSVs from the repository root, where this
repository keeps them under ``data/``. Symlinks are created for the run and
removed afterwards unless they already existed.

Usage:
    python3 tools/ab_evidence_backfill_rerank.py \\
      --python ../spec2mol_merged_fullrun_metabolon_only/.venv/bin/python \\
      --pool ../spec2mol_merged_fullrun_metabolon_only/models/\\
unified_candidate_pool_v14_lib209_negative/unified_candidates_per_control.csv \\
      --output-root /tmp/evidence_backfill_ab/lib209

Add ``--compare-only`` to re-read finished arms without retraining.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from structural_elucidation.utils.io import dataframe_to_markdown  # noqa: E402

# rerank.run hardcodes these at the repository root; this repo stores them in data/.
TRUTH_FILES: dict[str, str] = {
    "known_X.csv": "data/known_X.csv",
    "known_X_all_available_structures_simple.csv": "data/known_X_all_available_structures_simple.csv",
}

ARMS: dict[str, list[str]] = {
    "with_backfill": ["--evidence-backfill"],
    "no_backfill": ["--no-evidence-backfill"],
}

SUMMARY_FIELDS: tuple[str, ...] = (
    "n_candidates",
    "n_positive_labels",
    "selected_model",
    "final_rank1_ge_0p7",
    "final_top10_ge_0p7",
    "oof_rank1_ge_0p7",
    "oof_top10_ge_0p7",
    "g3_pass",
    "g3_rank1_delta_abs",
)


def link_truth_files() -> list[Path]:
    """Symlink the truth CSVs into the repository root; return what was created."""
    created: list[Path] = []
    for name, target in TRUTH_FILES.items():
        link = REPO_ROOT / name
        source = REPO_ROOT / target
        if link.exists() or link.is_symlink():
            continue
        if not source.exists():
            raise SystemExit(f"missing {source}; cannot satisfy {name} at the repo root")
        link.symlink_to(source)
        created.append(link)
    return created


def train_arm(python: Path, pool: Path, output_root: Path, arm: str, args: argparse.Namespace) -> float:
    """Run one rerank arm; return wall-clock seconds."""
    command = [
        str(python),
        "-m",
        "structural_elucidation.pipeline.rerank",
        "--unified-pool",
        str(pool),
        "--non-s0-pool",
        str(args.non_s0_pool) if args.non_s0_pool else str(pool.with_name("unified_candidates_per_control_non_s0.csv")),
        "--output-root",
        str(output_root / arm),
        "--lib-id",
        str(args.lib_id),
        "--ion-mode",
        args.ion_mode,
        "--top-k",
        str(args.top_k),
        *ARMS[arm],
    ]
    print(f"[ab] {arm}: {' '.join(command)}", flush=True)
    started = time.monotonic()
    subprocess.run(command, cwd=REPO_ROOT, check=True)
    elapsed = time.monotonic() - started
    print(f"[ab] {arm} finished in {elapsed / 60:.1f} min", flush=True)
    return elapsed


def read_summary(path: Path) -> dict[str, Any]:
    file = path / "unified_rerank_summary.json"
    if not file.exists():
        raise SystemExit(f"missing {file}; run the arm first")
    return json.loads(file.read_text())


def summary_table(summaries: dict[str, dict[str, Any]]) -> pd.DataFrame:
    rows = []
    for field in SUMMARY_FIELDS:
        on, off = summaries["with_backfill"].get(field), summaries["no_backfill"].get(field)
        delta = on - off if isinstance(on, (int, float)) and isinstance(off, (int, float)) else ""
        rows.append({"field": field, "with_backfill": on, "no_backfill": off, "delta": delta})
    return pd.DataFrame(rows)


def metrics_table(output_root: Path) -> pd.DataFrame:
    frames = []
    for arm in ARMS:
        file = output_root / arm / "unified_rerank_metrics.csv"
        if not file.exists():
            continue
        frame = pd.read_csv(file)
        frame.insert(0, "arm", arm)
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    keep = ["arm", "policy", "metric_scope", "n_controls_covered", "n_tanimoto_ge_0p7", "median_best_tanimoto"]
    out = pd.concat(frames, ignore_index=True)
    return out[[column for column in keep if column in out.columns]]


def per_control_movers(output_root: Path, threshold: float = 0.7) -> pd.DataFrame:
    """Controls whose rank-1 crossed the match threshold between the two arms."""
    frames = {}
    for arm in ARMS:
        file = output_root / arm / "unified_rerank_per_control.csv"
        if not file.exists():
            return pd.DataFrame()
        frame = pd.read_csv(file)
        if not {"control_id", "rank1_best_tanimoto"}.issubset(frame.columns):
            return pd.DataFrame()
        indexed = frame.drop_duplicates("control_id").set_index("control_id")
        frames[arm] = pd.to_numeric(indexed["rank1_best_tanimoto"], errors="coerce")
    merged = pd.DataFrame(frames).fillna(0.0)
    gained = merged["with_backfill"].ge(threshold) & merged["no_backfill"].lt(threshold)
    lost = merged["no_backfill"].ge(threshold) & merged["with_backfill"].lt(threshold)
    moved = merged[gained | lost].copy()
    moved["direction"] = gained[moved.index].map({True: "gained", False: "lost"})
    return moved.reset_index()


def evidence_importance(output_root: Path) -> pd.DataFrame:
    rows = []
    for arm in ARMS:
        file = output_root / arm / "feature_importances.csv"
        if not file.exists():
            continue
        frame = pd.read_csv(file)
        match = frame[frame["feature"].astype(str).eq("evidence_confidence")]
        rows.append(
            {
                "arm": arm,
                "evidence_confidence_importance": float(match["importance"].iloc[0]) if not match.empty else None,
                "rank_among_features": int(match.index[0]) + 1 if not match.empty else None,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--python", type=Path, required=True, help="interpreter with sklearn + rdkit")
    parser.add_argument("--pool", type=Path, required=True, help="pool CSV, before any backfill")
    parser.add_argument("--non-s0-pool", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--lib-id", type=int, default=209)
    parser.add_argument("--ion-mode", default="negative")
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--compare-only", action="store_true", help="read finished arms without retraining")
    args = parser.parse_args()

    pool = args.pool.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    if not args.compare_only:
        evidence = pd.read_csv(pool, low_memory=False, usecols=lambda column: column == "evidence_confidence")
        if "evidence_confidence" in evidence.columns:
            zeros = int(pd.to_numeric(evidence["evidence_confidence"], errors="coerce").fillna(0.0).le(0.0).sum())
            print(f"[ab] pool {pool.name}: {zeros:,} rows with evidence_confidence <= 0", flush=True)
            if not zeros:
                raise SystemExit("this pool has no zero-evidence rows; the two arms would be identical")
        else:
            print(f"[ab] pool {pool.name}: no evidence_confidence column at all", flush=True)
        created = link_truth_files()
        try:
            for arm in ARMS:
                train_arm(args.python.expanduser().resolve(), pool, output_root, arm, args)
        finally:
            for link in created:
                link.unlink(missing_ok=True)

    summaries = {arm: read_summary(output_root / arm) for arm in ARMS}
    sections = [
        "# evidence_confidence backfill A/B",
        "",
        f"Pool: `{pool}`",
        f"Scope: lib {args.lib_id} {args.ion_mode}",
        "",
        "## Recovery",
        dataframe_to_markdown(summary_table(summaries)),
        "",
        "## Metrics",
        dataframe_to_markdown(metrics_table(output_root)),
        "",
        "## Controls that crossed Tanimoto 0.70 at rank 1",
        dataframe_to_markdown(per_control_movers(output_root)),
        "",
        "## evidence_confidence as a model feature",
        dataframe_to_markdown(evidence_importance(output_root)),
        "",
    ]
    report = "\n".join(sections)
    (output_root / "ab_summary.md").write_text(report)
    print("\n" + report)
    print(f"[ab] written to {output_root / 'ab_summary.md'}")


if __name__ == "__main__":
    sys.exit(main())
