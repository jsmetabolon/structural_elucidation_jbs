#!/usr/bin/env python3
"""Top-K candidate ensemble for pathway assignment (T3).

The production predictor takes the pathway of the rank-1 candidate structure
only. The rerank emits ten scored candidates per control, and pathway
enrichment already labels every one of them — so when rank-1 is the wrong
structure but the pool is topically coherent, an ensemble over the top-K
candidate pathways can recover the right super-pathway. This tool builds those
ensemble predictions for scoring with ``tools/pathway_variant_benchmark.py``.

Strategies
----------
- ``rank1``: rank-1 candidate's pathway (reference behavior).
- ``uniform``: each of the top-K candidates votes its pathway with weight 1.
- ``rank_decay``: weight 1/rank.
- ``score_weighted``: weight from the rerank score column (shift-normalized
  within each control so negative scores cannot flip votes).
- ``agreement_gate``: rank-1 pathway, but abstain when fewer than
  ``--agreement-threshold`` of the top-K vote weight backs it (feeds the
  abstention policy work in T5).

Sub-pathways are voted only among candidates that backed the winning super.

Gated on the fetched rerank enrichment outputs
``models/pipeline_equalization/rerank_lib{lib}/unified_rerank_all_top10_with_pathways.csv``
(controls deduped across libraries first-lib-wins, matching the audit).

Usage::

    python tools/pathway_topk_ensemble.py run --out-prefix \
        docs/paper/data/pathway/topk_ensemble
    # then, per strategy:
    python tools/pathway_variant_benchmark.py score \
        --pred docs/paper/data/pathway/topk_ensemble_uniform.csv \
        --pred-sub-col predicted_sub --label topk_uniform
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_LIBS = (209, 305, 400, 402)
STRATEGIES = ("rank1", "uniform", "rank_decay", "score_weighted", "agreement_gate")
SCORE_COLUMN_CANDIDATES = (
    "rerank_score",
    "unified_score",
    "final_score",
    "gbm_score",
    "score",
    "verifier_score",
)


def clean(value: Any) -> str:
    text = str("" if value is None else value).strip()
    return "" if text.lower() in {"nan", "none", "null", "<na>"} else text


def pick_score_column(frame: pd.DataFrame, requested: str | None) -> str | None:
    if requested:
        if requested not in frame.columns:
            raise ValueError(f"score column {requested!r} not in input columns")
        return requested
    for candidate in SCORE_COLUMN_CANDIDATES:
        if candidate in frame.columns:
            return candidate
    return None


def candidate_weights(group: pd.DataFrame, strategy: str, score_col: str | None) -> list[float]:
    """Per-candidate vote weights for one control's top-K rows (rank order)."""

    n = len(group)
    if strategy in ("rank1", "agreement_gate", "uniform"):
        return [1.0] * n
    if strategy == "rank_decay":
        ranks = pd.to_numeric(group["rank"], errors="coerce").fillna(n).tolist()
        return [1.0 / max(1.0, float(r)) for r in ranks]
    if strategy == "score_weighted":
        if score_col is None:
            raise ValueError("score_weighted requires a score column; none found or given")
        scores = pd.to_numeric(group[score_col], errors="coerce")
        finite = scores.dropna()
        floor = float(finite.min()) if len(finite) else 0.0
        return [max(0.0, (float(s) if pd.notna(s) else floor) - floor) + 1e-6 for s in scores]
    raise ValueError(f"unknown strategy: {strategy}")


def ensemble_one(
    group: pd.DataFrame,
    strategy: str,
    *,
    k: int,
    score_col: str | None,
    agreement_threshold: float,
) -> dict[str, Any]:
    """Ensemble prediction for one control from its ranked candidate rows."""

    ranked = group.sort_values("rank", kind="mergesort").head(k)
    labeled = ranked[ranked["predicted_pathway"].map(clean).ne("")]
    result = {
        "predicted_super": "",
        "predicted_sub": "",
        "vote_share": 0.0,
        "n_candidates": int(len(ranked)),
        "n_labeled": int(len(labeled)),
    }
    if labeled.empty:
        return result

    weights = candidate_weights(labeled, strategy, score_col)
    super_tally: dict[str, float] = defaultdict(float)
    for weight, (_, row) in zip(weights, labeled.iterrows()):
        super_tally[clean(row["predicted_pathway"])] += weight
    total = sum(super_tally.values())
    rank1_super = clean(labeled.iloc[0]["predicted_pathway"])

    if strategy == "rank1":
        winner = rank1_super
    elif strategy == "agreement_gate":
        share = super_tally.get(rank1_super, 0.0) / total if total else 0.0
        if share < agreement_threshold:
            return result  # abstain
        winner = rank1_super
    else:
        winner = max(sorted(super_tally), key=lambda name: super_tally[name])

    sub_tally: dict[str, float] = defaultdict(float)
    for weight, (_, row) in zip(weights, labeled.iterrows()):
        if clean(row["predicted_pathway"]) == winner and clean(row.get("predicted_subpathway")):
            sub_tally[clean(row["predicted_subpathway"])] += weight
    winning_sub = max(sorted(sub_tally), key=lambda name: sub_tally[name]) if sub_tally else ""

    result.update(
        {
            "predicted_super": winner,
            "predicted_sub": winning_sub,
            "vote_share": round(super_tally.get(winner, 0.0) / total, 4) if total else 0.0,
        }
    )
    return result


def build_ensembles(
    frame: pd.DataFrame,
    strategies: tuple[str, ...] = STRATEGIES,
    *,
    k: int = 10,
    score_col: str | None = None,
    agreement_threshold: float = 0.6,
) -> dict[str, pd.DataFrame]:
    """Per-strategy prediction frames from labeled rerank rows (all controls)."""

    id_col = "control_id" if "control_id" in frame.columns else "x_compound_id"
    out: dict[str, list[dict[str, Any]]] = {s: [] for s in strategies}
    for control_id, group in frame.groupby(id_col, sort=False):
        for strategy in strategies:
            row = ensemble_one(
                group, strategy, k=k, score_col=score_col, agreement_threshold=agreement_threshold
            )
            out[strategy].append({"control_id": clean(control_id), **row})
    return {s: pd.DataFrame(rows) for s, rows in out.items()}


def load_rerank_rows(root: Path, libs: tuple[int, ...]) -> pd.DataFrame:
    """Labeled rerank rows, controls deduped across libs first-lib-wins."""

    frames = []
    seen: set[str] = set()
    for lib in libs:
        path = root / f"rerank_lib{lib}" / "unified_rerank_all_top10_with_pathways.csv"
        if not path.exists():
            print(f"[topk] WARN missing {path}")
            continue
        df = pd.read_csv(path, low_memory=False)
        id_col = "control_id" if "control_id" in df.columns else "x_compound_id"
        df = df.rename(columns={id_col: "control_id"})
        lib_controls = set(df["control_id"].map(clean)) - seen
        frames.append(df[df["control_id"].map(clean).isin(lib_controls)])
        seen |= lib_controls
    if not frames:
        raise FileNotFoundError(f"no *_with_pathways.csv under {root}")
    return pd.concat(frames, ignore_index=True)


def run(args: argparse.Namespace) -> int:
    root = args.equalization_root if args.equalization_root.is_absolute() else REPO_ROOT / args.equalization_root
    if not root.exists():
        print(
            f"[topk] SKIPPED: rerank enrichment outputs not present under {root}. "
            "Fetch models/pipeline_equalization and re-run."
        )
        return 0
    frame = load_rerank_rows(root, tuple(args.libs))
    score_col = pick_score_column(frame, args.score_col)
    if score_col is None:
        print("[topk] WARN no score column found; skipping score_weighted strategy")
        strategies = tuple(s for s in STRATEGIES if s != "score_weighted")
    else:
        strategies = STRATEGIES
    ensembles = build_ensembles(
        frame,
        strategies,
        k=args.k,
        score_col=score_col,
        agreement_threshold=args.agreement_threshold,
    )
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for strategy, pred in ensembles.items():
        out_path = args.out_prefix.parent / f"{args.out_prefix.name}_{strategy}.csv"
        pred.to_csv(out_path, index=False)
        n_called = int(pred["predicted_super"].map(clean).ne("").sum())
        print(f"[topk] {strategy}: {n_called}/{len(pred)} called -> {out_path}")
    print(
        "[topk] score each strategy with: python tools/pathway_variant_benchmark.py score "
        f"--pred {args.out_prefix}_<strategy>.csv --pred-sub-col predicted_sub --label topk_<strategy>"
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    pr = sub.add_parser("run", help="build ensemble predictions from fetched rerank outputs")
    pr.add_argument("--equalization-root", type=Path, default=Path("models/pipeline_equalization"))
    pr.add_argument("--libs", nargs="+", type=int, default=list(DEFAULT_LIBS))
    pr.add_argument("--k", type=int, default=10)
    pr.add_argument("--score-col", default=None)
    pr.add_argument("--agreement-threshold", type=float, default=0.6)
    pr.add_argument("--out-prefix", type=Path, default=REPO_ROOT / "docs/paper/data/pathway/topk_ensemble")
    args = parser.parse_args()
    if args.command == "run":
        return run(args)
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
