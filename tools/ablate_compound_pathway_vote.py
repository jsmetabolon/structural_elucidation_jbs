#!/usr/bin/env python3
"""Ablate compound-level vote policy against known_X (no new fingerprints)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd

TOOLS = Path(__file__).resolve().parent
ROOT = TOOLS.parent
for path in (ROOT, TOOLS):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from structural_elucidation.pipeline.enrichment.pathway import (
    VOTE_POLICIES,
    assign_compound_pathways,
    load_subpathway_name_priors,
)
from structural_elucidation.paths import REPO_ROOT, resolve

import audit_pathway_predictions as audit


def _first_lib_wins(rows: list[pd.DataFrame]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame()
    stacked = pd.concat(rows, ignore_index=True)
    stacked = stacked.sort_values(["control_id", "lib", "rank"], kind="mergesort")
    return stacked.drop_duplicates(subset=["control_id"], keep="first")


def rank1_from_assigned(assigned: pd.DataFrame, lib: str) -> pd.DataFrame:
    id_col = "control_id" if "control_id" in assigned.columns else "x_compound_id"
    work = assigned.sort_values("rank", kind="mergesort").drop_duplicates(subset=[id_col], keep="first")
    out = pd.DataFrame(
        {
            "control_id": work[id_col].map(audit.clean),
            "lib": lib,
            "rank": work["rank"] if "rank" in work.columns else pd.NA,
            "predicted_smiles": work.get("predicted_smiles", pd.Series([""] * len(work))).map(audit.clean),
            "predicted_pathway": work.get("predicted_pathway", pd.Series([""] * len(work))).map(audit.clean),
            "predicted_subpathway": work.get("predicted_subpathway", pd.Series([""] * len(work))).map(audit.clean),
            "compound_predicted_pathway": work.get("compound_predicted_pathway", pd.Series([""] * len(work))).map(
                audit.clean
            ),
            "compound_predicted_subpathway": work.get(
                "compound_predicted_subpathway", pd.Series([""] * len(work))
            ).map(audit.clean),
            "compound_predicted_pathway_raw": work.get(
                "compound_predicted_pathway_raw", pd.Series([""] * len(work))
            ).map(audit.clean),
            "compound_pathway_source": work.get("compound_pathway_source", pd.Series([""] * len(work))).map(audit.clean),
            "compound_pathway_status": work.get("compound_pathway_status", pd.Series([""] * len(work))).map(audit.clean),
            "compound_pathway_name_boosted": work.get("compound_pathway_name_boosted"),
            "source": work.get("pathway_source", pd.Series([""] * len(work))).map(audit.assignment_source),
            "tier": work.get("pathway_confidence_tier", pd.Series([""] * len(work))).map(audit.clean),
            "status": work.get("compound_pathway_status", work.get("pathway_assignment_status")).map(audit.clean)
            if "compound_pathway_status" in work.columns
            else work.get("pathway_assignment_status", pd.Series([""] * len(work))).map(audit.clean),
            "vote_share": pd.to_numeric(work.get("compound_pathway_vote_share"), errors="coerce"),
            "best_similarity": pd.to_numeric(work.get("pathway_best_reference_similarity"), errors="coerce"),
        }
    )
    return out


def score_variant(pred: pd.DataFrame, gt: pd.DataFrame) -> dict[str, Any]:
    compound = pred.copy()
    compound["predicted_pathway"] = compound["compound_predicted_pathway"]
    compound["predicted_subpathway"] = compound["compound_predicted_subpathway"]
    summary, _ = audit.build_summary(compound, gt)
    acc = summary["accuracy"]
    return {
        "n": acc["n_evaluable_super"],
        "super": acc["super_accuracy"],
        "super_aa_peptide": acc.get("super_accuracy_aa_peptide"),
        "sub": acc["sub_accuracy"],
        "sub_scoreable": acc.get("sub_accuracy_scoreable"),
        "n_scoreable_sub": acc.get("n_scoreable_sub"),
    }


def read_weighting_experiment(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"status": "missing", "path": str(path)}
    summaries = sorted(path.glob("*summary*.json")) + sorted(path.glob("*comparison*.json"))
    if not summaries:
        return {
            "status": "incomplete",
            "path": str(path),
            "note": "No comparison/summary JSON yet; weighting grid not adopted.",
            "present_files": sorted(p.name for p in path.iterdir() if p.is_file())[:20],
        }
    latest = max(summaries, key=lambda item: item.stat().st_mtime)
    try:
        payload = json.loads(latest.read_text())
    except (OSError, ValueError) as exc:
        return {"status": "unreadable", "path": str(latest), "error": str(exc)}
    return {"status": "ready", "path": str(latest), "payload": payload}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--equalization-root", type=Path, default=Path("models/pipeline_equalization"))
    parser.add_argument("--pred-name", default="unified_rerank_all_top10_with_pathways.csv")
    parser.add_argument("--libs", default="209,305,400,402")
    parser.add_argument("--known-x", type=Path, default=Path("data/known_X.csv"))
    parser.add_argument("--reference-db", type=Path, default=Path("vendor/metabolite_synonyms.sqlite"))
    parser.add_argument("--out", type=Path, default=Path("output/_pathway_vote_ablation.json"))
    parser.add_argument(
        "--experiment-root",
        type=Path,
        default=Path("../spec2mol_merged_fullrun_metabolon_only/models/weighted_pathway_assignment_experiment_v1"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = resolve(args.equalization_root, REPO_ROOT)
    db = resolve(args.reference_db, REPO_ROOT)
    known_x = resolve(args.known_x, REPO_ROOT)
    libs = [token.strip() for token in args.libs.split(",") if token.strip()]
    gt = audit.load_ground_truth(known_x)
    priors = load_subpathway_name_priors(db) if db.exists() else {}

    loaded: dict[str, pd.DataFrame] = {}
    for lib in libs:
        path = root / f"rerank_lib{lib}" / args.pred_name
        if not path.exists():
            print(f"[ablate] missing {path}")
            continue
        loaded[lib] = pd.read_csv(path, low_memory=False)
        print(f"[ablate] loaded lib {lib}: {len(loaded[lib])} rows", flush=True)

    variants: list[dict[str, Any]] = []
    baseline: dict[str, Any] | None = None
    for policy in VOTE_POLICIES:
        for collapse in (False, True):
            for name_boost in (False, True):
                ranked_rows = []
                for lib, frame in loaded.items():
                    assigned = assign_compound_pathways(
                        frame,
                        policy=policy,
                        collapse_amino_acid_peptide=collapse,
                        name_priors=priors if name_boost else None,
                        apply_name_boost=name_boost,
                        name_prior_holdout_frac=0.25 if name_boost else 0.0,
                    )
                    ranked_rows.append(rank1_from_assigned(assigned, lib))
                pred = _first_lib_wins(ranked_rows)
                metrics = score_variant(pred, gt)
                label = {
                    "policy": policy,
                    "collapse_amino_acid_peptide": collapse,
                    "name_boost": name_boost,
                    **metrics,
                }
                variants.append(label)
                print(json.dumps(label), flush=True)
                if policy == "exact_first" and not collapse and not name_boost:
                    baseline = metrics

    def _key(row: dict[str, Any]) -> tuple[float, float, float, int, int, int]:
        return (
            float(row.get("super_aa_peptide") or 0),
            float(row.get("sub_scoreable") or 0),
            float(row.get("super") or 0),
            int(row.get("policy") == "all_eligible"),
            int(not row.get("collapse_amino_acid_peptide")),
            int(not row.get("name_boost")),
        )

    winner = max(variants, key=_key) if variants else {}
    recommended = dict(winner)
    min_gain = 0.01
    if baseline:
        super_gain = float(winner.get("super") or 0) - float(baseline.get("super") or 0)
        aa_gain = float(winner.get("super_aa_peptide") or 0) - float(baseline.get("super_aa_peptide") or 0)
        sub_gain = float(winner.get("sub_scoreable") or 0) - float(baseline.get("sub_scoreable") or 0)
        helped = super_gain >= min_gain or aa_gain >= min_gain or sub_gain >= min_gain
        hurt = super_gain <= -min_gain or sub_gain <= -min_gain
        if not (helped and not hurt):
            recommended = {
                "policy": "exact_first",
                "collapse_amino_acid_peptide": False,
                "name_boost": False,
                **baseline,
            }

    experiment = read_weighting_experiment(resolve(args.experiment_root, REPO_ROOT))
    payload = {
        "baseline": baseline,
        "variants": variants,
        "winner": winner,
        "recommended": recommended,
        "weighting_experiment": {k: v for k, v in experiment.items() if k != "payload"} | (
            {"experiment_keys": list((experiment.get("payload") or {}).keys())[:30]}
            if experiment.get("status") == "ready"
            else {}
        ),
        "weighting_change_adopted": None,
        "note": (
            "No fingerprint recompute. Weighting-grid change adopted only if the experiment "
            "file is complete and beats this ablation on paired known_X super/sub."
        ),
    }
    if experiment.get("status") != "ready":
        payload["weighting_change_adopted"] = "none_experiment_incomplete"
    out = resolve(args.out, REPO_ROOT)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {"wrote": str(out), "winner": winner, "recommended": recommended, "experiment": experiment.get("status")},
            indent=2,
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
