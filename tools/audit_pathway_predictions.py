#!/usr/bin/env python3
"""Audit pathway/subpathway prediction results against known_X ground truth.

Assesses the neighbor-consensus pathway predictions that feed the compound
reports (``predicted_pathway`` / ``predicted_subpathway`` +
``pathway_confidence_tier`` / ``pathway_assignment_status`` /
``pathway_vote_share`` / ``pathway_best_reference_similarity``) carried on the
enriched rerank rows ``rerank_lib*/unified_rerank_all_top10_with_pathways.csv``.

Ground truth: ``data/known_X.csv`` columns ``SuperPathway`` / ``Draft Subpathway``
(same Metabolon vocabulary as the predictor), joined on ``control_id`` == ``Unknown``.

Reports coverage + confidence, super/sub accuracy (overall and by tier / status /
similarity bin), calibration, a confusion matrix, and predicted-vs-truth
distribution (bias). Emits ``audit_pathway_accuracy.csv``,
``audit_pathway_summary.json``, and ``audit_pathway_summary.md``.

Usage::

    python tools/audit_pathway_predictions.py --equalization-root output/pipeline_equalization
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

DEFAULT_LIBS = ("209", "305", "400", "402")
NULLS = {"", "nan", "none", "null", "<na>"}
# known_X data-entry sentinels that are not real pathways (exclude from accuracy).
EXCLUDE_TRUTH = {"unknown", "ask dts", "can add this", "do not include"}
SIM_BINS = [(0.0, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 1.0001)]


def norm(value: Any) -> str:
    """Normalize a label for comparison: lowercase, collapse whitespace."""

    text = re.sub(r"\s+", " ", str("" if value is None else value).strip()).lower()
    return "" if text in NULLS else text


def clean(value: Any) -> str:
    """Stripped display string, mapping null sentinels to empty."""

    text = str("" if value is None else value).strip()
    return "" if text.lower() in NULLS else text


def load_rank1_predictions(root: Path, libs: list[str], pred_name: str) -> pd.DataFrame:
    """Rank-1 prediction row per control, deduped across libs (first-lib-wins)."""

    seen: set[str] = set()
    rows: list[dict[str, Any]] = []
    for lib in libs:
        path = root / f"rerank_lib{lib}" / pred_name
        if not path.exists():
            print(f"[pathway-audit] WARN missing {path}")
            continue
        df = pd.read_csv(path, low_memory=False)
        id_col = "control_id" if "control_id" in df.columns else "x_compound_id"
        df = df.sort_values("rank", kind="mergesort")
        for _, r in df.drop_duplicates(subset=[id_col], keep="first").iterrows():
            cid = clean(r[id_col])
            if not cid or cid in seen:
                continue
            seen.add(cid)
            rows.append(
                {
                    "control_id": cid,
                    "lib": lib,
                    "predicted_pathway": clean(r.get("predicted_pathway")),
                    "predicted_subpathway": clean(r.get("predicted_subpathway")),
                    "tier": clean(r.get("pathway_confidence_tier")),
                    "status": clean(r.get("pathway_assignment_status")),
                    "vote_share": pd.to_numeric(r.get("pathway_vote_share"), errors="coerce"),
                    "best_similarity": pd.to_numeric(r.get("pathway_best_reference_similarity"), errors="coerce"),
                }
            )
    return pd.DataFrame(rows)


def load_ground_truth(known_x_path: Path) -> pd.DataFrame:
    """control_id -> (true super, true sub) from known_X, dropping empty/Unknown."""

    kx = pd.read_csv(known_x_path, low_memory=False)
    out: dict[str, dict[str, str]] = {}
    for _, r in kx.iterrows():
        cid = clean(r.get("Unknown"))
        sup = clean(r.get("SuperPathway"))
        if not cid or cid in out:
            continue
        out[cid] = {"true_super": sup, "true_sub": clean(r.get("Draft Subpathway"))}
    frame = pd.DataFrame([{"control_id": k, **v} for k, v in out.items()])
    return frame


def similarity_bin(value: float) -> str:
    """Label the similarity bin for a best-reference-similarity value."""

    if value is None or pd.isna(value):
        return "na"
    for lo, hi in SIM_BINS:
        if lo <= value < hi:
            return f"[{lo:.1f},{hi if hi <= 1 else 1.0:.1f})"
    return "na"


def accuracy_by(df: pd.DataFrame, key: str, correct_col: str) -> dict[str, dict[str, Any]]:
    """Accuracy of ``correct_col`` grouped by ``key``."""

    out: dict[str, dict[str, Any]] = {}
    for value, grp in df.groupby(key):
        n = len(grp)
        acc = float(grp[correct_col].mean()) if n else 0.0
        out[str(value)] = {"n": int(n), "accuracy": round(acc, 4)}
    return out


def build_summary(pred: pd.DataFrame, gt: pd.DataFrame) -> tuple[dict[str, Any], pd.DataFrame]:
    """Compute coverage, accuracy, calibration, confusion, and bias metrics."""

    n_pred = len(pred)
    coverage = {
        "n_controls_predicted": int(n_pred),
        "status_counts": {str(k): int(v) for k, v in pred["status"].value_counts().items()},
        "tier_counts": {str(k): int(v) for k, v in pred["tier"].value_counts().items()},
        "best_similarity_pctiles": {
            p: round(float(pred["best_similarity"].quantile(q)), 3)
            for p, q in (("p10", 0.1), ("p50", 0.5), ("p90", 0.9))
        },
        "vote_share_pctiles": {
            p: round(float(pred["vote_share"].quantile(q)), 3)
            for p, q in (("p10", 0.1), ("p50", 0.5), ("p90", 0.9))
        },
    }

    ev = pred.merge(gt, on="control_id", how="inner")
    true_norm = ev["true_super"].map(norm)
    ev = ev[true_norm.ne("") & ~true_norm.isin(EXCLUDE_TRUTH)].copy()
    ev["correct_super"] = ev.apply(lambda r: norm(r["predicted_pathway"]) == norm(r["true_super"]), axis=1)
    ev["correct_sub"] = ev.apply(
        lambda r: norm(r["predicted_subpathway"]) == norm(r["true_sub"]) and norm(r["true_sub"]) != "",
        axis=1,
    )
    ev["sim_bin"] = ev["best_similarity"].map(similarity_bin)

    n_ev = len(ev)
    sub_evaluable = ev[ev["true_sub"].map(norm).ne("")]
    accuracy = {
        "n_evaluable_super": int(n_ev),
        "super_accuracy": round(float(ev["correct_super"].mean()), 4) if n_ev else 0.0,
        "n_evaluable_sub": int(len(sub_evaluable)),
        "sub_accuracy": round(float(sub_evaluable["correct_sub"].mean()), 4) if len(sub_evaluable) else 0.0,
        "sub_accuracy_given_correct_super": (
            round(float(sub_evaluable.loc[sub_evaluable["correct_super"], "correct_sub"].mean()), 4)
            if len(sub_evaluable[sub_evaluable["correct_super"]])
            else 0.0
        ),
        "super_accuracy_by_tier": accuracy_by(ev, "tier", "correct_super"),
        "super_accuracy_by_status": accuracy_by(ev, "status", "correct_super"),
        "super_accuracy_by_similarity_bin": accuracy_by(ev, "sim_bin", "correct_super"),
    }

    # Vocabulary overlap (catch spelling drift that would understate accuracy).
    pred_vocab = {norm(v) for v in ev["predicted_pathway"] if norm(v)}
    true_vocab = {norm(v) for v in ev["true_super"] if norm(v)}
    vocab = {
        "predicted_supers_not_in_truth_vocab": sorted(pred_vocab - true_vocab)[:20],
        "truth_supers_not_predicted": sorted(true_vocab - pred_vocab)[:20],
    }

    # Confusion for the most common true super-pathways.
    confusion: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for _, r in ev.iterrows():
        confusion[clean(r["true_super"]) or "(blank)"][clean(r["predicted_pathway"]) or "(blank)"] += 1
    top_true = [t for t, _ in Counter(ev["true_super"].map(clean)).most_common(8)]
    confusion_top = {
        t: dict(sorted(confusion[t].items(), key=lambda kv: kv[1], reverse=True)[:5]) for t in top_true if t
    }

    # Bias: predicted vs truth super distribution on the evaluable set.
    pred_counts = Counter(ev["predicted_pathway"].map(clean)).most_common(12)
    true_counts = Counter(ev["true_super"].map(clean)).most_common(12)
    bias = {
        "predicted_super_distribution": {str(k): int(v) for k, v in pred_counts},
        "truth_super_distribution": {str(k): int(v) for k, v in true_counts},
    }

    summary = {
        "coverage": coverage,
        "accuracy": accuracy,
        "vocabulary": vocab,
        "confusion_top_true_super": confusion_top,
        "bias": bias,
    }
    return summary, ev


def leave_one_out_validation(db_path: Path, k: int, floor: float, power: float) -> dict[str, Any]:
    """Intrinsic accuracy of the neighbor-consensus mechanism via leave-one-out.

    For every metabolon-labeled reference metabolite, hide it and predict its
    super/sub pathway from a similarity-weighted vote of its k nearest neighbors
    (Tanimoto, Morgan r2). This measures the ceiling/behavior of the core
    neighbor-consensus itself (not the full production predictor with MW priors /
    SMARTS overrides - that is audited by the local build_v15 tool).
    """

    import sqlite3

    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.DataStructs import BulkTanimotoSimilarity

    con = sqlite3.connect(str(db_path))
    ref = pd.read_sql_query(
        "SELECT smiles, metabolon_super, metabolon_sub FROM metabolites "
        "WHERE TRIM(COALESCE(smiles,'')) <> '' AND TRIM(COALESCE(metabolon_super,'')) <> ''",
        con,
    )
    con.close()

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps: list[Any] = []
    supers: list[str] = []
    subs: list[str] = []
    for _, r in ref.iterrows():
        mol = Chem.MolFromSmiles(str(r["smiles"]).split("|")[0].strip())
        if mol is None:
            continue
        fps.append(gen.GetFingerprint(mol))
        supers.append(clean(r["metabolon_super"]))
        subs.append(clean(r["metabolon_sub"]))

    n = len(fps)
    correct_super = 0
    correct_sub = 0
    sub_evaluable = 0
    evaluated = 0
    by_bin: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for i in range(n):
        sims = BulkTanimotoSimilarity(fps[i], fps)
        cand = sorted(((sims[j], j) for j in range(n) if j != i and sims[j] >= floor), reverse=True)[:k]
        if not cand:
            continue
        evaluated += 1
        votes: dict[str, float] = defaultdict(float)
        for sim, j in cand:
            votes[supers[j]] += sim ** power
        pred_super = max(votes, key=votes.get)
        ok = norm(pred_super) == norm(supers[i])
        correct_super += int(ok)
        sub_votes: dict[str, float] = defaultdict(float)
        for sim, j in cand:
            if supers[j] == pred_super and subs[j]:
                sub_votes[subs[j]] += sim ** power
        pred_sub = max(sub_votes, key=sub_votes.get) if sub_votes else ""
        if subs[i]:
            sub_evaluable += 1
            correct_sub += int(norm(pred_sub) == norm(subs[i]))
        b = similarity_bin(cand[0][0])
        by_bin[b][0] += int(ok)
        by_bin[b][1] += 1

    return {
        "n_reference_labeled": int(n),
        "n_evaluated": int(evaluated),
        "k": k,
        "similarity_floor": floor,
        "similarity_power": power,
        "super_accuracy": round(correct_super / evaluated, 4) if evaluated else 0.0,
        "sub_accuracy": round(correct_sub / sub_evaluable, 4) if sub_evaluable else 0.0,
        "super_accuracy_by_similarity_bin": {
            b: {"n": c[1], "accuracy": round(c[0] / c[1], 4)} for b, c in sorted(by_bin.items())
        },
    }


def write_markdown(summary: dict[str, Any], out_path: Path) -> None:
    """Write a short human-readable audit summary."""

    acc = summary["accuracy"]
    cov = summary["coverage"]
    lines = [
        "# Pathway/subpathway prediction audit",
        "",
        f"- Controls predicted: **{cov['n_controls_predicted']}**",
        f"- Status: {cov['status_counts']}",
        f"- Tier: {cov['tier_counts']}",
        f"- best_reference_similarity p10/p50/p90: {cov['best_similarity_pctiles']}",
        "",
        "## Accuracy vs known_X ground truth",
        f"- Evaluable (super): **{acc['n_evaluable_super']}** controls",
        f"- **Super-pathway accuracy: {acc['super_accuracy']:.1%}**",
        f"- Sub-pathway accuracy: {acc['sub_accuracy']:.1%} "
        f"(given correct super: {acc['sub_accuracy_given_correct_super']:.1%}; n={acc['n_evaluable_sub']})",
        "",
        "### Super accuracy by confidence tier",
        "| tier | n | accuracy |",
        "|------|---|----------|",
    ]
    for tier in ("high", "medium", "low"):
        row = acc["super_accuracy_by_tier"].get(tier)
        if row:
            lines.append(f"| {tier} | {row['n']} | {row['accuracy']:.1%} |")
    lines += ["", "### Super accuracy by best-similarity bin", "| bin | n | accuracy |", "|-----|---|----------|"]
    for b, row in sorted(acc["super_accuracy_by_similarity_bin"].items()):
        lines.append(f"| {b} | {row['n']} | {row['accuracy']:.1%} |")
    lines += ["", "### Top true-super confusions (true -> predicted counts)"]
    for true_super, preds in summary["confusion_top_true_super"].items():
        lines.append(f"- **{true_super}**: {preds}")
    if summary["vocabulary"]["predicted_supers_not_in_truth_vocab"]:
        lines += ["", "### Predicted supers not in truth vocabulary (possible drift)",
                  f"- {summary['vocabulary']['predicted_supers_not_in_truth_vocab']}"]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """CLI entrypoint."""

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--equalization-root", default="output/pipeline_equalization")
    parser.add_argument("--pred-name", default="unified_rerank_all_top10_with_pathways.csv")
    parser.add_argument("--libs", default=",".join(DEFAULT_LIBS))
    parser.add_argument("--known-x", default="data/known_X.csv")
    parser.add_argument("--out-dir", default="output/_pathway_audit")
    parser.add_argument("--loo", action="store_true", help="also run leave-one-out validation over the reference DB")
    parser.add_argument("--db", default="vendor/metabolite_synonyms.sqlite")
    parser.add_argument("--loo-k", type=int, default=10)
    parser.add_argument("--loo-floor", type=float, default=0.3)
    parser.add_argument("--loo-power", type=float, default=2.0)
    args = parser.parse_args()

    libs = [x.strip() for x in args.libs.split(",") if x.strip()]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pred = load_rank1_predictions(Path(args.equalization_root), libs, args.pred_name)
    gt = load_ground_truth(Path(args.known_x))
    print(f"[pathway-audit] {len(pred)} predicted controls; {len(gt)} known_X ground-truth controls")

    summary, ev = build_summary(pred, gt)

    ev_out = ev[
        [
            "control_id", "lib", "predicted_pathway", "true_super", "correct_super",
            "predicted_subpathway", "true_sub", "correct_sub", "tier", "status", "best_similarity",
        ]
    ].sort_values(["correct_super", "tier"])
    ev_out.to_csv(out_dir / "audit_pathway_accuracy.csv", index=False)
    (out_dir / "audit_pathway_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_markdown(summary, out_dir / "audit_pathway_summary.md")

    a = summary["accuracy"]
    print(f"[pathway-audit] super accuracy: {a['super_accuracy']:.1%} on {a['n_evaluable_super']} evaluable")
    print(
        f"[pathway-audit] sub accuracy:   {a['sub_accuracy']:.1%} "
        f"(given correct super {a['sub_accuracy_given_correct_super']:.1%})"
    )
    print(f"[pathway-audit] by tier: {a['super_accuracy_by_tier']}")

    if args.loo:
        db_path = Path(args.db)
        if db_path.exists():
            print("[pathway-audit] running leave-one-out validation over reference DB ...")
            loo = leave_one_out_validation(db_path, args.loo_k, args.loo_floor, args.loo_power)
            (out_dir / "audit_pathway_loo.json").write_text(json.dumps(loo, indent=2), encoding="utf-8")
            print(
                f"[pathway-audit] LOO super accuracy: {loo['super_accuracy']:.1%} "
                f"(sub {loo['sub_accuracy']:.1%}) over {loo['n_evaluated']} reference metabolites"
            )
            print(f"[pathway-audit] LOO by similarity bin: {loo['super_accuracy_by_similarity_bin']}")
        else:
            print(f"[pathway-audit] --loo requested but DB not found at {db_path}")
    print(f"[pathway-audit] outputs in {out_dir}")


if __name__ == "__main__":
    main()
