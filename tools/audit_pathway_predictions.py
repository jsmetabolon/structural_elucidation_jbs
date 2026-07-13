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
import csv
import json
import math
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


def wilson(correct: int, n: int, z: float = 1.959963984540054) -> dict[str, Any]:
    """Binomial accuracy with a Wilson 95% confidence interval."""

    if n <= 0:
        return {"n": 0, "correct": 0, "accuracy": None, "ci95": [None, None]}
    p = correct / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denom
    half = z * math.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n) / denom
    return {
        "n": int(n),
        "correct": int(correct),
        "accuracy": round(p, 4),
        "ci95": [round(max(0.0, center - half), 4), round(min(1.0, center + half), 4)],
    }


def assignment_source(value: Any) -> str:
    """Collapse pathway provenance to direct lookup vs neighbor consensus."""

    return "exact_reference_lookup" if norm(value) == "exact" else "neighbor_consensus"


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
                    "predicted_smiles": clean(r.get("predicted_smiles")),
                    "predicted_pathway": clean(r.get("predicted_pathway")),
                    "predicted_subpathway": clean(r.get("predicted_subpathway")),
                    "source": assignment_source(r.get("pathway_source")),
                    "tier": clean(r.get("pathway_confidence_tier")),
                    "status": clean(r.get("pathway_assignment_status")),
                    "vote_share": pd.to_numeric(r.get("pathway_vote_share"), errors="coerce"),
                    "best_similarity": pd.to_numeric(r.get("pathway_best_reference_similarity"), errors="coerce"),
                }
            )
    return pd.DataFrame(rows)


def load_rank1_from_reports(root: Path) -> pd.DataFrame:
    """Load canonical rank-1 pathway + structure rows from delivered report sidecars."""

    rows: list[dict[str, Any]] = []
    for report_dir in sorted(root.glob("X-*")):
        if not report_dir.is_dir():
            continue
        slug = report_dir.name.lower()
        data_dir = report_dir / "data"
        pathway_path = data_dir / f"{slug}_pathway_prediction.csv"
        candidate_path = data_dir / f"{slug}_top_candidates.csv"
        if not pathway_path.exists() or not candidate_path.exists():
            continue
        with pathway_path.open(newline="") as handle:
            pathway = next(csv.DictReader(handle), {})
        with candidate_path.open(newline="") as handle:
            candidates = list(csv.DictReader(handle))
        candidates.sort(key=lambda row: pd.to_numeric(row.get("rank"), errors="coerce"))
        rank1 = candidates[0] if candidates else {}
        lib = clean(rank1.get("lib_id"))
        rows.append(
            {
                "control_id": report_dir.name,
                "lib": lib.removesuffix(".0"),
                "predicted_smiles": clean(rank1.get("predicted_smiles")),
                "predicted_pathway": clean(pathway.get("predicted_pathway")),
                "predicted_subpathway": clean(pathway.get("predicted_subpathway")),
                "source": assignment_source(pathway.get("pathway_source")),
                "tier": clean(pathway.get("pathway_confidence_tier")),
                "status": clean(pathway.get("pathway_assignment_status")),
                "vote_share": pd.to_numeric(pathway.get("pathway_vote_share"), errors="coerce"),
                "best_similarity": pd.to_numeric(
                    pathway.get("pathway_best_reference_similarity"), errors="coerce"
                ),
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


def build_db_inchikey_map(db_path: Path) -> dict[str, tuple[str, str]]:
    """InChIKey skeleton -> (metabolon_super, metabolon_sub) from the reference DB."""

    import sqlite3

    con = sqlite3.connect(str(db_path))
    db = pd.read_sql_query(
        "SELECT inchikeys, metabolon_super, metabolon_sub FROM metabolites "
        "WHERE TRIM(COALESCE(inchikeys,'')) <> '' AND TRIM(COALESCE(metabolon_super,'')) <> ''",
        con,
    )
    con.close()
    out: dict[str, tuple[str, str]] = {}
    for _, r in db.iterrows():
        for ik in str(r["inchikeys"]).replace("|", ";").split(";"):
            sk = ik.strip().split("-")[0]
            if len(sk) >= 14:
                out.setdefault(sk, (clean(r["metabolon_super"]), clean(r["metabolon_sub"])))
    return out


def inchikey_skeleton(smiles: str) -> str:
    """First InChIKey block for a SMILES (2D structure hash), or ''."""

    from rdkit import Chem

    mol = Chem.MolFromSmiles(str(smiles or "").split("|")[0].strip())
    if mol is None:
        return ""
    try:
        return Chem.MolToInchiKey(mol).split("-")[0]
    except Exception:
        return ""


def evaluate_vs_reference_db(
    pred: pd.DataFrame, gt: pd.DataFrame, db_map: dict[str, tuple[str, str]]
) -> dict[str, Any]:
    """Accuracy of predicted pathway vs the reference-DB Metabolon taxonomy that
    the predictor is built to reproduce (exact-structure subset via InChIKey), plus
    the DB-vs-known_X annotation agreement (the ceiling for a known_X-based score).
    """

    rows: list[dict[str, Any]] = []
    for _, r in pred.iterrows():
        sk = inchikey_skeleton(r.get("predicted_smiles"))
        if not sk or sk not in db_map:
            continue
        db_super, db_sub = db_map[sk]
        rows.append(
            {
                "control_id": r["control_id"],
                "predicted_pathway": r["predicted_pathway"],
                "predicted_subpathway": r["predicted_subpathway"],
                "db_super": db_super,
                "db_sub": db_sub,
                "source": r["source"],
                "tier": r["tier"],
                "status": r["status"],
                "best_similarity": r["best_similarity"],
            }
        )
    ex = pd.DataFrame(rows)
    if ex.empty:
        return {"n_exact_structure_in_db": 0}
    ex["correct_super"] = ex.apply(lambda x: norm(x["predicted_pathway"]) == norm(x["db_super"]), axis=1)
    ex["sub_evaluable"] = ex["db_sub"].map(norm).ne("")
    ex["correct_sub"] = ex.apply(
        lambda x: norm(x["predicted_subpathway"]) == norm(x["db_sub"]) and bool(x["sub_evaluable"]),
        axis=1,
    )
    ex["sim_bin"] = ex["best_similarity"].map(similarity_bin)
    metrics = super_sub_metrics(ex)
    result: dict[str, Any] = {
        "n_exact_structure_in_db": int(len(ex)),
        "super_accuracy_vs_db": metrics["super"]["accuracy"],
        "sub_accuracy_vs_db": metrics["sub"]["accuracy"],
        "overall": metrics,
        "by_assignment_source": {
            str(source): super_sub_metrics(group) for source, group in ex.groupby("source")
        },
        "by_pathway_tier": {
            str(tier): super_sub_metrics(group) for tier, group in ex.groupby("tier")
        },
        "by_assignment_status": {
            str(status): super_sub_metrics(group) for status, group in ex.groupby("status")
        },
        "by_similarity_bin": {
            str(bin_name): super_sub_metrics(group) for bin_name, group in ex.groupby("sim_bin")
        },
    }
    truth = gt.set_index("control_id")["true_super"].to_dict()
    agree = tot = 0
    for _, x in ex.iterrows():
        ks = norm(truth.get(x["control_id"], ""))
        if ks and ks not in EXCLUDE_TRUTH:
            tot += 1
            agree += int(ks == norm(x["db_super"]))
    if tot:
        result["db_vs_known_x_super_agreement"] = {"n": int(tot), "agreement": round(agree / tot, 4)}
    confusion: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for _, row in ex.iterrows():
        confusion[clean(row["db_super"])][clean(row["predicted_pathway"])] += 1
    result["confusion_super"] = {
        truth_label: dict(sorted(predicted.items(), key=lambda pair: pair[1], reverse=True))
        for truth_label, predicted in sorted(confusion.items())
    }
    result["distribution_super"] = {
        "predicted": {
            str(key): int(value) for key, value in Counter(ex["predicted_pathway"].map(clean)).most_common()
        },
        "truth": {str(key): int(value) for key, value in Counter(ex["db_super"].map(clean)).most_common()},
    }
    return result


def similarity_bin(value: float) -> str:
    """Label the similarity bin for a best-reference-similarity value."""

    if value is None or pd.isna(value):
        return "na"
    for lo, hi in SIM_BINS:
        if lo <= value < hi:
            return f"[{lo:.1f},{hi if hi <= 1 else 1.0:.1f})"
    return "na"


def accuracy_by(df: pd.DataFrame, key: str, correct_col: str) -> dict[str, dict[str, Any]]:
    """Accuracy + Wilson interval of ``correct_col`` grouped by ``key``."""

    out: dict[str, dict[str, Any]] = {}
    for value, grp in df.groupby(key):
        out[str(value)] = wilson(int(grp[correct_col].sum()), len(grp))
    return out


def accuracy_metrics(df: pd.DataFrame, correct_col: str) -> dict[str, Any]:
    """Accuracy + Wilson interval for one Boolean correctness column."""

    return wilson(int(df[correct_col].sum()), len(df))


def super_sub_metrics(df: pd.DataFrame) -> dict[str, Any]:
    """Super- and subpathway metrics for one evaluation cohort."""

    sub = df[df["sub_evaluable"]].copy() if "sub_evaluable" in df else df.iloc[0:0]
    return {
        "super": accuracy_metrics(df, "correct_super"),
        "sub": accuracy_metrics(sub, "correct_sub"),
    }


def build_summary(
    pred: pd.DataFrame,
    gt: pd.DataFrame,
    db_map: dict[str, tuple[str, str]] | None = None,
) -> tuple[dict[str, Any], pd.DataFrame]:
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
    ev["sub_evaluable"] = ev["true_sub"].map(norm).ne("")
    ev["sim_bin"] = ev["best_similarity"].map(similarity_bin)
    if db_map is not None:
        ev["in_reference_db"] = ev["predicted_smiles"].map(
            lambda smiles: inchikey_skeleton(smiles) in db_map
        )
    else:
        ev["in_reference_db"] = False

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
        "overall": super_sub_metrics(ev),
        "by_assignment_source": {
            str(source): super_sub_metrics(group) for source, group in ev.groupby("source")
        },
        "by_pathway_tier": {
            str(tier): super_sub_metrics(group) for tier, group in ev.groupby("tier")
        },
        "by_assignment_status": {
            str(status): super_sub_metrics(group) for status, group in ev.groupby("status")
        },
        "by_similarity_bin": {
            str(bin_name): super_sub_metrics(group) for bin_name, group in ev.groupby("sim_bin")
        },
        "by_reference_membership": {
            ("reference_structure" if bool(in_db) else "non_reference_structure"): super_sub_metrics(group)
            for in_db, group in ev.groupby("in_reference_db")
        },
    }
    neighbor = ev[ev["source"] == "neighbor_consensus"]
    accuracy["neighbor_consensus_by_reference_membership"] = {
        ("reference_structure" if bool(in_db) else "non_reference_structure"): super_sub_metrics(group)
        for in_db, group in neighbor.groupby("in_reference_db")
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
        "n_sub_evaluable": int(sub_evaluable),
        "k": k,
        "similarity_floor": floor,
        "similarity_power": power,
        "super_accuracy": round(correct_super / evaluated, 4) if evaluated else 0.0,
        "sub_accuracy": round(correct_sub / sub_evaluable, 4) if sub_evaluable else 0.0,
        "super": wilson(correct_super, evaluated),
        "sub": wilson(correct_sub, sub_evaluable),
        "super_accuracy_by_similarity_bin": {
            b: wilson(c[0], c[1]) for b, c in sorted(by_bin.items())
        },
        "scope_note": (
            "Intrinsic leave-one-out validation of the neighbor-consensus mechanism only; "
            "does not include production MW priors, SMARTS overrides, or other rules."
        ),
    }


def reference_db_md_lines(summary: dict[str, Any]) -> list[str]:
    """Markdown for reference-taxonomy consistency (not independent accuracy)."""

    ref = summary.get("reference_db_evaluation") or {}
    if not ref.get("n_exact_structure_in_db"):
        return []
    lines = [
        "## Reference-taxonomy consistency (not independent prediction accuracy)",
        f"- Exact-structure controls (InChIKey found in reference DB): **{ref['n_exact_structure_in_db']}**",
        f"- Overall super/sub consistency: **{ref['super_accuracy_vs_db']:.1%} / "
        f"{ref['sub_accuracy_vs_db']:.1%}**",
    ]
    for source, metrics in ref.get("by_assignment_source", {}).items():
        lines.append(
            f"- {source}: super {metrics['super']['accuracy']:.1%} "
            f"(n={metrics['super']['n']}), sub {metrics['sub']['accuracy']:.1%} "
            f"(n={metrics['sub']['n']})"
        )
    lines.append("- Direct exact-reference assignments copy the database label and are circular as accuracy tests.")
    agreement = ref.get("db_vs_known_x_super_agreement")
    if agreement:
        lines.append(
            f"- Caveat: the reference DB and known_X annotations agree only {agreement['agreement']:.1%} "
            f"(n={agreement['n']}) even on identical structures, so any known_X-based score is capped near that."
        )
    lines.append("")
    return lines


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
        *reference_db_md_lines(summary),
        "## Agreement with known_X annotation (independent labels, not the predictor's taxonomy)",
        f"- Evaluable (super): **{acc['n_evaluable_super']}** controls",
        f"- Super-pathway agreement with known_X: {acc['super_accuracy']:.1%} "
        "(understated: known_X uses a different Metabolon annotation than the reference DB the predictor targets)",
        f"- Sub-pathway agreement: {acc['sub_accuracy']:.1%} "
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
    parser.add_argument(
        "--reports-root",
        default=None,
        help="load canonical rank-1 rows from delivered report sidecars instead of equalization output",
    )
    parser.add_argument("--loo", action="store_true", help="also run leave-one-out validation over the reference DB")
    parser.add_argument("--db", default="vendor/metabolite_synonyms.sqlite")
    parser.add_argument("--loo-k", type=int, default=10)
    parser.add_argument("--loo-floor", type=float, default=0.3)
    parser.add_argument("--loo-power", type=float, default=2.0)
    args = parser.parse_args()

    libs = [x.strip() for x in args.libs.split(",") if x.strip()]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pred = (
        load_rank1_from_reports(Path(args.reports_root))
        if args.reports_root
        else load_rank1_predictions(Path(args.equalization_root), libs, args.pred_name)
    )
    gt = load_ground_truth(Path(args.known_x))
    print(f"[pathway-audit] {len(pred)} predicted controls; {len(gt)} known_X ground-truth controls")

    db_path = Path(args.db)
    if not db_path.exists() and Path("_bundle/metabolite_synonyms.sqlite").exists():
        db_path = Path("_bundle/metabolite_synonyms.sqlite")
    db_map = build_db_inchikey_map(db_path) if db_path.exists() else None
    summary, ev = build_summary(pred, gt, db_map)

    if db_path.exists():
        print("[pathway-audit] scoring vs reference-DB taxonomy (exact-structure subset) ...")
        summary["reference_db_evaluation"] = evaluate_vs_reference_db(pred, gt, db_map or {})

    ev_out = ev[
        [
            "control_id", "lib", "predicted_pathway", "true_super", "correct_super",
            "predicted_subpathway", "true_sub", "correct_sub", "source", "tier", "status",
            "best_similarity", "in_reference_db",
        ]
    ].sort_values(["correct_super", "tier"])
    ev_out.to_csv(out_dir / "audit_pathway_accuracy.csv", index=False)
    (out_dir / "audit_pathway_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_markdown(summary, out_dir / "audit_pathway_summary.md")

    ref = summary.get("reference_db_evaluation") or {}
    if ref.get("n_exact_structure_in_db"):
        print(
            f"[pathway-audit] reference-taxonomy consistency: super "
            f"{ref['super_accuracy_vs_db']:.1%} / sub {ref['sub_accuracy_vs_db']:.1%} "
            f"on {ref['n_exact_structure_in_db']} exact-structure controls"
        )
        g = ref.get("db_vs_known_x_super_agreement")
        if g:
            print(f"[pathway-audit] reference-DB vs known_X agreement: {g['agreement']:.1%} (n={g['n']})")
    a = summary["accuracy"]
    print(
        f"[pathway-audit] known_X agreement: super {a['super_accuracy']:.1%} on {a['n_evaluable_super']} "
        "evaluable (understated - known_X uses a different taxonomy than the predictor targets)"
    )
    print(f"[pathway-audit] known_X agreement by tier: {a['super_accuracy_by_tier']}")

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
