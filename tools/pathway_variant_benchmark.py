#!/usr/bin/env python3
"""Shared benchmark gate for pathway-predictor improvement variants.

Every candidate change to the pathway predictor (vendored neighbor vote, vote
tuning, top-K ensembles, evidence fusion) is scored here against the same
cohort and conventions before merge, so results are comparable across tracks.

Cohort and truth
----------------
The canonical evaluation cohort is the committed per-control audit
``docs/paper/data/pathway/audit_pathway_accuracy.csv`` (909 controls whose
known_X super label is real, produced by ``tools/audit_pathway_predictions.py``).
That file carries both the truth labels and the delivered structure-based
baseline, so any variant scored here uses byte-identical truth and an exact
head-to-head baseline.

Metrics
-------
- Raw and harmonized super-pathway accuracy with Wilson CIs. Harmonized merges
  the Amino Acid / Peptide structural continuum, mirroring
  ``tools/analyze_taxonomy_harmonization.py`` (see
  ``docs/paper/data/taxonomy_harmonization.json`` for why).
- Sub-pathway accuracy when the variant provides one.
- Coverage of the cohort (variants may abstain).
- Head-to-head vs the delivered baseline on the jointly-called set, with
  discordant-pair counts and an exact two-sided McNemar test.
- Leave-one-out of the neighbor-vote mechanism on the reference DB (gated on
  the fetched metabolite synonym SQLite; the tuning signal for vote changes).

Usage::

    python tools/pathway_variant_benchmark.py selfcheck
    python tools/pathway_variant_benchmark.py score \
        --pred my_variant_predictions.csv --label my_variant \
        --out-json results.json
    python tools/pathway_variant_benchmark.py abstention
    python tools/pathway_variant_benchmark.py loo --k 10 --floor 0.3 --power 2
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
BASELINE_AUDIT = REPO_ROOT / "docs" / "paper" / "data" / "pathway" / "audit_pathway_accuracy.csv"
AUDIT_SUMMARY = REPO_ROOT / "docs" / "paper" / "data" / "pathway" / "audit_pathway_summary.json"
REFERENCE_DB_CANDIDATES = (
    REPO_ROOT / "Metabolite Synonym Database" / "metabolite_synonyms_merged_with_mtx_giga.sqlite",
    REPO_ROOT / "vendor" / "metabolite_synonyms.sqlite",
)

# Structural continuum merged for harmonized scoring; mirrors SUPER_EQUIVALENCE
# in tools/analyze_taxonomy_harmonization.py (labels are compared normalized).
SUPER_EQUIVALENCE: tuple[frozenset[str], ...] = (frozenset({"amino acid", "peptide"}),)


def _load_tool(name: str):
    path = REPO_ROOT / "tools" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


_AUDIT = _load_tool("audit_pathway_predictions")
norm = _AUDIT.norm
clean = _AUDIT.clean
wilson = _AUDIT.wilson


def harmonize(label: str) -> str:
    """Normalized label collapsed to its equivalence class."""

    text = norm(label)
    for group in SUPER_EQUIVALENCE:
        if text in group:
            return "|".join(sorted(group))
    return text


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact McNemar p-value on discordant pair counts (binomial)."""

    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    p = sum(math.comb(n, i) for i in range(0, k + 1)) / 2.0**n
    return min(1.0, 2.0 * p)


def load_cohort(path: Path = BASELINE_AUDIT) -> list[dict[str, str]]:
    """The canonical 909-control cohort with truth and baseline predictions."""

    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"empty cohort file: {path}")
    return rows


def load_predictions(
    path: Path,
    *,
    id_col: str = "control_id",
    super_col: str = "predicted_super",
    sub_col: str | None = None,
) -> dict[str, dict[str, str]]:
    """control_id -> predicted labels; first occurrence wins; blanks = abstain."""

    out: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            cid = clean(row.get(id_col))
            if not cid or cid in out:
                continue
            out[cid] = {
                "predicted_super": clean(row.get(super_col)),
                "predicted_sub": clean(row.get(sub_col)) if sub_col else "",
            }
    return out


def score_rows(rows: list[dict[str, str]], pred_col: str, true_col: str) -> dict[str, Any]:
    """Raw + harmonized accuracy block for prediction/truth column pairs."""

    raw_correct = sum(1 for r in rows if norm(r[pred_col]) == norm(r[true_col]))
    harm_correct = sum(1 for r in rows if harmonize(r[pred_col]) == harmonize(r[true_col]))
    return {
        "raw": wilson(raw_correct, len(rows)),
        "harmonized": wilson(harm_correct, len(rows)),
    }


def evaluate_variant(
    predictions: dict[str, dict[str, str]],
    cohort: list[dict[str, str]],
    label: str,
) -> dict[str, Any]:
    """Score one variant against the cohort and head-to-head vs the baseline."""

    called: list[dict[str, str]] = []
    for row in cohort:
        pred = predictions.get(clean(row["control_id"]))
        if pred is None or not pred["predicted_super"]:
            continue
        called.append({**row, "variant_super": pred["predicted_super"], "variant_sub": pred["predicted_sub"]})

    n_cohort = len(cohort)
    result: dict[str, Any] = {
        "label": label,
        "n_cohort": n_cohort,
        "n_called": len(called),
        "coverage": round(len(called) / n_cohort, 4) if n_cohort else 0.0,
        "n_predictions_outside_cohort": len(set(predictions) - {clean(r["control_id"]) for r in cohort}),
    }
    if not called:
        return result

    result["super_accuracy_on_called"] = score_rows(called, "variant_super", "true_super")
    sub_rows = [r for r in called if norm(r["true_sub"]) and r["variant_sub"]]
    if sub_rows:
        sub_correct = sum(1 for r in sub_rows if norm(r["variant_sub"]) == norm(r["true_sub"]))
        result["sub_accuracy_on_called"] = wilson(sub_correct, len(sub_rows))

    by_true: dict[str, Any] = {}
    groups: dict[str, list[dict[str, str]]] = {}
    for r in called:
        groups.setdefault(clean(r["true_super"]), []).append(r)
    for true_super, rows in sorted(groups.items()):
        if len(rows) >= 5:
            by_true[true_super] = score_rows(rows, "variant_super", "true_super")
    result["accuracy_by_true_super"] = by_true

    # Head-to-head vs the delivered baseline on the jointly-called set.
    variant_raw = [norm(r["variant_super"]) == norm(r["true_super"]) for r in called]
    baseline_raw = [norm(r["predicted_pathway"]) == norm(r["true_super"]) for r in called]
    variant_harm = [harmonize(r["variant_super"]) == harmonize(r["true_super"]) for r in called]
    baseline_harm = [harmonize(r["predicted_pathway"]) == harmonize(r["true_super"]) for r in called]

    def head_to_head(variant_ok: list[bool], baseline_ok: list[bool]) -> dict[str, Any]:
        b = sum(1 for v, s in zip(variant_ok, baseline_ok) if v and not s)
        c = sum(1 for v, s in zip(variant_ok, baseline_ok) if not v and s)
        return {
            "n_joint": len(variant_ok),
            "variant_accuracy": round(sum(variant_ok) / len(variant_ok), 4),
            "baseline_accuracy": round(sum(baseline_ok) / len(baseline_ok), 4),
            "variant_only_correct": b,
            "baseline_only_correct": c,
            "both_correct": sum(1 for v, s in zip(variant_ok, baseline_ok) if v and s),
            "neither_correct": sum(1 for v, s in zip(variant_ok, baseline_ok) if not v and not s),
            "mcnemar_p_two_sided": float(f"{mcnemar_exact(b, c):.3g}"),
        }

    result["head_to_head_raw"] = head_to_head(variant_raw, baseline_raw)
    result["head_to_head_harmonized"] = head_to_head(variant_harm, baseline_harm)

    by_tier: dict[str, Any] = {}
    tier_groups: dict[str, list[int]] = {}
    for idx, r in enumerate(called):
        tier_groups.setdefault(clean(r["tier"]) or "(none)", []).append(idx)
    for tier, idxs in sorted(tier_groups.items()):
        by_tier[tier] = {
            "variant_raw": wilson(sum(variant_raw[i] for i in idxs), len(idxs)),
            "baseline_raw": wilson(sum(baseline_raw[i] for i in idxs), len(idxs)),
        }
    result["by_baseline_tier"] = by_tier
    return result


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def run_selfcheck() -> int:
    """Score the delivered baseline against itself and the committed summary.

    Verifies harness conventions: accuracy must equal the committed audit
    figure exactly and the self head-to-head must be perfectly concordant.
    """

    cohort = load_cohort()
    predictions = {
        clean(r["control_id"]): {
            "predicted_super": clean(r["predicted_pathway"]),
            "predicted_sub": clean(r["predicted_subpathway"]),
        }
        for r in cohort
    }
    result = evaluate_variant(predictions, cohort, label="baseline_selfcheck")
    committed = json.loads(AUDIT_SUMMARY.read_text())["accuracy"]

    checks = {
        "coverage_is_full": result["coverage"] == 1.0,
        "raw_accuracy_matches_committed_audit": (
            result["super_accuracy_on_called"]["raw"]["accuracy"] == committed["super_accuracy"]
        ),
        "self_head_to_head_concordant": (
            result["head_to_head_raw"]["variant_only_correct"] == 0
            and result["head_to_head_raw"]["baseline_only_correct"] == 0
            and result["head_to_head_raw"]["mcnemar_p_two_sided"] == 1.0
        ),
    }
    print(json.dumps({"checks": checks, "result": result}, indent=2))
    ok = all(checks.values())
    print(f"[selfcheck] {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def run_score(args: argparse.Namespace) -> int:
    cohort = load_cohort(args.cohort)
    predictions = load_predictions(
        args.pred, id_col=args.id_col, super_col=args.pred_super_col, sub_col=args.pred_sub_col
    )
    result = evaluate_variant(predictions, cohort, label=args.label)
    payload = json.dumps(result, indent=2)
    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(payload + "\n")
        print(f"[score] wrote {args.out_json}")
    print(payload)
    return 0


ABSTENTION_POLICIES: tuple[tuple[str, str], ...] = (
    ("all_assignments", "no abstention (delivered behavior)"),
    ("drop_low_tier", "abstain on low-confidence tier"),
    ("high_tier_only", "abstain unless high tier"),
    ("exact_status_only", "abstain unless exact reference match"),
    ("similarity_ge_0.5", "abstain when best reference similarity < 0.5"),
    ("similarity_ge_0.7", "abstain when best reference similarity < 0.7"),
)


def _policy_keep(policy: str, row: dict[str, str]) -> bool:
    tier = clean(row["tier"]).lower()
    status = clean(row["status"]).lower()
    try:
        similarity = float(row["best_similarity"])
    except (TypeError, ValueError):
        similarity = float("nan")
    if policy == "all_assignments":
        return True
    if policy == "drop_low_tier":
        return tier != "low"
    if policy == "high_tier_only":
        return tier == "high"
    if policy == "exact_status_only":
        return status == "exact_reference_match"
    if policy.startswith("similarity_ge_"):
        threshold = float(policy.rsplit("_", 1)[1])
        return similarity == similarity and similarity >= threshold
    raise ValueError(f"unknown policy: {policy}")


def run_abstention(args: argparse.Namespace) -> int:
    """Precision-coverage curve for abstention policies over the delivered audit."""

    cohort = load_cohort(args.cohort)
    n_cohort = len(cohort)
    points: list[dict[str, Any]] = []
    for policy, description in ABSTENTION_POLICIES:
        kept = [r for r in cohort if _policy_keep(policy, r)]
        if not kept:
            continue
        scores = score_rows(kept, "predicted_pathway", "true_super")
        points.append(
            {
                "policy": policy,
                "description": description,
                "n_called": len(kept),
                "coverage": round(len(kept) / n_cohort, 4),
                "super_accuracy_raw": scores["raw"],
                "super_accuracy_harmonized": scores["harmonized"],
            }
        )

    # Recommend the highest-coverage policy whose harmonized accuracy clears the
    # bar; the bar defaults to 0.70, near the measured annotation ceiling share.
    eligible = [p for p in points if (p["super_accuracy_harmonized"]["accuracy"] or 0) >= args.min_harmonized_accuracy]
    recommended = max(eligible, key=lambda p: p["coverage"]) if eligible else None
    summary = {
        "source": str(BASELINE_AUDIT.relative_to(REPO_ROOT)),
        "n_cohort": n_cohort,
        "min_harmonized_accuracy_bar": args.min_harmonized_accuracy,
        "points": points,
        "recommended_policy": recommended["policy"] if recommended else None,
        "note": (
            "Coverage is the share of the evaluable cohort still given a pathway; "
            "abstained compounds would be reported as insufficient evidence instead "
            "of carrying a low-tier label that is wrong most of the time."
        ),
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2) + "\n")

    lines = [
        "# Pathway abstention precision-coverage analysis",
        "",
        f"Cohort: {n_cohort} evaluable controls from `{summary['source']}`.",
        "",
        "| policy | n called | coverage | raw acc | harmonized acc |",
        "|---|---|---|---|---|",
    ]
    for p in points:
        lines.append(
            f"| {p['policy']} | {p['n_called']} | {p['coverage']:.1%}"
            f" | {p['super_accuracy_raw']['accuracy']:.1%}"
            f" | {p['super_accuracy_harmonized']['accuracy']:.1%} |"
        )
    lines += [
        "",
        (
            f"Recommended operating point (harmonized accuracy >= {args.min_harmonized_accuracy:.0%},"
            f" maximize coverage): **{summary['recommended_policy']}**."
            if recommended
            else "No policy clears the accuracy bar."
        ),
        "",
    ]
    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.write_text("\n".join(lines))
    print(f"[abstention] wrote {args.out_json} and {args.out_md}")
    print(json.dumps(summary, indent=2))
    return 0


def find_reference_db() -> Path | None:
    for candidate in REFERENCE_DB_CANDIDATES:
        if candidate.exists():
            return candidate
    return None


def run_loo(args: argparse.Namespace) -> int:
    """Leave-one-out of the neighbor vote on the reference DB (gated on asset)."""

    db_path = args.reference_db or find_reference_db()
    if db_path is None or not Path(db_path).exists():
        print(
            "[loo] SKIPPED: metabolite synonym SQLite not present "
            f"(looked for {[str(p) for p in REFERENCE_DB_CANDIDATES]}). "
            "Fetch it with ./fetch_assets.sh and re-run."
        )
        return 0
    result = _AUDIT.leave_one_out_validation(Path(db_path), k=args.k, floor=args.floor, power=args.power)
    payload = json.dumps(result, indent=2)
    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(payload + "\n")
    print(payload)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("selfcheck", help="verify harness conventions against the committed audit")

    ps = sub.add_parser("score", help="score a variant predictions CSV against the cohort")
    ps.add_argument("--pred", type=Path, required=True)
    ps.add_argument("--label", default="variant")
    ps.add_argument("--id-col", default="control_id")
    ps.add_argument("--pred-super-col", default="predicted_super")
    ps.add_argument("--pred-sub-col", default=None)
    ps.add_argument("--cohort", type=Path, default=BASELINE_AUDIT)
    ps.add_argument("--out-json", type=Path, default=None)

    curve_dir = REPO_ROOT / "docs" / "paper" / "data" / "pathway"
    pa = sub.add_parser("abstention", help="precision-coverage analysis of abstention policies")
    pa.add_argument("--cohort", type=Path, default=BASELINE_AUDIT)
    pa.add_argument("--min-harmonized-accuracy", type=float, default=0.70)
    pa.add_argument("--out-json", type=Path, default=curve_dir / "pathway_abstention_curve.json")
    pa.add_argument("--out-md", type=Path, default=curve_dir / "pathway_abstention_curve.md")

    pl = sub.add_parser("loo", help="leave-one-out of the neighbor vote (needs reference DB)")
    pl.add_argument("--reference-db", type=Path, default=None)
    pl.add_argument("--k", type=int, default=10)
    pl.add_argument("--floor", type=float, default=0.3)
    pl.add_argument("--power", type=float, default=2.0)
    pl.add_argument("--out-json", type=Path, default=None)

    args = parser.parse_args()
    if args.command == "selfcheck":
        return run_selfcheck()
    if args.command == "score":
        return run_score(args)
    if args.command == "abstention":
        return run_abstention(args)
    if args.command == "loo":
        return run_loo(args)
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
