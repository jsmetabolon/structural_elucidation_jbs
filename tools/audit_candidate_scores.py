#!/usr/bin/env python3
"""Audit per-candidate score behavior across reports (before vs after the rework).

The shipped top-10 widget leads with the pipeline's ``evidence_confidence`` for
each candidate. That number (a) saturates near 1.0 and is shared across
mass-consistent isomers, so every candidate looks identical, and (b) is a
reranker *input feature*, not the ranking objective, so it is frequently higher
for a lower-ranked candidate than for a higher-ranked one. This tool quantifies
both failure modes and confirms the new ``candidate_confidence`` fixes them.

Per report it reads ``data/`` sidecars:
  ``<slug>_candidate_score_summary.csv`` -> rank, candidate_report_score
      (= unified_final_score, the ordering objective), evidence_confidence,
      mass_match_score, predicted_formula, predicted_smiles
  ``<slug>_confidence_inputs.json``      -> confidence fields (overall anchor)
      and, once patched, the baked per-rank ``candidate_evidence``.

For the top-N candidates it computes, for each score series:
  - all_equal: max-min <= EPS (a non-discriminative column),
  - inversions: count of rank pairs (i<j) whose score rises with rank
      (score[i] < score[j]) -> a higher-ranked candidate scoring lower,
  - Kendall tau vs rank (1.0 = perfectly monotone non-increasing).

It reports these for the currently-featured ``evidence_confidence``, for the
ordering ``candidate_report_score``, and for the new rank-consistent
``candidate_confidence`` (``confidence.candidate_confidences``), which must show
zero inversions and only labelled tie groups.

Pure standard library.

Usage:
    python3 tools/audit_candidate_scores.py <reports_root> [--top-n 10] \
        [--json out.json] [--examples 12]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

csv.field_size_limit(10 ** 7)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
import confidence as C  # noqa: E402

EPS = 1e-6


def _num(value):
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return None if f != f else f


def _read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _rank_order(rows: list[dict]) -> list[dict]:
    return sorted(rows, key=lambda r: (_num(r.get("rank")) if _num(r.get("rank")) is not None else 1e9))


def all_equal(values: list[float]) -> bool:
    vals = [v for v in values if v is not None]
    return len(vals) >= 2 and (max(vals) - min(vals)) <= EPS


def inversion_count(values: list[float]) -> tuple[int, int]:
    """Pairs (i<j in rank order) where a later rank has a strictly higher score."""
    seq = [v for v in values if v is not None]
    inv = pairs = 0
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            pairs += 1
            if seq[j] > seq[i] + EPS:
                inv += 1
    return inv, pairs


def kendall_tau_vs_rank(values: list[float]) -> float | None:
    """tau in [-1,1]; 1.0 means the score is perfectly non-increasing with rank."""
    seq = [v for v in values if v is not None]
    n = len(seq)
    if n < 2:
        return None
    concordant = discordant = 0
    for i in range(n):
        for j in range(i + 1, n):
            diff = seq[i] - seq[j]  # expect >= 0 (non-increasing) for concordance
            if diff > EPS:
                concordant += 1
            elif diff < -EPS:
                discordant += 1
    total = concordant + discordant
    return None if total == 0 else round((concordant - discordant) / total, 4)


def _load_ci(data_dir: Path, slug: str) -> dict:
    p = data_dir / f"{slug}_confidence_inputs.json"
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text())
    except (ValueError, OSError):
        return {}


def _overall_anchor(data_dir: Path, slug: str, ci: dict) -> float | None:
    # confidence_inputs.json IS the confidence fields dict; extra keys are ignored.
    if ci:
        try:
            return C.compute_confidence(ci).overall
        except Exception:
            pass
    fields = C.read_fields_from_sidecars(data_dir, slug)
    if not fields:
        return None
    try:
        return C.compute_confidence(fields).overall
    except Exception:
        return None


def audit_report(report_dir: Path, top_n: int) -> dict | None:
    data = report_dir / "data"
    if not data.is_dir():
        return None
    slug = report_dir.name.lower()
    css = _read_csv(data / f"{slug}_candidate_score_summary.csv")
    if not css:
        return None
    rows = _rank_order(css)[:top_n]
    if len(rows) < 2:
        return None

    # The delivered widget features evidence_confidence from candidate_evidence in
    # confidence_inputs.json (injected by the v1 patch), NOT from the CSV (which has
    # no evidence_confidence column). Read the featured value from there.
    ci = _load_ci(data, slug)
    ev_by_rank: dict[int, dict] = {}
    for r in ci.get("candidate_evidence") or []:
        try:
            ev_by_rank[int(r["rank"])] = r
        except (KeyError, TypeError, ValueError):
            continue

    def _rk(row: dict, i: int) -> int:
        v = _num(row.get("rank"))
        return int(v) if v is not None else i + 1

    featured, mass, ordering, shipped_cc = [], [], [], []
    for i, r in enumerate(rows):
        ev = ev_by_rank.get(_rk(r, i), {})
        fe = ev.get("evidence_confidence")
        featured.append(_num(fe if fe is not None else r.get("evidence_confidence")))
        mm = ev.get("mass_match_score")
        mass.append(_num(mm if mm is not None else r.get("mass_match_score")))
        ordering.append(_num(r.get("candidate_report_score")))
        shipped_cc.append(_num(ev.get("candidate_confidence")))

    overall = _overall_anchor(data, slug, ci)
    after = C.candidate_confidences(rows, overall)
    after_scores = [a["candidate_confidence"] for a in after]
    tie_groups = {a["tie_group"] for a in after}
    tie_group_members = sum(1 for a in after if a["tied"])
    rank1_tied = bool(after and after[0]["tied"])
    # parity: does any already-refreshed candidate_confidence match the computed value?
    refreshed = all(v is not None for v in shipped_cc)
    shipped_matches = refreshed and all(
        abs(s - a) <= 1e-6 for s, a in zip(shipped_cc, after_scores)
    )

    f_inv, f_pairs = inversion_count(featured)
    o_inv, _ = inversion_count(ordering)
    a_inv, _ = inversion_count(after_scores)

    def spread(values: list[float]) -> float | None:
        vals = [v for v in values if v is not None]
        return round(max(vals) - min(vals), 4) if len(vals) >= 2 else None

    return {
        "control_id": report_dir.name,
        "n_candidates": len(rows),
        "overall_anchor": overall,
        "refreshed": refreshed,
        "shipped_matches_computed": shipped_matches,
        "featured_evidence": {
            "all_equal": all_equal(featured),
            "inversions": f_inv,
            "inversion_pairs": f_pairs,
            "kendall_tau_vs_rank": kendall_tau_vs_rank(featured),
            "spread": spread(featured),
        },
        "ordering_score": {
            "all_equal": all_equal(ordering),
            "inversions": o_inv,
            "kendall_tau_vs_rank": kendall_tau_vs_rank(ordering),
            "spread": spread(ordering),
        },
        "candidate_confidence": {
            "all_equal": all_equal(after_scores),
            "inversions": a_inv,
            "kendall_tau_vs_rank": kendall_tau_vs_rank(after_scores),
            "spread": spread(after_scores),
            "n_tie_groups": len(tie_groups),
            "tied_members": tie_group_members,
            "rank1": after_scores[0] if after_scores else None,
            "rank1_tied": rank1_tied,
        },
        "mass_all_equal": all_equal(mass),
    }


def summarize(reports: list[dict]) -> dict:
    n = len(reports)
    if not n:
        return {"n_reports": 0}

    def frac(pred) -> float:
        return round(sum(1 for r in reports if pred(r)) / n, 4)

    def median(values: list[float]) -> float | None:
        vals = sorted(v for v in values if v is not None)
        if not vals:
            return None
        mid = len(vals) // 2
        return round(vals[mid] if len(vals) % 2 else (vals[mid - 1] + vals[mid]) / 2, 4)

    n_refreshed = sum(1 for r in reports if r.get("refreshed"))
    return {
        "n_reports": n,
        "n_refreshed": n_refreshed,
        "refreshed_parity_ok": (
            all(r.get("shipped_matches_computed") for r in reports if r.get("refreshed"))
            if n_refreshed else None
        ),
        "featured_evidence": {
            "frac_all_equal": frac(lambda r: r["featured_evidence"]["all_equal"]),
            "frac_with_inversion": frac(lambda r: r["featured_evidence"]["inversions"] > 0),
            "total_inversions": sum(r["featured_evidence"]["inversions"] for r in reports),
            "median_spread": median([r["featured_evidence"]["spread"] for r in reports]),
            "median_tau": median([r["featured_evidence"]["kendall_tau_vs_rank"] for r in reports]),
        },
        "ordering_score": {
            "frac_all_equal": frac(lambda r: r["ordering_score"]["all_equal"]),
            "frac_with_inversion": frac(lambda r: r["ordering_score"]["inversions"] > 0),
            "median_tau": median([r["ordering_score"]["kendall_tau_vs_rank"] for r in reports]),
        },
        "candidate_confidence": {
            "frac_all_equal": frac(lambda r: r["candidate_confidence"]["all_equal"]),
            "frac_with_inversion": frac(lambda r: r["candidate_confidence"]["inversions"] > 0),
            "total_inversions": sum(r["candidate_confidence"]["inversions"] for r in reports),
            "median_spread": median([r["candidate_confidence"]["spread"] for r in reports]),
            "median_tau": median([r["candidate_confidence"]["kendall_tau_vs_rank"] for r in reports]),
            "frac_with_tie_group": frac(lambda r: r["candidate_confidence"]["n_tie_groups"] < r["n_candidates"]),
            "frac_rank1_tied": frac(lambda r: r["candidate_confidence"].get("rank1_tied", False)),
            "median_rank1": median([r["candidate_confidence"]["rank1"] for r in reports]),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reports_root", type=Path, nargs="?", default=Path("_audit_pull/reports"))
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--json", type=Path, default=None)
    parser.add_argument("--examples", type=int, default=12, help="worst featured-inversion examples to list")
    args = parser.parse_args()

    if not args.reports_root.is_dir():
        raise SystemExit(f"reports root not found: {args.reports_root}")

    reports: list[dict] = []
    skipped = 0
    for d in sorted(p for p in args.reports_root.iterdir() if p.is_dir() and p.name.upper().startswith("X-")):
        rec = audit_report(d, args.top_n)
        if rec is None:
            skipped += 1
            continue
        reports.append(rec)

    summary = summarize(reports)
    worst = sorted(
        reports, key=lambda r: (r["featured_evidence"]["inversions"], r["featured_evidence"]["all_equal"]), reverse=True
    )[: args.examples]

    payload = {"summary": summary, "skipped": skipped, "worst_featured_inversions": worst, "top_n": args.top_n}
    if args.json:
        args.json.write_text(json.dumps(payload, indent=2))

    s = summary
    print(f"[audit-candidate-scores] reports audited: {s.get('n_reports', 0)} (skipped {skipped})")
    if not s.get("n_reports"):
        return
    if s.get("n_refreshed"):
        print(f"  refreshed reports: {s['n_refreshed']}  shipped==computed parity: {s['refreshed_parity_ok']}")
    fe, os_, cc = s["featured_evidence"], s["ordering_score"], s["candidate_confidence"]
    print("  CURRENTLY FEATURED (evidence_confidence):")
    print(f"    all-equal across top-N : {fe['frac_all_equal']:.1%}")
    print(f"    >=1 rank/score inversion: {fe['frac_with_inversion']:.1%}  (total {fe['total_inversions']})")
    print(f"    median spread / tau     : {fe['median_spread']} / {fe['median_tau']}")
    print("  ORDERING (candidate_report_score = unified_final_score):")
    print(f"    >=1 inversion           : {os_['frac_with_inversion']:.1%}   median tau {os_['median_tau']}")
    print("  NEW candidate_confidence:")
    print(f"    all-equal across top-N  : {cc['frac_all_equal']:.1%}  (honest tie groups)")
    print(f"    >=1 inversion           : {cc['frac_with_inversion']:.1%}   (total {cc['total_inversions']})")
    print(f"    median spread / tau     : {cc['median_spread']} / {cc['median_tau']}")
    print(f"    median rank-1 conf      : {cc['median_rank1']}")
    print(f"    rank-1 tied (ambiguous) : {cc['frac_rank1_tied']:.1%}")
    if worst:
        print("  worst featured-inversion examples:", ", ".join(w["control_id"] for w in worst[:args.examples]))


if __name__ == "__main__":
    main()
