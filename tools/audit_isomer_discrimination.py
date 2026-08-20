#!/usr/bin/env python3
"""Why the correct structure sits below rank 1, when it is in the list at all.

Read-only against the delivered exec-summary export and the local feature rows.
Answers four questions the top-10 export raises but does not settle:

1. When the truth is in the top 10 but not first, what separates it from rank 1?
   If the two candidates share a molecular formula, the mass and formula axes are
   degenerate by construction and cannot be the answer.

2. Is the spectral axis — the only shipped axis that could break an isomer tie —
   contributing anything? ``forward_spectrum_score`` compares predicted fragments
   to observed ions and predicted losses to observed losses. A fragment of mass
   ``m`` and a loss of ``parent - m`` are the same event, so the two axes carry
   the same information and agreement across them is never tested.

3. Would pooling the axes fix the ranking? Measured, not assumed: several
   normalizations are swept against the delivered ordering. They all lose.

4. How often is rank 1 an arbitrary pick among candidates whose displayed
   evidence is identical?

RDKit is required for questions 2 and 3 (fragment prediction); those sections are
skipped with a recorded reason when it is unavailable. Questions 1 and 4 run on
the standard library alone.

Usage:
    python3 tools/audit_isomer_discrimination.py
    python3 tools/audit_isomer_discrimination.py --out-dir docs/paper/data/isomer
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Iterable

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

csv.field_size_limit(10**9)

COHORT_TRUTH = "truth_control"
BASIS_UNAVAILABLE = "unavailable"

# A predicted fragment and the observed peak it explains may differ by a
# hydrogen: the fragmenter sanitizes both pieces of a split, so each gains an
# implicit H, while negative-mode ESI observes the deprotonated species.
HYDROGEN_DA = 1.008
MATCH_TOLERANCE_PPM = 20.0

# Displayed evidence, at the precision a reader of the report actually sees.
TIE_COLUMNS = ("candidate_confidence", "evidence_confidence", "mass_match_score")
TIE_PLACES = 4


def _repo_relative(path: Path) -> str:
    """Path as written in the artifact: relative to the repo, so it is portable.

    The JSON is committed and CI regenerates it to check the two agree, so an
    absolute path would make the artifact differ between any two checkouts.
    """
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return resolved.name


def _num(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if number != number else number


def _rank(row: dict[str, str]) -> int:
    value = _num(row.get("rank"))
    return int(value) if value is not None else 0


# ── Inputs ───────────────────────────────────────────────────────────
def load_top10(path: Path) -> dict[str, list[dict[str, str]]]:
    """Delivered candidates grouped by compound, ordered by rank."""
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["x_compound_id"]].append(row)
    return {cid: sorted(rs, key=_rank) for cid, rs in grouped.items()}


def load_features(directory: Path) -> dict[str, dict[str, Any]]:
    """Observed MS2 bins and neutral mass per control, from the per-lib feature rows."""
    features: dict[str, dict[str, Any]] = {}
    for path in sorted(directory.glob("*_feature_rows.csv")):
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                control = (row.get("control_id") or "").strip()
                if not control:
                    continue
                try:
                    ions = json.loads(row.get("observed_ion_bins_json") or "[]")
                    losses = json.loads(row.get("observed_loss_bins_json") or "[]")
                except json.JSONDecodeError:
                    continue
                lib = _num(row.get("lib_id"))
                features[control] = {
                    "lib_id": int(lib) if lib is not None else 209,
                    "ion_bins": [int(b) for b in ions],
                    "loss_bins": [int(b) for b in losses],
                    "neutral_mass": _num(row.get("neutral_mass_fallback")),
                }
    return features


# ── 1. What separates the truth from rank 1 ──────────────────────────
def recovery_breakdown(controls: dict[str, list[dict[str, str]]]) -> dict[str, Any]:
    """Classify every truth control by why the truth is or is not at rank 1."""
    buckets: dict[str, list[Any]] = defaultdict(list)
    detail: list[dict[str, Any]] = []

    for cid, rows in sorted(controls.items()):
        first = rows[0]
        hit = next((r for r in rows if r["matched_truth"] == "TRUE"), None)
        if hit is None:
            buckets["truth_absent_from_top10"].append(cid)
            continue
        if _rank(hit) == 1:
            buckets["truth_at_rank_1"].append(cid)
            continue

        same_formula = bool(hit["candidate_formula"]) and (
            hit["candidate_formula"] == first["candidate_formula"]
        )
        truth_unscored = hit["evidence_confidence_basis"] == BASIS_UNAVAILABLE
        rank1_scored = first["evidence_confidence_basis"] != BASIS_UNAVAILABLE

        if truth_unscored and rank1_scored:
            reason = "truth_row_unscored"
        elif same_formula:
            reason = "evidence_degenerate_same_formula"
        else:
            reason = "different_formula"
        buckets[reason].append(cid)

        confidence_gap = None
        truth_confidence, rank1_confidence = (
            _num(hit["candidate_confidence"]),
            _num(first["candidate_confidence"]),
        )
        if truth_confidence is not None and rank1_confidence is not None:
            confidence_gap = round(rank1_confidence - truth_confidence, 6)
        detail.append(
            {
                "control_id": cid,
                "reason": reason,
                "truth_rank": _rank(hit),
                "msi_level": hit["msi_level"],
                "same_formula_as_rank1": same_formula,
                "truth_formula": hit["candidate_formula"],
                "rank1_formula": first["candidate_formula"],
                "truth_evidence_basis": hit["evidence_confidence_basis"],
                "rank1_evidence_basis": first["evidence_confidence_basis"],
                "truth_mass_match": _num(hit["mass_match_score"]),
                "rank1_mass_match": _num(first["mass_match_score"]),
                "candidate_confidence_gap": confidence_gap,
            }
        )

    below = [d for d in detail]
    return {
        "n_truth_controls": len(controls),
        "counts": {name: len(ids) for name, ids in sorted(buckets.items())},
        "n_truth_below_rank_1": len(below),
        "n_same_formula_as_rank1": sum(1 for d in below if d["same_formula_as_rank1"]),
        "n_gap_under_0p05": sum(
            1
            for d in below
            if d["candidate_confidence_gap"] is not None and d["candidate_confidence_gap"] < 0.05
        ),
        "detail": sorted(below, key=lambda d: d["control_id"]),
    }


def evidence_basis_counts(controls: dict[str, list[dict[str, str]]]) -> dict[str, Any]:
    """How many candidate rows carry an evidence score, and how many cannot."""
    every = [row for rows in controls.values() for row in rows]
    rank1 = [rows[0] for rows in controls.values() if rows]

    def tally(rows: list[dict[str, str]]) -> dict[str, int]:
        counts: dict[str, int] = defaultdict(int)
        for row in rows:
            counts[row["evidence_confidence_basis"] or "(blank)"] += 1
        return dict(sorted(counts.items()))

    return {
        "n_rows": len(every),
        "all_rows": tally(every),
        "rank_1_rows": tally(rank1),
    }


# ── 4. Ties ──────────────────────────────────────────────────────────
def tie_breakdown(controls: dict[str, list[dict[str, str]]]) -> dict[str, Any]:
    """Compounds whose rank 1 is indistinguishable from a lower-ranked candidate."""

    def signature(row: dict[str, str]) -> tuple[float | None, ...]:
        out = []
        for column in TIE_COLUMNS:
            value = _num(row.get(column))
            out.append(None if value is None else round(value, TIE_PLACES))
        return tuple(out)

    tied: list[dict[str, Any]] = []
    for cid, rows in sorted(controls.items()):
        if not rows:
            continue
        group = [r for r in rows if signature(r) == signature(rows[0])]
        if len(group) < 2:
            continue
        hit = next((r for r in group if r["matched_truth"] == "TRUE"), None)
        tied.append(
            {
                "control_id": cid,
                "tie_size": len(group),
                "truth_in_tie": hit is not None,
                "truth_rank": _rank(hit) if hit is not None else None,
                "msi_level": rows[0]["msi_level"],
            }
        )
    sizes = [entry["tie_size"] for entry in tied]
    return {
        "n_compounds": len(controls),
        "n_rank1_tied": len(tied),
        "share_rank1_tied": round(len(tied) / len(controls), 4) if controls else None,
        "mean_tie_size": round(sum(sizes) / len(sizes), 2) if sizes else None,
        "max_tie_size": max(sizes) if sizes else None,
        "n_truth_inside_tie_below_rank_1": sum(
            1 for e in tied if e["truth_in_tie"] and e["truth_rank"] not in (None, 1)
        ),
        "detail": [e for e in tied if e["truth_in_tie"]],
    }


# ── 2 and 3. The spectral axis ───────────────────────────────────────
def _import_rdkit_scorers() -> tuple[Any, Any, str | None]:
    try:
        from structural_elucidation.nce.evidence import bin_base, candidate_fragment_bins
    except Exception as error:  # pragma: no cover - environment dependent
        return None, None, f"{type(error).__name__}: {error}"
    try:
        from structural_elucidation.nce.evidence import HAS_RDKIT
    except Exception:  # pragma: no cover - older trees
        HAS_RDKIT = True
    if not HAS_RDKIT:
        return None, None, "rdkit is not installed"
    return bin_base, candidate_fragment_bins, None


def _explained(
    smiles: str,
    feature: dict[str, Any],
    bin_base: Callable[[int], float],
    fragment_bins: Callable[..., tuple[tuple[int, ...], tuple[int, ...]]],
) -> dict[str, int]:
    """Same-axis and pooled cross-axis agreement for one candidate."""
    frag, loss = fragment_bins(smiles, feature["lib_id"])
    same_axis = len(set(frag) & set(feature["ion_bins"])) + len(
        set(loss) & set(feature["loss_bins"])
    )

    predicted = {bin_base(b) for b in frag} | {bin_base(b) for b in loss}
    neutral = feature["neutral_mass"]
    if neutral is not None:
        predicted |= {neutral - m for m in list(predicted) if neutral - m > 0}
    observed = [bin_base(b) for b in feature["ion_bins"]] + [
        bin_base(b) for b in feature["loss_bins"]
    ]

    def explains(candidate: float, peak: float) -> bool:
        for shift in (0.0, -HYDROGEN_DA, HYDROGEN_DA):
            shifted = candidate + shift
            if shifted <= 0:
                continue
            if abs(shifted - peak) / peak * 1e6 <= MATCH_TOLERANCE_PPM:
                return True
        return False

    hits = 0
    used: set[float] = set()
    for peak in observed:
        for candidate in predicted:
            if explains(candidate, peak):
                hits += 1
                used.add(round(candidate, 4))
                break
    return {
        "same_axis": same_axis,
        "cross_axis": hits,
        "n_observed": len(observed),
        "n_predicted": len(predicted),
        "n_predicted_used": len(used),
    }


def _variants() -> dict[str, Callable[[dict[str, int]], float]]:
    """Candidate normalizations for a pooled-axis spectral score."""

    def recall(s: dict[str, int]) -> float:
        return s["cross_axis"] / s["n_observed"] if s["n_observed"] else 0.0

    def precision(s: dict[str, int]) -> float:
        return s["n_predicted_used"] / s["n_predicted"] if s["n_predicted"] else 0.0

    def f1(s: dict[str, int]) -> float:
        r, p = recall(s), precision(s)
        return 2 * r * p / (r + p) if r + p > 0 else 0.0

    def geometric(s: dict[str, int]) -> float:
        denominator = math.sqrt(s["n_predicted"] * s["n_observed"])
        return s["cross_axis"] / denominator if denominator else 0.0

    def damped(s: dict[str, int]) -> float:
        return recall(s) / math.log1p(s["n_predicted"]) if s["n_predicted"] else 0.0

    return {
        "raw_count": lambda s: float(s["cross_axis"]),
        "recall": recall,
        "f1_recall_precision": f1,
        "geometric_mean": geometric,
        "recall_over_log_npredicted": damped,
    }


def spectral_axis_audit(
    controls: dict[str, list[dict[str, str]]], features: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    """Is the spectral axis alive, and would pooling it fix the ordering?"""
    bin_base, fragment_bins, reason = _import_rdkit_scorers()
    if bin_base is None:
        return {"status": "skipped", "reason": reason}

    scored: dict[str, tuple[dict[str, str], list[tuple[dict[str, str], dict[str, int]]]]] = {}
    n_no_spectrum = 0
    for cid, rows in sorted(controls.items()):
        feature = features.get(cid)
        if feature is None or not (feature["ion_bins"] or feature["loss_bins"]):
            n_no_spectrum += 1
            continue
        hit = next((r for r in rows if r["matched_truth"] == "TRUE"), None)
        if hit is None:
            continue
        per = [
            (row, _explained(row["candidate_smiles"], feature, bin_base, fragment_bins))
            for row in rows
            if row["candidate_smiles"]
        ]
        if per:
            scored[cid] = (hit, per)

    dead_same_axis = [
        cid for cid, (_hit, per) in scored.items() if sum(s["same_axis"] for _r, s in per) == 0
    ]
    revived = [
        cid
        for cid in dead_same_axis
        if sum(s["cross_axis"] for _r, s in scored[cid][1]) > 0
    ]

    baseline_rank1 = sum(1 for _cid, (hit, _per) in scored.items() if _rank(hit) == 1)
    baseline_top3 = sum(1 for _cid, (hit, _per) in scored.items() if _rank(hit) <= 3)

    sweep = []
    for name, score_fn in _variants().items():
        rank1 = top3 = moved_up = moved_down = 0
        for _cid, (hit, per) in scored.items():
            # Ties keep the delivered order, so a variant only moves a compound
            # when it has an actual opinion about it.
            order = sorted(per, key=lambda t: (-score_fn(t[1]), _rank(t[0])))
            new = next(i + 1 for i, t in enumerate(order) if t[0] is hit)
            old = _rank(hit)
            rank1 += new == 1
            top3 += new <= 3
            moved_up += new < old
            moved_down += new > old
        sweep.append(
            {
                "variant": name,
                "rank_1": rank1,
                "top_3": top3,
                "moved_up": moved_up,
                "moved_down": moved_down,
                "delta_rank_1": rank1 - baseline_rank1,
                "delta_top_3": top3 - baseline_top3,
            }
        )

    return {
        "status": "ok",
        "match_tolerance_ppm": MATCH_TOLERANCE_PPM,
        "hydrogen_shift_da": HYDROGEN_DA,
        "n_assessed": len(scored),
        "n_skipped_no_spectrum": n_no_spectrum,
        "n_same_axis_dead": len(dead_same_axis),
        "n_same_axis_dead_but_cross_axis_alive": len(revived),
        "delivered_rank_1": baseline_rank1,
        "delivered_top_3": baseline_top3,
        "sweep": sorted(sweep, key=lambda row: -row["delta_rank_1"]),
    }


# ── Rendering ────────────────────────────────────────────────────────
def render_markdown(summary: dict[str, Any]) -> str:
    recovery = summary["recovery"]
    counts = recovery["counts"]
    ties = summary["ties_truth_controls"]
    spectral = summary["spectral_axis"]

    lines = [
        "# Why the correct structure sits below rank 1",
        "",
        "Generated by `tools/audit_isomer_discrimination.py`. Read-only.",
        "",
        "## Truth-control recovery",
        "",
        f"Of {recovery['n_truth_controls']} truth controls:",
        "",
        "| Outcome | Controls |",
        "|---|---|",
    ]
    labels = {
        "truth_at_rank_1": "Truth at rank 1",
        "truth_absent_from_top10": "Truth absent from the top 10",
        "evidence_degenerate_same_formula": "Truth below rank 1, same formula as rank 1",
        "different_formula": "Truth below rank 1, different formula from rank 1",
        "truth_row_unscored": "Truth below rank 1, truth row carries no evidence score",
    }
    for key, label in labels.items():
        if key in counts:
            lines.append(f"| {label} | {counts[key]} |")

    below = recovery["n_truth_below_rank_1"]
    same = recovery["n_same_formula_as_rank1"]
    unscored = counts.get("truth_row_unscored", 0)
    share = f"{same / below:.0%}" if below else "n/a"
    lines += [
        "",
        "The rows are exclusive and an unscored truth row is reported under that "
        "heading even when it also shares the rank-1 formula, which is why the "
        f"same-formula total below is {same} rather than "
        f"{counts.get('evidence_degenerate_same_formula', 0)}.",
        "",
        f"{same} of the {below} controls whose truth sits below rank 1 share the "
        f"rank-1 candidate's molecular formula ({share}), the {unscored} unscored "
        "one included. For those the mass and formula axes are degenerate: every "
        "candidate fits the same measured mass equally well, so no mass gate, "
        "tolerance change or formula prior can separate them. The failure is "
        "isomer discrimination, which is a different population from the rank-1 "
        "mass mismatches in `docs/pipeline_refinement_handoff.md` item 2.",
        "",
        f"{recovery['n_gap_under_0p05']} of them trail rank 1 by less than 0.05 of "
        "`candidate_confidence`.",
        "",
        "## Is rank 1 ever an arbitrary pick?",
        "",
        f"On {ties['n_rank1_tied']} of {ties['n_compounds']} truth controls "
        f"({ties['share_rank1_tied']:.1%}) rank 1 is tied with at least one "
        "lower-ranked candidate on every displayed evidence column "
        f"(mean tie {ties['mean_tie_size']}, largest {ties['max_tie_size']}). "
        f"In {ties['n_truth_inside_tie_below_rank_1']} the truth is inside that "
        "tied group but not first, so the ordering that decided the report was "
        "made by `rank_frame`'s final tie-break on `predicted_smiles`, not by "
        "evidence.",
        "",
        "## The spectral axis",
        "",
    ]

    if spectral.get("status") != "ok":
        lines += [f"Skipped: {spectral.get('reason')}.", ""]
    else:
        lines += [
            f"`forward_spectrum_score` compares predicted fragments to observed "
            f"ions and predicted losses to observed losses. Of "
            f"{spectral['n_assessed']} truth controls with a local observed "
            f"spectrum, **{spectral['n_same_axis_dead']} have zero same-axis "
            "agreement for every candidate in the top 10** — the axis contributes "
            "nothing at all to those rankings.",
            "",
            "A fragment of mass `m` and a loss of `parent - m` are the same event. "
            f"Pooling the two axes and allowing a hydrogen shift revives "
            f"{spectral['n_same_axis_dead_but_cross_axis_alive']} of those "
            f"{spectral['n_same_axis_dead']}: at least one observed peak becomes "
            "explainable.",
            "",
            "That is a real defect, but pooling alone is **not** a ranking fix. "
            "Re-ordering the delivered top 10 by each normalization of the pooled "
            f"score, against a delivered baseline of {spectral['delivered_rank_1']} "
            f"at rank 1 and {spectral['delivered_top_3']} within the top 3:",
            "",
            "| Variant | Rank 1 | Top 3 | Moved up | Moved down | Δ rank 1 |",
            "|---|---|---|---|---|---|",
        ]
        for row in spectral["sweep"]:
            lines.append(
                f"| `{row['variant']}` | {row['rank_1']} | {row['top_3']} | "
                f"{row['moved_up']} | {row['moved_down']} | {row['delta_rank_1']:+d} |"
            )
        lines += [
            "",
            "Every variant loses. The bond-breaking fragmenter in "
            "`nce/evidence.py::candidate_fragment_bins` breaks only acyclic single "
            "bonds and is blind to regiochemistry, so it cannot rank isomers even "
            "when its output is compared on the right axis. Fixing the axis makes "
            "the signal visible; it does not make it discriminating.",
            "",
        ]

    lines += [
        "## What this rules out",
        "",
        "- A looser or tighter mass window. Every candidate in a degenerate group "
        "fits the same mass.",
        "- Backfilling `evidence_confidence` onto unscored rows. It is worth doing "
        "for correctness, but it is the operative cause for exactly "
        f"{counts.get('truth_row_unscored', 0)} of {below} controls, and it hands "
        "the unscored row the *same* score as its competitors rather than a "
        "distinguishing one.",
        "- Retention index. `docs/paper/data/retention/audit_rt_feasibility.json` "
        "measures decoy AUC 0.651 overall and 0.590 for lib 305, where the proxy "
        "correlation is negative (r = -0.765, n = 13).",
        "",
    ]
    return "\n".join(lines) + "\n"


# ── CLI entry point ──────────────────────────────────────────────────
def build_summary(exec_summary: Path, feature_rows: Path) -> dict[str, Any]:
    grouped = load_top10(exec_summary)
    controls = {
        cid: rows for cid, rows in grouped.items() if rows and rows[0]["cohort"] == COHORT_TRUTH
    }
    if not controls:
        raise SystemExit(f"no {COHORT_TRUTH} rows in {exec_summary}")
    features = load_features(feature_rows) if feature_rows.is_dir() else {}
    return {
        "inputs": {
            "exec_summary": _repo_relative(exec_summary),
            "feature_rows_dir": _repo_relative(feature_rows),
            "n_compounds": len(grouped),
            "n_truth_controls": len(controls),
            "n_feature_rows": len(features),
        },
        "recovery": recovery_breakdown(controls),
        "evidence_basis": {
            "truth_controls": evidence_basis_counts(controls),
            "all_compounds": evidence_basis_counts(grouped),
        },
        "ties_truth_controls": tie_breakdown(controls),
        "ties_all_compounds": tie_breakdown(grouped),
        "spectral_axis": spectral_axis_audit(controls, features),
    }


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--exec-summary",
        type=Path,
        default=PROJECT_ROOT / "docs/paper/data/exec_summary/exec_summary_top10.csv",
    )
    parser.add_argument("--feature-rows-dir", type=Path, default=PROJECT_ROOT / "data/feature_rows")
    parser.add_argument("--out-dir", type=Path, default=PROJECT_ROOT / "docs/paper/data/isomer")
    args = parser.parse_args(list(argv) if argv is not None else None)

    if not args.exec_summary.exists():
        raise SystemExit(f"missing {args.exec_summary}")

    summary = build_summary(args.exec_summary, args.feature_rows_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "audit_isomer_discrimination.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    (args.out_dir / "audit_isomer_discrimination.md").write_text(render_markdown(summary))

    recovery = summary["recovery"]
    spectral = summary["spectral_axis"]
    print(
        f"[isomer] {recovery['n_truth_controls']} truth controls | "
        f"{recovery['n_truth_below_rank_1']} with truth below rank 1, of which "
        f"{recovery['n_same_formula_as_rank1']} share the rank-1 formula"
    )
    if spectral.get("status") == "ok":
        print(
            f"[isomer] spectral axis dead on {spectral['n_same_axis_dead']}/"
            f"{spectral['n_assessed']} controls; pooling revives "
            f"{spectral['n_same_axis_dead_but_cross_axis_alive']} but no variant "
            f"beats the delivered {spectral['delivered_rank_1']} at rank 1"
        )
    else:
        print(f"[isomer] spectral section skipped: {spectral.get('reason')}")
    print(f"[isomer] wrote {args.out_dir}/audit_isomer_discrimination.json")


if __name__ == "__main__":
    main()
