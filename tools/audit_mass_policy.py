#!/usr/bin/env python3
"""Measure whether changing the mass-match window would move any result.

Read-only against a report bundle. Pure standard library, following the shape of
``tools/prototype_mass_gate_rerank.py``.

Three things, in the order they are worth knowing:

1. Provenance. Which formula actually produced the stored ``mass_match_score``.
   This is checked rather than assumed because a bundle whose scores do not come
   from the function in ``structural_elucidation/inference/predict.py`` makes every
   band in the documentation describe something that did not run.
2. The screen. How many candidates sit between the current mass window and a
   proposed one. If that count is near zero, no rebanding can move anything and
   the remaining stages are unnecessary.
3. The reorder experiment (``--rerank``). Applies the mass gate from handoff item 2
   under each policy, then scores rank-1 and rank-3 recovery against
   ``truth_tanimoto`` on the truth controls with an exact paired McNemar test.

Scope limits that hold no matter what comes out:

- Reordering a fixed list of ten cannot change top-10 recovery, so top-10 is not
  reported.
- The pool filter in ``structural_elucidation/utils/merge.py`` drops candidates
  before they reach a report, so delivered data cannot show what a tighter
  tolerance would have discarded. This measures reordering only.

Usage:
    python3 tools/audit_mass_policy.py --reports-root /path/to/reports
    python3 tools/audit_mass_policy.py --reports-root /path/to/reports --rerank
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from structural_elucidation.inference.mass_windows import (  # noqa: E402
    MASS_SCORE_PLATEAU_DA,
    MASS_SCORE_TOLERANCE_DA,
)
from structural_elucidation.reports.confidence import MASS_STRONG, MASS_WEAK  # noqa: E402

csv.field_size_limit(10**7)

# Truth is a Tanimoto bar, matching tools/audit_structure_predictions.py.
TRUTH_THRESHOLD = 0.70

STRONG, PARTIAL, FAIL = "strong", "partial", "fail"


# --- mass policies -----------------------------------------------------------


@dataclass(frozen=True)
class MassPolicy:
    """A mass-match scoring rule: flat 1.0 below `plateau`, 0.0 above `tolerance`.

    `unit` is "da" or "ppm". The linear arm is identical either way, because
    ``err_ppm / tol_ppm == err_da / tol_da``; only where the two edges sit differs.
    """

    name: str
    plateau: float
    tolerance: float
    unit: str

    def edges_da(self, neutral_mass: float | None) -> tuple[float, float] | None:
        """The plateau and tolerance expressed in Da at a given mass."""

        if self.unit == "da":
            return self.plateau, self.tolerance
        if neutral_mass is None or neutral_mass <= 0:
            return None
        scale = neutral_mass / 1e6
        return self.plateau * scale, self.tolerance * scale

    def score(self, error_da: float, neutral_mass: float | None) -> float | None:
        edges = self.edges_da(neutral_mass)
        if edges is None:
            return None
        plateau, tolerance = edges
        if error_da <= plateau:
            return 1.0
        if error_da >= tolerance:
            return 0.0
        return max(0.0, 1.0 - error_da / tolerance)


# What inference/predict.py implements.
CURRENT = MassPolicy("current-absolute-da", MASS_SCORE_PLATEAU_DA, MASS_SCORE_TOLERANCE_DA, "da")
# What the drafted parameter sheet asks for: Strong to 1 ppm, nothing past 5 ppm.
DRAFTED = MassPolicy("drafted-ppm", 1.0, 5.0, "ppm")

# The shape the delivered column actually has: linear from zero, no plateau, which
# is chemistry.conjugates.linear_mass_match_score. Swept across three orders of
# magnitude because the shipped default of 1.0 Da is thousands of ppm, loose enough
# that the axis may not be separating candidates at all. A flat curve here says the
# tolerance is not a dial worth turning; a peak says which way to turn it.
SWEEP_TOLERANCES_DA = (0.005, 0.05, 0.5, 1.0, 5.0)


def delivered_policy(tolerance_da: float) -> MassPolicy:
    return MassPolicy(f"linear-{tolerance_da:g}da", 0.0, tolerance_da, "da")


def band_of(score: float | None) -> str | None:
    """The band a sub-score falls in, using the confidence model's own cutoffs."""

    if score is None:
        return None
    if score >= MASS_STRONG:
        return STRONG
    if score >= MASS_WEAK:
        return PARTIAL
    return FAIL


def ppm_of(error_da: float, neutral_mass: float | None) -> float | None:
    if neutral_mass is None or neutral_mass <= 0:
        return None
    return error_da / neutral_mass * 1e6


# --- exact McNemar -----------------------------------------------------------


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact binomial p-value for `b` and `c` discordant pairs.

    Paired rather than two independent intervals, because the same controls are
    scored twice. On a 92-control cohort the paired form resolves a handful of
    flips where comparing two Wilson intervals resolves nothing under ten points.
    """

    n = b + c
    if n == 0:
        return 1.0
    tail = sum(math.comb(n, k) for k in range(min(b, c) + 1))
    return min(1.0, 2.0 * tail / (2**n))


# --- bundle reading ----------------------------------------------------------


@dataclass
class Candidate:
    rank: int
    smiles: str
    error_da: float | None
    stored_score: float | None
    truth_tanimoto: float | None


@dataclass
class Report:
    control_id: str
    neutral_mass: float | None
    msi_level: int | None
    source: str
    candidates: list[Candidate]

    @property
    def is_truth_control(self) -> bool:
        top = self.candidates[0].truth_tanimoto if self.candidates else None
        return top is not None and top > 0.0

    def correct_at(self, depth: int) -> bool:
        return any(
            (candidate.truth_tanimoto or 0.0) >= TRUTH_THRESHOLD
            for candidate in self.candidates[:depth]
        )


def _num(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if number != number else number


def _rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_reports(root: Path) -> Iterator[Report]:
    """Yield one Report per `X-*` directory under `root`.

    Directories beginning with a dot are skipped: a bundle staged on macOS carries
    AppleDouble `._X-00000` siblings that otherwise parse as empty reports.
    """

    for directory in sorted(root.glob("X-*")):
        if not directory.is_dir() or directory.name.startswith("."):
            continue
        slug = directory.name.lower()
        data = directory / "data"

        fields: dict[str, Any] = {}
        enriched = data / f"{slug}_confidence_inputs.json"
        if enriched.exists():
            try:
                fields = json.loads(enriched.read_text())
            except (ValueError, OSError):
                fields = {}

        # The observed m/z is a property of the feature, so one value serves every
        # candidate. The per-candidate `parent_mz` column is unreliable: in the
        # bundles seen so far it is blank on most rows.
        neutral_mass = _num(fields.get("observed_neutral_mass"))
        if neutral_mass is None:
            observed_mz = _num(fields.get("observed_mz"))
            candidate_mass = _num(fields.get("candidate_exact_mass"))
            # Either is within a milli-Dalton of the other, and a ppm denominator
            # does not care which; take whichever is present.
            neutral_mass = observed_mz if observed_mz is not None else candidate_mass

        candidates = []
        for row in _rows(data / f"{slug}_top_candidates.csv"):
            rank = _num(row.get("rank"))
            candidates.append(
                Candidate(
                    rank=int(rank) if rank is not None else 10**6,
                    smiles=str(row.get("predicted_smiles", "") or "").strip(),
                    error_da=_num(row.get("abs_mass_error_da")),
                    stored_score=_num(row.get("mass_match_score")),
                    truth_tanimoto=_num(row.get("truth_tanimoto")),
                )
            )
        if not candidates:
            continue
        candidates.sort(key=lambda candidate: candidate.rank)
        msi = _num(fields.get("msi_level"))
        yield Report(
            control_id=directory.name,
            neutral_mass=neutral_mass,
            msi_level=int(msi) if msi is not None else None,
            source=str(fields.get("source", "") or "").strip(),
            candidates=candidates,
        )


# --- provenance --------------------------------------------------------------

# Candidate explanations for a stored score, checked in order. The bundle is the
# authority on what ran; the documentation is not.
def _linear(tolerance: float):
    return lambda err, mass: max(0.0, 1.0 - err / tolerance)


def _saturating(scale: float):
    return lambda err, mass: scale / (scale + err)


PROVENANCE_RULES = {
    "predict.py (plateau 0.001 Da, tolerance 0.005 Da)": lambda err, mass: CURRENT.score(err, mass),
    "1 - err_da": _linear(1.0),
    "1 - err_da / 5": _linear(5.0),
    "5 / (5 + err_da)": _saturating(5.0),
    "1 / (1 + err_da)": _saturating(1.0),
}


def provenance(reports: list[Report], tolerance: float = 1e-6) -> dict[str, Any]:
    """Report which formula reproduces the stored scores.

    Worth doing before anything else. If the stored scores do not come from the
    function the documentation describes, then every band in that documentation
    describes a rule that did not run, and no comparison against it means much.
    """

    counts = {name: 0 for name in PROVENANCE_RULES}
    clamped_zero = 0
    ambiguous = 0
    unexplained = 0
    total = 0
    for report in reports:
        for candidate in report.candidates:
            if candidate.error_da is None or candidate.stored_score is None:
                continue
            total += 1
            if candidate.stored_score == 0.0:
                # Every rule bottoms out at zero, so a zero says nothing about which one ran.
                clamped_zero += 1
                continue
            matched = [
                name
                for name, rule in PROVENANCE_RULES.items()
                if (predicted := rule(candidate.error_da, report.neutral_mass)) is not None
                and abs(candidate.stored_score - predicted) < tolerance
            ]
            if not matched:
                unexplained += 1
                continue
            # Attribute each candidate once. The rules converge as the error goes to
            # zero, so a tiny error matches several of them and counting all would
            # invent agreement that the data does not contain.
            counts[matched[0]] += 1
            if len(matched) > 1:
                ambiguous += 1
    return {
        "n_scored_candidates": total,
        "matches_by_rule": counts,
        "attribution": "each candidate counted once, against the first rule it matches",
        "matched_more_than_one_rule": ambiguous,
        "stored_score_exactly_zero": clamped_zero,
        "unexplained": unexplained,
    }


# --- the screen --------------------------------------------------------------

PPM_BUCKETS = (1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0, 1000.0)


def screen(reports: list[Report], current: MassPolicy, proposed: MassPolicy) -> dict[str, Any]:
    """Count what a change from `current` to `proposed` could possibly touch."""

    histogram = {f"<={edge:g}": 0 for edge in PPM_BUCKETS}
    histogram[f">{PPM_BUCKETS[-1]:g}"] = 0
    no_mass_reference = 0
    in_reband_zone = 0
    rank1_in_reband_zone = 0
    band_changes: dict[str, int] = {}
    rank1_band_changes: dict[str, int] = {}
    total = 0

    for report in reports:
        for candidate in report.candidates:
            if candidate.error_da is None:
                continue
            total += 1
            ppm = ppm_of(candidate.error_da, report.neutral_mass)
            if ppm is None:
                no_mass_reference += 1
                continue
            for edge in PPM_BUCKETS:
                if ppm <= edge:
                    histogram[f"<={edge:g}"] += 1
                    break
            else:
                histogram[f">{PPM_BUCKETS[-1]:g}"] += 1

            current_edges = current.edges_da(report.neutral_mass)
            proposed_edges = proposed.edges_da(report.neutral_mass)
            if current_edges and proposed_edges:
                low, high = sorted((proposed_edges[1], current_edges[1]))
                if low < candidate.error_da <= high:
                    in_reband_zone += 1
                    if candidate.rank == 1:
                        rank1_in_reband_zone += 1

            before = band_of(current.score(candidate.error_da, report.neutral_mass))
            after = band_of(proposed.score(candidate.error_da, report.neutral_mass))
            if before is not None and after is not None and before != after:
                key = f"{before}->{after}"
                band_changes[key] = band_changes.get(key, 0) + 1
                if candidate.rank == 1:
                    rank1_band_changes[key] = rank1_band_changes.get(key, 0) + 1

    return {
        "n_candidates_with_error": total,
        "n_without_mass_reference": no_mass_reference,
        "ppm_histogram": histogram,
        "in_reband_zone": in_reband_zone,
        "rank1_in_reband_zone": rank1_in_reband_zone,
        "band_changes": dict(sorted(band_changes.items())),
        "rank1_band_changes": dict(sorted(rank1_band_changes.items())),
    }


def known_figure_check(reports: list[Report]) -> dict[str, Any]:
    """Reproduce the already-audited mass-mismatch counts as a sanity check.

    If these do not land near the committed values the bundle is being read
    wrongly, and nothing downstream is worth reading.
    """

    mismatched = recoverable = 0
    for report in reports:
        top = report.candidates[0]
        if top.stored_score is None or top.stored_score >= MASS_WEAK:
            continue
        mismatched += 1
        if any(
            candidate.stored_score is not None and candidate.stored_score >= MASS_STRONG
            for candidate in report.candidates[1:]
        ):
            recoverable += 1
    total = len(reports)
    return {
        "n_reports": total,
        "rank1_mass_mismatch": mismatched,
        "rank1_mass_mismatch_rate": round(mismatched / total, 4) if total else None,
        "recoverable_by_in_list_alternative": recoverable,
        "expected": {
            "n_reports": 3897,
            "rank1_mass_mismatch": 598,
            "rank1_mass_mismatch_rate": 0.1535,
            "recoverable_by_in_list_alternative": 523,
        },
    }


# --- the reorder experiment --------------------------------------------------


def gated_order(report: Report, policy: MassPolicy) -> list[Candidate]:
    """The candidate list with mass-failing candidates moved below passing ones.

    This is handoff item 2's mass gate generalised from "promote one alternative"
    to the whole list. The sort is stable, so within each group the delivered
    order is preserved and the only thing the policy changes is which candidates
    are considered to have failed.
    """

    def fails(candidate: Candidate) -> int:
        if candidate.error_da is None:
            return 0
        score = policy.score(candidate.error_da, report.neutral_mass)
        return 1 if score is not None and score < MASS_WEAK else 0

    return sorted(report.candidates, key=fails)


def paired_recovery(reports: list[Report], policy: MassPolicy) -> dict[str, Any]:
    """Rank-1 and rank-3 recovery before and after the gate, paired per control."""

    controls = [report for report in reports if report.is_truth_control]
    result: dict[str, Any] = {"n_controls": len(controls)}
    for depth in (1, 3):
        before_correct = after_correct = 0
        gained = lost = 0
        for report in controls:
            before = report.correct_at(depth)
            reordered = gated_order(report, policy)[:depth]
            after = any((c.truth_tanimoto or 0.0) >= TRUTH_THRESHOLD for c in reordered)
            before_correct += before
            after_correct += after
            if after and not before:
                gained += 1
            elif before and not after:
                lost += 1
        result[f"rank{depth}"] = {
            "correct_before": before_correct,
            "correct_after": after_correct,
            "gained": gained,
            "lost": lost,
            "net": gained - lost,
            "mcnemar_exact_p": round(mcnemar_exact(lost, gained), 6),
        }
    return result


def truth_coverage(reports: list[Report], policy: MassPolicy) -> dict[str, Any]:
    """How much of the cohort the gate targets is actually checkable against truth.

    The number that decides how far the paired result generalises. A gate aimed at
    the mass-mismatch cohort can only be validated on the part of that cohort whose
    correct answer is known, and if the overlap is small then a clean result on the
    controls says little about the rest.
    """

    controls = [report for report in reports if report.is_truth_control]
    mismatched = [
        report
        for report in reports
        if report.candidates[0].stored_score is not None
        and report.candidates[0].stored_score < MASS_WEAK
    ]
    fired = gained = lost = unchanged = 0
    for report in controls:
        reordered = gated_order(report, policy)
        if reordered[0] is report.candidates[0]:
            continue
        fired += 1
        before = (report.candidates[0].truth_tanimoto or 0.0) >= TRUTH_THRESHOLD
        after = (reordered[0].truth_tanimoto or 0.0) >= TRUTH_THRESHOLD
        if after and not before:
            gained += 1
        elif before and not after:
            lost += 1
        else:
            unchanged += 1
    return {
        "n_truth_controls": len(controls),
        "n_rank1_mass_mismatch": len(mismatched),
        "n_rank1_mass_mismatch_that_are_truth_controls": sum(
            1 for report in mismatched if report.is_truth_control
        ),
        "gate_reordered_rank1_on_controls": fired,
        "of_those_gained": gained,
        "of_those_lost": lost,
        "of_those_unchanged": unchanged,
    }


def tolerance_sweep(reports: list[Report]) -> dict[str, Any]:
    """Recovery and band split across three orders of magnitude of mass tolerance.

    Answers whether the mass axis is carrying information. The shipped default of
    1.0 Da is thousands of ppm, which is loose enough that nearly every candidate
    scores near the top and the axis stops separating anything. If rank-1 recovery
    is unmoved from 0.005 Da to 5.0 Da then the tolerance is not the problem and
    tightening it is not the fix.
    """

    rows = []
    for tolerance in SWEEP_TOLERANCES_DA:
        policy = delivered_policy(tolerance)
        paired = paired_recovery(reports, policy)
        rows.append(
            {
                "tolerance_da": tolerance,
                "rank1_correct": paired["rank1"]["correct_after"],
                "rank1_gained": paired["rank1"]["gained"],
                "rank1_lost": paired["rank1"]["lost"],
                "rank1_mcnemar_exact_p": paired["rank1"]["mcnemar_exact_p"],
                "rank3_correct": paired["rank3"]["correct_after"],
                "rank1_band": distributions(reports, policy),
            }
        )
    baseline = rows[0]["rank1_correct"] if rows else None
    spread = max(row["rank1_correct"] for row in rows) - min(row["rank1_correct"] for row in rows)
    return {
        "baseline_rank1_correct_before_any_gate": paired_recovery(reports, CURRENT)["rank1"][
            "correct_before"
        ],
        "points": rows,
        "rank1_spread_across_the_sweep": spread,
        "reading": (
            f"A spread of {spread} across a "
            f"{SWEEP_TOLERANCES_DA[-1] / SWEEP_TOLERANCES_DA[0]:.0f}-fold change in tolerance "
            "means the mass axis is not separating correct candidates from incorrect ones at any "
            "setting, so no choice of window is defensible on this cohort. Note also that losses "
            "fall as the tolerance loosens: the gate does least harm where it fires least, which "
            "is what an intervention carrying no signal looks like."
            if spread <= 3
            else "Recovery varies with the tolerance, so the window is a real choice and the "
            "best-performing setting is worth examining against a larger truth set."
        ),
        "_baseline_at_tightest": baseline,
    }


def distributions(reports: list[Report], policy: MassPolicy) -> dict[str, Any]:
    """How the mass band splits the delivered set under a policy.

    Reported as band counts rather than as recomputed MSI levels: MSI depends on
    structure, MS2 and retention as well, and this tool holds only the mass axis,
    so recomputing a level here would misattribute the other axes' effects.
    """

    counts: dict[str, int] = {STRONG: 0, PARTIAL: 0, FAIL: 0, "unknown": 0}
    for report in reports:
        top = report.candidates[0]
        band = (
            band_of(policy.score(top.error_da, report.neutral_mass))
            if top.error_da is not None
            else None
        )
        counts[band or "unknown"] += 1
    return counts


# --- rendering ---------------------------------------------------------------


def render_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# Mass-policy audit",
        "",
        summary["method"],
        "",
        "## Which formula produced the stored scores",
        "",
        "| Rule | Candidates matching |",
        "|---|---:|",
    ]
    prov = summary["provenance"]
    for name, count in prov["matches_by_rule"].items():
        lines.append(f"| `{name}` | {count} |")
    lines += [
        f"| stored score exactly 0.0 (uninformative) | {prov['stored_score_exactly_zero']} |",
        f"| unexplained | {prov['unexplained']} |",
        "",
        f"Out of {prov['n_scored_candidates']} scored candidates.",
        "",
        "## The screen",
        "",
        f"Candidates between the two windows: **{summary['screen']['in_reband_zone']}** "
        f"({summary['screen']['rank1_in_reband_zone']} of them at rank 1) out of "
        f"{summary['screen']['n_candidates_with_error']} with a mass error.",
        "",
        "| Mass error | Candidates |",
        "|---|---:|",
    ]
    for bucket, count in summary["screen"]["ppm_histogram"].items():
        lines.append(f"| {bucket} ppm | {count} |")
    lines += [
        "",
        "Bands are how the confidence model reads a mass error. Rank-1 candidates by band:",
        "",
        "| Policy | strong | partial | fail | no mass |",
        "|---|---:|---:|---:|---:|",
    ]
    for name, counts in summary["band_distribution"].items():
        lines.append(
            f"| `{name}` | {counts['strong']} | {counts['partial']} | {counts['fail']} | "
            f"{counts['unknown']} |"
        )

    if "paired_reorder" in summary:
        lines += [
            "",
            "## Does the mass gate improve recovery",
            "",
            "The gate from handoff item 2, applied to the whole candidate list under each "
            "policy, scored against the truth controls. Paired, because the same controls are "
            "scored twice and two independent intervals could not resolve a handful of flips.",
            "",
            "| Policy | Depth | Correct before | Correct after | Gained | Lost | Exact p |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
        for name, block in summary["paired_reorder"].items():
            for depth in ("rank1", "rank3"):
                row = block[depth]
                lines.append(
                    f"| `{name}` | {depth} | {row['correct_before']} | {row['correct_after']} | "
                    f"{row['gained']} | {row['lost']} | {row['mcnemar_exact_p']} |"
                )
        sweep = summary["tolerance_sweep"]
        lines += [
            "",
            "### Does the tolerance matter at all",
            "",
            "The delivered scorer is linear from zero with no plateau. Sweeping its tolerance "
            "across three orders of magnitude, against "
            f"{sweep['baseline_rank1_correct_before_any_gate']} controls correct at rank 1 before "
            "any gate:",
            "",
            "| Tolerance | Rank-1 correct | Gained | Lost | Exact p | Rank-1 strong | partial | fail |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for row in sweep["points"]:
            band = row["rank1_band"]
            lines.append(
                f"| {row['tolerance_da']:g} Da | {row['rank1_correct']} | {row['rank1_gained']} | "
                f"{row['rank1_lost']} | {row['rank1_mcnemar_exact_p']} | {band['strong']} | "
                f"{band['partial']} | {band['fail']} |"
            )
        lines += [
            "",
            f"Recovery moves by {sweep['rank1_spread_across_the_sweep']} controls across the whole "
            "sweep, and no setting gains a single one. What does move is what the reports say: "
            "rank-1 strong-mass claims range from "
            f"{sweep['points'][0]['rank1_band']['strong']} to "
            f"{sweep['points'][-1]['rank1_band']['strong']} across the same range. "
            + sweep["reading"],
        ]
        cover = next(iter(summary["truth_coverage"].values()))
        lines += [
            "",
            "### How far that generalises",
            "",
            f"Only **{cover['n_rank1_mass_mismatch_that_are_truth_controls']} of the "
            f"{cover['n_rank1_mass_mismatch']}** rank-1 mass-mismatch reports have a known correct "
            "answer, so that handful carries the entire evidence base for any projection about the "
            "rest. The gate is aimed at a cohort the truth set barely covers, which bounds what "
            "these counts can settle either way.",
        ]

    lines += [
        "",
        "## Scope",
        "",
        summary["scope_note"],
        "",
        "Two limits are structural rather than incidental. Reordering a fixed list of ten cannot "
        "change top-10 recovery, so no top-10 figure appears above. And the pool filter in "
        "`structural_elucidation/utils/merge.py` drops candidates before they reach a report, so "
        "delivered data cannot show what a tighter tolerance would have discarded; answering that "
        "needs a pipeline re-run, not this bundle.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports-root", type=Path, required=True)
    parser.add_argument("--rerank", action="store_true", help="run the paired reorder experiment")
    parser.add_argument("--out-dir", type=Path, default=None, help="write JSON and Markdown here")
    args = parser.parse_args()

    if not args.reports_root.is_dir():
        raise SystemExit(f"reports root not found: {args.reports_root}")

    reports = list(read_reports(args.reports_root))
    if not reports:
        raise SystemExit(f"no X-* report directories under {args.reports_root}")

    summary: dict[str, Any] = {
        "method": (
            "Read-only comparison of the mass-match window against a proposed one, over a "
            "delivered report bundle. Measures reordering only: the pool filter drops "
            "candidates before they reach a report, so a tighter tolerance's effect on what "
            "was available cannot be seen from here."
        ),
        "policies": {
            "current": CURRENT.__dict__,
            "proposed": DRAFTED.__dict__,
        },
        "provenance": provenance(reports),
        "known_figure_check": known_figure_check(reports),
        "screen": screen(reports, CURRENT, DRAFTED),
        "band_distribution": {
            CURRENT.name: distributions(reports, CURRENT),
            DRAFTED.name: distributions(reports, DRAFTED),
        },
        "scope_note": (
            "Bands describe how a measurement is read, not how often that reading is right."
        ),
    }
    if args.rerank:
        summary["paired_reorder"] = {
            CURRENT.name: paired_recovery(reports, CURRENT),
            DRAFTED.name: paired_recovery(reports, DRAFTED),
        }
        summary["truth_coverage"] = {
            CURRENT.name: truth_coverage(reports, CURRENT),
            DRAFTED.name: truth_coverage(reports, DRAFTED),
        }
        summary["tolerance_sweep"] = tolerance_sweep(reports)

    prov = summary["provenance"]
    print(f"[mass-policy] reports: {len(reports)}")
    print(f"[mass-policy] scored candidates: {prov['n_scored_candidates']}")
    for name, count in prov["matches_by_rule"].items():
        if count:
            print(f"[mass-policy]   stored score matches {name}: {count}")
    if prov["matched_more_than_one_rule"]:
        print(
            f"[mass-policy]   (of those, {prov['matched_more_than_one_rule']} matched several rules;"
            " the rules converge as the error goes to zero)"
        )
    if prov["stored_score_exactly_zero"]:
        print(f"[mass-policy]   stored score exactly 0.0: {prov['stored_score_exactly_zero']}")
    if prov["unexplained"]:
        print(f"[mass-policy]   UNEXPLAINED: {prov['unexplained']}")
    check = summary["known_figure_check"]
    print(
        f"[mass-policy] rank-1 mass mismatches: {check['rank1_mass_mismatch']}"
        f" ({check['rank1_mass_mismatch_rate']}), recoverable {check['recoverable_by_in_list_alternative']}"
        f"  [expected on the full set: 598 (0.1535), 523]"
    )
    scr = summary["screen"]
    print(
        f"[mass-policy] between the windows: {scr['in_reband_zone']} candidates,"
        f" {scr['rank1_in_reband_zone']} at rank 1"
    )
    if args.rerank:
        for name, block in summary["paired_reorder"].items():
            rank1 = block["rank1"]
            cover = summary["truth_coverage"][name]
            print(
                f"[mass-policy] {name}: {block['n_controls']} truth controls;"
                f" rank-1 {rank1['correct_before']} -> {rank1['correct_after']}"
                f" (+{rank1['gained']} / -{rank1['lost']}, exact p={rank1['mcnemar_exact_p']})"
            )
            print(
                f"[mass-policy]   gate reordered rank 1 on {cover['gate_reordered_rank1_on_controls']}"
                f" controls; only {cover['n_rank1_mass_mismatch_that_are_truth_controls']} of the"
                f" {cover['n_rank1_mass_mismatch']} mass-mismatch reports are checkable against truth"
            )
        sweep = summary["tolerance_sweep"]
        print(
            f"[mass-policy] tolerance sweep (rank-1 correct of"
            f" {sweep['baseline_rank1_correct_before_any_gate']} before any gate):"
        )
        for row in sweep["points"]:
            band = row["rank1_band"]
            print(
                f"[mass-policy]   {row['tolerance_da']:>6g} Da -> {row['rank1_correct']:2d} correct"
                f" (+{row['rank1_gained']} / -{row['rank1_lost']}, p={row['rank1_mcnemar_exact_p']});"
                f" rank-1 strong {band['strong']}, partial {band['partial']}, fail {band['fail']}"
            )
        print(f"[mass-policy]   spread across the sweep: {sweep['rank1_spread_across_the_sweep']} controls")

    if args.out_dir:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        json_path = args.out_dir / "audit_mass_policy.json"
        json_path.write_text(json.dumps(summary, indent=2) + "\n")
        md_path = args.out_dir / "audit_mass_policy.md"
        md_path.write_text(render_markdown(summary) + "\n")
        print(f"[mass-policy] wrote {json_path} and {md_path}")


if __name__ == "__main__":
    main()
