#!/usr/bin/env python3
"""Measure what each candidate cutoff would actually buy, before anyone picks one.

Three policy numbers were unresolved: how many candidate structures a public
report should show, what evidence justifies asserting a name rather than
demoting an MSI level, and whether a threshold on the evidence score can carry
either decision. This computes the operating curve for each of them so the
numbers are chosen from measured precision and coverage rather than by hand.

Three sweeps:

1. ``score_cutoff`` - accuracy of the calls kept at each threshold on the
   leave-one-out identification-confidence probability, plus the rank AUC. This
   is the sweep that tests the "cut at 0.1" idea directly.
2. ``candidate_depth`` - how often the correct structure is among the top N, and
   how many incorrect structures a report must display to buy that. This is the
   public-candidate-count decision.
3. ``msi_rule`` - MSI boundary parameters swept against naming precision. Needs
   per-axis evidence, so it runs only when the controls file carries the axis
   columns or a report bundle is supplied.

Sweeps 1 and 2 read the committed per-control audit output and run in any
checkout. Sweep 3 needs axis-level evidence, which
``tools/audit_structure_predictions.py`` writes into the controls CSV when it is
run against a report bundle.

Usage:
    python3 tools/audit_thresholds.py
    python3 tools/audit_thresholds.py --reports-root _audit_pull/reports
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
import confidence as C  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
# The Platt fit is imported rather than reimplemented so the measurement below
# cannot drift away from the calibration it is measuring.
import calibrate_confidence as CAL  # noqa: E402

# Thresholds on the calibrated evidence score. 0.10 is the value that prompted
# this audit. The tier cutoffs are unioned into the grid rather than written into
# it by hand, because the named points below are looked up by exact membership: a
# retuned cutoff that no longer appeared in a hand-maintained grid would be
# dropped from the artifact silently rather than failing.
BASE_SCORE_GRID = (0.10, 0.20, 0.30, 0.45, 0.55, 0.60, 0.70)
# Keys name the tier, not its current value. They used to embed the number
# ("tier_high_0.65"), so a retune left three keys in the committed JSON - and in
# the test that compares against it - asserting cutoffs that no longer existed.
NAMED_SCORE_POINTS = {
    "proposed_0.10": 0.10,
    "tier_low": C.TIER_LOW,
    "tier_moderate": C.TIER_MODERATE,
    "tier_high": C.TIER_HIGH,
}
SCORE_GRID = tuple(sorted(set(BASE_SCORE_GRID) | {round(v, 4) for v in NAMED_SCORE_POINTS.values()}))


def _assert_named_points_are_on_the_grid() -> None:
    """Fail at import if a named cutoff could drop out of the artifact.

    The sweep reports named operating points by looking them up in the grid, so a
    cutoff missing from it is omitted rather than reported. The union above makes
    that impossible; this is what notices if someone unpicks the union.
    """
    missing = {
        name: value
        for name, value in NAMED_SCORE_POINTS.items()
        if round(float(value), 4) not in SCORE_GRID
    }
    if missing:
        raise AssertionError(
            f"named operating points are off the sweep grid and would be dropped: {missing}"
        )


_assert_named_points_are_on_the_grid()
# Report depths the committed audit records recovery for.
DEPTH_GRID = (1, 3, 10)
# MSI L2 boundary candidates. The L2 rule is the one that decides whether a
# report asserts a reference name, so it is the rule worth sweeping.
STRUCT_GRID = (0.60, 0.70, 0.80, 0.85)
MS2_GRID = (0.30, 0.40, 0.50, 0.60)
MS2_EXACT_GRID = (0.30, 0.40, 0.50, 0.60, 0.70)
MASS_GRID = (0.90, 0.95)

BOOTSTRAP_SEED = 20260805
BOOTSTRAP_N = 2000


def wilson(correct: int, n: int, z: float = 1.959963984540054) -> dict[str, Any]:
    if n <= 0:
        return {"n": 0, "correct": 0, "accuracy": None, "ci95": [None, None]}
    p = correct / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denom
    half = z * math.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n) / denom
    return {
        "n": n,
        "correct": correct,
        "accuracy": round(p, 4),
        "ci95": [round(max(0.0, center - half), 4), round(min(1.0, center + half), 4)],
    }


def _float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def auc(scores: Sequence[float], labels: Sequence[bool]) -> float | None:
    """Rank AUC: the chance a correct call outscores an incorrect one, ties at 0.5.

    0.5 means the score carries no information about correctness.
    """

    n_pos = sum(1 for label in labels if label)
    n_neg = len(labels) - n_pos
    if not n_pos or not n_neg:
        return None
    order = sorted(range(len(scores)), key=lambda index: scores[index])
    ranks = [0.0] * len(scores)
    start = 0
    while start < len(order):
        stop = start
        while stop + 1 < len(order) and scores[order[stop + 1]] == scores[order[start]]:
            stop += 1
        average = (start + stop) / 2.0 + 1.0
        for index in range(start, stop + 1):
            ranks[order[index]] = average
        start = stop + 1
    positive_ranks = sum(rank for rank, label in zip(ranks, labels) if label)
    return (positive_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def bootstrap_auc_ci(scores: Sequence[float], labels: Sequence[bool]) -> list[float | None]:
    if not scores:
        return [None, None]
    rng = random.Random(BOOTSTRAP_SEED)
    size = len(scores)
    values: list[float] = []
    for _ in range(BOOTSTRAP_N):
        picks = [rng.randrange(size) for _ in range(size)]
        sample_scores = [scores[index] for index in picks]
        sample_labels = [labels[index] for index in picks]
        value = auc(sample_scores, sample_labels)
        if value is not None:
            values.append(value)
    if not values:
        return [None, None]
    values.sort()
    low = values[max(0, int(0.025 * len(values)) - 1)]
    high = values[min(len(values) - 1, int(0.975 * len(values)))]
    return [round(low, 4), round(high, 4)]


def score_cutoff_sweep(
    scores: Sequence[float],
    labels: Sequence[bool],
    grid: Iterable[float] = SCORE_GRID,
) -> dict[str, Any]:
    """Accuracy and coverage of the calls kept at each threshold on the score."""

    total = len(scores)
    base_rate = (sum(1 for label in labels if label) / total) if total else None
    points = []
    for threshold in grid:
        kept = [label for score, label in zip(scores, labels) if score >= threshold]
        point = wilson(sum(1 for label in kept if label), len(kept))
        point["threshold"] = round(float(threshold), 4)
        point["coverage"] = round(len(kept) / total, 4) if total else None
        point["lift_over_base_rate"] = (
            round(point["accuracy"] - base_rate, 4)
            if point["accuracy"] is not None and base_rate is not None
            else None
        )
        points.append(point)
    area = auc(scores, labels)
    return {
        "ranked_by": "leave-one-out identification-confidence probability",
        "n": total,
        "base_rate": round(base_rate, 4) if base_rate is not None else None,
        "auc": round(area, 4) if area is not None else None,
        "auc_ci95": bootstrap_auc_ci(scores, labels),
        "grid": points,
        "named_operating_points": _named_points(points, grid),
    }


def _named_points(points: list[dict[str, Any]], grid: Sequence[float]) -> dict[str, Any]:
    """Pull the named cutoffs out of the sweep.

    A named point absent from the grid is omitted, which is right for a caller
    that deliberately passes a narrow grid. What must never happen is the
    *default* grid losing a tier cutoff after a retune, and that is prevented by
    construction: :data:`SCORE_GRID` unions the constants in, and
    :func:`_assert_named_points_are_on_the_grid` checks it at import.
    """
    by_threshold = {point["threshold"]: point for point in points}
    return {
        name: by_threshold[round(float(value), 4)]
        for name, value in NAMED_SCORE_POINTS.items()
        if round(float(value), 4) in by_threshold
    }


def load_msi_levels(path: Path) -> dict[str, str]:
    """Per-compound MSI level, keyed by compound id. Missing file is not an error."""
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        return {
            str(row.get("x_compound_id", "")).strip(): str(row.get("msi_level", "")).strip()
            for row in csv.DictReader(handle)
            if str(row.get("x_compound_id", "")).strip()
        }


def source_stratification(
    rows: list[dict[str, str]], levels: dict[str, str] | None = None
) -> dict[str, Any]:
    """Test whether pooling the two candidate sources is what flattens the score.

    Exact catalog hits and de novo predictions have opposite base rates, and the
    evidence score pins both the structure and the pathway axis to 1.0 for a
    catalog hit. If the pooled curve is flat while the within-source curves are
    not, the score is ordering provenance rather than correctness, and the fix is
    to calibrate per source rather than to move a threshold.
    """

    def measure(subset: list[dict[str, str]], column: str) -> dict[str, Any]:
        scores = [float(row[column]) for row in subset if _float(row.get(column)) is not None]
        labels = [
            _truthy(row.get("top1_correct_0p70", ""))
            for row in subset
            if _float(row.get(column)) is not None
        ]
        area = auc(scores, labels)
        ordered = sorted(scores)
        middle = len(ordered) // 2
        return {
            "n": len(scores),
            "base_rate": round(sum(labels) / len(labels), 4) if labels else None,
            "auc": round(area, 4) if area is not None else None,
            "median_score": round(
                ordered[middle] if len(ordered) % 2 else (ordered[middle - 1] + ordered[middle]) / 2,
                4,
            )
            if ordered
            else None,
        }

    strata = {"pooled": rows}
    for source in sorted({str(row.get("source", "")).strip() for row in rows} - {""}):
        strata[source] = [row for row in rows if str(row.get("source", "")).strip() == source]
    out = {
        stratum: {column: measure(subset, column) for column in ("raw_full", "loo_probability_full")}
        for stratum, subset in strata.items()
    }
    within = [
        block["raw_full"]["auc"]
        for stratum, block in out.items()
        if stratum != "pooled" and block["raw_full"]["auc"] is not None
    ]
    pooled = out["pooled"]["raw_full"]["auc"]
    return {
        "method": (
            "Rank AUC of the raw evidence score and of its calibrated probability, "
            "pooled and within each candidate source."
        ),
        "strata": out,
        "pooled_auc": pooled,
        "within_source_auc_range": [min(within), max(within)] if within else None,
        "confounded": bool(within and pooled is not None and pooled < min(within)),
        "per_source_calibration": per_source_calibration(rows),
        "positive_class_degeneracy": positive_class_degeneracy(rows, levels),
    }


def positive_class_degeneracy(
    rows: list[dict[str, str]], levels: dict[str, str] | None
) -> dict[str, Any]:
    """Does this cohort's "correct" label coincide with MSI Level 1?

    It does, and it invalidates every AUC on this page as evidence about the
    reports a confidence tier is for.

    The 92 controls are Metabolon Level-1 compounds - that is what makes them a
    truth set - and ``assign_msi_level`` returns Level 1 exactly when a Level-1
    standard's top candidate matches it. So "the pipeline got this right" and
    "this report is Level 1" are the same event here, and a score that separates
    the positive class is separating the already-confirmed.

    That is fine for asking whether the pipeline recovers known compounds. It
    cannot answer whether a score orders correctness among the reports nobody has
    a standard for, which is the entire population a High tier would apply to.
    Measured rather than reasoned, because the promising within-source AUCs above
    invite exactly the opposite conclusion.
    """
    if not levels:
        return {
            "status": "not_evaluated",
            "reason": "no per-compound MSI levels supplied; pass --exec-summary to evaluate",
        }
    matched = [row for row in rows if str(row.get("control_id", "")).strip() in levels]
    if not matched:
        return {"status": "not_evaluated", "reason": "no control ids joined to the exec summary"}
    positives = [row for row in matched if _truthy(row.get("top1_correct_0p70", ""))]
    level_1 = [row for row in matched if levels[str(row["control_id"]).strip()].strip() == "1"]
    positive_ids = {str(row["control_id"]).strip() for row in positives}
    level_1_ids = {str(row["control_id"]).strip() for row in level_1}
    non_pinned_positives = len(positive_ids - level_1_ids)
    return {
        "n_joined": len(matched),
        "n_correct": len(positive_ids),
        "n_msi_level_1": len(level_1_ids),
        "correct_and_level_1": len(positive_ids & level_1_ids),
        "correct_outside_level_1": non_pinned_positives,
        "positive_class_is_level_1": bool(positive_ids and positive_ids == level_1_ids),
        "usable_for_non_level_1_reports": bool(non_pinned_positives),
        "note": (
            "Every correct call in this cohort is MSI Level 1, so any AUC here measures "
            "separation of the already-confirmed rather than ordering of correctness. It "
            "cannot support a confidence tier for reports without a standard."
            if positive_ids and positive_ids == level_1_ids
            else f"{non_pinned_positives} correct calls sit outside Level 1, so the cohort "
                 "carries some signal about unconfirmed reports."
        ),
    }


def per_source_calibration(rows: list[dict[str, str]]) -> dict[str, Any]:
    """Does fitting the calibration per source actually recover the lost ordering?

    ``source_stratification`` above reports that pooling flattens the score and
    recommends calibrating per source. This measures whether that works, because
    a recommendation nobody checked is how the pooled fit got shipped.

    Three orderings are compared inside each stratum:

    - ``raw`` - the evidence score, uncalibrated.
    - ``pooled_loo`` - the shipped approach, one calibration over both sources.
    - ``per_source_loo`` - one calibration fitted inside the stratum.

    All three are leave-one-out where they are calibrated at all. That is not a
    nicety: a single Platt fit is monotone, so applying one in-sample would leave
    every rank statistic identical to ``raw`` and the comparison would be
    vacuous. Only a per-point refit produces a genuinely different ordering, and
    only a held-out one says anything about a report the fit has not seen.

    The number that decides the question is ``per_source_beats_pooled``. If
    calibration cannot beat the raw score either, the honest reading is that the
    Platt step is costing signal rather than adding calibration, and the bands
    belong on the raw score.
    """
    usable = [
        row for row in rows
        if _float(row.get("raw_full")) is not None and str(row.get("source", "")).strip()
    ]
    if not usable:
        return {"status": "not_evaluated", "reason": "no rows carry both raw_full and a source"}

    def label(row: dict[str, str]) -> bool:
        return _truthy(row.get("top1_correct_0p70", ""))

    sources = sorted({str(row.get("source", "")).strip() for row in usable})
    pooled_loo = CAL.loo_probabilities(
        [float(row["raw_full"]) for row in usable], [float(label(row)) for row in usable]
    )
    pooled_by_id = {id(row): value for row, value in zip(usable, pooled_loo)}

    strata: dict[str, Any] = {}
    for source in sources:
        subset = [row for row in usable if str(row.get("source", "")).strip() == source]
        raw = [float(row["raw_full"]) for row in subset]
        labels = [label(row) for row in subset]
        own = CAL.loo_probabilities(raw, [float(v) for v in labels])
        fit = CAL.fit_platt(raw, [float(v) for v in labels])

        def scored(values: Sequence[float | None]) -> float | None:
            pairs = [(v, lab) for v, lab in zip(values, labels) if v is not None]
            if not pairs:
                return None
            area = auc([v for v, _ in pairs], [lab for _, lab in pairs])
            return round(area, 4) if area is not None else None

        strata[source] = {
            "n": len(subset),
            "base_rate": round(sum(labels) / len(labels), 4) if labels else None,
            "auc_raw": scored(raw),
            "auc_pooled_loo": scored([pooled_by_id[id(row)] for row in subset]),
            "auc_per_source_loo": scored(own),
            "in_sample_fit": fit.as_dict(),
        }

    def gathered(key: str) -> list[float]:
        return [block[key] for block in strata.values() if block.get(key) is not None]

    per_source = gathered("auc_per_source_loo")
    pooled_auc = gathered("auc_pooled_loo")
    raw_auc = gathered("auc_raw")
    beats_pooled = bool(
        per_source and pooled_auc and len(per_source) == len(pooled_auc)
        and all(
            block["auc_per_source_loo"] > block["auc_pooled_loo"]
            for block in strata.values()
            if block.get("auc_per_source_loo") is not None and block.get("auc_pooled_loo") is not None
        )
    )
    beats_raw = bool(
        per_source and raw_auc
        and all(
            block["auc_per_source_loo"] > block["auc_raw"]
            for block in strata.values()
            if block.get("auc_per_source_loo") is not None and block.get("auc_raw") is not None
        )
    )
    return {
        "method": (
            "Leave-one-out rank AUC inside each candidate source under three orderings: the raw "
            "evidence score, the pooled calibration, and a calibration fitted within the source. "
            "Leave-one-out throughout, because a single Platt fit is monotone and would leave "
            "every rank statistic equal to the raw score."
        ),
        "strata": strata,
        "per_source_beats_pooled": beats_pooled,
        "per_source_beats_raw": beats_raw,
        "recommendation": (
            "calibrate per source" if beats_pooled and beats_raw else
            "band the raw evidence score per source; calibration is not adding ordering"
            if beats_pooled else
            "neither calibration helps on this cohort"
        ),
    }


def candidate_depth_sweep(rows: list[dict[str, str]]) -> dict[str, Any]:
    """How often the correct structure is within the top N, and what N costs.

    Each control has at most one correct structure, so displaying N of them puts
    ``N - hit_rate`` incorrect structures in front of the reader on average. The
    decision-relevant quantity is the exchange rate between the two: how many
    extra incorrect structures a deeper list displays per additional compound
    whose correct structure it surfaces.
    """

    total = len(rows)
    points = []
    for depth in DEPTH_GRID:
        column = f"top{depth}_correct_0p70"
        point = wilson(sum(1 for row in rows if _truthy(row.get(column, ""))), total)
        point["n_shown"] = depth
        point["hit_rate"] = point.pop("accuracy")
        point["incorrect_shown_per_compound"] = (
            round(depth - point["hit_rate"], 4) if point["hit_rate"] is not None else None
        )
        points.append(point)
    for previous, current in zip(points, points[1:]):
        gain = (current["hit_rate"] or 0.0) - (previous["hit_rate"] or 0.0)
        cost = (current["incorrect_shown_per_compound"] or 0.0) - (
            previous["incorrect_shown_per_compound"] or 0.0
        )
        current["hit_rate_gain_vs_previous"] = round(gain, 4)
        current["incorrect_shown_per_additional_hit"] = round(cost / gain, 1) if gain > 0 else None
    return {
        "correctness": "correct structure within the displayed list at truth Tanimoto >= 0.70",
        "n": total,
        "points": points,
    }


def _axis_record(fields: dict) -> dict[str, Any]:
    axes = {axis.key: axis for axis in C.build_axes(fields)}
    return {
        "structure": axes["structure"].score or 0.0,
        "mass": axes["mass"].score if axes["mass"].available else None,
        "ms2": axes["ms2"].score if axes["ms2"].available else None,
        "observed_mz_available": bool(fields.get("observed_mz_available")),
        "source": str(fields.get("source", "")).strip(),
    }


def _msi_level(record: dict[str, Any], thresholds: C.MsiThresholds) -> int:
    """Replay :func:`confidence.assign_msi_level` for one control's axes.

    Builds the axis objects the real rule consumes so the sweep and the delivered
    reports cannot disagree about what a level means.
    """

    fields = {
        "observed_mz_available": record["observed_mz_available"],
        "source": record.get("source", ""),
    }
    axes = [
        C.Axis("structure", "", record["structure"], C.WEIGHTS["structure"], True),
        C.Axis("mass", "", record["mass"], C.WEIGHTS["mass"], record["mass"] is not None),
        C.Axis("ms2", "", record["ms2"], C.WEIGHTS["ms2"], record["ms2"] is not None),
        C.Axis("pathway", "", None, C.WEIGHTS["pathway"], False),
    ]
    conf = C.Confidence(overall=0.0, tier="Low", tier_note="", completeness=0.0, descriptor="", axes=axes)
    level, _label, _rationale = C.assign_msi_level(fields, conf, thresholds)
    return level


def msi_rule_sweep(records: list[dict[str, Any]]) -> dict[str, Any]:
    """Sweep the MSI L2 boundary and score the naming decision it implies.

    A report asserts a reference name at L1-L2 and must hedge below it, so the
    quality of the L2 boundary is the precision of the names it lets through and
    the share of the cohort it lets through at all.
    """

    combinations = []
    for struct in STRUCT_GRID:
        for ms2 in MS2_GRID:
            for mass in MASS_GRID:
                for ms2_exact in MS2_EXACT_GRID:
                    if ms2_exact < ms2:
                        # A catalog hit is never allowed an easier bar than a de
                        # novo prediction, so these combinations are not options.
                        continue
                    thresholds = C.MsiThresholds(
                        struct_moderate=struct,
                        ms2_present=ms2,
                        mass_strong=mass,
                        ms2_probable_exact=ms2_exact,
                    )
                    named = [record for record in records if _msi_level(record, thresholds) <= 2]
                    point = wilson(sum(1 for record in named if record["correct"]), len(named))
                    point.update(
                        {
                            "struct_moderate": struct,
                            "ms2_present": ms2,
                            "ms2_probable_exact": ms2_exact,
                            "mass_strong": mass,
                            "named_share": round(len(named) / len(records), 4) if records else None,
                            "incorrect_names": len(named) - point["correct"],
                        }
                    )
                    combinations.append(point)
    default = C.DEFAULT_MSI_THRESHOLDS
    level_counts: dict[str, int] = {}
    for record in records:
        key = f"L{_msi_level(record, default)}"
        level_counts[key] = level_counts.get(key, 0) + 1
    return {
        "cohort": "truth-verified controls with axis-level evidence",
        "n": len(records),
        "naming_rule": "a report asserts a reference name at MSI L1-L2 and hedges below it",
        "current_rule": {
            "struct_moderate": default.struct_moderate,
            "ms2_present": default.ms2_present,
            "ms2_probable_exact": default.ms2_probable_exact,
            "mass_strong": default.mass_strong,
            "level_distribution": dict(sorted(level_counts.items())),
        },
        "sweep": sorted(combinations, key=lambda point: (-(point["accuracy"] or 0.0), -(point["n"] or 0))),
    }


def load_controls(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def axis_records_from_controls(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    """Build axis records from the controls CSV when it carries the axis columns."""

    records = []
    for row in rows:
        if "axis_structure" not in row:
            return []
        structure = _float(row.get("axis_structure"))
        if structure is None:
            continue
        records.append(
            {
                "control_id": row.get("control_id", ""),
                "structure": structure,
                "mass": _float(row.get("axis_mass")),
                "ms2": _float(row.get("axis_ms2")),
                "observed_mz_available": _truthy(row.get("observed_mz_available", "")),
                "correct": _truthy(row.get("top1_correct_0p70", "")),
            }
        )
    return records


def axis_records_from_reports(root: Path) -> list[dict[str, Any]]:
    """Build axis records straight from a report bundle."""

    records = []
    for report_dir in sorted(root.glob("X-*")):
        if not report_dir.is_dir():
            continue
        slug = report_dir.name.lower()
        data = report_dir / "data"
        fields = C.read_fields_from_sidecars(data, slug)
        if not fields:
            continue
        candidates_path = data / f"{slug}_top_candidates.csv"
        if not candidates_path.exists():
            continue
        with candidates_path.open(newline="") as handle:
            candidates = sorted(csv.DictReader(handle), key=lambda row: _float(row.get("rank")) or 1e9)
        truth = _float(candidates[0].get("truth_tanimoto")) if candidates else None
        if truth is None or truth <= 0.0:
            continue
        record = _axis_record(fields)
        record["control_id"] = report_dir.name
        record["correct"] = truth >= 0.70
        records.append(record)
    return records


def render_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# Threshold audit",
        "",
        "Generated by `tools/audit_thresholds.py`. Every cutoff proposed for the",
        "reports should be traceable to a row in one of these tables.",
        "",
        "## Cutting on the evidence score",
        "",
    ]
    score = summary["score_cutoff"]
    lines += [
        f"Cohort: {score['n']} truth-verified controls, base rate {score['base_rate']:.1%}. "
        f"Rank AUC {score['auc']:.3f} (95% CI {score['auc_ci95'][0]:.3f}-{score['auc_ci95'][1]:.3f}); "
        "0.5 means the score carries no information about correctness.",
        "",
        "| Threshold | Kept | Coverage | Accuracy | 95% CI | Lift over base rate |",
        "|---|---|---|---|---|---|",
    ]
    for point in score["grid"]:
        accuracy = "n/a" if point["accuracy"] is None else f"{point['accuracy']:.1%}"
        interval = (
            "n/a"
            if point["ci95"][0] is None
            else f"{point['ci95'][0]:.1%}-{point['ci95'][1]:.1%}"
        )
        lift = "n/a" if point["lift_over_base_rate"] is None else f"{point['lift_over_base_rate']:+.1%}"
        lines.append(
            f"| {point['threshold']:.2f} | {point['n']} | {point['coverage']:.0%} | "
            f"{accuracy} | {interval} | {lift} |"
        )
    lines += ["", "## Why the pooled score is flat", ""]
    strat = summary["source_stratification"]
    if strat["confounded"]:
        lines.append(
            "The pooled curve is flatter than either source on its own, so the score "
            "is ordering provenance rather than correctness. Calibrating per source is "
            "the fix; moving a threshold on the pooled score is not."
        )
    else:
        lines.append("Pooling does not cost discrimination on this cohort.")
    lines += [
        "",
        "| Cohort | n | Base rate | Median raw score | Raw AUC | Calibrated AUC |",
        "|---|---|---|---|---|---|",
    ]
    for stratum, block in strat["strata"].items():
        raw, calibrated = block["raw_full"], block["loo_probability_full"]
        lines.append(
            f"| {stratum} | {raw['n']} | {raw['base_rate']:.1%} | {raw['median_score']:.3f} | "
            f"{raw['auc']:.3f} | {calibrated['auc']:.3f} |"
        )

    per_source = strat.get("per_source_calibration") or {}
    if per_source.get("strata"):
        lines += [
            "",
            "### Does calibrating per source recover the ordering",
            "",
            "Leave-one-out throughout: a single Platt fit is monotone, so an in-sample "
            "\"calibrated AUC\" would restate the raw one rather than test anything.",
            "",
            "| Source | n | Base rate | Raw AUC | Pooled-calibration AUC | Per-source-calibration AUC |",
            "|---|---:|---:|---:|---:|---:|",
        ]

        def cell(value: Any) -> str:
            return "n/a" if value is None else f"{value:.3f}"

        for source, block in per_source["strata"].items():
            lines.append(
                f"| {source} | {block['n']} | {block['base_rate']:.1%} | "
                f"{cell(block['auc_raw'])} | {cell(block['auc_pooled_loo'])} | "
                f"{cell(block['auc_per_source_loo'])} |"
            )
        lines += ["", f"**{per_source['recommendation'].capitalize()}.**"]
        if per_source.get("per_source_beats_pooled") and not per_source.get("per_source_beats_raw"):
            lines.append(
                " Splitting the calibration does help - it beats the pooled fit in every "
                "stratum - but the uncalibrated evidence score beats both, in every stratum. "
                "So the Platt step is spending ordering rather than adding calibration, and a "
                "tier boundary belongs on the raw score within a source."
            )
        lines += [
            "",
            "The in-sample fit per source is recorded in the JSON, including the score range each "
            "one can reach. A cutoff above that range names a tier no report can enter.",
        ]

    degeneracy = strat.get("positive_class_degeneracy") or {}
    if degeneracy.get("positive_class_is_level_1"):
        lines += [
            "",
            "### What none of these AUCs can tell you",
            "",
            f"All {degeneracy['n_correct']} correct calls in this cohort are MSI Level 1, and so are "
            f"exactly the same {degeneracy['n_msi_level_1']} reports - the sets coincide, with "
            f"{degeneracy['correct_outside_level_1']} correct calls outside Level 1.",
            "",
            "That is structural rather than unlucky. These controls are Metabolon Level-1 compounds, "
            "which is what makes them a truth set, and a report reaches MSI Level 1 exactly when a "
            "Level-1 standard's top candidate matches it. So \"the pipeline got this right\" and "
            "\"this report is already confirmed\" are the same event, and every AUC above is "
            "measuring how well a score separates the already-confirmed.",
            "",
            "It follows that this cohort cannot establish a confidence tier for reports without a "
            "standard, which is the whole population such a tier is for. The within-source numbers "
            "invite the opposite reading, which is why this is measured and stated here rather than "
            "left to the reader.",
        ]

    lines += ["", "## How many candidates to display", ""]
    depth = summary["candidate_depth"]
    lines += [
        "Each control has one correct structure, so a list of N puts `N - hit rate`",
        "incorrect structures in front of the reader.",
        "",
        "| Shown | Correct structure on the page | 95% CI | Incorrect shown per compound "
        "| Incorrect shown per additional hit |",
        "|---|---|---|---|---|",
    ]
    for point in depth["points"]:
        interval = (
            "n/a"
            if point["ci95"][0] is None
            else f"{point['ci95'][0]:.1%}-{point['ci95'][1]:.1%}"
        )
        exchange = point.get("incorrect_shown_per_additional_hit")
        lines.append(
            f"| {point['n_shown']} | {point['hit_rate']:.1%} | {interval} | "
            f"{point['incorrect_shown_per_compound']:.2f} | "
            f"{'-' if exchange is None else f'{exchange:.1f}'} |"
        )
    lines += ["", "## MSI boundary and the naming decision", ""]
    msi = summary["msi_rule"]
    if msi.get("status") == "not_evaluated":
        lines += [f"Not evaluated: {msi['reason']}", ""]
    else:
        current = msi["current_rule"]
        lines += [
            f"Cohort: {msi['n']} controls. Current rule: structure >= "
            f"{current['struct_moderate']}, MS2 >= {current['ms2_present']} "
            f"({current['ms2_probable_exact']} for a catalog hit), mass >= "
            f"{current['mass_strong']}. Level distribution: {current['level_distribution']}.",
            "",
            "| structure | MS2 | MS2 (catalog hit) | mass | Names asserted | Share "
            "| Name accuracy | Incorrect names |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for point in msi["sweep"]:
            accuracy = "n/a" if point["accuracy"] is None else f"{point['accuracy']:.1%}"
            lines.append(
                f"| {point['struct_moderate']} | {point['ms2_present']} | "
                f"{point['ms2_probable_exact']} | {point['mass_strong']} | "
                f"{point['n']} | {point['named_share']:.0%} | {accuracy} | {point['incorrect_names']} |"
            )
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--structure-controls",
        type=Path,
        default=Path("docs/paper/data/structure/audit_structure_controls.csv"),
    )
    parser.add_argument(
        "--reports-root",
        type=Path,
        default=None,
        help="Report bundle to read axis-level evidence from when the controls CSV lacks it.",
    )
    parser.add_argument(
        "--exec-summary",
        type=Path,
        default=Path("docs/paper/data/exec_summary/exec_summary_by_compound.csv"),
        help=(
            "Per-compound exec summary, joined to the controls to check whether this cohort's "
            "correct calls are all MSI Level 1. They are, and it bounds what the AUCs can mean."
        ),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("docs/paper/data/thresholds"))
    args = parser.parse_args()

    rows = load_controls(args.structure_controls)
    levels = load_msi_levels(args.exec_summary)
    scored = [row for row in rows if _float(row.get("loo_probability_full")) is not None]
    scores = [float(row["loo_probability_full"]) for row in scored]
    labels = [_truthy(row.get("top1_correct_0p70", "")) for row in scored]

    if args.reports_root is not None:
        records = axis_records_from_reports(args.reports_root)
        axis_source = str(args.reports_root)
    else:
        records = axis_records_from_controls(rows)
        axis_source = str(args.structure_controls)

    summary = {
        "method": (
            "Precision and coverage for each candidate cutoff, measured on the "
            "truth-verified controls, so a threshold can be chosen from an "
            "operating point rather than by hand."
        ),
        "inputs": {
            "structure_controls": str(args.structure_controls),
            "axis_source": axis_source if records else None,
        },
        "score_cutoff": score_cutoff_sweep(scores, labels),
        "source_stratification": source_stratification(rows, levels),
        "candidate_depth": candidate_depth_sweep(rows),
        "msi_rule": msi_rule_sweep(records)
        if records
        else {
            "status": "not_evaluated",
            "reason": (
                "Axis-level evidence is unavailable. Run "
                "tools/audit_structure_predictions.py against a report bundle to add the "
                "axis_* columns to the controls CSV, or pass --reports-root."
            ),
        },
        "scope_note": (
            "Retrospective operating curves on a 92-control cohort. Intervals are wide; "
            "they bound how confidently a cutoff can be placed, and a flat curve means no "
            "cutoff on that quantity is defensible at any value."
        ),
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.out_dir / "audit_thresholds_summary.json"
    json_path.write_text(json.dumps(summary, indent=2) + "\n")
    markdown_path = args.out_dir / "audit_thresholds_summary.md"
    markdown_path.write_text(render_markdown(summary) + "\n")

    score = summary["score_cutoff"]
    print(
        f"[thresholds] score AUC {score['auc']:.3f} "
        f"(95% CI {score['auc_ci95'][0]:.3f}-{score['auc_ci95'][1]:.3f}), "
        f"base rate {score['base_rate']:.1%}"
    )
    for point in summary["candidate_depth"]["points"]:
        print(
            f"[thresholds] top-{point['n_shown']}: correct structure on the page "
            f"{point['hit_rate']:.1%}, incorrect shown/compound "
            f"{point['incorrect_shown_per_compound']:.2f}"
        )
    print(f"[thresholds] wrote {json_path} and {markdown_path}")


if __name__ == "__main__":
    main()
