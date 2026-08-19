#!/usr/bin/env python3
"""Can an MSI sub-level (2a/2b/2c, 3a/3b/3c) be cut from a confidence score?

Read-only, pure standard library, and it runs against inputs committed to this
repository rather than a delivered report bundle, so it can join the CI
regenerate-and-diff set.

A sub-level scheme that bands a score into letters makes two claims, and this
audit checks both against the delivered 3,897-report set:

1. **The cutoffs are attainable.** A band nobody can enter is not a band. The
   identification tiers (``TIER_HIGH``/``TIER_MODERATE``/``TIER_LOW``) are
   defined on the Platt-scaled evidence score, and Platt scaling with a fixed
   ``(A, B)`` maps the raw evidence average onto a bounded interval - so whether
   a cutoff is reachable is arithmetic, not opinion. Section 1 inverts the Platt
   map at each cutoff and reports the raw score it would take to get there.

2. **The letters are ordered by correctness.** "2b" has to beat "2c" or the
   scheme misinforms the reader. Section 3 measures P(rank-1 correct) inside
   each band on the truth-verified cohort.

The second measurement has a trap in it that the audit is built around. The 44
authentic-standard-confirmed reports are correct *by construction* - truth
matching is what makes them Level 1 - and their displayed score is pinned by
policy (``CONFIRMED_OVERALL``), not produced by the model. Leaving them in the
denominator measures the pin. Every discrimination statistic here is therefore
reported twice: over the whole truth-verified cohort, and over the model-scored
reports alone.

Usage:
    python3 tools/audit_sublevel_bands.py
    python3 tools/audit_sublevel_bands.py --exec-summary path/to/by_compound.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "structural_elucidation" / "reports"))

import confidence as C  # noqa: E402  (path insert above is deliberate)

csv.field_size_limit(10**7)

DEFAULT_EXEC_SUMMARY = PROJECT_ROOT / "docs" / "paper" / "data" / "exec_summary" / "exec_summary_by_compound.csv"
DEFAULT_TOP10 = PROJECT_ROOT / "docs" / "paper" / "data" / "exec_summary" / "exec_summary_top10.csv"
DEFAULT_OUT_DIR = PROJECT_ROOT / "docs" / "paper" / "data" / "identification"

# The per-compound scores a sub-level scheme could be cut from. ``overall`` is the
# report headline (the exec summary carries it at rank 1 as
# ``rank1_candidate_confidence``); ``evidence_confidence`` is the pipeline's own
# candidate score, which the matrix workbook surfaces as "Confidence Score";
# ``raw_evidence`` is the pre-Platt weighted axis average, which the candidate
# widget on the candidate-confidence rework re-anchors rank 1 to. The third is
# reconstructed rather than read - see :func:`raw_evidence_rows`.
RAW_EVIDENCE_COLUMN = "_raw_evidence"

SCORES = {
    "candidate_confidence": {
        "column": "rank1_candidate_confidence",
        "label": "candidate_confidence (report overall, Platt-scaled)",
    },
    "evidence_confidence": {
        "column": "rank1_evidence_confidence",
        "label": "evidence_confidence (pipeline candidate score / matrix \"Confidence Score\")",
    },
    "raw_evidence": {
        "column": RAW_EVIDENCE_COLUMN,
        "label": "raw_evidence (pre-Platt axis average; the candidate-widget re-anchor)",
    },
}

# Band schemes to evaluate. The first is the one the tier constants imply; the
# others exist to show that the failure is a property of the scores rather than
# of one unlucky choice of cutoff.
BAND_SCHEMES = {
    "identification_tiers": {
        "cuts": (C.TIER_HIGH, C.TIER_MODERATE, C.TIER_LOW),
        "note": "TIER_HIGH / TIER_MODERATE / TIER_LOW, the cutoffs the report tiers already use",
    },
    "reachable_ceiling": {
        "cuts": (0.585, 0.48, 0.40),
        "note": "top band lowered to sit inside the Platt-reachable range",
    },
    "empirical_terciles": {
        "cuts": (0.55, 0.46, 0.40),
        "note": "roughly equal-sized bands over the model-scored population",
    },
}

BAND_LETTERS = ("a", "b", "c")
BELOW = "below"

# A sub-level scheme is only worth shipping if each letter is meaningfully more
# accurate than the next. Anything under this is noise on a cohort this size.
MONOTONE_MARGIN = 0.02

# The curated sub-levels that are actually a *letter within a level* - the ones a
# derived 2a/2b/2c scheme would have to reproduce. The rest of the SubLevel
# vocabulary (`4`, `4?`, `5`, `KNOWN`, `Released`, `Exclude`, `Artifact`) is
# either a bare level or a curation state, and mixing those in inflates the
# apparent spread with the Level-1-adjacent entries.
FINE_SUBLEVELS = ("2C", "3A", "3B", "3C", "3D")


def _num(value) -> float | None:
    text = str("" if value is None else value).strip()
    if not text or text.lower() in {"nan", "none", "null", "<na>"}:
        return None
    try:
        result = float(text)
    except ValueError:
        return None
    return None if math.isnan(result) else result


def _truthy(value) -> bool:
    return str("" if value is None else value).strip().lower() in {"true", "1", "yes"}


def _quantiles(values: list[float], points=(0.0, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1.0)) -> dict[str, float]:
    if not values:
        return {}
    ordered = sorted(values)
    out = {}
    for q in points:
        index = min(len(ordered) - 1, max(0, int(round(q * (len(ordered) - 1)))))
        out[f"p{int(round(q * 100))}"] = round(ordered[index], 4)
    return out


def _auc(positive: list[float], negative: list[float]) -> float | None:
    """Mann-Whitney AUC: P(score(correct) > score(incorrect)), ties counted half.

    0.50 means the score carries no information about correctness; below 0.50
    means it carries information with the sign reversed.
    """
    if not positive or not negative:
        return None
    merged = sorted([(v, 1) for v in positive] + [(v, 0) for v in negative])
    ranks: dict[int, float] = {}
    i = 0
    while i < len(merged):
        j = i
        while j + 1 < len(merged) and merged[j + 1][0] == merged[i][0]:
            j += 1
        shared = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = shared
        i = j + 1
    rank_sum = sum(ranks[k] for k in range(len(merged)) if merged[k][1] == 1)
    n_pos, n_neg = len(positive), len(negative)
    return round((rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg), 4)


def _band_of(value: float | None, cuts) -> str:
    if value is None:
        return "missing"
    for letter, cut in zip(BAND_LETTERS, cuts):
        if value >= cut:
            return letter
    return BELOW


# --- recovering the pre-Platt score ------------------------------------------


def raw_evidence_rows(rows: list[dict]) -> dict[str, object]:
    """Recover ``raw_evidence`` per compound by inverting the Platt map, in place.

    The candidate-confidence rework re-anchors the report's rank-1 value from the
    Platt ``overall`` to ``raw_evidence``, on the grounds that ``overall`` caps
    near 0.59 and deflates strong best matches. That is a real display defect and
    the re-anchor fixes it. What it cannot do is change the ordering:
    ``overall = sigmoid(A * raw + B)`` with ``A > 0`` is strictly increasing, so
    the two scores rank every report identically and every rank statistic - AUC
    included - is the same number for both.

    Rather than assert that, this recovers ``raw = (logit(overall) - B) / A`` and
    measures it, so the claim is checked rather than argued. The reconstruction is
    exact wherever ``overall`` is the model's own output.

    It is *not* recoverable where a policy override replaced ``overall``: an
    authentic-standard confirmation pins the headline, and inverting a pinned
    value yields a raw score above 1.0, which no weighted average of 0..1 axes can
    produce. Those are left unset rather than clamped, because clamping would
    invent a value for exactly the reports whose influence is in question.
    """
    ceiling = C._sigmoid(C.PLATT_A + C.PLATT_B)
    floor = C._sigmoid(C.PLATT_B)
    recovered = 0
    unrecoverable = 0
    for row in rows:
        overall = _num(row.get("rank1_candidate_confidence"))
        if overall is None or overall <= 0.0 or overall >= 1.0:
            row[RAW_EVIDENCE_COLUMN] = ""
            unrecoverable += 1
            continue
        raw = (math.log(overall / (1.0 - overall)) - C.PLATT_B) / C.PLATT_A
        if raw > 1.0 or raw < 0.0:
            row[RAW_EVIDENCE_COLUMN] = ""
            unrecoverable += 1
            continue
        row[RAW_EVIDENCE_COLUMN] = f"{raw:.6f}"
        recovered += 1
    values = [v for v in (_num(r.get(RAW_EVIDENCE_COLUMN)) for r in rows) if v is not None]
    return {
        "recovered": recovered,
        "unrecoverable_policy_pinned": unrecoverable,
        "platt_floor": round(floor, 4),
        "platt_ceiling": round(ceiling, 4),
        "distribution": _quantiles(values),
        "note": (
            "raw = (logit(overall) - B) / A. Strictly monotone in overall, so every "
            "rank statistic is identical between the two."
        ),
    }


# --- section 1: is the cutoff reachable at all -------------------------------


def platt_reachability() -> dict[str, object]:
    """Invert the Platt map at each tier cutoff.

    The raw evidence score is a weighted average of 0..1 axes, so it cannot
    exceed 1.0. Any cutoff whose Platt pre-image is above 1.0 is unreachable by
    the evidence model, and can only ever be entered by a policy override.
    """
    lo = C._sigmoid(C.PLATT_B)
    hi = C._sigmoid(C.PLATT_A + C.PLATT_B)
    cutoffs = {}
    for name, value in (("TIER_HIGH", C.TIER_HIGH), ("TIER_MODERATE", C.TIER_MODERATE), ("TIER_LOW", C.TIER_LOW)):
        raw_needed = (math.log(value / (1.0 - value)) - C.PLATT_B) / C.PLATT_A
        cutoffs[name] = {
            "cutoff": value,
            "raw_evidence_required": round(raw_needed, 4),
            "reachable": bool(raw_needed <= 1.0),
        }
    return {
        "platt_a": C.PLATT_A,
        "platt_b": C.PLATT_B,
        "evidence_score_min": round(lo, 4),
        "evidence_score_max": round(hi, 4),
        "cutoffs": cutoffs,
        "unreachable_cutoffs": sorted(k for k, v in cutoffs.items() if not v["reachable"]),
    }


# --- section 2: how the population spreads over the bands --------------------


def band_occupancy(rows: list[dict], column: str, cuts) -> dict[str, object]:
    counts = Counter(_band_of(_num(r.get(column)), cuts) for r in rows)
    total = len(rows)
    return {
        "counts": {b: counts.get(b, 0) for b in (*BAND_LETTERS, BELOW, "missing")},
        "shares": {
            b: round(counts.get(b, 0) / total, 4) if total else 0.0
            for b in (*BAND_LETTERS, BELOW, "missing")
        },
        "empty_bands": sorted(b for b in BAND_LETTERS if not counts.get(b, 0)),
    }


def by_msi_level(rows: list[dict], column: str, cuts) -> dict[str, dict[str, int]]:
    """Cross-tab of MSI level against band - the cells a sub-level table would carry."""
    table: dict[str, Counter] = defaultdict(Counter)
    for row in rows:
        level = str(row.get("msi_level") or "").strip() or "unknown"
        table[level][_band_of(_num(row.get(column)), cuts)] += 1
    return {
        level: {b: table[level].get(b, 0) for b in (*BAND_LETTERS, BELOW, "missing")}
        for level in sorted(table)
    }


# --- section 3: are the bands ordered by correctness --------------------------


def _precision(subset: list[dict], column: str, cuts) -> dict[str, object]:
    total: Counter = Counter()
    correct: Counter = Counter()
    for row in subset:
        band = _band_of(_num(row.get(column)), cuts)
        total[band] += 1
        if _truthy(row.get("matched_at_rank_1")):
            correct[band] += 1
    per_band = {}
    for band in (*BAND_LETTERS, BELOW):
        n = total.get(band, 0)
        per_band[band] = {
            "n": n,
            "correct": correct.get(band, 0),
            "p_correct": round(correct.get(band, 0) / n, 4) if n else None,
        }
    # Monotone means every letter beats the next by a margin worth printing.
    observed = [(b, per_band[b]["p_correct"]) for b in (*BAND_LETTERS, BELOW) if per_band[b]["p_correct"] is not None]
    monotone = all(
        observed[i][1] >= observed[i + 1][1] + MONOTONE_MARGIN for i in range(len(observed) - 1)
    )
    inversions = [
        f"{observed[i][0]}<{observed[i + 1][0]}"
        for i in range(len(observed) - 1)
        if observed[i][1] < observed[i + 1][1]
    ]
    return {"bands": per_band, "monotone": monotone, "inversions": inversions}


def _median(values: list[float]) -> float | None:
    return round(sorted(values)[len(values) // 2], 4) if values else None


def discrimination(rows: list[dict], column: str) -> dict[str, object]:
    positive = [_num(r.get(column)) for r in rows if _truthy(r.get("matched_at_rank_1"))]
    negative = [_num(r.get(column)) for r in rows if not _truthy(r.get("matched_at_rank_1"))]
    positive = [v for v in positive if v is not None]
    negative = [v for v in negative if v is not None]
    return {
        "n_correct": len(positive),
        "n_incorrect": len(negative),
        "median_correct": _median(positive),
        "median_incorrect": _median(negative),
        "auc": _auc(positive, negative),
    }


def deciles(rows: list[dict], column: str, k: int = 10) -> list[dict[str, object]]:
    """P(correct) by score decile - the finest ordering check the cohort supports."""
    scored = sorted(
        (r for r in rows if _num(r.get(column)) is not None),
        key=lambda r: (_num(r.get(column)), str(r.get("x_compound_id") or "")),
    )
    if len(scored) < k:
        return []
    size = len(scored) / k
    out = []
    for i in range(k):
        chunk = scored[int(i * size):int((i + 1) * size)]
        if not chunk:
            continue
        correct = sum(_truthy(r.get("matched_at_rank_1")) for r in chunk)
        out.append({
            "decile": i + 1,
            "score_min": round(_num(chunk[0][column]), 4),
            "score_max": round(_num(chunk[-1][column]), 4),
            "n": len(chunk),
            "correct": correct,
            "p_correct": round(correct / len(chunk), 4),
        })
    return out


# --- section 4: do the curated sub-levels separate on either score -----------


def curated_sublevel_separation(rows: list[dict], column: str, min_n: int = 25) -> dict[str, object]:
    """Median score per curated Metabolon SubLevel.

    If a score could be banded into sub-levels, the curated letters - which are
    the only ground truth for what a sub-level means here - would separate on it.
    """
    grouped: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        value = _num(row.get(column))
        raw = str(row.get("metabolon_sublevel") or "").strip()
        if value is None or not raw:
            continue
        for piece in raw.split(";"):
            piece = piece.strip()
            if piece:
                grouped[piece].append(value)
    per_sublevel = {}
    for sublevel, values in grouped.items():
        if len(values) < min_n:
            continue
        quart = _quantiles(values, (0.25, 0.5, 0.75))
        per_sublevel[sublevel] = {
            "n": len(values),
            "median": quart.get("p50"),
            "p25": quart.get("p25"),
            "p75": quart.get("p75"),
        }
    medians = [v["median"] for v in per_sublevel.values() if v["median"] is not None]
    fine = [
        v["median"] for k, v in per_sublevel.items()
        if k in FINE_SUBLEVELS and v["median"] is not None
    ]
    return {
        "per_sublevel": dict(sorted(per_sublevel.items())),
        "median_spread": round(max(medians) - min(medians), 4) if medians else None,
        "fine_sublevels": sorted(k for k in per_sublevel if k in FINE_SUBLEVELS),
        "fine_median_spread": round(max(fine) - min(fine), 4) if fine else None,
    }


# --- section 6: the one derived sub-letter ------------------------------------

# The candidate-confidence rework derives exactly one fine sub-letter, on the
# grounds that same-formula multiplicity among the top candidates is a real
# semantic signal rather than a score bin - which it is. The question this
# section asks is narrower: the derived label reuses the curated letter ``3A``,
# so does it mean the same thing.
DERIVED_3A_MIN_ISOMERS = 2


def same_formula_multiplicity(rows: list[dict]) -> int:
    """Top-10 candidates sharing the rank-1 predicted formula."""
    top = sorted(rows, key=lambda r: _num(r.get("rank")) or 1e9)[:10]
    if not top:
        return 0
    formula = str(top[0].get("candidate_formula") or "").strip()
    if not formula:
        return 0
    return sum(1 for r in top if str(r.get("candidate_formula") or "").strip() == formula)


def derived_3a_agreement(rows: list[dict], top10: list[dict]) -> dict[str, object]:
    """Does the derived ``3A`` agree with the curated ``3A`` whose name it takes?

    A derived label and a curated label that share a letter will be read as the
    same claim. They are only the same claim if they select the same compounds.
    """
    by_compound: dict[str, list[dict]] = defaultdict(list)
    for row in top10:
        by_compound[str(row.get("x_compound_id") or "").strip()].append(row)
    multiplicity = {cid: same_formula_multiplicity(group) for cid, group in by_compound.items()}

    level_3 = [r for r in rows if str(r.get("msi_level") or "").strip() == "3"]
    derived = [
        r for r in level_3
        if multiplicity.get(str(r.get("x_compound_id") or "").strip(), 0) >= DERIVED_3A_MIN_ISOMERS
    ]

    def has_curated_3a(row: dict) -> bool:
        return any(p.strip() == "3A" for p in str(row.get("metabolon_sublevel") or "").split(";"))

    def has_any_curated(row: dict) -> bool:
        return bool(str(row.get("metabolon_sublevel") or "").strip())

    curated_3a = [r for r in rows if has_curated_3a(r)]
    derived_ids = {str(r.get("x_compound_id") or "").strip() for r in derived}
    # Precision is only meaningful over compounds a curator actually labeled;
    # an uncurated compound is not a disagreement, it is an absence.
    derived_curated = [r for r in derived if has_any_curated(r)]
    derived_and_curated_3a = [r for r in derived_curated if has_curated_3a(r)]
    curated_3a_recovered = [r for r in curated_3a if str(r.get("x_compound_id") or "").strip() in derived_ids]

    letters: Counter = Counter()
    for r in derived_curated:
        for piece in str(r.get("metabolon_sublevel") or "").split(";"):
            piece = piece.strip()
            if piece:
                letters[piece] += 1

    return {
        "rule": f"msi_level == 3 and n_isomers >= {DERIVED_3A_MIN_ISOMERS}",
        "derived_3a": len(derived),
        "level_3_total": len(level_3),
        "share_of_level_3": round(len(derived) / len(level_3), 4) if level_3 else None,
        "curated_3a_total": len(curated_3a),
        "derived_with_any_curated_letter": len(derived_curated),
        "derived_that_are_curated_3a": len(derived_and_curated_3a),
        "precision_vs_curated": (
            round(len(derived_and_curated_3a) / len(derived_curated), 4) if derived_curated else None
        ),
        "curated_3a_also_derived": len(curated_3a_recovered),
        "recall_vs_curated": (
            round(len(curated_3a_recovered) / len(curated_3a), 4) if curated_3a else None
        ),
        "curated_letters_inside_derived_3a": dict(letters.most_common()),
    }


# --- assembly ----------------------------------------------------------------


def _reanchor_equivalence(per_score: dict) -> dict[str, object]:
    """Does re-anchoring rank 1 to ``raw_evidence`` change what the score can do.

    Two separate questions, and the answers differ, which is the whole point of
    reporting them together. Re-anchoring moves the score onto a scale where the
    tier cutoffs are reachable, so it fixes the empty top band. It leaves every
    ordering statistic untouched, so it does not fix the inverted bands.
    """
    platt = per_score["candidate_confidence"]
    raw = per_score["raw_evidence"]
    platt_auc = platt["discrimination_model_scored"]["auc"]
    raw_auc = raw["discrimination_model_scored"]["auc"]
    platt_top = platt["band_schemes"]["identification_tiers"]["occupancy"]["counts"]["a"]
    raw_top = raw["band_schemes"]["identification_tiers"]["occupancy"]["counts"]["a"]
    return {
        "auc_model_scored_platt": platt_auc,
        "auc_model_scored_raw": raw_auc,
        "auc_unchanged": platt_auc == raw_auc,
        "top_band_platt": platt_top,
        "top_band_raw": raw_top,
        "median_platt": platt["distribution"].get("p50"),
        "median_raw": raw["distribution"].get("p50"),
        "fixes_reachability": raw_top > platt_top,
        "fixes_ordering": any(
            scheme["precision_model_scored"]["monotone"] for scheme in raw["band_schemes"].values()
        ),
    }


def audit(rows: list[dict], top10: list[dict] | None = None) -> dict[str, object]:
    raw_recovery = raw_evidence_rows(rows)
    truth = [
        r for r in rows
        if str(r.get("truth_basis") or "").strip() and str(r.get("reference_answer") or "").strip()
    ]
    # The pinned reports are correct by construction and their score is a display
    # value, so they are held out of every discrimination statistic as well as
    # reported inside it.
    pinned = [r for r in truth if str(r.get("msi_level") or "").strip() == "1"]
    model_scored = [r for r in truth if str(r.get("msi_level") or "").strip() != "1"]

    per_score: dict[str, object] = {}
    for key, meta in SCORES.items():
        column = meta["column"]
        values = [v for v in (_num(r.get(column)) for r in rows) if v is not None]
        saturated = sum(1 for v in values if v >= 0.98)
        schemes = {}
        for name, scheme in BAND_SCHEMES.items():
            cuts = tuple(scheme["cuts"])
            schemes[name] = {
                "cuts": list(cuts),
                "note": scheme["note"],
                "occupancy": band_occupancy(rows, column, cuts),
                "by_msi_level": by_msi_level(rows, column, cuts),
                "precision_all_truth": _precision(truth, column, cuts),
                "precision_model_scored": _precision(model_scored, column, cuts),
            }
        per_score[key] = {
            "label": meta["label"],
            "column": column,
            "n_present": len(values),
            "distribution": _quantiles(values),
            "share_at_or_above_0_98": round(saturated / len(values), 4) if values else None,
            "discrimination_all_truth": discrimination(truth, column),
            "discrimination_model_scored": discrimination(model_scored, column),
            "deciles_model_scored": deciles(model_scored, column),
            "curated_sublevel_separation": curated_sublevel_separation(rows, column),
            "band_schemes": schemes,
        }

    accuracy_by_level: dict[str, dict[str, object]] = {}
    level_total: Counter = Counter()
    level_correct: Counter = Counter()
    for row in truth:
        level = str(row.get("msi_level") or "").strip() or "unknown"
        level_total[level] += 1
        if _truthy(row.get("matched_at_rank_1")):
            level_correct[level] += 1
    for level in sorted(level_total):
        accuracy_by_level[f"L{level}"] = {
            "n": level_total[level],
            "correct": level_correct[level],
            "p_correct": round(level_correct[level] / level_total[level], 4),
        }

    return {
        "cohort": {
            "compounds": len(rows),
            "truth_verified": len(truth),
            "truth_verified_correct": sum(_truthy(r.get("matched_at_rank_1")) for r in truth),
            "policy_pinned_level_1": len(pinned),
            "model_scored": len(model_scored),
            "model_scored_correct": sum(_truthy(r.get("matched_at_rank_1")) for r in model_scored),
        },
        "platt_reachability": platt_reachability(),
        "raw_evidence_recovery": raw_recovery,
        "reanchor_equivalence": _reanchor_equivalence(per_score),
        "scores": per_score,
        "accuracy_by_msi_level": accuracy_by_level,
        "derived_3a": derived_3a_agreement(rows, top10) if top10 else None,
        "monotone_margin": MONOTONE_MARGIN,
    }


# --- rendering ---------------------------------------------------------------


def _pct(value) -> str:
    return "n/a" if value is None else f"{value:.1%}"


def render_markdown(result: dict) -> str:
    cohort = result["cohort"]
    reach = result["platt_reachability"]
    lines: list[str] = []
    add = lines.append

    add("# MSI sub-level confidence-band audit")
    add("")
    add(
        "Asks whether an MSI sub-level (`2a`/`2b`/`2c`, `3a`/`3b`/`3c`) can be cut from a "
        "confidence score, by checking the two things such a scheme claims: that its cutoffs are "
        "reachable, and that its letters are ordered by correctness. Generated by "
        "`tools/audit_sublevel_bands.py`; read-only, and reproducible from committed inputs."
    )
    add("")
    add(
        f"{cohort['compounds']:,} delivered compounds, of which {cohort['truth_verified']:,} are "
        f"truth-verified ({cohort['truth_verified_correct']:,} correct at rank 1). "
        f"{cohort['policy_pinned_level_1']:,} of those are authentic-standard-confirmed Level 1 "
        "reports, which are correct by construction and carry a policy-pinned display score; the "
        f"remaining {cohort['model_scored']:,} are model-scored "
        f"({cohort['model_scored_correct']:,} correct)."
    )
    add("")

    add("## 1. Are the cutoffs reachable")
    add("")
    add(
        f"The evidence score is `sigmoid(A * raw + B)` with A={reach['platt_a']}, B={reach['platt_b']}, "
        "and `raw` is a weighted average of 0..1 axes, so it cannot exceed 1.0. That bounds the "
        f"score to **[{reach['evidence_score_min']}, {reach['evidence_score_max']}]** whatever the "
        "evidence says. Inverting the map at each tier cutoff gives the raw score it would take to "
        "get there."
    )
    add("")
    add("| Cutoff | Value | Raw evidence required | Reachable |")
    add("|---|---:|---:|---|")
    for name, block in reach["cutoffs"].items():
        mark = "yes" if block["reachable"] else "**no**"
        add(f"| `{name}` | {block['cutoff']} | {block['raw_evidence_required']} | {mark} |")
    add("")
    if reach["unreachable_cutoffs"]:
        add(
            f"**{', '.join('`' + c + '`' for c in reach['unreachable_cutoffs'])} cannot be reached by "
            "the evidence model at all.** A report only lands above it through a policy override - "
            "`CONFIRMED_OVERALL` for an authentic standard, or the curator-identity promotion - so a "
            "top band cut there holds exactly the reports that are already Level 1."
        )
        add("")

    recovery = result["raw_evidence_recovery"]
    eq = result["reanchor_equivalence"]
    add("## 1b. What re-anchoring rank 1 to `raw_evidence` does, and does not, fix")
    add("")
    add(
        "The candidate widget can re-anchor its rank-1 value from the Platt `overall` to the "
        "pre-Platt `raw_evidence`, because `overall` caps near the ceiling above and deflates "
        "strong best matches. `raw_evidence` is recovered here by inverting the Platt map, "
        f"`raw = (logit(overall) - B) / A`, for the {recovery['recovered']:,} reports whose "
        f"headline is the model's own output. The other {recovery['unrecoverable_policy_pinned']:,} "
        "carry a policy-pinned headline that inverts to a raw score above 1.0, which no weighted "
        "average of 0..1 axes can produce, so they are left unset rather than clamped."
    )
    add("")
    add("| | Platt `overall` | Re-anchored `raw_evidence` |")
    add("|---|---:|---:|")
    add(f"| Median | {eq['median_platt']} | {eq['median_raw']} |")
    add(
        f"| Reports in the top band at {C.TIER_HIGH} | {eq['top_band_platt']:,} "
        f"| {eq['top_band_raw']:,} |"
    )
    add(f"| AUC on model-scored reports | {eq['auc_model_scored_platt']} | {eq['auc_model_scored_raw']} |")
    add("")
    verdict = []
    verdict.append(
        "**Reachability: fixed.**" if eq["fixes_reachability"]
        else "**Reachability: not fixed.**"
    )
    verdict.append(
        "The tier cutoffs sit inside the raw scale, so the top band stops being empty."
        if eq["fixes_reachability"] else
        "The top band is no better populated on the raw scale."
    )
    verdict.append(
        "**Ordering: unchanged.**" if eq["auc_unchanged"]
        else "**Ordering: changed.**"
    )
    verdict.append(
        "`sigmoid(A * raw + B)` is strictly increasing, so the two scores rank every report "
        "identically and AUC is the same number for both. Re-anchoring is a rescale, and "
        "rescaling a score that is inverted leaves it inverted."
        if eq["auc_unchanged"] else
        "The two scores do not rank reports identically, which needs explaining before either is used."
    )
    add(" ".join(verdict))
    add("")
    if not eq["fixes_ordering"]:
        add(
            "So a validation figure quoted for `raw_evidence` is a figure for the Platt score under "
            "another name. Any such figure that looks discriminative is being carried by the "
            "policy-pinned reports, which are correct by construction."
        )
        add("")

    add("## 2. Where the population lands")
    add("")
    for key, block in result["scores"].items():
        dist = block["distribution"]
        add(f"### `{key}`")
        add("")
        add(f"{block['label']}. Column `{block['column']}`, {block['n_present']:,} values present.")
        add("")
        add("| p0 | p10 | p25 | p50 | p75 | p90 | p100 | >=0.98 |")
        add("|---:|---:|---:|---:|---:|---:|---:|---:|")
        add(
            f"| {dist.get('p0')} | {dist.get('p10')} | {dist.get('p25')} | {dist.get('p50')} | "
            f"{dist.get('p75')} | {dist.get('p90')} | {dist.get('p100')} | "
            f"{_pct(block['share_at_or_above_0_98'])} |"
        )
        add("")
        add("| Band scheme | Cuts | a | b | c | below | Empty bands |")
        add("|---|---|---:|---:|---:|---:|---|")
        for name, scheme in block["band_schemes"].items():
            counts = scheme["occupancy"]["counts"]
            empty = ", ".join(scheme["occupancy"]["empty_bands"]) or "none"
            cuts = " / ".join(str(c) for c in scheme["cuts"])
            add(
                f"| {name} | {cuts} | {counts['a']:,} | {counts['b']:,} | {counts['c']:,} | "
                f"{counts['below']:,} | {empty} |"
            )
        add("")

    add("## 3. Are the bands ordered by correctness")
    add("")
    add(
        "P(rank-1 candidate correct) inside each band. A band scheme that reads as confidence has "
        "to fall from `a` to `below`; the model-scored column is the honest one, because the "
        "pinned Level 1 reports are correct by construction."
    )
    add("")
    for key, block in result["scores"].items():
        add(f"### `{key}`")
        add("")
        disc_all = block["discrimination_all_truth"]
        disc_model = block["discrimination_model_scored"]
        add(
            f"AUC over all truth-verified reports **{disc_all['auc']}**, over model-scored reports "
            f"alone **{disc_model['auc']}** (0.50 is no discrimination; below 0.50 is discrimination "
            "with the sign reversed)."
        )
        add("")
        add("| Band scheme | a | b | c | below | Ordered |")
        add("|---|---:|---:|---:|---:|---|")
        for name, scheme in block["band_schemes"].items():
            bands = scheme["precision_model_scored"]["bands"]
            ordered = "yes" if scheme["precision_model_scored"]["monotone"] else "**no**"
            cells = " | ".join(_pct(bands[b]["p_correct"]) for b in (*BAND_LETTERS, BELOW))
            add(f"| {name} | {cells} | {ordered} |")
        add("")
        rows_dec = block["deciles_model_scored"]
        if rows_dec:
            add("Model-scored reports by score decile:")
            add("")
            add("| Decile | Score range | n | Correct | P(correct) |")
            add("|---:|---|---:|---:|---:|")
            for entry in rows_dec:
                add(
                    f"| {entry['decile']} | {entry['score_min']}-{entry['score_max']} | {entry['n']} | "
                    f"{entry['correct']} | {_pct(entry['p_correct'])} |"
                )
            add("")

    add("## 4. Do the curated sub-levels separate on either score")
    add("")
    add(
        "Metabolon's curated `SubLevel` in `data/known_X.csv` is the only ground truth for what a "
        "sub-level means here. If a score could be banded into sub-levels, the curated letters "
        "would separate on it."
    )
    add("")
    for key, block in result["scores"].items():
        sep = block["curated_sublevel_separation"]
        add(f"### `{key}`")
        add("")
        fine = ", ".join(f"`{s}`" for s in sep["fine_sublevels"]) or "none"
        add(
            f"Spread between the highest and lowest median across the fine sub-level letters "
            f"({fine}): **{sep['fine_median_spread']}**. Across every curated value including the "
            f"bare levels and curation states: {sep['median_spread']}."
        )
        add("")
        add("| Curated SubLevel | n | p25 | Median | p75 |")
        add("|---|---:|---:|---:|---:|")
        for sublevel, stats in sep["per_sublevel"].items():
            add(f"| `{sublevel}` | {stats['n']:,} | {stats['p25']} | {stats['median']} | {stats['p75']} |")
        add("")

    derived = result.get("derived_3a")
    if derived:
        add("## 4b. The one derived sub-letter, against the curated letter it reuses")
        add("")
        add(
            "Same-formula multiplicity among the top candidates is a real semantic signal and not a "
            "score bin, so deriving a positional-isomer label from it is a different kind of claim "
            "from banding a score. The question here is narrower: the derived label takes the name "
            f"`3A`, which a curator already uses, so does it select the same compounds. Rule: "
            f"`{derived['rule']}`."
        )
        add("")
        add("| | Count |")
        add("|---|---:|")
        add(f"| Derived `3A` | {derived['derived_3a']:,} |")
        add(f"| Share of all Level 3 | {_pct(derived['share_of_level_3'])} |")
        add(f"| Curated `3A` | {derived['curated_3a_total']:,} |")
        add(f"| Derived `3A` that a curator also labeled at all | {derived['derived_with_any_curated_letter']:,} |")
        add(f"| ...of those, curated `3A` | {derived['derived_that_are_curated_3a']:,} |")
        add("")
        add(
            f"Agreement with the curator, over the compounds a curator labeled: "
            f"**{_pct(derived['precision_vs_curated'])}** of derived `3A` carry curated `3A`, and "
            f"**{_pct(derived['recall_vs_curated'])}** of curated `3A` is recovered."
        )
        add("")
        add("Curated letters actually carried by the compounds the rule labels `3A`:")
        add("")
        add("| Curated letter | n |")
        add("|---|---:|")
        for letter, count in derived["curated_letters_inside_derived_3a"].items():
            add(f"| `{letter}` | {count:,} |")
        add("")

    add("## 5. What is ordered")
    add("")
    add(
        "For contrast, the coarse MSI level on the same truth-verified cohort. Level 1 is correct by "
        "construction; the rest is what the evidence rule achieves on its own."
    )
    add("")
    add("| MSI level | Truth-verified n | Correct | P(correct) |")
    add("|---|---:|---:|---:|")
    for level, stats in result["accuracy_by_msi_level"].items():
        add(f"| {level} | {stats['n']:,} | {stats['correct']:,} | {_pct(stats['p_correct'])} |")
    add("")

    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--exec-summary", type=Path, default=DEFAULT_EXEC_SUMMARY,
                        help="per-compound exec summary CSV (default: the committed one)")
    parser.add_argument("--exec-summary-top10", type=Path, default=DEFAULT_TOP10,
                        help="per-candidate exec summary CSV, for the derived-3A section")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                        help="directory for audit_sublevel_bands.json / .md")
    args = parser.parse_args(argv)

    if not args.exec_summary.exists():
        parser.error(f"exec summary not found: {args.exec_summary}")
    with args.exec_summary.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        parser.error(f"exec summary is empty: {args.exec_summary}")

    top10: list[dict] = []
    if args.exec_summary_top10.exists():
        with args.exec_summary_top10.open(newline="") as handle:
            top10 = list(csv.DictReader(handle))

    result = audit(rows, top10)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "audit_sublevel_bands.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    (args.out_dir / "audit_sublevel_bands.md").write_text(render_markdown(result))

    reach = result["platt_reachability"]
    print(f"compounds {result['cohort']['compounds']:,}; model-scored truth {result['cohort']['model_scored']:,}")
    if reach["unreachable_cutoffs"]:
        print(f"unreachable cutoffs: {', '.join(reach['unreachable_cutoffs'])}")
    for key, block in result["scores"].items():
        print(f"{key}: AUC(all)={block['discrimination_all_truth']['auc']} "
              f"AUC(model-scored)={block['discrimination_model_scored']['auc']}")
    print(f"wrote {args.out_dir / 'audit_sublevel_bands.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
