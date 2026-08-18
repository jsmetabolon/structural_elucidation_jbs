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
DEFAULT_OUT_DIR = PROJECT_ROOT / "docs" / "paper" / "data" / "identification"

# The two per-compound scores a sub-level scheme could be cut from. ``overall``
# is the report headline (the exec summary carries it re-anchored at rank 1 as
# ``rank1_candidate_confidence``); ``evidence_confidence`` is the pipeline's own
# candidate score, which the matrix workbook surfaces as "Confidence Score".
SCORES = {
    "candidate_confidence": {
        "column": "rank1_candidate_confidence",
        "label": "candidate_confidence (report overall, Platt-scaled)",
    },
    "evidence_confidence": {
        "column": "rank1_evidence_confidence",
        "label": "evidence_confidence (pipeline candidate score / matrix \"Confidence Score\")",
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


# --- assembly ----------------------------------------------------------------


def audit(rows: list[dict]) -> dict[str, object]:
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
        "scores": per_score,
        "accuracy_by_msi_level": accuracy_by_level,
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
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                        help="directory for audit_sublevel_bands.json / .md")
    args = parser.parse_args(argv)

    if not args.exec_summary.exists():
        parser.error(f"exec summary not found: {args.exec_summary}")
    with args.exec_summary.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        parser.error(f"exec summary is empty: {args.exec_summary}")

    result = audit(rows)
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
