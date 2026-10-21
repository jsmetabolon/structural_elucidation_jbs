#!/usr/bin/env python3
"""Calibrate candidate_identity_confidence against the truth controls.

The score means confidence in the structure, not confidence in whatever the
level happens to assert, so one bar applies to every tier: is the rank-1
structure actually the true molecule, by the same Tanimoto match already used
for ``matched_truth`` elsewhere (``TRUTH_THRESHOLD``, 0.70).

Every compound with a reference structure is already reported at Level 1, so
no Level 2-5 call carries a structural label. The replay therefore masks each
structure control's naming and curator fields, reruns the naming cascade on
its own top-10 candidates, and records which call it would have received and
whether rank 1 is the true structure. ``formula_only`` controls have no
reference structure, so they fall back to a formula match, which is a weaker
label flagged as such in ``label_basis``.

Subcommands:
  replay  write calibration_rows.csv and calibration_table.csv
  fit     replay, then write the calibration JSON (anchors, bands, slope)
  report  write calibration_report.md (spread, overlap, reliability, counts)
  apply   copy a delivery, add pure_nce_score, recompute the confidence column

Examples:
  python3 tools/calibrate_identity_confidence.py fit \\
      --delivery output/pipeline_equalization_reported_v2_delivery \\
      --pure-nce output/pure_nce_scores.csv --out-dir output/calibration \\
      --json tools/identity_confidence_calibration.json
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from pathlib import Path
from typing import Any

import pandas as pd

_TOOLS = Path(__file__).resolve().parent
PROJECT_ROOT = _TOOLS.parent
for path in (str(PROJECT_ROOT), str(_TOOLS)):
    if path not in sys.path:
        sys.path.insert(0, path)

import reported_identification as R  # noqa: E402
from build_exec_summary import TRUTH_THRESHOLD, normalise_formula  # noqa: E402

NAMING_FIELDS = (
    "feature_name",
    "report_display_name",
    "metabolon_level",
    "msi_level",
    "exception_basis",
    "identity_exclusion_reason",
    "excluded_library_reference",
    "library_collision",
    "known_x_tentative_information",
    "known_x_cleaned_external_info",
)
TEXT_COLUMNS = (*NAMING_FIELDS, "candidate_smiles", "candidate_formula", "reference_formula", "truth_basis")
PRIOR_STRENGTH = 10.0
MIN_LABELS_TO_CALIBRATE = 5
MIN_LABELS_FOR_SLOPE = 20
BASIS_STRUCTURE = "structure"
BASIS_FORMULA = "formula_proxy"


# --- loading -----------------------------------------------------------------


def read_delivery(delivery: Path) -> pd.DataFrame:
    frame = pd.read_csv(delivery / "exec_summary_top10.csv", low_memory=False)
    for column in TEXT_COLUMNS:
        if column in frame.columns:
            frame[column] = frame[column].fillna("").astype(str)
    return frame


def read_pure_nce(path: Path | None) -> dict[tuple[str, str], float]:
    if path is None or not Path(path).exists():
        return {}
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    frame = frame[~frame["nce_status"].isin(["no_observed_bins", "error"])]
    return {
        (row.x_compound_id, row.candidate_smiles): float(row.pure_nce_score)
        for row in frame.itertuples(index=False)
        if row.pure_nce_score != ""
    }


def _value(value: Any) -> Any:
    return None if isinstance(value, float) and math.isnan(value) else value


def compounds(frame: pd.DataFrame, pure_nce: dict[tuple[str, str], float]) -> dict[str, dict[str, Any]]:
    """One record per compound: naming fields, truth basis, and its candidates."""

    out: dict[str, dict[str, Any]] = {}
    for compound, group in frame.sort_values(["x_compound_id", "rank"]).groupby("x_compound_id", sort=False):
        first = group.iloc[0]
        candidates = []
        for row in group.itertuples(index=False):
            if not row.candidate_smiles and not row.candidate_formula:
                continue
            candidates.append(
                {
                    "rank": _value(row.rank),
                    "smiles": row.candidate_smiles,
                    "formula": row.candidate_formula,
                    "evidence_confidence": _value(row.evidence_confidence),
                    "mass_match_score": _value(row.mass_match_score),
                    "truth_tanimoto": _value(row.truth_tanimoto),
                    "pure_nce_score": pure_nce.get((compound, row.candidate_smiles)),
                }
            )
        out[compound] = {
            "fields": {key: first.get(key, "") for key in NAMING_FIELDS},
            "truth_basis": first.get("truth_basis", ""),
            "reference_formula": first.get("reference_formula", ""),
            "candidates": candidates,
        }
    return out


# --- replay ------------------------------------------------------------------


def masked(fields: dict[str, Any]) -> dict[str, Any]:
    """Naming and curator fields cleared, so only the candidates decide the call."""

    return {key: "" for key in fields}


def _same_formula(left: str, right: str) -> bool:
    left, right = normalise_formula(left or ""), normalise_formula(right or "")
    return bool(left) and left == right


def call_is_correct(
    source: str, candidates: list[dict[str, Any]], reference_formula: str
) -> tuple[bool | None, str]:
    """Whether rank 1 is the true structure, and on what basis (None when unanswerable).

    For a structure control (``source == "masked_level1"``) the bar is the
    rank-1 candidate's own precomputed ``truth_tanimoto`` against
    ``TRUTH_THRESHOLD`` (0.70) -- the same test ``matched_truth`` uses
    everywhere else in this delivery, applied uniformly whatever the
    replayed tier claims. A ``formula_only`` control has no reference
    structure, so the only question it can answer is whether rank 1's
    formula matches, a weaker label recorded as ``BASIS_FORMULA``.
    """

    ranked = sorted(candidates, key=R._rank)
    rank1 = ranked[0] if ranked else {}
    if source == "masked_level1":
        tanimoto = R._num(rank1.get("truth_tanimoto"))
        if tanimoto is None:
            return None, ""
        return tanimoto >= TRUTH_THRESHOLD, BASIS_STRUCTURE
    references = [part for part in reference_formula.split(";") if part.strip()]
    if references and rank1:
        return any(_same_formula(R._clean(rank1.get("formula")), ref) for ref in references), BASIS_FORMULA
    return None, ""


def replay(records: dict[str, dict[str, Any]]) -> pd.DataFrame:
    rows = []
    for compound, record in records.items():
        basis = record["truth_basis"]
        if basis == "structure_tanimoto":
            fields = masked(record["fields"])
            source = "masked_level1"
        elif basis == "formula_only":
            fields = record["fields"]
            source = "formula_only"
        else:
            continue
        candidates = record["candidates"]
        call = R.reported_identification(compound, candidates, fields)
        tier = R._confidence_tier(call)
        correct, label_basis = call_is_correct(source, candidates, record["reference_formula"])
        rows.append(
            {
                "x_compound_id": compound,
                "source": source,
                "replayed_call": call,
                "tier": tier,
                "margin": R._cluster_identity_confidence(candidates),
                "margin_basis": "pure_nce" if any(c.get("pure_nce_score") is not None for c in candidates) else "evidence",
                "n_candidates": len(candidates),
                "correct": correct,
                "label_basis": label_basis,
            }
        )
    return pd.DataFrame(rows)


# --- statistics --------------------------------------------------------------


def wilson(successes: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n <= 0:
        return 0.0, 1.0
    p = successes / n
    denominator = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denominator
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def beta_shrink(successes: int, n: int, prior: float, strength: float = PRIOR_STRENGTH) -> float:
    """Posterior mean of a Beta(strength*prior, strength*(1-prior)) prior."""

    return (successes + strength * prior) / (n + strength)


def pav_nonincreasing(values: list[float], weights: list[float]) -> list[float]:
    """Weighted pooled-adjacent-violators fit, non-increasing along the list."""

    blocks: list[list[float]] = []  # [weighted mean, weight, count]
    for value, weight in zip(values, weights):
        blocks.append([value, weight, 1])
        while len(blocks) > 1 and blocks[-2][0] < blocks[-1][0]:
            mean_b, weight_b, count_b = blocks.pop()
            mean_a, weight_a, count_a = blocks.pop()
            total = weight_a + weight_b
            blocks.append([(mean_a * weight_a + mean_b * weight_b) / total, total, count_a + count_b])
    out: list[float] = []
    for mean, _weight, count in blocks:
        out.extend([mean] * count)
    return out


def calibration_table(rows: pd.DataFrame) -> pd.DataFrame:
    labelled = rows[rows["correct"].notna()] if len(rows) else rows
    out = []
    for tier in R.CALIBRATED_TIER_ORDER:
        for basis in (BASIS_STRUCTURE, BASIS_FORMULA):
            group = labelled[(labelled["tier"] == tier) & (labelled["label_basis"] == basis)] if len(labelled) else labelled
            n = int(len(group))
            hits = int(group["correct"].astype(bool).sum()) if n else 0
            lo, hi = wilson(hits, n)
            out.append(
                {
                    "tier": tier,
                    "label_basis": basis,
                    "n": n,
                    "n_correct": hits,
                    "rate": round(hits / n, 4) if n else None,
                    "lo": round(lo, 4) if n else None,
                    "hi": round(hi, 4) if n else None,
                }
            )
    return pd.DataFrame(out)


def _structure_counts(rows: pd.DataFrame, tier: str) -> tuple[int, int]:
    if not len(rows):
        return 0, 0
    group = rows[(rows["tier"] == tier) & (rows["label_basis"] == BASIS_STRUCTURE) & rows["correct"].notna()]
    return int(len(group)), int(group["correct"].astype(bool).sum())


def fit_anchors(rows: pd.DataFrame) -> dict[str, dict[str, Any]]:
    """Per-tier anchor (percent), shrunk toward the hand-set band midpoint, PAV-monotone.

    Only structure-labelled tiers with at least ``MIN_LABELS_TO_CALIBRATE``
    controls are calibrated; the rest keep the hand-set midpoint and enter
    PAV with the prior weight alone.
    """

    raw, weights, info = [], [], {}
    for tier in R.CALIBRATED_TIER_ORDER:
        floor, ceiling = R._CONFIDENCE_BANDS[tier]
        prior = (floor + ceiling) / 200.0
        n, hits = _structure_counts(rows, tier)
        calibrated = n >= MIN_LABELS_TO_CALIBRATE
        posterior = beta_shrink(hits, n, prior) if calibrated else prior
        raw.append(posterior)
        weights.append((n if calibrated else 0) + PRIOR_STRENGTH)
        info[tier] = {"n": n, "n_correct": hits, "prior": round(100 * prior, 2), "calibrated": calibrated}
    monotone = pav_nonincreasing(raw, weights)
    for tier, posterior, anchor in zip(R.CALIBRATED_TIER_ORDER, raw, monotone):
        info[tier]["posterior"] = round(100 * posterior, 2)
        info[tier]["anchor"] = round(100 * anchor, 2)
    return info


def fit_bands(info: dict[str, dict[str, Any]]) -> None:
    """Band = anchor +/- Wilson half-width, clipped halfway to neighbouring anchors.

    Uncalibrated tiers keep the hand-set half-width. The clip makes bands of
    adjacent tiers disjoint, so a lower call can never outscore a higher one.
    """

    order = list(R.CALIBRATED_TIER_ORDER)
    for index, tier in enumerate(order):
        spec = info[tier]
        anchor = spec["anchor"]
        if spec["calibrated"]:
            lo, hi = wilson(spec["n_correct"], spec["n"])
            half = 50.0 * (hi - lo)
        else:
            floor, ceiling = R._CONFIDENCE_BANDS[tier]
            half = (ceiling - floor) / 2.0
        upper = R.LEVEL_2_STAR_FLOOR - 0.1 if index == 0 else (anchor + info[order[index - 1]]["anchor"]) / 2.0
        lower = 0.0 if index == len(order) - 1 else (anchor + info[order[index + 1]]["anchor"]) / 2.0
        spec["ceiling"] = round(min(anchor + half, upper), 2)
        spec["floor"] = round(max(anchor - half, lower), 2)
        spec["floor"] = min(spec["floor"], anchor)
        spec["ceiling"] = max(spec["ceiling"], anchor)


def _sigmoid(value: float) -> float:
    return 1.0 / (1.0 + math.exp(-value))


def _logit(p: float) -> float:
    p = min(max(p, 1e-4), 1 - 1e-4)
    return math.log(p / (1 - p))


def fit_slope(rows: pd.DataFrame, info: dict[str, dict[str, Any]], centers: dict[str, float]) -> float | None:
    """One logistic slope on the tier-centred margin, by grid MLE; None if not supported."""

    labelled = rows[(rows["label_basis"] == BASIS_STRUCTURE) & rows["correct"].notna() & rows["margin"].notna()]
    labelled = labelled[labelled["tier"].isin(list(info))]
    if len(labelled) < MIN_LABELS_FOR_SLOPE:
        return None
    best, best_ll = 0.0, None
    for step in range(0, 201):
        slope = step * 0.1
        ll = 0.0
        for row in labelled.itertuples(index=False):
            offset = (float(row.margin) - centers.get(row.tier, 50.0)) / 100.0
            p = _sigmoid(_logit(info[row.tier]["anchor"] / 100.0) + slope * offset)
            p = min(max(p, 1e-6), 1 - 1e-6)
            ll += math.log(p) if bool(row.correct) else math.log(1 - p)
        if best_ll is None or ll > best_ll:
            best, best_ll = slope, ll
    return round(best, 2) if best > 0 else None


def tier_centers(records: dict[str, dict[str, Any]]) -> dict[str, float]:
    """Median margin of the delivered (unmasked) compounds in each tier."""

    margins: dict[str, list[float]] = {tier: [] for tier in R.CALIBRATED_TIER_ORDER}
    for compound, record in records.items():
        if R._is_level1(record["fields"]):
            continue
        tier = R._confidence_tier(R.reported_identification(compound, record["candidates"], record["fields"]))
        margin = R._cluster_identity_confidence(record["candidates"])
        if tier in margins and margin is not None:
            margins[tier].append(margin)
    out = {}
    for tier, values in margins.items():
        values = sorted(values)
        out[tier] = round(values[len(values) // 2], 1) if values else 50.0
    return out


def build_calibration(
    rows: pd.DataFrame, records: dict[str, dict[str, Any]], normalization: str = "relative"
) -> dict[str, Any]:
    info = fit_anchors(rows)
    fit_bands(info)
    centers = tier_centers(records)
    for tier in info:
        info[tier]["center"] = centers[tier]
    labelled = rows[(rows["label_basis"] == BASIS_STRUCTURE) & rows["correct"].notna()] if len(rows) else rows
    return {
        "version": 1,
        "nce_normalization": normalization,
        "prior_strength": PRIOR_STRENGTH,
        "n_labelled": int(len(labelled)),
        "slope": fit_slope(rows, info, centers),
        "tiers": info,
    }


def loo_brier(rows: pd.DataFrame) -> dict[str, float | None]:
    """Leave-one-out Brier score of the calibrated anchors versus the hand-set midpoints."""

    labelled = rows[(rows["label_basis"] == BASIS_STRUCTURE) & rows["correct"].notna()]
    labelled = labelled[labelled["tier"].isin(list(R.CALIBRATED_TIER_ORDER))].reset_index(drop=True)
    if not len(labelled):
        return {"n": 0, "calibrated": None, "hand_set": None}
    calibrated, hand = [], []
    for index, row in labelled.iterrows():
        info = fit_anchors(labelled.drop(index))
        y = 1.0 if bool(row["correct"]) else 0.0
        calibrated.append((info[row["tier"]]["anchor"] / 100.0 - y) ** 2)
        floor, ceiling = R._CONFIDENCE_BANDS[row["tier"]]
        hand.append(((floor + ceiling) / 200.0 - y) ** 2)
    return {
        "n": int(len(labelled)),
        "calibrated": round(sum(calibrated) / len(calibrated), 4),
        "hand_set": round(sum(hand) / len(hand), 4),
    }


# --- scoring a delivery ------------------------------------------------------


def score_delivery(records: dict[str, dict[str, Any]]) -> pd.DataFrame:
    """Reported level, tier, margin, and confidence for every compound."""

    out = []
    for compound, record in records.items():
        call = R.reported_identification(compound, record["candidates"], record["fields"])
        out.append(
            {
                "x_compound_id": compound,
                "reported_level": R.reported_level(call),
                "reported_identification": call,
                "tier": R._confidence_tier(call),
                "margin": R._cluster_identity_confidence(record["candidates"]),
                "candidate_identity_confidence": R.candidate_identity_confidence(record["candidates"], record["fields"]),
            }
        )
    return pd.DataFrame(out)


def spread_and_overlap(scored: pd.DataFrame, bands: dict[str, tuple[float, float]]) -> dict[str, Any]:
    """Share near each band's midpoint, and how often a lower tier outscores a higher one."""

    spread = {}
    for tier, (floor, ceiling) in bands.items():
        values = scored.loc[scored["tier"] == tier, "candidate_identity_confidence"].dropna().astype(float)
        if not len(values):
            spread[tier] = {"n": 0, "near_midpoint": None, "min": None, "max": None, "iqr_share": None}
            continue
        mid = (floor + ceiling) / 2.0
        width = ceiling - floor
        iqr = float(values.quantile(0.75) - values.quantile(0.25))
        spread[tier] = {
            "n": int(len(values)),
            "near_midpoint": round(float((values.sub(mid).abs() <= 2.0).mean()), 3),
            "min": round(float(values.min()), 1),
            "max": round(float(values.max()), 1),
            "iqr_share": round(iqr / width, 3) if width > 0 else None,
        }

    def above(lower: str, upper: str) -> float | None:
        low = scored.loc[scored["tier"] == lower, "candidate_identity_confidence"].dropna().astype(float)
        high = scored.loc[scored["tier"] == upper, "candidate_identity_confidence"].dropna().astype(float)
        if not len(low) or not len(high):
            return None
        return round(float((low > high.min()).mean()), 3)

    return {
        "spread": spread,
        "level4_above_weakest_group_isomer": above("4", "3_probable_isomer_group"),
        "level5_above_lowest_level4": above("5", "4"),
    }


# --- commands ----------------------------------------------------------------


def _use_calibration(path: Path | None, normalization: str = "relative") -> None:
    """Load ``path``, or none; with none, margins still honour ``normalization``."""

    R.CALIBRATION_PATH = Path(path) if path else Path("/nonexistent/identity_confidence_calibration.json")
    R._CALIBRATION_CACHE.clear()
    if path is None and normalization != "relative":
        R._CALIBRATION_CACHE[str(R.CALIBRATION_PATH)] = {"nce_normalization": normalization, "tiers": {}}


def _records(args: argparse.Namespace) -> dict[str, dict[str, Any]]:
    return compounds(read_delivery(args.delivery), read_pure_nce(args.pure_nce))


def cmd_replay(args: argparse.Namespace) -> pd.DataFrame:
    _use_calibration(None, args.normalization)
    rows = replay(_records(args))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows.to_csv(args.out_dir / "calibration_rows.csv", index=False)
    table = calibration_table(rows)
    table.to_csv(args.out_dir / "calibration_table.csv", index=False)
    print(table.to_string(index=False))
    return rows


def cmd_fit(args: argparse.Namespace) -> dict[str, Any]:
    rows = cmd_replay(args)
    records = _records(args)
    calibration = build_calibration(rows, records, args.normalization)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(calibration, indent=2) + "\n", encoding="utf-8")
    R.load_calibration(args.json)
    _use_calibration(None)
    print(json.dumps({tier: {k: spec[k] for k in ("n", "anchor", "floor", "ceiling", "calibrated")}
                      for tier, spec in calibration["tiers"].items()}, indent=1))
    print(f"[calibrate] slope={calibration['slope']} wrote {args.json}")
    return calibration


def _bands_of(calibration: dict[str, Any] | None) -> dict[str, tuple[float, float]]:
    if calibration is None:
        return dict(R._CONFIDENCE_BANDS)
    return {tier: (spec["floor"], spec["ceiling"]) for tier, spec in calibration["tiers"].items()}


def cmd_report(args: argparse.Namespace) -> str:
    delivery = read_delivery(args.delivery)
    nce = read_pure_nce(args.pure_nce)
    _use_calibration(None)
    before = score_delivery(compounds(delivery, {}))
    nce_hand = score_delivery(compounds(delivery, nce))
    _use_calibration(args.json)
    calibration = R.active_calibration()
    after = score_delivery(compounds(delivery, nce))
    rows = pd.read_csv(args.out_dir / "calibration_rows.csv")
    rows["correct"] = rows["correct"].map({True: True, False: False, "True": True, "False": False})
    table = pd.read_csv(args.out_dir / "calibration_table.csv")
    brier = loo_brier(rows)
    stages = {
        "v2 (evidence margin, hand-set bands)": spread_and_overlap(before, dict(R._CONFIDENCE_BANDS)),
        "pure NCE margin, hand-set bands": spread_and_overlap(nce_hand, dict(R._CONFIDENCE_BANDS)),
        "pure NCE margin, calibrated": spread_and_overlap(after, _bands_of(calibration)),
    }
    level_diffs = int((before["reported_identification"] != after["reported_identification"]).sum())
    text = render_report(args, rows, table, calibration, stages, brier, level_diffs, after, nce)
    path = args.out_dir / "calibration_report.md"
    path.write_text(text, encoding="utf-8")
    print(f"[calibrate] wrote {path}")
    return text


def _pct(value: float | None) -> str:
    return "n/a" if value is None else f"{100 * value:.1f}%"


def render_report(args, rows, table, calibration, stages, brier, level_diffs, after, nce) -> str:
    lines = ["# candidate_identity_confidence calibration", ""]
    n_pairs = len(nce)
    normalization = (calibration or {}).get("nce_normalization", "relative")
    lines += [
        f"- Delivery: `{args.delivery}`; pure NCE scores: `{args.pure_nce}` ({n_pairs} scored compound-candidate pairs).",
        f"- Compounds with candidates, rescored: {after['x_compound_id'].nunique()}. Compounds without candidates keep "
        f"their earlier value. Reported calls changed by this work: {level_diffs} (must be 0).",
        f"- Margin: pure NCE, `{normalization}` normalization "
        + ("(rank-1 and rival scores on the matcher's own scale)." if normalization == "raw"
           else "(scores divided by the compound's best NCE score)."),
        f"- Replay: {int((rows['source'] == 'masked_level1').sum())} masked Level 1 structure controls and "
        f"{int((rows['source'] == 'formula_only').sum())} formula-only controls; "
        f"{int(rows['correct'].notna().sum())} carry a usable label.",
        "",
        "## Replay labels per tier",
        "",
        "| tier | label basis | n | correct | rate | Wilson 95% |",
        "|---|---|---|---|---|---|",
    ]
    for row in table.itertuples(index=False):
        if not row.n:
            continue
        lines.append(f"| {row.tier} | {row.label_basis} | {row.n} | {row.n_correct} | {_pct(row.rate)} | "
                     f"{_pct(row.lo)} - {_pct(row.hi)} |")
    lines += ["", "## Calibrated anchors", ""]
    if calibration is None:
        lines.append("No calibration file was loaded; scores use the hand-set bands.")
    else:
        lines += ["| tier | n labels | prior | posterior | anchor (PAV) | band | center margin | status |",
                  "|---|---|---|---|---|---|---|---|"]
        for tier, spec in calibration["tiers"].items():
            if spec["calibrated"]:
                status = "calibrated" if spec["n"] >= 10 else "calibrated, prior-dominated (fewer than 10 labels)"
            elif spec["n"]:
                status = f"uncalibrated: {spec['n']} labels, below {MIN_LABELS_TO_CALIBRATE}; hand-set midpoint"
            else:
                status = "uncalibrated: no labels; hand-set midpoint"
            if spec["anchor"] != spec["posterior"]:
                status += "; moved by PAV to keep level order"
            lines.append(f"| {tier} | {spec['n']} | {spec['prior']} | {spec['posterior']} | {spec['anchor']} | "
                         f"{spec['floor']}-{spec['ceiling']} | {spec['center']} | {status} |")
        lines += ["", f"Within-tier logistic slope on margin: {calibration['slope']} "
                  f"(fit only with at least {MIN_LABELS_FOR_SLOPE} pooled labels and a positive slope; "
                  "otherwise linear placement inside the band)."]
    lines += ["", "## Reliability (leave-one-out on the structure-labelled replay)", "",
              f"- Labelled controls: {brier['n']}",
              f"- Brier score, calibrated anchors: {brier['calibrated']}",
              f"- Brier score, hand-set band midpoints: {brier['hand_set']}", ""]
    lines += ["## Spread and overlap", "",
              "Near midpoint = share of a tier's compounds within 2 points of its band midpoint (target under 25%). "
              "It is not meaningful for a band narrower than about 8 points, where every value is near the "
              "midpoint; IQR share (interquartile range over band width) is the width-free spread measure.",
              ""]
    for stage, result in stages.items():
        lines += [f"### {stage}", "", "| tier | compounds | near midpoint | IQR share | min | max |",
                  "|---|---|---|---|---|---|"]
        for tier, spec in result["spread"].items():
            lines.append(f"| {tier} | {spec['n']} | {_pct(spec['near_midpoint'])} | {_pct(spec['iqr_share'])} | "
                         f"{spec['min']} | {spec['max']} |")
        lines += ["",
                  f"- Level 4 compounds above the weakest Level 3 group-isomer score: "
                  f"{_pct(result['level4_above_weakest_group_isomer'])} (target under 5%)",
                  f"- Level 5 compounds above the lowest Level 4 score: {_pct(result['level5_above_lowest_level4'])} "
                  "(target under 5%)", ""]
    lines += ["## Caveats", "",
              "- 'Correct' is the same bar for every structure-labelled tier: is rank 1 the true structure, by "
              "truth_tanimoto >= the published threshold (the same test behind matched_truth elsewhere). A Level "
              "4 or 5 call that is right on this bar is a right structure, not just a right formula, so its rate "
              "is a direct measure of the confidence this score reports. Where raw rates still come out of level "
              "order, pooling (PAV) shares an anchor across the affected tiers and the band clip keeps the "
              "higher level's band above the lower level's.",
              "- Uncalibrated tiers enter PAV with the prior weight only, so a calibrated neighbour can move their "
              "anchor; the report lists the prior, posterior, and pooled anchor for each.",
              "- No Level 2 or Level 3 name-based call (Probable name, isomer of, class phrase from curator text) can "
              "be replayed, because masking removes the name that produces it; those tiers start from the hand-set "
              "band midpoints.",
              "- The replay set is library compounds with stored-mass candidate lists; novel X compounds may have "
              "weaker lists, so calibrated rates may be optimistic for them.",
              ""]
    return "\n".join(lines)


def cmd_apply(args: argparse.Namespace) -> None:
    """Copy a delivery, add pure_nce_score, and recompute candidate_identity_confidence."""

    _use_calibration(args.json)
    if (args.dst / "exec_summary_top10.csv").exists():
        raise SystemExit(f"{args.dst} already holds a delivery; pick a new delivery directory")
    shutil.copytree(args.delivery, args.dst, dirs_exist_ok=True)
    nce = read_pure_nce(args.pure_nce)
    delivery = read_delivery(args.delivery)
    scored = score_delivery(compounds(delivery, nce)).set_index("x_compound_id")
    confidence = scored["candidate_identity_confidence"]
    for name in ("exec_summary_top10.csv", "exec_summary_top10_truncated.csv", "exec_summary_by_compound.csv"):
        path = args.dst / name
        frame = pd.read_csv(path, low_memory=False, dtype=str, keep_default_na=False)
        # Compounds with no candidate rows (unpublished) are not rescored; keep their value.
        rescored = frame["x_compound_id"].isin(confidence.index)
        frame.loc[rescored, "candidate_identity_confidence"] = frame.loc[rescored, "x_compound_id"].map(confidence).map(
            lambda value: "" if value is None or (isinstance(value, float) and math.isnan(value)) else f"{value}"
        )
        if "candidate_smiles" in frame.columns:
            frame["pure_nce_score"] = [
                "" if (cid, smiles) not in nce else f"{round(nce[(cid, smiles)], 6)}"
                for cid, smiles in zip(frame["x_compound_id"], frame["candidate_smiles"])
            ]
        frame.to_csv(path, index=False)
    rebuild_workbooks(args.dst, confidence, nce)
    print(f"[calibrate] wrote {args.dst} ({len(scored)} compounds)")


REVIEW_WORKBOOK = "executive_summary_review_columns.xlsx"


def rebuild_workbooks(dst: Path, confidence: pd.Series, nce: dict[tuple[str, str], float]) -> None:
    """Rebuild executive_summary_top10.xlsx and refresh the review-columns workbook in place."""

    from build_exec_summary_workbook import WORKBOOK_NAME, build_workbook

    build_workbook(dst, dst / WORKBOOK_NAME)
    review_path = dst / REVIEW_WORKBOOK
    if not review_path.exists():
        return
    review = pd.read_excel(review_path)
    top10 = pd.read_csv(dst / "exec_summary_top10.csv", dtype=str, keep_default_na=False,
                        usecols=["x_compound_id", "rank", "candidate_smiles"])
    smiles = {(cid, str(int(float(rank)))): s for cid, rank, s in top10.itertuples(index=False) if rank}
    rescored = review["x_compound_id"].isin(confidence.index)
    review.loc[rescored, "candidate_identity_confidence"] = pd.to_numeric(
        review.loc[rescored, "x_compound_id"].map(confidence), errors="coerce"
    )
    pure = []
    for cid, rank in zip(review["x_compound_id"], review["rank"]):
        key = (cid, smiles.get((cid, "" if pd.isna(rank) else str(int(rank))), ""))
        pure.append(nce.get(key))
    review = review.drop(columns=["pure_nce_score"], errors="ignore")
    at = review.columns.get_loc("fallback_score") + 1 if "fallback_score" in review.columns else len(review.columns)
    review.insert(at, "pure_nce_score", pd.to_numeric(pd.Series(pure, index=review.index), errors="coerce"))
    with pd.ExcelWriter(review_path, engine="openpyxl") as writer:
        review.to_excel(writer, sheet_name="Sheet1", index=False)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, func in (("replay", cmd_replay), ("fit", cmd_fit), ("report", cmd_report), ("apply", cmd_apply)):
        command = sub.add_parser(name)
        command.add_argument("--delivery", type=Path, required=True)
        command.add_argument("--pure-nce", type=Path, default=None)
        command.add_argument("--out-dir", type=Path, default=Path("output/identity_confidence_calibration"))
        command.add_argument("--json", type=Path, default=R.CALIBRATION_PATH)
        command.add_argument("--normalization", choices=("relative", "raw"), default="relative")
        if name == "apply":
            command.add_argument("--dst", type=Path, required=True)
        command.set_defaults(func=func)
    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
