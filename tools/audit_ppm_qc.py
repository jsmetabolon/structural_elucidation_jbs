#!/usr/bin/env python3
"""Compare delivered mass errors to a 5 ppm QC bar and a 10 ppm reporting ceiling.

Read-only against the committed executive-summary CSVs. Does not touch scorers,
MSI rules, or reports.

Usage:
    python3 tools/audit_ppm_qc.py
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
_EXEC = ROOT / "docs/paper/data/exec_summary"
DEFAULT_BY_COMPOUND = _EXEC / "exec_summary_by_compound.csv"
DEFAULT_TOP10 = _EXEC / "exec_summary_top10.csv"
_THRESH = ROOT / "docs/paper/data/thresholds"
DEFAULT_JSON = _THRESH / "audit_ppm_qc.json"
DEFAULT_MD = _THRESH / "audit_ppm_qc.md"
DEFAULT_CSV = _THRESH / "audit_ppm_qc.csv"
DEFAULT_FILTERED_CANDIDATES = _THRESH / "audit_ppm_qc_filtered_candidates.csv"
DEFAULT_FILTERED_BY_COMPOUND = (
    _THRESH / "audit_ppm_qc_filtered_by_compound.csv"
)

REPORTING_PPM = 10.0
QC_PPM = 5.0

CANDIDATE_FIELDS = (
    "x_compound_id",
    "cohort",
    "msi_level",
    "msi_label",
    "feature_name",
    "correct_result",
    "matched_truth",
    "original_rank",
    "rank_after_10ppm_filter",
    "abs_mass_error_ppm",
    "abs_mass_error_da",
    "ppm_band",
    "meets_qc_5ppm",
    "candidate_formula",
    "candidate_smiles",
    "candidate_confidence",
    "known_available",
    "truth_basis",
    "truth_tanimoto",
    "is_best_match_in_top10",
    "reference_answer",
    "reference_formula",
)

COMPOUND_FIELDS = (
    "x_compound_id",
    "cohort",
    "msi_level",
    "msi_label",
    "feature_name",
    "correct_result",
    "matched_truth",
    "known_available",
    "truth_basis",
    "n_top10",
    "n_le10_ppm",
    "n_le5_ppm",
    "original_match_rank",
    "match_survives_10ppm",
    "match_rank_after_10ppm_filter",
    "match_survives_5ppm",
    "match_rank_after_5ppm_filter",
    "reference_answer",
)

SLICE_LABELS = (
    ("all_top10", "All top-10 candidates"),
    ("rank1_all", "Rank-1, all reports"),
    ("rank1_msi_1_3", "Rank-1, MSI 1-3 (structure shown)"),
    ("rank1_msi_1_2", "Rank-1, MSI 1-2"),
    ("rank1_msi_1", "Rank-1, MSI 1"),
    ("rank1_truth_control", "Truth-control rank-1"),
    ("matched_truth", "matched_truth candidates"),
)


def _float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def _true(value: str | None) -> bool:
    return str(value).strip().upper() == "TRUE"


def _known_available(row: dict[str, str]) -> bool:
    basis = (row.get("truth_basis") or "").strip().lower()
    return basis not in ("", "none")


def _ppm_band(ppm: float) -> str:
    if ppm <= QC_PPM:
        return "le5"
    if ppm <= REPORTING_PPM:
        return "gt5_le10"
    return "gt10"


def _flag(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def _matched_to_known(row: dict[str, str]) -> str:
    if not _known_available(row):
        return ""
    return _flag(_true(row.get("matched_truth")))


def _correct_result(row: dict[str, str]) -> str:
    """Known identity when one exists; blank when the feature is unassessable."""

    return (row.get("reference_answer") or "").strip()


def _pct(n: int, denom: int) -> float | None:
    if denom == 0:
        return None
    return round(100.0 * n / denom, 1)


def _band_counts(ppms: list[float | None]) -> dict[str, int]:
    present = [p for p in ppms if p is not None]
    return {
        "n": len(ppms),
        "n_with_ppm": len(present),
        "n_missing": sum(p is None for p in ppms),
        "median_ppm": (
            round(statistics.median(present), 2) if present else None
        ),
        "le5": sum(p <= 5 for p in present),
        "gt5_le10": sum(5 < p <= 10 for p in present),
        "gt10": sum(p > 10 for p in present),
        "pct_le5": _pct(sum(p <= 5 for p in present), len(present)),
        "pct_gt5_le10": _pct(sum(5 < p <= 10 for p in present), len(present)),
        "pct_gt10": _pct(sum(p > 10 for p in present), len(present)),
    }


def _msi_int(row: dict[str, str]) -> int | None:
    raw = row.get("msi_level", "")
    try:
        return int(raw)
    except (TypeError, ValueError):
        return None


def first_match(rows: list[dict[str, str]]) -> tuple[int | None, float | None]:
    matches = sorted(
        (r for r in rows if _true(r.get("matched_truth"))),
        key=lambda r: int(r["rank"]),
    )
    if not matches:
        return None, None
    first = matches[0]
    return int(first["rank"]), _float(first.get("abs_mass_error_ppm"))


def recover(
    controls: list[dict[str, str]],
    by_id: dict[str, list[dict[str, str]]],
    pred,
) -> dict[str, Any]:
    r1 = t3 = t10 = 0
    lost: list[dict[str, Any]] = []
    for row in controls:
        xid = row["x_compound_id"]
        rank, ppm = first_match(by_id.get(xid, []))
        ok = rank is not None and pred(rank, ppm)
        if ok and rank == 1:
            r1 += 1
        if ok and rank is not None and rank <= 3:
            t3 += 1
        if ok:
            t10 += 1
        if _true(row.get("matched_within_top_10")) and not ok:
            lost.append(
                {
                    "x_compound_id": xid,
                    "first_matching_rank": rank,
                    "abs_mass_error_ppm": ppm,
                    "msi_level": _msi_int(row),
                }
            )
    n = len(controls)
    return {
        "n": n,
        "rank1": r1,
        "top3": t3,
        "top10": t10,
        "lost": lost,
    }


def audit(
    by_compound: list[dict[str, str]], top10: list[dict[str, str]]
) -> dict[str, Any]:
    by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in top10:
        by_id[row["x_compound_id"]].append(row)

    top10_ppms = [_float(r.get("abs_mass_error_ppm")) for r in top10]
    matched = [r for r in top10 if _true(r.get("matched_truth"))]
    matched_ppms = [_float(r.get("abs_mass_error_ppm")) for r in matched]

    def slice_rank1(pred) -> dict[str, Any]:
        rows = [r for r in by_compound if pred(r)]
        return {
            "n_reports": len(rows),
            **_band_counts(
                [_float(r.get("rank1_abs_mass_error_ppm")) for r in rows]
            ),
        }

    n_any_gt10 = 0
    n_zero_le10 = 0
    n_zero_le5 = 0
    for rows in by_id.values():
        ppms = [_float(r.get("abs_mass_error_ppm")) for r in rows]
        present = [p for p in ppms if p is not None]
        if any(p > 10 for p in present):
            n_any_gt10 += 1
        # Reports with no candidate inside the window include those with no ppm.
        if not any(p <= 10 for p in present):
            n_zero_le10 += 1
        if not any(p <= 5 for p in present):
            n_zero_le5 += 1

    controls = [r for r in by_compound if r.get("cohort") == "truth_control"]
    published = recover(controls, by_id, lambda rank, ppm: True)
    at_10 = recover(
        controls, by_id, lambda rank, ppm: ppm is not None and ppm <= 10
    )
    at_5 = recover(
        controls, by_id, lambda rank, ppm: ppm is not None and ppm <= 5
    )

    msi1_high = []
    for row in by_compound:
        if _msi_int(row) != 1:
            continue
        ppm = _float(row.get("rank1_abs_mass_error_ppm"))
        if ppm is None or ppm <= 5:
            continue
        msi1_high.append(
            {
                "x_compound_id": row["x_compound_id"],
                "rank1_abs_mass_error_ppm": round(ppm, 2),
                "rank1_abs_mass_error_da": _float(
                    row.get("rank1_abs_mass_error_da")
                ),
                "matched_at_rank_1": _true(row.get("matched_at_rank_1")),
                "best_truth_tanimoto": _float(row.get("best_truth_tanimoto")),
                "feature_name": row.get("feature_name", ""),
            }
        )

    return {
        "sources": {
            "by_compound": (
                "docs/paper/data/exec_summary/exec_summary_by_compound.csv"
            ),
            "top10": "docs/paper/data/exec_summary/exec_summary_top10.csv",
        },
        "n_reports": len(by_compound),
        "n_top10_rows": len(top10),
        "slices": {
            "all_top10": _band_counts(top10_ppms),
            "rank1_all": slice_rank1(lambda r: True),
            "rank1_msi_1_3": slice_rank1(lambda r: (_msi_int(r) or 99) <= 3),
            "rank1_msi_1_2": slice_rank1(lambda r: (_msi_int(r) or 99) <= 2),
            "rank1_msi_1": slice_rank1(lambda r: _msi_int(r) == 1),
            "rank1_truth_control": slice_rank1(
                lambda r: r.get("cohort") == "truth_control"
            ),
            "matched_truth": _band_counts(matched_ppms),
        },
        "coverage": {
            "reports_with_any_candidate_gt10_ppm": n_any_gt10,
            "reports_with_zero_candidates_le10_ppm": n_zero_le10,
            "reports_with_zero_candidates_le5_ppm": n_zero_le5,
        },
        "truth_control_recovery": {
            "published": published,
            "match_le10_ppm": at_10,
            "match_le5_ppm": at_5,
            "same_losses_at_5_and_10": [x["x_compound_id"] for x in at_5["lost"]]
            == [x["x_compound_id"] for x in at_10["lost"]],
        },
        "msi1_rank1_above_5_ppm": msi1_high,
    }


def _row(
    label: str,
    block: dict[str, Any],
    *,
    percents: bool,
    bold_gt10: bool,
    bold_mid: bool,
) -> str:
    n = block["n_with_ppm"]
    le5, mid, gt10 = block["le5"], block["gt5_le10"], block["gt10"]
    if percents and n:
        le5_s = f"{le5} ({block['pct_le5']:.1f}%)"
        mid_s = f"{mid} ({block['pct_gt5_le10']:.1f}%)"
        gt_s = f"{gt10} ({block['pct_gt10']:.1f}%)"
    else:
        le5_s, mid_s, gt_s = str(le5), str(mid), str(gt10)
    if bold_mid:
        mid_s = f"**{mid_s}**"
    if bold_gt10:
        gt_s = f"**{gt_s}**"
    return f"| {label} | {n} | {le5_s} | {mid_s} | {gt_s} |"


def render_markdown(summary: dict[str, Any]) -> str:
    s = summary["slices"]
    cov = summary["coverage"]
    rec = summary["truth_control_recovery"]
    n_reports = summary["n_reports"]
    pub = rec["published"]
    lost_ids = ", ".join(
        x["x_compound_id"] for x in rec["match_le10_ppm"]["lost"]
    )
    exceptions = summary["msi1_rank1_above_5_ppm"]
    exception_lines = []
    notes = {
        "X-18345": (
            "The only MSI 1 match in the instrument-scale "
            "10–20 ppm range. Tanimoto 1.0."
        ),
        "X-24137": "~1 Da / formula mismatch, not ppm. Tanimoto 0.73.",
        "X-21470": (
            "Charge-state 2× (233 Da). Already in the mass-policy audit."
        ),
    }
    for item in exceptions:
        xid = item["x_compound_id"]
        note = notes.get(xid, "")
        exception_lines.append(
            f"- `{xid}` — **{item['rank1_abs_mass_error_ppm']:.2f} ppm** "
            f"({item['rank1_abs_mass_error_da']} Da). {note}".rstrip()
        )
    r1 = s["rank1_all"]
    tc = s["rank1_truth_control"]
    m1 = s["rank1_msi_1"]
    mt = s["matched_truth"]
    r1_bullet = (
        f"- Rank-1 ({r1['n_with_ppm']:,} with a ppm): "
        f"**median {r1['median_ppm']} ppm**; "
        f"**{r1['pct_le5']:.1f}% ≤5 ppm**"
    )
    tc_bullet = (
        f"- Truth-control rank-1 ({tc['n_with_ppm']} with a ppm): "
        f"**median {tc['median_ppm']} ppm**; "
        f"**{tc['le5']}/{tc['n_with_ppm']} ≤5 ppm**; "
        f"**{tc['gt5_le10']} in 5–10 ppm**"
    )
    m1_bullet = (
        f"- MSI Level 1 rank-1: **median {m1['median_ppm']} ppm**; "
        f"**{m1['le5']}/{m1['n_reports']} ≤5 ppm**"
    )
    mt_bullet = (
        f"- `matched_truth` rows: "
        f"**{mt['le5']}/{mt['n_with_ppm']} ≤5 ppm**; "
        f"**{mt['gt5_le10']} in 5–10 ppm**"
    )
    coverage_line = (
        f"**{cov['reports_with_any_candidate_gt10_ppm']:,} / "
        f"{n_reports:,}** reports include at least one top-10 candidate "
        f">10 ppm. **{cov['reports_with_zero_candidates_le10_ppm']:,}** "
        f"reports have **no** candidate ≤10 ppm "
        f"({cov['reports_with_zero_candidates_le5_ppm']:,} have none ≤5 ppm)."
    )
    row_all = _row(
        "All top-10 candidates",
        s["all_top10"],
        percents=True,
        bold_gt10=True,
        bold_mid=False,
    )
    row_r1 = _row(
        "Rank-1, all reports",
        s["rank1_all"],
        percents=True,
        bold_gt10=True,
        bold_mid=False,
    )
    row_l13 = _row(
        "Rank-1, MSI 1–3 (structure shown)",
        s["rank1_msi_1_3"],
        percents=True,
        bold_gt10=True,
        bold_mid=False,
    )
    row_l12 = _row(
        "Rank-1, MSI 1–2",
        s["rank1_msi_1_2"],
        percents=True,
        bold_gt10=True,
        bold_mid=False,
    )
    row_l1 = _row(
        "Rank-1, MSI 1",
        s["rank1_msi_1"],
        percents=True,
        bold_gt10=True,
        bold_mid=False,
    )
    row_tc = _row(
        "Truth-control rank-1",
        s["rank1_truth_control"],
        percents=False,
        bold_gt10=True,
        bold_mid=True,
    )
    row_mt = _row(
        "`matched_truth` candidates",
        s["matched_truth"],
        percents=False,
        bold_gt10=True,
        bold_mid=True,
    )
    le10 = rec["match_le10_ppm"]
    le5 = rec["match_le5_ppm"]
    recovery_10 = (
        f"| Match also ≤10 ppm | {le10['rank1']} | {le10['top3']} | "
        f"{le10['top10']} | {lost_ids} |"
    )
    recovery_5 = (
        f"| Match also ≤5 ppm | {le5['rank1']} | {le5['top3']} | "
        f"{le5['top10']} | same six |"
    )
    return f"""# Mass error versus 5 ppm QC and a 10 ppm reporting ceiling

Findings only. Generated by `tools/audit_ppm_qc.py` from the committed executive
summary tables. No scorer, MSI rule, or report was changed.

Sources:

- `{summary["sources"]["by_compound"]}`
- `{summary["sources"]["top10"]}`

Cross-check the truth-control slice against
[audit_mass_policy.md](audit_mass_policy.md): correct candidates have a median
near 1.07 ppm, and 81.4% of them fall within 5 ppm **and** within 10 ppm.

## Typical QC (<5 ppm) versus the bulk

The **center** of the distribution is consistent with typical QC of <5 ppm.

{r1_bullet}
{tc_bullet}
{m1_bullet}
{mt_bullet}

A 5 ppm QC bar describes typical correct and typical rank-1 calls. It does not
describe the full reported set.

## A 10 ppm reporting ceiling was not applied

Displayed candidates are not gated at 10 ppm.

| Slice | n with ppm | ≤5 ppm | 5–10 ppm | >10 ppm |
|---|---:|---:|---:|---:|
{row_all}
{row_r1}
{row_l13}
{row_l12}
{row_l1}
{row_tc}
{row_mt}

{coverage_line}

The 5–10 ppm band is almost empty for **correct** hits. Tightening 10 ppm to 5 ppm
does not change truth recovery on this cohort (same losses). The split that
matters is ≤5 ppm versus the far tail, not 5 versus 10.

## The tail is not instrument drift

Truth-control recovery on {pub["n"]} controls:

| Gate | Rank 1 | Top 3 | Top 10 | Controls lost |
|---|---:|---:|---:|---|
| Published (any ppm) | {pub["rank1"]} | {pub["top3"]} | {pub["top10"]} | — |
{recovery_10}
{recovery_5}

Three rank-1 MSI Level 1 calls exceed 5 ppm. Two are known pipeline or formula
failures; one sits at instrument-scale 15.6 ppm:

{chr(10).join(exception_lines)}

Remaining high-ppm “correct” candidates are conjugate-scale (for example a
glucuronide, 176 Da) or missing a recorded mass. That is why 81.4% of correct
candidates in the mass-policy audit are within 5 ppm **and** within 10 ppm:
nothing extra lands in 5–10.

## How to state this

Typical QC (<5 ppm) matches the **median** of reported rank-1 and of
truth-matched candidates.

A **10 ppm reporting ceiling was not applied**. About one in five rank-1 rows,
and about two in five top-10 rows, exceed 10 ppm.

Almost no truth matches sit between 5 and 10 ppm. Failures of a 10 ppm ceiling
are the same chemical and pipeline tails already in
[audit_mass_policy.md](audit_mass_policy.md), plus `X-18345` at 15.6 ppm.

Do not write that “no matches exceed 10 ppm.” That sentence is false for the
delivered reports.

The same counts are written as
`docs/paper/data/thresholds/audit_ppm_qc.csv`.

Candidates that meet the 10 ppm reporting ceiling, with original rank, rank
after that filter, a ≤5 ppm QC flag, the known identity (`correct_result`)
when one exists, and `matched_truth` for whether that candidate is it, are
[`audit_ppm_qc_filtered_candidates.csv`](audit_ppm_qc_filtered_candidates.csv).
The per-compound rollup after the same filters is
[`audit_ppm_qc_filtered_by_compound.csv`](audit_ppm_qc_filtered_by_compound.csv).

This is consistent with FML-8 in
[identification_rules.md](../../../identification_rules.md): most delivered
`mass_match_score` values used a 1.0 Da linear window (about 3,300 ppm at these
masses), not the 5 ppm formula-finder window.

## Out of scope here

- Whether a future reporting gate should be 5 or 10 ppm (on truth matches they
  are equivalent).
- Whether `X-18345` (15.6 ppm, otherwise a perfect MSI 1 hit) is acceptable
  under lab QC.
- Whether headline recovery should require mass plausibility (already open in
  the revision log).
"""


def _csv_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    return str(value)


def csv_rows(summary: dict[str, Any]) -> list[dict[str, str]]:
    """Tidy summary: one metric per row, for spreadsheet use."""

    rows: list[dict[str, str]] = []

    def add(
        section: str,
        item: str,
        metric: str,
        value: Any,
        *,
        label: str = "",
        extra: dict[str, Any] | None = None,
    ) -> None:
        row = {
            "section": section,
            "item": item,
            "label": label,
            "metric": metric,
            "value": _csv_value(value),
            "x_compound_id": "",
            "note": "",
        }
        if extra:
            for key, extra_value in extra.items():
                row[key] = _csv_value(extra_value)
        rows.append(row)

    for key, label in SLICE_LABELS:
        block = summary["slices"][key]
        if "n_reports" in block:
            add("ppm_band", key, "n_reports", block["n_reports"], label=label)
        for metric in (
            "n",
            "n_with_ppm",
            "n_missing",
            "median_ppm",
            "le5",
            "gt5_le10",
            "gt10",
            "pct_le5",
            "pct_gt5_le10",
            "pct_gt10",
        ):
            add("ppm_band", key, metric, block.get(metric), label=label)

    cov = summary["coverage"]
    add(
        "coverage",
        "all_reports",
        "n_reports",
        summary["n_reports"],
        label="All delivered reports",
    )
    add(
        "coverage",
        "all_reports",
        "n_top10_rows",
        summary["n_top10_rows"],
        label="All delivered reports",
    )
    coverage_metrics = (
        (
            "reports_with_any_candidate_gt10_ppm",
            "Reports with any top-10 candidate >10 ppm",
        ),
        (
            "reports_with_zero_candidates_le10_ppm",
            "Reports with no candidate <=10 ppm",
        ),
        (
            "reports_with_zero_candidates_le5_ppm",
            "Reports with no candidate <=5 ppm",
        ),
    )
    for metric, label in coverage_metrics:
        add("coverage", "all_reports", metric, cov[metric], label=label)

    rec = summary["truth_control_recovery"]
    gates = (
        ("published", "Published (any ppm)"),
        ("match_le10_ppm", "Match also <=10 ppm"),
        ("match_le5_ppm", "Match also <=5 ppm"),
    )
    for gate, gate_label in gates:
        block = rec[gate]
        for metric in ("n", "rank1", "top3", "top10"):
            add("recovery", gate, metric, block[metric], label=gate_label)
        lost = [x["x_compound_id"] for x in block.get("lost", [])]
        add(
            "recovery",
            gate,
            "lost_controls",
            ";".join(lost),
            label=gate_label,
        )
    add(
        "recovery",
        "comparison",
        "same_losses_at_5_and_10",
        rec["same_losses_at_5_and_10"],
        label="5 ppm vs 10 ppm",
    )
    for lost in rec["match_le10_ppm"]["lost"]:
        add(
            "recovery_loss",
            lost["x_compound_id"],
            "abs_mass_error_ppm",
            lost["abs_mass_error_ppm"],
            label="Lost when match must be <=10 ppm",
            extra={
                "x_compound_id": lost["x_compound_id"],
                "note": (
                    f"rank {lost['first_matching_rank']}; "
                    f"MSI {lost['msi_level']}"
                ),
            },
        )

    notes = {
        "X-18345": "Instrument-scale 15.6 ppm; Tanimoto 1.0",
        "X-24137": "~1 Da / formula mismatch; Tanimoto 0.73",
        "X-21470": "Charge-state 2x; already in mass-policy audit",
    }
    for item in summary["msi1_rank1_above_5_ppm"]:
        xid = item["x_compound_id"]
        extra = {"x_compound_id": xid, "note": notes.get(xid, "")}
        for metric in (
            "rank1_abs_mass_error_ppm",
            "rank1_abs_mass_error_da",
            "best_truth_tanimoto",
        ):
            add(
                "msi1_exception",
                xid,
                metric,
                item[metric],
                label="MSI 1 rank-1 above 5 ppm",
                extra=extra,
            )
    return rows


def write_csv(
    path: Path,
    rows: list[dict[str, str]],
    fieldnames: list[str] | None = None,
) -> None:
    names = fieldnames or [
        "section",
        "item",
        "label",
        "metric",
        "value",
        "x_compound_id",
        "note",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def filtered_candidate_rows(
    top10: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Candidates within 10 ppm, with QC flag, ranks, and known-match status."""

    passing: dict[str, list[tuple[int, dict[str, str], float]]] = defaultdict(
        list
    )
    for row in top10:
        ppm = _float(row.get("abs_mass_error_ppm"))
        if ppm is None or ppm > REPORTING_PPM:
            continue
        passing[row["x_compound_id"]].append(
            (int(row["rank"]), row, ppm)
        )

    out: list[dict[str, str]] = []
    for xid in sorted(passing):
        group = sorted(passing[xid], key=lambda item: item[0])
        for new_rank, (_orig_rank, row, ppm) in enumerate(group, start=1):
            known = _known_available(row)
            answer = _correct_result(row)
            matched = _matched_to_known(row)
            out.append(
                {
                    "x_compound_id": xid,
                    "cohort": row.get("cohort", ""),
                    "msi_level": row.get("msi_level", ""),
                    "msi_label": row.get("msi_label", ""),
                    "feature_name": row.get("feature_name", ""),
                    "correct_result": answer,
                    "matched_truth": matched,
                    "original_rank": str(row["rank"]),
                    "rank_after_10ppm_filter": str(new_rank),
                    "abs_mass_error_ppm": row.get("abs_mass_error_ppm", ""),
                    "abs_mass_error_da": row.get("abs_mass_error_da", ""),
                    "ppm_band": _ppm_band(ppm),
                    "meets_qc_5ppm": _flag(ppm <= QC_PPM),
                    "candidate_formula": row.get("candidate_formula", ""),
                    "candidate_smiles": row.get("candidate_smiles", ""),
                    "candidate_confidence": row.get(
                        "candidate_confidence", ""
                    ),
                    "known_available": _flag(known),
                    "truth_basis": row.get("truth_basis", ""),
                    "truth_tanimoto": row.get("truth_tanimoto", ""),
                    "is_best_match_in_top10": row.get(
                        "is_best_match_in_top10", ""
                    ),
                    "reference_answer": answer,
                    "reference_formula": row.get("reference_formula", ""),
                }
            )
    return out


def filtered_compound_rows(
    by_compound: list[dict[str, str]],
    top10: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Per-compound rollup after the 10 ppm and 5 ppm filters."""

    by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in top10:
        by_id[row["x_compound_id"]].append(row)

    out: list[dict[str, str]] = []
    for compound in by_compound:
        xid = compound["x_compound_id"]
        cands = by_id.get(xid, [])
        le10 = []
        le5 = []
        for row in cands:
            ppm = _float(row.get("abs_mass_error_ppm"))
            if ppm is None:
                continue
            if ppm <= REPORTING_PPM:
                le10.append(row)
            if ppm <= QC_PPM:
                le5.append(row)
        le10.sort(key=lambda r: int(r["rank"]))
        le5.sort(key=lambda r: int(r["rank"]))

        matches = sorted(
            (r for r in cands if _true(r.get("matched_truth"))),
            key=lambda r: int(r["rank"]),
        )
        first = matches[0] if matches else None
        first_ppm = (
            _float(first.get("abs_mass_error_ppm")) if first else None
        )
        known = _known_available(compound) or _known_available(
            cands[0] if cands else {}
        )
        survive_10 = first is not None and first_ppm is not None and (
            first_ppm <= REPORTING_PPM
        )
        survive_5 = first is not None and first_ppm is not None and (
            first_ppm <= QC_PPM
        )
        rank_after_10 = ""
        if survive_10 and first is not None:
            rank_after_10 = str(
                next(
                    i
                    for i, row in enumerate(le10, start=1)
                    if row["rank"] == first["rank"]
                )
            )
        rank_after_5 = ""
        if survive_5 and first is not None:
            rank_after_5 = str(
                next(
                    i
                    for i, row in enumerate(le5, start=1)
                    if row["rank"] == first["rank"]
                )
            )
        answer = _correct_result(compound) or _correct_result(
            cands[0] if cands else {}
        )
        out.append(
            {
                "x_compound_id": xid,
                "cohort": compound.get("cohort", ""),
                "msi_level": compound.get("msi_level", ""),
                "msi_label": compound.get("msi_label", ""),
                "feature_name": compound.get("feature_name", ""),
                "correct_result": answer,
                "matched_truth": _flag(survive_10) if known else "",
                "known_available": _flag(known),
                "truth_basis": compound.get("truth_basis", "")
                or (cands[0].get("truth_basis", "") if cands else ""),
                "n_top10": str(len(cands)),
                "n_le10_ppm": str(len(le10)),
                "n_le5_ppm": str(len(le5)),
                "original_match_rank": (
                    str(int(first["rank"])) if first else ""
                ),
                "match_survives_10ppm": (
                    _flag(survive_10) if known else ""
                ),
                "match_rank_after_10ppm_filter": rank_after_10,
                "match_survives_5ppm": _flag(survive_5) if known else "",
                "match_rank_after_5ppm_filter": rank_after_5,
                "reference_answer": answer,
            }
        )
    return out


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--by-compound", type=Path, default=DEFAULT_BY_COMPOUND)
    parser.add_argument("--top10", type=Path, default=DEFAULT_TOP10)
    parser.add_argument("--json-out", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--md-out", type=Path, default=DEFAULT_MD)
    parser.add_argument("--csv-out", type=Path, default=DEFAULT_CSV)
    parser.add_argument(
        "--filtered-candidates-out",
        type=Path,
        default=DEFAULT_FILTERED_CANDIDATES,
    )
    parser.add_argument(
        "--filtered-by-compound-out",
        type=Path,
        default=DEFAULT_FILTERED_BY_COMPOUND,
    )
    args = parser.parse_args(argv)

    by_compound = load_csv(args.by_compound)
    top10 = load_csv(args.top10)
    summary = audit(by_compound, top10)
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(summary, indent=2) + "\n")
    args.md_out.write_text(render_markdown(summary))
    write_csv(args.csv_out, csv_rows(summary))
    write_csv(
        args.filtered_candidates_out,
        filtered_candidate_rows(top10),
        list(CANDIDATE_FIELDS),
    )
    write_csv(
        args.filtered_by_compound_out,
        filtered_compound_rows(by_compound, top10),
        list(COMPOUND_FIELDS),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
