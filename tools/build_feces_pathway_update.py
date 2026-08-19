#!/usr/bin/env python3
"""Re-weight the chemical-pathway composition of human feces with X-compounds.

The panel's own view of what feces contains is built from named metabolites.
Every X-compound is a reproducibly detected feature that the panel cannot name,
so it is absent from that view by construction. This adds the pathway
assignments the pipeline made for those features and reports how the
composition changes.

Three inputs, each supplying one part:

- The Microbiome-vs-Global workbook, sheet ``GLOBAL 20260129``. Detection is
  read from the **GDP** (Global Discovery Panel) column, ``GDP Bucket - hFecal``,
  and pathway from ``Orig SuperPathway`` / ``Orig SubPathway``. The MFP columns
  beside them describe the narrower microbiome panel and are not used.
- The human-feces X-compound detection history, which says which X-compounds are
  seen in feces at all and how often.
- The delivered report bundle, which supplies each X-compound's predicted
  pathway.

What the result is not
----------------------

These are predicted pathways, not measured ones. Super-pathway agreement
against independent labels was measured at 52.2% against a ceiling of 62.5%
(see ``docs/paper/data/pathway/audit_pathway_summary.json``), so roughly half of
the added assignments are expected to be wrong at the class level. The
``--stratify`` output splits by the pipeline's own confidence tier, which is the
honest way to read the shift. Treat the updated share as an indication of
direction and rough size, not as a measurement.

Two structural caveats that no confidence tier fixes:

- The vocabularies are not the same size. The pipeline can only emit pathways it
  was trained to emit, so classes it never predicts are diluted purely by the
  denominator growing. Those classes are flagged in the output.
- Compounds are counted once each, regardless of how often they are detected.
  ``--weight-by-detection`` reports the alternative: the expected number of
  metabolites detected in a pathway, the sum of each compound's detection rate
  (knowns plus pathway-annotated X-compounds).

Usage:
    python3 tools/build_feces_pathway_update.py \
        --global-workbook "20260130 Microbiome vs. Global Detection (9).xlsx" \
        --detection-history "20250728 hFeces X-compound (Unknown) Detection History.xlsx" \
        --reports-root /path/to/reports
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Iterable

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

csv.field_size_limit(10**9)

PROJECT_DATA = PROJECT_ROOT / "data"
# Curated pathway labels for X-compounds live here, alongside their Metabolon
# identification level.
KNOWN_X_CSV = PROJECT_DATA / "known_X.csv"

GLOBAL_SHEET = "GLOBAL 20260129"
DETECTION_SHEET = "unk_feces_detection_rates"

# The Global sheet carries a merged group header above the real one.
GLOBAL_HEADER_ROW = 2

SUPER_COLUMN = "Orig SuperPathway"
SUB_COLUMN = "Orig SubPathway"
# GDP is the Global Discovery Panel. The MFP columns beside it describe the
# narrower microbiome-focused panel; "detected on global" means this one.
GDP_BUCKET_COLUMN = "GDP Bucket - hFecal"
GDP_RATE_COLUMN = "hFecal HD4"

NEVER = "Never"
COMMON_BUCKETS = ("Often (50-90%)", "Very Often (>90%)")
# The rate that puts an X-compound in the same company as those buckets.
COMMON_DETECTION_FLOOR = 0.50

# The sheet spells the same class two ways. Collapsing them is a presentation
# fix, not a taxonomy decision: left apart, one class reports as two.
CANONICAL_PATHWAY = {"secondary metabolism": "Secondary metabolism"}


def canonical_pathway(name: Any) -> str:
    text = "" if name is None else str(name).strip()
    if not text or text.lower() in {"nan", "none"}:
        return "(unassigned)"
    return CANONICAL_PATHWAY.get(text.lower(), text)


def _clean(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    return "" if text.lower() in {"nan", "none"} else text


def _num(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if number != number else number


def read_global_panel(path: Path) -> list[dict[str, Any]]:
    """Named metabolites with their pathway and GDP fecal detection."""

    import openpyxl

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    if GLOBAL_SHEET not in workbook.sheetnames:
        raise SystemExit(f"{path.name} has no sheet {GLOBAL_SHEET!r}")
    sheet = workbook[GLOBAL_SHEET]
    rows = sheet.iter_rows(values_only=True)
    for _ in range(GLOBAL_HEADER_ROW - 1):
        next(rows)
    header = [_clean(cell) for cell in next(rows)]
    index = {name: position for position, name in enumerate(header) if name}
    for required in (SUPER_COLUMN, SUB_COLUMN, GDP_BUCKET_COLUMN, GDP_RATE_COLUMN):
        if required not in index:
            raise SystemExit(f"{path.name}/{GLOBAL_SHEET} has no column {required!r}")

    records = []
    for row in rows:
        if not row or row[0] is None:
            continue
        bucket = _clean(row[index[GDP_BUCKET_COLUMN]]) or NEVER
        records.append(
            {
                "chemical_id": _clean(row[0]),
                "metabolite": _clean(row[1]),
                "super_pathway": canonical_pathway(row[index[SUPER_COLUMN]]),
                "sub_pathway": canonical_pathway(row[index[SUB_COLUMN]]),
                "gdp_bucket": bucket,
                "gdp_rate": _num(row[index[GDP_RATE_COLUMN]]),
            }
        )
    workbook.close()
    return records


def read_inclusion_list(path: Path, sheet: str = "Fecal Inclusion List") -> set[str]:
    """CHEMICAL_IDs on the curated fecal panel.

    A narrower and more purposeful denominator than every metabolite the Global
    panel has ever seen in feces: it is what the panel is actually meant to
    carry.
    """

    import openpyxl

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    if sheet not in workbook.sheetnames:
        raise SystemExit(f"{path.name} has no sheet {sheet!r}")
    rows = workbook[sheet].iter_rows(values_only=True)
    identifiers = set()
    for row in rows:
        identifier = _clean(row[0] if row else None)
        # The sheet opens with a spacer row and a header row; both fail this.
        if identifier and identifier.upper() != "CHEMICAL_ID":
            identifiers.add(identifier)
    workbook.close()
    return identifiers


def read_detection_history(path: Path) -> dict[str, float]:
    """X-compound to fecal detection rate, as a fraction."""

    import openpyxl

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    sheet = workbook[DETECTION_SHEET] if DETECTION_SHEET in workbook.sheetnames else workbook.active
    rows = sheet.iter_rows(values_only=True)
    next(rows)
    detected: dict[str, float] = {}
    for row in rows:
        compound = _clean(row[0] if row else None)
        if not compound:
            continue
        percent = _num(row[2]) if len(row) > 2 else None
        detected[compound] = (percent or 0.0) / 100.0
    workbook.close()
    return detected


def read_curated_pathways(path: Path, vocabulary: set[str]) -> dict[str, dict[str, str]]:
    """Curated pathway labels for X-compounds, where a curator has recorded one.

    These are the *known* assignments. Where one exists there is no reason to use
    a prediction instead, and good reason not to: the two agree less than half
    the time, which is the predictor's measured accuracy rather than a surprise.

    A fifth of the column is not taxonomy at all but review notes -- ``Unknown``,
    ``DO NOT INCLUDE``, ``Ask DTS`` -- which is the same defect the manuscript
    already records for the subpathway labels. Anything outside the reference
    vocabulary is rejected here and the compound falls back to its prediction,
    so a curator's open question cannot become a chemical class.
    """

    if not path.exists():
        return {}
    curated: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            compound = _clean(row.get("Unknown"))
            super_pathway = _clean(row.get("SuperPathway"))
            if not compound or not super_pathway:
                continue
            if canonical_pathway(super_pathway) not in vocabulary:
                continue
            # First non-empty wins, matching how the delivered set resolves a
            # compound that appears once per library.
            curated.setdefault(
                compound,
                {
                    "super_pathway": canonical_pathway(super_pathway),
                    "sub_pathway": canonical_pathway(row.get("Draft Subpathway")),
                },
            )
    return curated


def read_pathway_assignments(
    root: Path, wanted: Iterable[str], curated: dict[str, dict[str, str]] | None = None
) -> dict[str, dict[str, str]]:
    """Each X-compound's pathway: curated where one exists, predicted otherwise."""

    curated = curated or {}
    assignments: dict[str, dict[str, str]] = {}
    for compound in sorted(set(wanted)):
        candidates = root / compound / "data" / f"{compound.lower()}_top_candidates.csv"
        if not candidates.exists():
            continue
        with candidates.open(newline="") as handle:
            row = next(csv.DictReader(handle), None)
        if row is None:
            continue
        known = curated.get(compound)
        assignments[compound] = {
            "super_pathway": (known or {}).get("super_pathway")
            or canonical_pathway(row.get("predicted_pathway")),
            "sub_pathway": (known or {}).get("sub_pathway")
            or canonical_pathway(row.get("predicted_subpathway")),
            "pathway_source": "curated" if known else "predicted",
            "predicted_super_pathway": canonical_pathway(row.get("predicted_pathway")),
            "confidence_tier": _clean(row.get("pathway_confidence_tier")) or "(none)",
            "assignment_status": _clean(row.get("pathway_assignment_status")) or "(none)",
        }
    return assignments


def compose(
    baseline: list[dict[str, Any]],
    added: list[dict[str, Any]],
    level: str,
    weight_by_detection: bool,
) -> list[dict[str, Any]]:
    """Per-pathway baseline, addition and updated totals.

    ``baseline`` and ``added`` are already filtered to what counts as detected.
    Each record contributes 1, or its detection rate when weighting.

    Under ``weight_by_detection``, the primary reported figure is
    ``avg_metabolites_detected``: the sum of detection rates for knowns and
    pathway-annotated X-compounds in that pathway (five metabolites each seen
    20% of the time contribute 1.0), rounded to the nearest integer.
    Compositional ``updated_pct`` is omitted in that mode; by-compound counting
    still reports percentage shares.
    """

    key = "super_pathway" if level == "super" else "sub_pathway"

    def tally(records: list[dict[str, Any]]) -> dict[str, float]:
        totals: dict[str, float] = {}
        for record in records:
            weight = (record.get("detection_rate") or 0.0) if weight_by_detection else 1.0
            totals[record[key]] = totals.get(record[key], 0.0) + weight
        return totals

    if weight_by_detection:
        # A record with no rate contributes nothing when weighting, which would
        # quietly shrink a class rather than fail. Counted so it is visible.
        missing = sum(
            1 for record in baseline + added if record.get("detection_rate") is None
        )
        if missing:
            print(
                f"[feces-pathway] warning: {missing} records have no detection rate and "
                "contribute zero to the weighted view",
                file=sys.stderr,
            )

    base_totals, add_totals = tally(baseline), tally(added)
    base_sum, add_sum = sum(base_totals.values()), sum(add_totals.values())
    updated_sum = base_sum + add_sum
    predicted_vocabulary = {record[key] for record in added}

    rows = []
    for pathway in sorted(set(base_totals) | set(add_totals)):
        base = base_totals.get(pathway, 0.0)
        extra = add_totals.get(pathway, 0.0)
        combined = base + extra
        base_pct = (base / base_sum * 100) if base_sum else 0.0
        updated_pct = (combined / updated_sum * 100) if updated_sum else 0.0
        if weight_by_detection:
            known_avg = int(round(base))
            x_avg = int(round(extra))
            total_avg = int(round(combined))
            row = {
                "chemical_pathway": pathway,
                "identified_detected": known_avg,
                "baseline_pct": round(base_pct, 2),
                "x_compounds_added": x_avg,
                "updated_total": total_avg,
                # A class the pipeline never predicts can only lose share here,
                # so its movement says nothing about feces.
                "outside_predicted_vocabulary": str(
                    pathway not in predicted_vocabulary
                ).upper(),
                # Primary reported average: expected metabolites detected = Σ rates,
                # rounded to the nearest integer. No updated_pct: absolute averages
                # always rise or stay flat when X are added; compositional shares
                # do not, and were misleading here.
                "avg_metabolites_detected_known": known_avg,
                "avg_metabolites_detected_x": x_avg,
                "avg_metabolites_detected": total_avg,
            }
        else:
            row = {
                "chemical_pathway": pathway,
                "identified_detected": int(base),
                "baseline_pct": round(base_pct, 2),
                "x_compounds_added": int(extra),
                "updated_total": int(combined),
                "updated_pct": round(updated_pct, 2),
                "delta_pct_points": round(updated_pct - base_pct, 2),
                # A class the pipeline never predicts can only lose share here,
                # so its movement says nothing about feces.
                "outside_predicted_vocabulary": str(
                    pathway not in predicted_vocabulary
                ).upper(),
            }
        rows.append(row)
    sort_key = "avg_metabolites_detected" if weight_by_detection else "updated_total"
    return sorted(rows, key=lambda row: -float(row[sort_key]))


def build(
    panel: list[dict[str, Any]],
    detection: dict[str, float],
    assignments: dict[str, dict[str, str]],
    *,
    common_only: bool = False,
    tier: str | None = None,
    weight_by_detection: bool = False,
    restrict_to: set[str] | None = None,
) -> dict[str, Any]:
    if restrict_to is not None:
        panel = [row for row in panel if row["chemical_id"] in restrict_to]
    if common_only:
        baseline = [row for row in panel if row["gdp_bucket"] in COMMON_BUCKETS]
    else:
        baseline = [row for row in panel if row["gdp_bucket"] != NEVER]
    for row in baseline:
        row["detection_rate"] = row.get("gdp_rate")

    added = []
    for compound, assignment in assignments.items():
        if compound not in detection:
            continue
        # The same detection bar has to apply to both sides. Restricting the
        # named metabolites to the frequently seen ones while adding every
        # X-compound would inflate the shift by comparing unlike cohorts.
        if common_only and detection[compound] < COMMON_DETECTION_FLOOR:
            continue
        if tier and assignment["confidence_tier"] != tier:
            continue
        added.append(
            {
                "x_compound_id": compound,
                "super_pathway": assignment["super_pathway"],
                "sub_pathway": assignment["sub_pathway"],
                "detection_rate": detection[compound],
            }
        )

    return {
        "n_identified_detected": len(baseline),
        "n_x_compounds_added": len(added),
        "super_pathway": compose(baseline, added, "super", weight_by_detection),
        "sub_pathway": compose(baseline, added, "sub", weight_by_detection),
    }


def pathway_source_split(assignments: dict[str, dict[str, str]]) -> dict[str, Any]:
    """How many added pathways are known versus predicted, and how far they differ.

    The disagreement rate is the point: it is what the table would have been
    built on had the curated labels been ignored.
    """

    curated = [row for row in assignments.values() if row.get("pathway_source") == "curated"]
    disagree = [
        row for row in curated if row["super_pathway"] != row.get("predicted_super_pathway")
    ]
    return {
        "curated": len(curated),
        "predicted": len(assignments) - len(curated),
        "curated_disagreeing_with_prediction": len(disagree),
        "curated_agreement_rate": round(1 - len(disagree) / len(curated), 4) if curated else None,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_workbook(
    path: Path,
    pathways: list[dict[str, Any]],
    subpathways: list[dict[str, Any]],
) -> None:
    """Write pathway and sub-pathway results as tabs in one workbook."""

    import openpyxl
    from openpyxl.styles import Font
    from openpyxl.utils import get_column_letter

    path.parent.mkdir(parents=True, exist_ok=True)
    workbook = openpyxl.Workbook()
    workbook.remove(workbook.active)

    for title, rows in (("Pathway", pathways), ("Subpathway", subpathways)):
        sheet = workbook.create_sheet(title)
        fields = list(rows[0])
        sheet.append(fields)
        for row in rows:
            sheet.append([row[field] for field in fields])

        for cell in sheet[1]:
            cell.font = Font(bold=True)
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions
        for column_index, field in enumerate(fields, start=1):
            longest = max(
                len(str(field)),
                max((len(str(row[field])) for row in rows), default=0),
            )
            sheet.column_dimensions[get_column_letter(column_index)].width = min(
                longest + 2, 48
            )

    workbook.save(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--global-workbook", type=Path, required=True)
    parser.add_argument("--detection-history", type=Path, required=True)
    parser.add_argument("--reports-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("docs/paper/data/feces_pathway"))
    parser.add_argument(
        "--stratify", action="store_true", help="also report by pathway-confidence tier"
    )
    parser.add_argument(
        "--weight-by-detection",
        action="store_true",
        help=(
            "report expected metabolites detected per pathway (sum of detection "
            "rates for knowns and pathway-annotated X-compounds) rather than "
            "counting each compound once"
        ),
    )
    parser.add_argument(
        "--inclusion-list",
        type=Path,
        default=None,
        help="restrict the known side to the curated fecal panel in this workbook",
    )
    parser.add_argument(
        "--no-curated-pathways",
        action="store_true",
        help="predict every X-compound's pathway even where a curated one exists",
    )
    args = parser.parse_args()

    for path in (args.global_workbook, args.detection_history):
        if not path.is_file():
            raise SystemExit(f"missing input: {path}")
    if not args.reports_root.is_dir():
        raise SystemExit(f"reports root not found: {args.reports_root}")

    panel = read_global_panel(args.global_workbook)
    detection = read_detection_history(args.detection_history)
    # Global is the reference taxonomy, so it also decides which curated labels
    # are chemical classes rather than curation notes.
    vocabulary = {row["super_pathway"] for row in panel if row["super_pathway"] != "(unassigned)"}
    curated = {} if args.no_curated_pathways else read_curated_pathways(KNOWN_X_CSV, vocabulary)
    assignments = read_pathway_assignments(args.reports_root, detection, curated)

    restrict_to = None
    if args.inclusion_list is not None:
        if not args.inclusion_list.is_file():
            raise SystemExit(f"missing input: {args.inclusion_list}")
        restrict_to = read_inclusion_list(args.inclusion_list)

    primary = build(
        panel,
        detection,
        assignments,
        weight_by_detection=args.weight_by_detection,
        restrict_to=restrict_to,
    )
    summary: dict[str, Any] = {
        "method": (
            "Chemical-pathway composition of human feces, recomputed with X-compound "
            "pathway predictions added to the named metabolites the Global Discovery Panel "
            "detects. Detection is GDP, not MFP. Predicted pathways carry 52.2% measured "
            "super-pathway agreement against independent labels, so this shows direction "
            "and rough size rather than a measurement."
        ),
        "detection_basis": {
            "identified": f"{GDP_BUCKET_COLUMN} != {NEVER!r} in {GLOBAL_SHEET!r}",
            "known_side": "curated fecal inclusion list"
            if restrict_to is not None
            else "every metabolite the panel detects in feces",
            "x_compounds": (
                "present in the human-feces detection history and carrying a "
                "pathway annotation (curated-in-vocabulary or predicted)"
            ),
            "counting": (
                "expected metabolites detected per pathway: sum of detection "
                "rates for knowns plus pathway-annotated X-compounds, each "
                "average rounded to the nearest integer "
                "(no compositional updated_pct in this mode)"
                if args.weight_by_detection
                else "one compound, one count"
            ),
        },
        "x_compound_pathway_source": pathway_source_split(assignments),
        "coverage": {
            "global_panel_rows": len(panel),
            "known_side_size": len(restrict_to) if restrict_to is not None else len(panel),
            "identified_detected_in_feces": primary["n_identified_detected"],
            "x_compounds_in_detection_history": len(detection),
            "x_compounds_with_an_assignment": len(assignments),
            "x_compounds_added": primary["n_x_compounds_added"],
            "x_compounds_without_a_report": sorted(set(detection) - set(assignments)),
        },
        "super_pathway": primary["super_pathway"],
        "sub_pathway": primary["sub_pathway"],
    }

    if args.stratify:
        tiers = sorted({row["confidence_tier"] for row in assignments.values()})
        summary["by_confidence_tier"] = {
            tier: build(panel, detection, assignments, tier=tier)["super_pathway"] for tier in tiers
        }
        summary["common_detection_only"] = build(panel, detection, assignments, common_only=True)[
            "super_pathway"
        ]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    # The two counting modes answer different questions and are both worth
    # keeping, so they are named apart rather than overwriting each other.
    stem = "feces_pathway_update_weighted" if args.weight_by_detection else "feces_pathway_update"
    write_csv(args.out_dir / f"{stem}.csv", primary["super_pathway"])
    write_csv(args.out_dir / f"{stem}_subpathway.csv", primary["sub_pathway"])
    write_workbook(
        args.out_dir / f"{stem}.xlsx",
        primary["super_pathway"],
        primary["sub_pathway"],
    )
    (args.out_dir / f"{stem}.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(
        f"[feces-pathway] known detected in feces (GDP): "
        f"{primary['n_identified_detected']}  +  X-compounds added: "
        f"{primary['n_x_compounds_added']}"
    )
    split = summary["x_compound_pathway_source"]
    if split["curated"]:
        print(
            f"[feces-pathway] X pathways: {split['curated']} curated, {split['predicted']} "
            f"predicted; the curated ones agree with the prediction "
            f"{split['curated_agreement_rate']:.1%} of the time"
        )
    print()
    if args.weight_by_detection:
        header = (
            f"{'Chemical pathway':<40s} {'known':>8s} {'+X':>8s} {'avg det':>8s}"
        )
        print(header)
        print("-" * len(header))
        for row in primary["super_pathway"]:
            flag = " *" if row["outside_predicted_vocabulary"] == "TRUE" else ""
            print(
                f"{row['chemical_pathway']:<40s} "
                f"{row['avg_metabolites_detected_known']:>8d} "
                f"{row['avg_metabolites_detected_x']:>8d} "
                f"{row['avg_metabolites_detected']:>8d}{flag}"
            )
        print()
        print(
            "avg det = expected metabolites detected (sum of detection rates for "
            "knowns + pathway-annotated X-compounds)"
        )
    else:
        header = (
            f"{'Chemical pathway':<40s} {'ident':>6s} {'base%':>7s} {'+X':>5s} {'new%':>7s} "
            f"{'delta':>7s}"
        )
        print(header)
        print("-" * len(header))
        for row in primary["super_pathway"]:
            flag = " *" if row["outside_predicted_vocabulary"] == "TRUE" else ""
            print(
                f"{row['chemical_pathway']:<40s} {row['identified_detected']:>6} "
                f"{row['baseline_pct']:>6.2f}% {row['x_compounds_added']:>5} "
                f"{row['updated_pct']:>6.2f}% {row['delta_pct_points']:>+6.2f}{flag}"
            )
        print()
    print(
        "* the pipeline never predicts this class, so its share falls only because the "
        "denominator grew"
    )
    print(f"[feces-pathway] wrote {args.out_dir}/{stem}.csv and .xlsx")


if __name__ == "__main__":
    main()
