#!/usr/bin/env python3
"""Executive summary of the top ten candidates proposed for every X-compound.

Read-only against a delivered report bundle. Pure standard library.

Emits two CSVs:

- ``exec_summary_top10.csv`` — one row per candidate, ten per compound at most.
- ``exec_summary_by_compound.csv`` — one row per compound, rolling the ten up
  into whether and where the truth was recovered.

The point of interest is the ``matched_truth`` column, and it is only honestly
answerable for some compounds. Three cohorts, three different situations:

``truth_control``
    92 compounds carrying a blinded reference structure. Every candidate has a
    precomputed ``truth_tanimoto`` against it, so a match is a structural claim
    at the published threshold. These are almost all Metabolon Level 1.

``level_2c``
    163 delivered compounds annotated at Metabolon Level 2C. There is no
    structure to compare against: the annotation is a tentative *name*, no 2C
    compound has a single named candidate in its top ten, and none carries a
    ``truth_tanimoto``. Where the annotation records a molecular formula, and
    only 26 of the 163 do, formula agreement is reported instead. That is a
    weaker claim than a structural match against a weaker reference, and the
    ``truth_basis`` column says so on every row.

``other``
    No reference of either kind. ``matched_truth`` is left blank rather than
    false, because unassessable and wrong are different things.

Usage:
    python3 tools/build_exec_summary.py --reports-root /path/to/reports
    python3 tools/build_exec_summary.py --reports-root /path/to/reports --out-dir out/
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Any, Iterator

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

csv.field_size_limit(10**9)

# The published definition of a correct structure, matching audit_mass_policy.py
# and audit_structure_predictions.py.
TRUTH_THRESHOLD = 0.70

KNOWN_X_CSV = PROJECT_ROOT / "data" / "known_X.csv"
STRUCTURES_CSV = PROJECT_ROOT / "data" / "known_X_all_available_structures_simple.csv"

TOP_N = 10

COHORT_TRUTH = "truth_control"
COHORT_2C = "level_2c"
COHORT_OTHER = "other"

BASIS_STRUCTURE = "structure_tanimoto"
BASIS_FORMULA = "formula_only"
BASIS_NONE = "none"

LONG_COLUMNS = [
    "x_compound_id",
    "cohort",
    "metabolon_level",
    "metabolon_sublevel",
    "msi_level",
    "msi_label",
    "feature_name",
    "reference_answer",
    "reference_formula",
    "rank",
    "candidate_smiles",
    "candidate_formula",
    "candidate_confidence",
    "rank_score_relative",
    "evidence_confidence",
    "mass_match_score",
    "abs_mass_error_da",
    "truth_basis",
    "truth_tanimoto",
    "matched_truth",
    "is_best_match_in_top10",
]

ROLLUP_COLUMNS = [
    "x_compound_id",
    "cohort",
    "metabolon_level",
    "metabolon_sublevel",
    "msi_level",
    "msi_label",
    "feature_name",
    "reference_answer",
    "truth_basis",
    "n_candidates",
    "matched_at_rank_1",
    "matched_within_top_3",
    "matched_within_top_10",
    "first_matching_rank",
    "best_truth_tanimoto",
    "rank1_candidate_smiles",
    "rank1_candidate_formula",
    "rank1_candidate_confidence",
    "rank1_rank_score_relative",
    "rank1_evidence_confidence",
    "rank1_mass_match_score",
]


def _num(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if number != number else number


def _clean(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    return "" if text.lower() in {"nan", "none"} else text


def _round(value: float | None, places: int) -> str:
    return "" if value is None else f"{round(value, places)}"


FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def normalise_formula(formula: str) -> str:
    """Canonicalise a molecular formula so element order stops mattering.

    The annotations and the candidates do not agree on ordering: a reference may
    read ``C6H8SO7`` where the pipeline writes ``C6H8O7S``. Comparing the strings
    would call an exact match a mismatch, so the formula is parsed into element
    counts and re-emitted in a fixed order. Anything unparseable falls back to
    the stripped string rather than silently matching.
    """

    text = "".join(formula.split())
    if not text:
        return ""
    counts: dict[str, int] = {}
    position = 0
    for token in FORMULA_TOKEN.finditer(text):
        if token.start() != position:
            return text.upper()
        position = token.end()
        element, digits = token.groups()
        counts[element] = counts.get(element, 0) + (int(digits) if digits else 1)
    if position != len(text):
        return text.upper()
    return "".join(f"{element}{counts[element]}" for element in sorted(counts))


def load_metabolon_levels() -> dict[str, dict[str, Any]]:
    """Level, sublevel, tentative names and formulae, keyed by X-compound.

    A compound can appear several times in ``known_X.csv``, once per library, so
    values are collected as sets and reported joined.
    """

    if not KNOWN_X_CSV.exists():
        raise SystemExit(f"missing {KNOWN_X_CSV}")
    table: dict[str, dict[str, Any]] = {}
    with KNOWN_X_CSV.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            compound = _clean(row.get("Unknown"))
            if not compound:
                continue
            entry = table.setdefault(
                compound, {"levels": set(), "sublevels": set(), "names": set(), "formulae": set()}
            )
            if level := _clean(row.get("Level")):
                entry["levels"].add(level)
            if sublevel := _clean(row.get("SubLevel")):
                entry["sublevels"].add(sublevel)
            if name := _clean(row.get("Tentative information")):
                entry["names"].add(name)
            if formula := _clean(row.get("Formula")):
                entry["formulae"].add(formula)
    return table


def load_reference_structures() -> dict[str, str]:
    """Known names for the blinded controls, keyed by X-compound."""

    if not STRUCTURES_CSV.exists():
        return {}
    with STRUCTURES_CSV.open(newline="", encoding="utf-8-sig") as handle:
        return {
            _clean(row.get("x_compound_id")): _clean(row.get("known_name"))
            for row in csv.DictReader(handle)
            if _clean(row.get("x_compound_id"))
        }


def read_report(directory: Path) -> dict[str, Any] | None:
    """One compound's top ten, joined across the two files that describe it.

    ``candidate_evidence`` in the JSON and the rows of the candidates CSV are
    both ordered by rank and always the same length, so they join on rank.
    """

    stem = directory.name.lower()
    candidates_csv = directory / "data" / f"{stem}_top_candidates.csv"
    inputs_json = directory / "data" / f"{stem}_confidence_inputs.json"
    if not candidates_csv.exists() or not inputs_json.exists():
        return None

    fields = json.loads(inputs_json.read_text())
    evidence = {
        int(item["rank"]): item
        for item in (fields.get("candidate_evidence") or [])
        if item.get("rank") is not None
    }

    with candidates_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))[:TOP_N]

    candidates = []
    for row in rows:
        rank_value = _num(row.get("rank"))
        rank = int(rank_value) if rank_value is not None else len(candidates) + 1
        extra = evidence.get(rank, {})
        candidates.append(
            {
                "rank": rank,
                "smiles": _clean(row.get("predicted_smiles")) or _clean(row.get("smiles")),
                "formula": _clean(row.get("predicted_formula"))
                or _clean(row.get("candidate_formula")),
                "candidate_confidence": _num(extra.get("candidate_confidence")),
                "rank_score_relative": _num(row.get("candidate_report_score")),
                "evidence_confidence": _num(row.get("evidence_confidence")),
                "mass_match_score": _num(row.get("mass_match_score")),
                "abs_mass_error_da": _num(row.get("abs_mass_error_da")),
                "truth_tanimoto": _num(row.get("truth_tanimoto")),
            }
        )

    return {
        "x_compound_id": directory.name,
        "display_name": _clean(fields.get("display_name")),
        "msi_level": fields.get("msi_level"),
        "msi_label": _clean(fields.get("msi_label")),
        "candidates": candidates,
    }


def read_reports(root: Path) -> Iterator[dict[str, Any]]:
    for directory in sorted(root.iterdir()):
        if not directory.is_dir() or not directory.name.upper().startswith("X-"):
            continue
        if (report := read_report(directory)) is not None:
            yield report


def classify(
    report: dict[str, Any], levels: dict[str, dict[str, Any]], structures: dict[str, str]
) -> dict[str, Any]:
    """Assign the cohort and pick the strongest reference available for it.

    Structural truth wins wherever it exists, because it is the only basis that
    supports a claim about the structure. Level 2C falls back to the tentative
    annotation, and everything else gets no basis at all.
    """

    compound = report["x_compound_id"]
    entry = levels.get(compound, {})
    sublevels = entry.get("sublevels") or set()
    has_structure = any(
        candidate["truth_tanimoto"] is not None and candidate["truth_tanimoto"] > 0
        for candidate in report["candidates"]
    )

    if has_structure:
        return {
            "cohort": COHORT_TRUTH,
            "truth_basis": BASIS_STRUCTURE,
            "reference_answer": structures.get(compound, ""),
            "reference_formula": "",
        }
    if "2C" in sublevels:
        formulae = sorted(entry.get("formulae") or set())
        return {
            "cohort": COHORT_2C,
            "truth_basis": BASIS_FORMULA if formulae else BASIS_NONE,
            "reference_answer": "; ".join(sorted(entry.get("names") or set())),
            "reference_formula": "; ".join(formulae),
        }
    return {
        "cohort": COHORT_OTHER,
        "truth_basis": BASIS_NONE,
        "reference_answer": "",
        "reference_formula": "",
    }


def match_flag(candidate: dict[str, Any], context: dict[str, Any]) -> bool | None:
    """Whether one candidate matches the reference, on whatever basis exists.

    Returns ``None`` when the question cannot be asked, which the CSV writes as
    an empty cell so it never reads as a failed match.
    """

    basis = context["truth_basis"]
    if basis == BASIS_STRUCTURE:
        tanimoto = candidate["truth_tanimoto"]
        return tanimoto is not None and tanimoto >= TRUTH_THRESHOLD
    if basis == BASIS_FORMULA:
        reference = {
            normalise_formula(part) for part in context["reference_formula"].split(";") if part
        }
        return bool(reference) and normalise_formula(candidate["formula"]) in reference
    return None


def build_rows(
    report: dict[str, Any], context: dict[str, Any], entry: dict[str, Any]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    matches = [match_flag(candidate, context) for candidate in report["candidates"]]
    first_match = next(
        (candidate["rank"] for candidate, hit in zip(report["candidates"], matches) if hit), None
    )
    tanimotos = [
        candidate["truth_tanimoto"]
        for candidate in report["candidates"]
        if candidate["truth_tanimoto"] is not None
    ]
    shared = {
        "x_compound_id": report["x_compound_id"],
        "cohort": context["cohort"],
        "metabolon_level": "; ".join(sorted(entry.get("levels") or set())),
        "metabolon_sublevel": "; ".join(sorted(entry.get("sublevels") or set())),
        "msi_level": "" if report["msi_level"] is None else report["msi_level"],
        "msi_label": report["msi_label"],
        "feature_name": report["display_name"],
        "reference_answer": context["reference_answer"],
    }

    rows = []
    for candidate, hit in zip(report["candidates"], matches):
        rows.append(
            {
                **shared,
                "reference_formula": context["reference_formula"],
                "rank": candidate["rank"],
                "candidate_smiles": candidate["smiles"],
                "candidate_formula": candidate["formula"],
                "candidate_confidence": _round(candidate["candidate_confidence"], 4),
                "rank_score_relative": _round(candidate["rank_score_relative"], 6),
                "evidence_confidence": _round(candidate["evidence_confidence"], 6),
                "mass_match_score": _round(candidate["mass_match_score"], 6),
                "abs_mass_error_da": _round(candidate["abs_mass_error_da"], 6),
                "truth_basis": context["truth_basis"],
                "truth_tanimoto": _round(candidate["truth_tanimoto"], 4),
                "matched_truth": "" if hit is None else str(hit).upper(),
                "is_best_match_in_top10": str(
                    first_match is not None and candidate["rank"] == first_match
                ).upper(),
            }
        )

    top = report["candidates"][0] if report["candidates"] else {}
    assessable = context["truth_basis"] != BASIS_NONE
    rollup = {
        **shared,
        "truth_basis": context["truth_basis"],
        "n_candidates": len(report["candidates"]),
        "matched_at_rank_1": str(first_match == 1).upper() if assessable else "",
        "matched_within_top_3": str(
            first_match is not None and first_match <= 3
        ).upper()
        if assessable
        else "",
        "matched_within_top_10": str(first_match is not None).upper() if assessable else "",
        "first_matching_rank": "" if first_match is None else first_match,
        "best_truth_tanimoto": _round(max(tanimotos), 4) if tanimotos else "",
        "rank1_candidate_smiles": top.get("smiles", ""),
        "rank1_candidate_formula": top.get("formula", ""),
        "rank1_candidate_confidence": _round(top.get("candidate_confidence"), 4),
        "rank1_rank_score_relative": _round(top.get("rank_score_relative"), 6),
        "rank1_evidence_confidence": _round(top.get("evidence_confidence"), 6),
        "rank1_mass_match_score": _round(top.get("mass_match_score"), 6),
    }
    return rows, rollup


def summarise(rollups: list[dict[str, Any]]) -> dict[str, Any]:
    def cohort(name: str) -> list[dict[str, Any]]:
        return [row for row in rollups if row["cohort"] == name]

    def rate(rows: list[dict[str, Any]], column: str) -> dict[str, Any]:
        assessable = [row for row in rows if row[column] != ""]
        hits = sum(1 for row in assessable if row[column] == "TRUE")
        return {
            "n_assessable": len(assessable),
            "n_matched": hits,
            "rate": round(hits / len(assessable), 4) if assessable else None,
        }

    controls = cohort(COHORT_TRUTH)
    two_c = cohort(COHORT_2C)
    return {
        "n_compounds": len(rollups),
        "truth_control": {
            "n": len(controls),
            "basis": BASIS_STRUCTURE,
            "truth_threshold": TRUTH_THRESHOLD,
            "rank_1": rate(controls, "matched_at_rank_1"),
            "top_3": rate(controls, "matched_within_top_3"),
            "top_10": rate(controls, "matched_within_top_10"),
        },
        "level_2c": {
            "n": len(two_c),
            "basis": BASIS_FORMULA,
            "n_with_a_reference_formula": sum(
                1 for row in two_c if row["truth_basis"] == BASIS_FORMULA
            ),
            "rank_1": rate(two_c, "matched_at_rank_1"),
            "top_3": rate(two_c, "matched_within_top_3"),
            "top_10": rate(two_c, "matched_within_top_10"),
            "caveat": (
                "Formula agreement against a tentative annotation, not a structural match. "
                "No Level 2C compound has a reference structure or a named candidate, and "
                "only a minority record a formula at all; the rest are not assessable."
            ),
        },
        "other": {"n": len(cohort(COHORT_OTHER)), "basis": BASIS_NONE},
    }


def write_csv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports-root", type=Path, required=True)
    parser.add_argument(
        "--out-dir", type=Path, default=Path("docs/paper/data/exec_summary"), help="write CSVs here"
    )
    args = parser.parse_args()

    if not args.reports_root.is_dir():
        raise SystemExit(f"reports root not found: {args.reports_root}")

    levels = load_metabolon_levels()
    structures = load_reference_structures()

    long_rows: list[dict[str, Any]] = []
    rollups: list[dict[str, Any]] = []
    for report in read_reports(args.reports_root):
        context = classify(report, levels, structures)
        entry = levels.get(report["x_compound_id"], {})
        rows, rollup = build_rows(report, context, entry)
        long_rows.extend(rows)
        rollups.append(rollup)

    if not rollups:
        raise SystemExit(f"no X-* report directories under {args.reports_root}")

    summary = summarise(rollups)
    write_csv(args.out_dir / "exec_summary_top10.csv", LONG_COLUMNS, long_rows)
    write_csv(args.out_dir / "exec_summary_by_compound.csv", ROLLUP_COLUMNS, rollups)
    (args.out_dir / "exec_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(f"[exec-summary] {summary['n_compounds']} compounds, {len(long_rows)} candidate rows")
    controls = summary["truth_control"]
    print(
        f"[exec-summary] truth controls: {controls['n']} | structural match at rank 1 "
        f"{controls['rank_1']['n_matched']}/{controls['rank_1']['n_assessable']}"
        f" ({controls['rank_1']['rate']:.1%}), top 3 {controls['top_3']['rate']:.1%},"
        f" top 10 {controls['top_10']['rate']:.1%}"
    )
    two_c = summary["level_2c"]
    print(
        f"[exec-summary] level 2C: {two_c['n']} compounds, only "
        f"{two_c['n_with_a_reference_formula']} carry a reference formula; of those, formula "
        f"agreement in the top 10 is {two_c['top_10']['n_matched']}/"
        f"{two_c['top_10']['n_assessable']}"
    )
    print(f"[exec-summary] not assessable: {summary['other']['n']} compounds")
    print(f"[exec-summary] wrote {args.out_dir}/exec_summary_top10.csv")


if __name__ == "__main__":
    main()
