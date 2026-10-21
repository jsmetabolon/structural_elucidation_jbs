#!/usr/bin/env python3
"""Turn an executive-summary export directory into ``executive_summary_top10.xlsx``.

Read-only over the CSVs that ``tools/build_exec_summary.py`` writes. Five sheets:

- ``exec_summary_top10`` — the reviewed top-10 candidate sheet (one row per
  candidate) with the JSON provenance blobs left out; the readable
  ``excluded_library_reference`` column carries what they say about names.
- ``by_compound`` — one row per compound.
- ``naming_conventions`` — what each ``feature_name`` pattern means, so a
  reader of the sheet does not have to know the naming policy.
- ``summary`` — compound counts, match rates and how many features fall under
  each naming pattern.
- ``known_library_msi`` — every ``data/known_X.csv`` record with its Metabolon
  level translated to an MSI level and label. ``4?`` becomes MSI 4.

Text is read and written as UTF-8 throughout, so em dashes and primes survive
(the earlier hand-built workbook carried ``‚Äî`` where the CSV had ``—``).

Usage:
    python3 tools/build_exec_summary_workbook.py --exec-dir output/.../exec_summary
    python3 tools/build_exec_summary_workbook.py --exec-dir DIR --out path/to/workbook.xlsx
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd  # noqa: E402
from openpyxl.utils import get_column_letter  # noqa: E402

_TOOLS = Path(__file__).resolve().parent
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))
from reported_identification import (  # noqa: E402
    CALIBRATED_TIER_ORDER,
    LEVEL_2_STAR_FLOOR,
    _CONFIDENCE_BANDS,
    active_calibration,
    known_library_msi_rows,
)

KNOWN_X_CSV = PROJECT_ROOT / "data" / "known_X.csv"

csv.field_size_limit(10**9)

WORKBOOK_NAME = "executive_summary_top10.xlsx"

# Machine-readable provenance kept in the CSVs but too noisy for a spreadsheet.
JSON_COLUMNS = {"excluded_library_matches", "raw_rank1_candidate", "identity_candidate"}

NUMERIC_COLUMNS = {
    "msi_level", "rank", "identity_candidate_rank", "candidate_confidence",
    "rank_score_relative", "evidence_confidence", "mass_match_score",
    "abs_mass_error_da", "abs_mass_error_ppm", "predicted_spectrum_score_v2", "truth_tanimoto", "matching_ranks",
    "candidate_identity_confidence", "pure_nce_score", "nce_score", "fallback_score",
    "n_candidates", "first_matching_rank", "best_truth_tanimoto",
    "rank1_candidate_confidence", "rank1_rank_score_relative", "rank1_evidence_confidence",
    "rank1_mass_match_score", "rank1_abs_mass_error_da", "rank1_abs_mass_error_ppm",
}

NAMING_CONVENTIONS: list[dict[str, str]] = [
    {
        "feature_name pattern": "<name>",
        "msi_level": "1",
        "meaning": "Identity confirmed against an authentic standard; the bare name is asserted.",
        "name_source": "authentic_standard",
    },
    {
        "feature_name pattern": "Probable <name>",
        "msi_level": "2",
        "meaning": "Probable structure from the highest-ranked identity-eligible candidate.",
        "name_source": "ranked_candidate",
    },
    {
        "feature_name pattern": "Possible <name>",
        "msi_level": "3",
        "meaning": (
            "Compound-specific tentative annotation carried on the X-compound's own known_X row "
            "(curated hypothesis); supports at most Level 3."
        ),
        "name_source": "known_x_tentative",
    },
    {
        "feature_name pattern": "Isomer of <library metabolite>",
        "msi_level": "4 or 5",
        "meaning": (
            "The rank-1 predicted structure is already a catalogued global-library metabolite. "
            "An X-compound cannot be that metabolite, so the feature is reported as an isomer of it: "
            "same mass/formula, different structure. The report itself shows the formula or exact-mass "
            "placeholder (see report_display_name); every excluded library metabolite is listed in "
            "excluded_library_reference."
        ),
        "name_source": "ranked_candidate / formula_or_mass",
    },
    {
        "feature_name pattern": "<formula>; structure unresolved",
        "msi_level": "4",
        "meaning": "Molecular formula supported; no structure is asserted at this level.",
        "name_source": "ranked_candidate / formula_or_mass",
    },
    {
        "feature_name pattern": "Exact mass <m>; unidentified",
        "msi_level": "5",
        "meaning": "Only the exact mass is reported; no formula or structure is asserted.",
        "name_source": "formula_or_mass",
    },
    {
        "feature_name pattern": "Tentative Structure",
        "msi_level": "3",
        "meaning": "A de novo structure with no library name (formerly 'De novo structure (<formula>)').",
        "name_source": "ranked_candidate",
    },
    {
        "feature_name pattern": "<class> conjugate of <formula> / Putative analog of <name> / class descriptors",
        "msi_level": "3",
        "meaning": "Descriptor synthesized from the structure; not a library metabolite's name and never prefixed.",
        "name_source": "known_x_tentative / ranked_candidate",
    },
    {
        "feature_name pattern": "[column] known_x_tentative_information / known_x_cleaned_external_info",
        "msi_level": "",
        "meaning": (
            "Curator annotation text from known_X.csv, verbatim: the Level 2-3 'Tentative information' "
            "and the 'Cleaned External Info' name. Joined with '; ' when a compound is in several libraries."
        ),
        "name_source": "",
    },
    {
        "feature_name pattern": "[column] reported_level",
        "msi_level": "",
        "meaning": (
            "White-paper level alone: 1, 2*, 2, 3, 4, or 5. The publication name for that "
            "level is in reported_identification."
        ),
        "name_source": "reported_level",
    },
    {
        "feature_name pattern": "[column] reported_identification",
        "msi_level": "",
        "meaning": (
            "White-paper publication call, the same value on every candidate row of a compound. "
            "Level 1 is an authentic-standard name. Level 2* is a named-library entry that is not "
            "one structure, such as a glucuronide of a formula. Level 2 is "
            "'[X-######], Probable {name}' for a specific structure name (a name is enough; "
            "the candidates do not have to agree). Level 3 is 'isomer of {library name}' for an "
            "exact match to a library entry besides retention time, 'Probable isomer of {top "
            "candidate}' for a related family of candidates that clears unrelated candidates by a "
            "distinctive Match Score, a class phrase such as 'Glucoside', or a conjugate, sugar, "
            "or peptide shared by most of the candidates. "
            "Level 4 is the molecular formula when the compound has no name and no shared class. "
            "Level 5 is 'exact mass of interest' when no formula is supported."
        ),
        "name_source": "reported_identification",
    },
    {
        "feature_name pattern": "[column] feature_name_on_prod_library_panel / prod_library_panel_matches",
        "msi_level": "",
        "meaning": (
            "TRUE/FALSE: whether the library metabolite the feature is named after (or is an isomer of) is on "
            "the production global-library panel (library_global_discovery_panel_v5, 5,413 entries); blank for "
            "formula / mass / tentative placeholders. prod_library_panel_matches lists every panel metabolite "
            "the report references as 'name (CHEMICAL_ID)'."
        ),
        "name_source": "",
    },
    {
        "feature_name pattern": "(retired) De novo structure (<formula>)",
        "msi_level": "",
        "meaning": "Replaced by 'Tentative Structure'.",
        "name_source": "",
    },
    {
        "feature_name pattern": "(retired) <name> (proposed structure)",
        "msi_level": "",
        "meaning": "Level 3 hedge replaced by 'Possible <name>' or 'Isomer of <name>'.",
        "name_source": "",
    },
]

TIER_LABELS = {
    "2_probable": "Level 2 Probable name",
    "3_isomer_exact": "Level 3 isomer of",
    "3_class": "Level 3 class/motif",
    "3_probable_isomer_group": "Level 3 Probable isomer of",
    "4": "Level 4 formula",
    "5": "Level 5 exact mass",
}


def identity_confidence_conventions() -> list[dict[str, str]]:
    """Glossary rows for candidate_identity_confidence and pure_nce_score, from the active calibration."""

    calibration = active_calibration()
    if calibration is None:
        tiers = "; ".join(
            f"{TIER_LABELS[tier]} {lo:g}-{hi:g}" for tier, (lo, hi) in _CONFIDENCE_BANDS.items()
        )
        meaning = (
            "0-100 confidence that the reported identification, not a rival, is right. Level 1 is 100; "
            f"Level 2* is at least {LEVEL_2_STAR_FLOOR:g}. Other levels use hand-set, uncalibrated bands "
            f"({tiers}); the candidate margin places a compound inside its band."
        )
    else:
        parts = []
        for tier in CALIBRATED_TIER_ORDER:
            spec = calibration["tiers"][tier]
            if not spec.get("calibrated"):
                status = "uncalibrated"
            elif spec["n"] < 10:
                status = f"n={spec['n']}, prior-dominated"
            else:
                status = f"n={spec['n']}"
            parts.append(f"{TIER_LABELS[tier]} {spec['anchor']:g} [{spec['floor']:g}-{spec['ceiling']:g}] ({status})")
        meaning = (
            "0-100 estimate of the chance the reported structure, not a rival, is right. Level 1 is 100; "
            f"Level 2* is at least {LEVEL_2_STAR_FLOOR:g} (curator review). Other levels are calibrated on "
            f"{calibration.get('n_labelled')} truth controls replayed with their names hidden: for each, is rank 1 "
            "the true structure, by the same Tanimoto test behind matched_truth elsewhere. Rates are shrunk toward "
            "the earlier hand-set bands and kept in level order: "
            + "; ".join(parts)
            + ". 'Uncalibrated' tiers had no replayable labels: they start from the hand-set band midpoint and "
            "move only as far as needed to keep level order. A formula_only control has no reference structure, so "
            "it can only confirm the formula, a weaker label flagged in the calibration audit. The pure NCE "
            "margin of rank 1 over unrelated rivals places a compound inside its band."
        )
    return [
        {
            "feature_name pattern": "[column] candidate_identity_confidence",
            "msi_level": "",
            "meaning": meaning,
            "name_source": "candidate_identity_confidence",
        },
        {
            "feature_name pattern": "[column] pure_nce_score",
            "msi_level": "",
            "meaning": (
                "NCE predicted-spectrum match for this candidate with no fallback score (0 when the matcher finds "
                "no overlap; blank when the compound has no observed MS2 bins). Unlike evidence_confidence it does "
                "not saturate near 1, so it separates rival candidates."
            ),
            "name_source": "pure_nce_score",
        },
    ]


NAME_PATTERNS: list[tuple[str, Any]] = [
    ("Isomer of <library metabolite>", lambda s: s.lower().startswith("isomer of ")),
    ("Possible <name>", lambda s: s.lower().startswith("possible ")),
    ("Probable <name>", lambda s: s.lower().startswith("probable ")),
    ("<formula>; structure unresolved", lambda s: "structure unresolved" in s.lower()),
    ("Exact mass <m>; unidentified", lambda s: s.lower().startswith("exact mass ")),
    ("Tentative Structure", lambda s: s.lower().startswith("tentative structure")),
    ("(retired) De novo structure", lambda s: s.lower().startswith("de novo structure")),
    ("(retired) (proposed structure)", lambda s: "(proposed structure)" in s.lower()),
]


def name_pattern(label: str) -> str:
    text = str(label or "").strip()
    if not text:
        return "<blank>"
    for pattern, test in NAME_PATTERNS:
        if test(text):
            return pattern
    return "<name> / descriptor"


def read_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8")
    frame = frame[[column for column in frame.columns if column not in JSON_COLUMNS]]
    for column in frame.columns:
        if column in NUMERIC_COLUMNS:
            frame[column] = pd.to_numeric(frame[column].replace("", pd.NA), errors="coerce")
    return frame


def summary_rows(exec_dir: Path, by_compound: pd.DataFrame, top10: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    summary_path = exec_dir / "exec_summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        rows.append({"section": "scope", "item": "compounds", "value": summary.get("n_compounds")})
        rows.append({"section": "scope", "item": "candidate rows", "value": len(top10)})
        controls = summary.get("truth_control") or {}
        rows.append({"section": "truth controls", "item": "n", "value": controls.get("n")})
        threshold = controls.get("truth_threshold")
        for key in ("rank_1", "top_3", "top_10"):
            rate = controls.get(key) or {}
            rows.append({
                "section": "truth controls",
                "item": f"structural match within {key.replace('_', ' ')} (Tanimoto >= {threshold})",
                "value": f"{rate.get('n_matched')}/{rate.get('n_assessable')} ({rate.get('rate')})",
            })
        two_c = summary.get("level_2c") or {}
        rows.append({"section": "level 2C", "item": "n", "value": two_c.get("n")})
        rows.append({
            "section": "level 2C", "item": "with a reference formula",
            "value": two_c.get("n_with_a_reference_formula"),
        })
        rows.append({"section": "not assessable", "item": "n", "value": (summary.get("other") or {}).get("n")})
    policies = Counter(by_compound.get("naming_policy_version", pd.Series(dtype=str)).fillna(""))
    for policy, count in sorted(policies.items()):
        rows.append({"section": "naming policy", "item": policy or "<unset>", "value": count})
    if "msi_level" in by_compound:
        for level, count in sorted(Counter(by_compound["msi_level"].fillna(-1)).items()):
            item = "" if level == -1 else int(level)
            rows.append({"section": "compounds by MSI level", "item": item, "value": count})
    for pattern, count in Counter(name_pattern(v) for v in by_compound["feature_name"]).most_common():
        rows.append({"section": "compounds by feature_name pattern", "item": pattern, "value": count})
    if "known_x_tentative_information" in by_compound:
        rows.append({
            "section": "known_X annotations", "item": "compounds with Tentative information",
            "value": int((by_compound["known_x_tentative_information"] != "").sum()),
        })
        rows.append({
            "section": "known_X annotations", "item": "compounds with Cleaned External Info",
            "value": int((by_compound["known_x_cleaned_external_info"] != "").sum()),
        })
    if "feature_name_on_prod_library_panel" in by_compound:
        flags = Counter(by_compound["feature_name_on_prod_library_panel"])
        for flag, label in (("TRUE", "feature named after a panel metabolite"),
                            ("FALSE", "feature named after a metabolite NOT on the panel"),
                            ("", "placeholder label (no metabolite to check)")):
            rows.append({"section": "production library panel", "item": label, "value": flags.get(flag, 0)})
        rows.append({
            "section": "production library panel", "item": "compounds referencing >=1 panel metabolite anywhere",
            "value": int((by_compound["prod_library_panel_matches"] != "").sum()),
        })
    return rows


def autosize(sheet, frame: pd.DataFrame, *, max_width: int = 60) -> None:
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions
    for index, column in enumerate(frame.columns, start=1):
        sample = frame[column].head(500).map(lambda value: len(str(value)))
        width = max([len(str(column))] + sample.tolist()) if len(sample) else len(str(column))
        sheet.column_dimensions[get_column_letter(index)].width = min(max(10, width + 2), max_width)


def build_workbook(exec_dir: Path, out: Path) -> dict[str, int]:
    top10_path = exec_dir / "exec_summary_top10_truncated.csv"
    if not top10_path.exists():
        top10_path = exec_dir / "exec_summary_top10.csv"
    by_compound_path = exec_dir / "exec_summary_by_compound.csv"
    for path in (top10_path, by_compound_path):
        if not path.exists():
            raise SystemExit(f"missing {path}")

    top10 = read_table(top10_path)
    by_compound = read_table(by_compound_path)
    conventions = pd.DataFrame(NAMING_CONVENTIONS + identity_confidence_conventions())
    summary = pd.DataFrame(summary_rows(exec_dir, by_compound, top10))

    sheets = [
        ("exec_summary_top10", top10),
        ("by_compound", by_compound),
        ("naming_conventions", conventions),
        ("summary", summary),
    ]
    if KNOWN_X_CSV.exists():
        sheets.append(("known_library_msi", pd.DataFrame(known_library_msi_rows(KNOWN_X_CSV))))

    out.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        for name, frame in sheets:
            frame.to_excel(writer, sheet_name=name, index=False)
            autosize(writer.sheets[name], frame)
    return {"candidate_rows": len(top10), "compounds": len(by_compound)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--exec-dir", type=Path, required=True, help="directory holding the exec_summary CSVs")
    parser.add_argument("--out", type=Path, default=None, help=f"workbook path (default: <exec-dir>/{WORKBOOK_NAME})")
    args = parser.parse_args()
    if not args.exec_dir.is_dir():
        raise SystemExit(f"exec dir not found: {args.exec_dir}")
    out = args.out or args.exec_dir / WORKBOOK_NAME
    counts = build_workbook(args.exec_dir, out)
    print(f"[exec-summary-workbook] {counts['compounds']} compounds, {counts['candidate_rows']} candidate rows")
    print(f"[exec-summary-workbook] wrote {out}")


if __name__ == "__main__":
    main()
