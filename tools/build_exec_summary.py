#!/usr/bin/env python3
"""Executive summary of the top ten candidates proposed for every X-compound.

Read-only against a delivered report bundle. The displayed ``evidence_confidence``
is filled through the pipeline's own canonical definition
(``structural_elucidation/chemistry/evidence.py``) so an uncomputed score renders
blank rather than as a misleading ``0.0``, and a populated one matches the number
the reranker consumes as a feature.

Emits three CSVs:

- ``exec_summary_top10.csv`` — one row per candidate, ten per compound at most.
  ``matching_ranks`` is the highest-ranked match only (e.g. ``2``, never ``2;4``);
  ``matched_truth`` is still TRUE for every candidate at or above the threshold;
  ``is_best_match_in_top10`` marks that same first rank.
- ``exec_summary_top10_truncated.csv`` — the same rows restricted to the
  columns kept in the reviewed ``executive_summary_top10.xlsx`` top-ten sheet
  (identity, the best Metabolon level, the known_X annotation text, the
  displayed name, and the candidate's rank, structure, evidence, mass error,
  and match flags).
- ``exec_summary_by_compound.csv`` — one row per compound, rolling the ten up
  into whether and where the truth was recovered, including the same
  ``matching_ranks`` value.

``feature_name`` follows the report's naming policy. For Naming v2 records a
compound whose rank-1 structure is an existing global-library metabolite, and
which the report therefore left at a formula / exact-mass placeholder, is
written ``Isomer of <library metabolite>``; the report's own label is kept in
``report_display_name`` and every excluded library metabolite is listed in
``excluded_library_reference``. See :func:`isomer_of_feature_name`.

``reported_level`` is the white-paper level on its own (``1``, ``2*``, ``2``,
``3``, ``4``, ``5``). ``reported_identification`` is the publication name that
goes with it. See :mod:`reported_identification`.

Two further groups of columns come from outside the report bundle:

- ``known_x_tentative_information`` / ``known_x_cleaned_external_info`` carry
  the curator's Level 2-3 annotation text from ``data/known_X.csv`` verbatim
  (joined with ``; `` when the compound sits in several libraries).
- ``feature_name_on_prod_library_panel`` / ``prod_library_panel_matches`` mark
  whether the library metabolite a feature is named after, or collided with,
  is on the production global-library panel (``--prod-library``, default
  ``data/prod_library/library_global_discovery_panel_v5.xlsx``, 5,413 entries).

``tools/build_exec_summary_workbook.py`` turns these CSVs into the
``executive_summary_top10.xlsx`` workbook.

The point of interest is the ``matched_truth`` column (per candidate) and
``matching_ranks`` (the highest-ranked hit). Both are only honestly answerable
for some compounds. Three cohorts, three different situations:

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

import pandas as pd  # noqa: E402  (imported after the sys.path bootstrap above)

from structural_elucidation.chemistry.evidence import (  # noqa: E402
    BASIS_BACKFILLED,
    BASIS_EMITTED,
    BASIS_UNAVAILABLE,
    backfill_evidence_confidence,
)
from structural_elucidation.reports.tentative_naming import (  # noqa: E402
    FORMULA_ONLY_NOTE,
    UNNAMED_STRUCTURE_LABEL,
    isomer_of_label,
    read_known_x_name,
    sheet_feature_name,
)

_TOOLS = Path(__file__).resolve().parent
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))
from reported_identification import (  # noqa: E402
    candidate_identity_confidence,
    msi_from_library_level,
    reported_identification,
    reported_level,
)

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
    "known_x_tentative_information",
    "known_x_cleaned_external_info",
    "msi_level",
    "msi_label",
    "feature_name",
    "report_display_name",
    "reported_level",
    "reported_identification",
    "candidate_identity_confidence",
    "excluded_library_reference",
    "feature_name_on_prod_library_panel",
    "prod_library_panel_matches",
    "naming_policy_version", "name_source", "identity_candidate_rank", "library_collision",
    "identity_eligible", "identity_exclusion_reason", "exception_basis", "review_verdict",
    "excluded_library_matches", "raw_rank1_candidate", "identity_candidate",
    "reference_answer",
    "reference_formula",
    "rank",
    "candidate_smiles",
    "candidate_formula",
    "candidate_confidence",
    "rank_score_relative",
    "evidence_confidence",
    "evidence_confidence_basis",
    "mass_match_score",
    "abs_mass_error_da",
    "abs_mass_error_ppm",
    "truth_basis",
    "truth_tanimoto",
    "matched_truth",
    "is_best_match_in_top10",
    "matching_ranks",
    "pure_nce_score",
]

# Reviewed top-10 sheet: the columns kept in the hand-edited
# ``executive_summary_top10.xlsx`` (exec_summary_top10 tab). The by-compound
# sheet and the full CSV keep the rest.
TRUNCATED_COLUMNS = [
    "x_compound_id",
    "metabolon_level",
    "metabolon_sublevel",
    "known_x_tentative_information",
    "known_x_cleaned_external_info",
    "msi_level",
    "msi_label",
    "feature_name",
    "report_display_name",
    "reported_level",
    "reported_identification",
    "candidate_identity_confidence",
    "identity_candidate_rank",
    "rank",
    "candidate_smiles",
    "candidate_formula",
    "evidence_confidence",
    "mass_match_score",
    "abs_mass_error_ppm",
    "is_best_match_in_top10",
    "matching_ranks",
    "pure_nce_score",
]

ROLLUP_COLUMNS = [
    "x_compound_id",
    "cohort",
    "metabolon_level",
    "metabolon_sublevel",
    "known_x_tentative_information",
    "known_x_cleaned_external_info",
    "msi_level",
    "msi_label",
    "feature_name",
    "report_display_name",
    "reported_level",
    "reported_identification",
    "candidate_identity_confidence",
    "excluded_library_reference",
    "feature_name_on_prod_library_panel",
    "prod_library_panel_matches",
    "naming_policy_version", "name_source", "identity_candidate_rank", "library_collision",
    "identity_eligible", "identity_exclusion_reason", "exception_basis", "review_verdict",
    "excluded_library_matches", "raw_rank1_candidate", "identity_candidate",
    "reference_answer",
    "truth_basis",
    "n_candidates",
    "matched_at_rank_1",
    "matched_within_top_3",
    "matched_within_top_10",
    "first_matching_rank",
    "matching_ranks",
    "best_truth_tanimoto",
    "rank1_candidate_smiles",
    "rank1_candidate_formula",
    "rank1_candidate_confidence",
    "rank1_rank_score_relative",
    "rank1_evidence_confidence",
    "rank1_evidence_confidence_basis",
    "rank1_mass_match_score",
    "rank1_abs_mass_error_da",
    "rank1_abs_mass_error_ppm",
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


# Parts-per-million relative to the candidate exact mass (or, if that is
# missing, the observed precursor m/z). Same conversion as
# ``audit_mass_policy.ppm_of`` / ``adduct_neutral_mass_error_ppm``.
_MASS_REFERENCE_KEYS = (
    "predicted_exact_mass",
    "candidate_exact_mass",
    "theoretical_neutral_mass",
    "mass_hypothesis_neutral_mass",
    "target_neutral_mass",
    "parent_mz",
)


def _mass_reference(row: dict[str, Any], fallback: float | None = None) -> float | None:
    for key in _MASS_REFERENCE_KEYS:
        mass = _num(row.get(key))
        if mass is not None and mass > 0:
            return mass
    if fallback is not None and fallback > 0:
        return fallback
    return None


def abs_mass_error_ppm(error_da: float | None, reference_mass: float | None) -> float | None:
    """Convert an absolute Da error into ppm. Blank when either input is missing."""

    if error_da is None or reference_mass is None or reference_mass <= 0:
        return None
    return abs(error_da) / reference_mass * 1e6


# evidence_confidence is only emitted by some candidate generators; the rest reach
# the sidecar as 0.0 through the pool's missing-column fill, so an absent measurement
# reads as a measured zero. Fill it through the pipeline's canonical definition so the
# displayed number equals the reranker's own feature, and an uncomputed score renders
# blank instead of 0.0. See structural_elucidation/chemistry/evidence.py.
EVIDENCE_EMITTED = BASIS_EMITTED
EVIDENCE_BACKFILLED = BASIS_BACKFILLED
EVIDENCE_UNAVAILABLE = BASIS_UNAVAILABLE


def monotone_confidence(candidates: list[dict[str, Any]]) -> None:
    """Keep published confidence from rising as rank gets worse.

    Rank order is the rerank order and is not changed. Rank 1 keeps its
    confidence. A later candidate keeps its own value when that value is lower
    or equal, and is lowered to the previous published value when it is higher.
    Equal candidates stay equal. A missing confidence does not move the floor.
    """

    ordered = sorted(
        candidates,
        key=lambda candidate: candidate["rank"] if candidate.get("rank") is not None else 10**9,
    )
    floor: float | None = None
    for candidate in ordered:
        value = candidate.get("evidence_confidence")
        if value is None:
            continue
        if floor is not None and value > floor:
            candidate["evidence_confidence"] = floor
        else:
            floor = value


def _evidence_for_rows(rows: list[dict[str, Any]]) -> list[tuple[float | None, str]]:
    """Canonical evidence_confidence + basis for a report's candidate rows.

    ``value`` is ``None`` where basis is ``unavailable`` so the cell renders blank
    rather than a misleading ``0.0``.
    """
    if not rows:
        return []
    score, basis = backfill_evidence_confidence(pd.DataFrame(rows))
    out: list[tuple[float | None, str]] = []
    for position in range(len(rows)):
        label = str(basis.iloc[position])
        value = None if label == EVIDENCE_UNAVAILABLE else float(score.iloc[position])
        out.append((value, label))
    return out


def _display_evidence(row: dict[str, Any]) -> tuple[float | None, str]:
    """Single-row convenience wrapper over :func:`_evidence_for_rows`."""
    return _evidence_for_rows([row])[0]


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
    values are collected as sets. The displayed level is the best one across
    those libraries (see :func:`best_metabolon`); cohort assignment still sees
    every sublevel.
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
                compound,
                {
                    "levels": set(), "sublevels": set(), "names": set(), "formulae": set(),
                    "cleaned": set(), "external": set(), "known": set(), "level_rows": set(),
                },
            )
            level = _clean(row.get("Level"))
            sublevel = _clean(row.get("SubLevel"))
            if level or sublevel:
                entry["level_rows"].add((level, sublevel))
            if level:
                entry["levels"].add(level)
            if sublevel:
                entry["sublevels"].add(sublevel)
            if name := _clean(row.get("Tentative information")):
                entry["names"].add(name)
            if formula := _clean(row.get("Formula")):
                entry["formulae"].add(formula)
            if cleaned := _clean(row.get("Cleaned External Info")):
                entry["cleaned"].add(cleaned)
            if external := _clean(row.get("External Information (Class or modification)")):
                entry["external"].add(external)
            if known := _clean(row.get("KNOWN")):
                entry["known"].add(known)
    return table


def _joined(entry: dict[str, Any], key: str) -> str:
    return "; ".join(sorted(entry.get(key) or set()))


_LEVEL_NUMBER = re.compile(r"^(\d+)")


def _level_sort_key(level: str) -> tuple[int, int, int, str]:
    """Lower is better: a numbered level beats a label, a smaller number beats a
    larger one, and a definite level beats the same number with a ``?``."""

    match = _LEVEL_NUMBER.match(level.strip())
    if not match:
        return (1, 99, 0, level)
    return (0, int(match.group(1)), 1 if "?" in level else 0, level)


def best_metabolon(entry: dict[str, Any]) -> tuple[str, str]:
    """The best Metabolon level across libraries, and the sublevel(s) at that level.

    ``3`` in one library and ``5`` in another reports ``3``, not ``3; 5``.
    Sublevels are those recorded on the winning level (``3A; 3C`` when two
    libraries both sit at 3). Rows with no level do not count.
    """

    rows = [(level, sublevel) for level, sublevel in entry.get("level_rows") or set() if level]
    if not rows:
        return "", ""
    winning = min((level for level, _ in rows), key=_level_sort_key)
    sublevels = sorted({sublevel for level, sublevel in rows if level == winning and sublevel})
    return winning, "; ".join(sublevels)


# --- production library panel -------------------------------------------------
# The list of named metabolites on the production global library
# (CHEMICAL_ID / CHEMICAL_NAME / SUPER_PATHWAY / PATHWAY, 5,413 rows). The sheet
# marks whether the library metabolite a feature is named after, or collided
# with, is on it. Chemical IDs are the primary key; names are the fallback for
# labels that carry no ID (authentic standards, curated hypotheses).

PROD_LIBRARY_PANEL = PROJECT_ROOT / "data" / "prod_library" / "library_global_discovery_panel_v5.xlsx"
_LABEL_PREFIX = re.compile(r"^(?:isomer of|possible|probable)\s+", re.I)
_LABEL_MORE = re.compile(r"\s*\(\+\d+ more\)\s*$", re.I)


def load_prod_library(path: Path | None) -> dict[str, dict[str, str]]:
    """``{"ids": {chemical_id: name}, "names": {casefolded name: chemical_id}}``.

    Empty when no panel is supplied, in which case the panel columns are blank.
    """

    if path is None or not path.exists():
        return {"ids": {}, "names": {}}
    if path.suffix.lower() in {".xlsx", ".xlsm", ".xls"}:
        frame = pd.read_excel(path, dtype=str)
    else:
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    columns = {str(column).strip().upper(): column for column in frame.columns}
    id_column, name_column = columns.get("CHEMICAL_ID"), columns.get("CHEMICAL_NAME")
    if id_column is None or name_column is None:
        raise SystemExit(f"{path}: expected CHEMICAL_ID and CHEMICAL_NAME columns")
    ids: dict[str, str] = {}
    names: dict[str, str] = {}
    for chemical_id, name in zip(frame[id_column], frame[name_column]):
        chemical_id, name = _clean(chemical_id), _clean(name)
        if chemical_id.endswith(".0"):
            chemical_id = chemical_id[:-2]
        if not chemical_id:
            continue
        ids[chemical_id] = name
        if name:
            names.setdefault(name.casefold(), chemical_id)
    return {"ids": ids, "names": names}


def _label_names(label: str) -> list[str]:
    """Library names a sheet label refers to: ``Isomer of a / b (+2 more)`` -> [a, b]."""

    text = _LABEL_MORE.sub("", _LABEL_PREFIX.sub("", _clean(label)))
    if not text or _names_no_structure(text):
        return []
    return [part.strip() for part in text.split(ISOMER_NAME_SEPARATOR) if part.strip()]


def prod_library_marks(
    feature_name: str,
    fields: dict[str, Any],
    panel: dict[str, dict[str, str]],
    *,
    db_chemical_id: str = "",
) -> tuple[str, str]:
    """``(feature_name_on_prod_library_panel, prod_library_panel_matches)``.

    The first is TRUE/FALSE for a feature whose label names a library
    metabolite (by the rank-1 collision's chemical IDs for an ``Isomer of``
    label, else by the report's resolved library ``db_chemical_id``, else by
    name) and blank for formula / mass / tentative placeholders. The second
    lists every panel metabolite the report references anywhere (rank-1
    collision, deeper excluded matches, the identity candidate's library name,
    the label itself) as ``name (CHEMICAL_ID)``.
    """

    if not panel["ids"]:
        return "", ""
    ids, names = panel["ids"], panel["names"]

    def hits(chemical_ids: list[str], labels: list[str]) -> list[str]:
        found: list[str] = []
        for chemical_id in chemical_ids:
            if chemical_id in ids and chemical_id not in found:
                found.append(chemical_id)
        for label in labels:
            chemical_id = names.get(label.casefold())
            if chemical_id and chemical_id not in found:
                found.append(chemical_id)
        return found

    raw = fields.get("raw_rank1_candidate") or {}
    rank1_ids = [_clean(m.get("chemical_id")) for m in raw.get("matched_library") or []]
    label_names = _label_names(feature_name)
    is_isomer_label = _clean(feature_name).lower().startswith("isomer of ")
    db_ids = [_clean(db_chemical_id).removesuffix(".0")] if _clean(db_chemical_id) else []
    if label_names:
        label_hits = hits(rank1_ids if is_isomer_label else db_ids, label_names)
        flag = "TRUE" if label_hits else "FALSE"
    else:
        flag = ""

    all_ids = list(rank1_ids) + db_ids
    for item in fields.get("excluded_library_matches") or []:
        all_ids.extend(_clean(m.get("chemical_id")) for m in item.get("matches") or [])
    all_names = list(label_names)
    for candidate in (raw, fields.get("identity_candidate") or {}):
        all_names.extend(n for n in (_clean(candidate.get("name")), _clean(candidate.get("library_name"))) if n)
    matches = "; ".join(f"{ids[chemical_id]} ({chemical_id})" for chemical_id in hits(all_ids, all_names))
    return flag, matches


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


# --- "Isomer of" labelling for Naming v2 collisions ---------------------------
# Under Naming v2 an X-compound whose rank-1 structure is already a catalogued
# global-library metabolite cannot *be* that metabolite (it would not be an X),
# so the report withholds the name and records the collision in
# ``excluded_library_matches``. The X feature shares that metabolite's mass and
# formula, which is exactly what "Isomer of <metabolite>" says, so the sheet
# shows the relationship instead of a bare formula placeholder. The report's
# own label is kept alongside in ``report_display_name``.

NAMING_POLICY_V2 = "v2"
COLLISION_REASON = "existing_global_library_structure"
ISOMER_NAME_SEPARATOR = " / "
# Stereoisomers share a connectivity, so one collision can carry several library
# names. The label lists this many; ``excluded_library_reference`` has them all.
ISOMER_NAME_LIMIT = 3


def _names_no_structure(label: str) -> bool:
    """Whether a display label is a placeholder rather than a structural name."""

    text = _clean(label)
    lower = text.lower()
    return (
        not text
        or FORMULA_ONLY_NOTE in lower
        or lower.startswith("exact mass ")
        or lower.startswith(UNNAMED_STRUCTURE_LABEL.lower())
    )


def _library_names(matches: Any) -> list[str]:
    names: list[str] = []
    seen: set[str] = set()
    for match in matches or []:
        name = _clean((match or {}).get("name"))
        if name and name.casefold() not in seen:
            seen.add(name.casefold())
            names.append(name)
    return names


def excluded_library_reference(fields: dict[str, Any]) -> str:
    """Global-library metabolites the report excluded as identity collisions."""

    names: list[str] = []
    seen: set[str] = set()
    for item in fields.get("excluded_library_matches") or []:
        if _clean((item or {}).get("reason")) != COLLISION_REASON:
            continue
        for name in _library_names(item.get("matches")):
            if name.casefold() not in seen:
                seen.add(name.casefold())
                names.append(name)
    return "; ".join(names)


def isomer_of_feature_name(display_name: str, fields: dict[str, Any]) -> str:
    """Sheet label for a Naming v2 record.

    Applies only when the report itself named no structure (formula, exact
    mass, or ``Tentative Structure``) *and* its rank-1 candidate was excluded
    because that structure is already a global-library metabolite. A curated
    ``known_X`` hypothesis (``Possible *``), an authentic-standard identity, or
    any synthesized descriptor is the report's settled label and passes
    through unchanged.
    """

    text = _clean(display_name)
    if _clean(fields.get("naming_policy_version")) != NAMING_POLICY_V2:
        return text
    if _clean(fields.get("exception_basis")) or not _names_no_structure(text):
        return text
    raw = fields.get("raw_rank1_candidate") or {}
    if _clean(raw.get("identity_exclusion_reason")) != COLLISION_REASON:
        return text
    names = _library_names(raw.get("matched_library"))
    if not names:
        return text
    shown = ISOMER_NAME_SEPARATOR.join(names[:ISOMER_NAME_LIMIT])
    if len(names) > ISOMER_NAME_LIMIT:
        shown = f"{shown} (+{len(names) - ISOMER_NAME_LIMIT} more)"
    return isomer_of_label(shown)


def read_report(
    directory: Path, panel: dict[str, dict[str, str]] | None = None
) -> dict[str, Any] | None:
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

    evidence_by_position = _evidence_for_rows(rows)
    observed_mz = _num(fields.get("observed_mz"))
    candidates = []
    for position, row in enumerate(rows):
        rank_value = _num(row.get("rank"))
        rank = int(rank_value) if rank_value is not None else len(candidates) + 1
        extra = evidence.get(rank, {})
        evidence_value, evidence_basis = evidence_by_position[position]
        error_da = _num(row.get("abs_mass_error_da"))
        candidates.append(
            {
                "rank": rank,
                "smiles": _clean(row.get("predicted_smiles")) or _clean(row.get("smiles")),
                "formula": _clean(row.get("predicted_formula"))
                or _clean(row.get("candidate_formula")),
                "candidate_confidence": _num(extra.get("candidate_confidence")),
                "rank_score_relative": _num(row.get("candidate_report_score")),
                "evidence_confidence": evidence_value,
                "evidence_confidence_basis": evidence_basis,
                "mass_match_score": _num(row.get("mass_match_score")),
                "abs_mass_error_da": error_da,
                "abs_mass_error_ppm": abs_mass_error_ppm(error_da, _mass_reference(row, observed_mz)),
                "truth_tanimoto": _num(row.get("truth_tanimoto")),
            }
        )

    display_name = _clean(fields.get("display_name"))
    report_display_name = display_name
    if fields.get("naming_policy_version") == NAMING_POLICY_V2:
        display_name = isomer_of_feature_name(display_name, fields)
    else:
        display_name = sheet_feature_name(
            display_name,
            name_resolved=bool(fields.get("name_resolved")),
            source=_clean(fields.get("source")),
            msi_level=fields.get("msi_level"),
            identity_status=_clean(fields.get("identity_status")),
            truth_match=bool(fields.get("truth_match")),
            descriptor_kind=_clean(fields.get("descriptor_kind")),
            from_known_x=_clean(fields.get("descriptor_source")) == "known_x_row"
            or bool(read_known_x_name(directory)),
        )
    evidence_json = directory / "data" / f"{stem}_evidence_summary.json"
    db_chemical_id = ""
    if evidence_json.exists():
        try:
            evidence_summary = json.loads(evidence_json.read_text())
        except ValueError:
            evidence_summary = {}
        db_chemical_id = _clean((evidence_summary.get("db_detail") or {}).get("chemical_id"))
    on_panel, panel_matches = prod_library_marks(
        display_name, fields, panel or {"ids": {}, "names": {}}, db_chemical_id=db_chemical_id
    )
    return {
        "x_compound_id": directory.name,
        "display_name": display_name,
        "report_display_name": report_display_name,
        "excluded_library_reference": excluded_library_reference(fields),
        "feature_name_on_prod_library_panel": on_panel,
        "prod_library_panel_matches": panel_matches,
        "msi_level": fields.get("msi_level"),
        "msi_label": _clean(fields.get("msi_label")),
        "candidates": candidates,
        "naming_policy_version": _clean(fields.get("naming_policy_version")),
        "name_source": _clean(fields.get("name_source")),
        "identity_candidate_rank": fields.get("identity_candidate_rank"),
        "library_collision": fields.get("library_collision"),
        "identity_eligible": fields.get("identity_eligible"),
        "identity_exclusion_reason": _clean(fields.get("identity_exclusion_reason")),
        "exception_basis": _clean(fields.get("exception_basis")),
        "review_verdict": _clean(fields.get("review_verdict")),
        "excluded_library_matches": json.dumps(fields.get("excluded_library_matches") or [], sort_keys=True),
        "raw_rank1_candidate": json.dumps(fields.get("raw_rank1_candidate") or {}, sort_keys=True),
        "identity_candidate": json.dumps(fields.get("identity_candidate") or {}, sort_keys=True),
    }


def read_reports(
    root: Path, panel: dict[str, dict[str, str]] | None = None
) -> Iterator[dict[str, Any]]:
    for directory in sorted(root.iterdir()):
        if not directory.is_dir() or not directory.name.upper().startswith("X-"):
            continue
        if (report := read_report(directory, panel)) is not None:
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
    monotone_confidence(report["candidates"])
    matches = [match_flag(candidate, context) for candidate in report["candidates"]]
    first_match = next(
        (candidate["rank"] for candidate, hit in zip(report["candidates"], matches) if hit),
        None,
    )
    matching_ranks = "" if first_match is None else str(first_match)
    tanimotos = [
        candidate["truth_tanimoto"]
        for candidate in report["candidates"]
        if candidate["truth_tanimoto"] is not None
    ]
    level, sublevel = best_metabolon(entry)
    naming_fields = {
        "feature_name": report["display_name"],
        "report_display_name": report.get("report_display_name", ""),
        "metabolon_level": level,
        "msi_level": report.get("msi_level"),
        "exception_basis": report.get("exception_basis", ""),
        "identity_exclusion_reason": report.get("identity_exclusion_reason", ""),
        "excluded_library_reference": report.get("excluded_library_reference", ""),
        "library_collision": report.get("library_collision", ""),
        "known_x_tentative_information": _joined(entry, "names"),
        "known_x_cleaned_external_info": _joined(entry, "cleaned"),
    }
    identification = reported_identification(
        report["x_compound_id"],
        report["candidates"],
        naming_fields,
    )
    candidate_identity_confidence_value = candidate_identity_confidence(
        report["candidates"],
        naming_fields,
    )
    shared = {
        "x_compound_id": report["x_compound_id"],
        "cohort": context["cohort"],
        "metabolon_level": level,
        "metabolon_sublevel": sublevel,
        "known_x_tentative_information": _joined(entry, "names"),
        "known_x_cleaned_external_info": _joined(entry, "cleaned"),
        "msi_level": "" if report["msi_level"] is None else report["msi_level"],
        "msi_label": report["msi_label"],
        "feature_name": report["display_name"],
        "report_display_name": report.get("report_display_name", ""),
        "reported_level": reported_level(identification),
        "reported_identification": identification,
        "candidate_identity_confidence": "" if candidate_identity_confidence_value is None else candidate_identity_confidence_value,
        "excluded_library_reference": report.get("excluded_library_reference", ""),
        "feature_name_on_prod_library_panel": report.get("feature_name_on_prod_library_panel", ""),
        "prod_library_panel_matches": report.get("prod_library_panel_matches", ""),
        "naming_policy_version": report.get("naming_policy_version", ""),
        "name_source": report.get("name_source", ""),
        "identity_candidate_rank": report.get("identity_candidate_rank") or "",
        "library_collision": report.get("library_collision", ""),
        "identity_eligible": report.get("identity_eligible", ""),
        "identity_exclusion_reason": report.get("identity_exclusion_reason", ""),
        "exception_basis": report.get("exception_basis", ""),
        "review_verdict": report.get("review_verdict", ""),
        "excluded_library_matches": report.get("excluded_library_matches", ""),
        "raw_rank1_candidate": report.get("raw_rank1_candidate", ""),
        "identity_candidate": report.get("identity_candidate", ""),
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
                "evidence_confidence_basis": candidate.get("evidence_confidence_basis", ""),
                "mass_match_score": _round(candidate["mass_match_score"], 6),
                "abs_mass_error_da": _round(candidate["abs_mass_error_da"], 6),
                "abs_mass_error_ppm": _round(candidate["abs_mass_error_ppm"], 2),
                "truth_basis": context["truth_basis"],
                "truth_tanimoto": _round(candidate["truth_tanimoto"], 4),
                "matched_truth": "" if hit is None else str(hit).upper(),
                "is_best_match_in_top10": str(
                    first_match is not None and candidate["rank"] == first_match
                ).upper(),
                "matching_ranks": matching_ranks,
                "pure_nce_score": _round(candidate.get("pure_nce_score"), 6),
            }
        )

    top = report["candidates"][0] if report["candidates"] else {}
    assessable = context["truth_basis"] != BASIS_NONE
    rollup = {
        **shared,
        "truth_basis": context["truth_basis"],
        "n_candidates": len(report["candidates"]),
        "matched_at_rank_1": str(first_match == 1).upper() if assessable else "",
        "matched_within_top_3": str(first_match is not None and first_match <= 3).upper()
        if assessable
        else "",
        "matched_within_top_10": str(first_match is not None).upper() if assessable else "",
        "first_matching_rank": "" if first_match is None else first_match,
        "matching_ranks": matching_ranks,
        "best_truth_tanimoto": _round(max(tanimotos), 4) if tanimotos else "",
        "rank1_candidate_smiles": top.get("smiles", ""),
        "rank1_candidate_formula": top.get("formula", ""),
        "rank1_candidate_confidence": _round(top.get("candidate_confidence"), 4),
        "rank1_rank_score_relative": _round(top.get("rank_score_relative"), 6),
        "rank1_evidence_confidence": _round(top.get("evidence_confidence"), 6),
        "rank1_evidence_confidence_basis": top.get("evidence_confidence_basis", ""),
        "rank1_mass_match_score": _round(top.get("mass_match_score"), 6),
        "rank1_abs_mass_error_da": _round(top.get("abs_mass_error_da"), 6),
        "rank1_abs_mass_error_ppm": _round(top.get("abs_mass_error_ppm"), 2),
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
        writer = csv.DictWriter(
            handle, fieldnames=columns, lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def load_pure_nce_scores(path: Path) -> dict[tuple[str, str], float]:
    """(x_compound_id, candidate_smiles) -> pure NCE score; unscorable pairs are left out."""

    scores: dict[tuple[str, str], float] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("nce_status") in {"no_observed_bins", "error"}:
                continue
            value = _num(row.get("pure_nce_score"))
            if value is not None:
                scores[(_clean(row.get("x_compound_id")), _clean(row.get("candidate_smiles")))] = value
    return scores


def attach_pure_nce_scores(report: dict[str, Any], scores: dict[tuple[str, str], float]) -> None:
    for candidate in report["candidates"]:
        candidate["pure_nce_score"] = scores.get((report["x_compound_id"], candidate["smiles"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports-root", type=Path, required=True)
    parser.add_argument(
        "--out-dir", type=Path, default=Path("docs/paper/data/exec_summary"), help="write CSVs here"
    )
    parser.add_argument(
        "--prod-library",
        type=Path,
        default=PROD_LIBRARY_PANEL,
        help="production library panel (CHEMICAL_ID / CHEMICAL_NAME); pass a missing path to skip",
    )
    parser.add_argument(
        "--pure-nce-scores",
        type=Path,
        default=None,
        help="tools/rescore_pure_nce.py output; drives the candidate_identity_confidence margin",
    )
    args = parser.parse_args()

    if not args.reports_root.is_dir():
        raise SystemExit(f"reports root not found: {args.reports_root}")

    pure_nce = load_pure_nce_scores(args.pure_nce_scores) if args.pure_nce_scores else {}
    levels = load_metabolon_levels()
    structures = load_reference_structures()
    panel = load_prod_library(args.prod_library)
    if panel["ids"]:
        print(f"[exec-summary] production library panel: {len(panel['ids'])} entries from {args.prod_library}")
    else:
        print(f"[exec-summary] no production library panel at {args.prod_library}; panel columns left blank")

    long_rows: list[dict[str, Any]] = []
    rollups: list[dict[str, Any]] = []
    for report in read_reports(args.reports_root, panel):
        attach_pure_nce_scores(report, pure_nce)
        context = classify(report, levels, structures)
        entry = levels.get(report["x_compound_id"], {})
        rows, rollup = build_rows(report, context, entry)
        long_rows.extend(rows)
        rollups.append(rollup)

    if not rollups:
        raise SystemExit(f"no X-* report directories under {args.reports_root}")

    summary = summarise(rollups)
    write_csv(args.out_dir / "exec_summary_top10.csv", LONG_COLUMNS, long_rows)
    write_csv(args.out_dir / "exec_summary_top10_truncated.csv", TRUNCATED_COLUMNS, long_rows)
    write_csv(args.out_dir / "exec_summary_by_compound.csv", ROLLUP_COLUMNS, rollups)
    (args.out_dir / "exec_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(f"[exec-summary] {summary['n_compounds']} compounds, {len(long_rows)} candidate rows")
    controls = summary["truth_control"]
    def pct(value: float | None) -> str:
        return f"{value:.1%}" if value is not None else "n/a"
    print(
        f"[exec-summary] truth controls: {controls['n']} | structural match at rank 1 "
        f"{controls['rank_1']['n_matched']}/{controls['rank_1']['n_assessable']}"
        f" ({pct(controls['rank_1']['rate'])}), top 3 {pct(controls['top_3']['rate'])},"
        f" top 10 {pct(controls['top_10']['rate'])}"
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
    print(f"[exec-summary] wrote {args.out_dir}/exec_summary_top10_truncated.csv")


if __name__ == "__main__":
    main()
