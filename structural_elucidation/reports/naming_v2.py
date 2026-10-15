"""Identity-safe naming policy for X-compound reports.

An InChIKey connectivity match is a *candidate structure* lookup.  It never
establishes that an X feature is the already catalogued global-library peak.
The raw rank is retained independently from the candidate allowed to drive the
report.  Only a recorded authentic-standard confirmation can waive the gate.
"""
from __future__ import annotations

import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from structural_elucidation.utils.denovo import inchikey_from_smiles


POLICY_VERSION = "v2"
_GENERIC = re.compile(
    r"\b(unknown|unidentified|compound|molecule|derivative|metabolite|" 
    r"amino acid derivative|no info|not available)\b", re.I
)
_HEDGE = re.compile(r"^(?:possible|possibly|putative|tentative|likely|probably)\s+", re.I)


def _clean(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text.lower() in {"nan", "none", "null", "<na>"} else text


def _num(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if number == number else None


def _true(value: Any) -> bool:
    return "true" in _clean(value).lower() or _clean(value).lower() in {"yes", "1"}


def connectivity(smiles: str) -> str:
    if not _clean(smiles):
        return ""
    try:
        return inchikey_from_smiles(smiles).split("-")[0]
    except Exception:
        return ""


@dataclass(frozen=True)
class LibraryRecord:
    chemical_id: str
    name: str
    connectivity: str

    def as_dict(self) -> dict[str, str]:
        return {"chemical_id": self.chemical_id, "name": self.name,
                "connectivity": self.connectivity}


class CollisionIndex:
    """Global-library structures and names from the report reference SQLite."""

    def __init__(self, db_path: Path | str, *, supplied_names: Sequence[str] = ()):
        self.by_connectivity: dict[str, list[LibraryRecord]] = {}
        self.by_name: dict[str, list[LibraryRecord]] = {}
        self.by_chemical_id: dict[str, list[LibraryRecord]] = {}
        uri = Path(db_path).resolve().as_uri() + "?mode=ro"
        with sqlite3.connect(uri, uri=True) as connection:
            rows = connection.execute(
                "SELECT chemical_id, canonical_name, inchikeys, smiles, in_global_library "
                "FROM metabolites ORDER BY chemical_id"
            )
            for chemical_id, name, keys, smiles, global_flag in rows:
                if not _true(global_flag):
                    continue
                blocks = sorted({
                    token.strip().split("-")[0]
                    for token in re.split(r"[|;,]", _clean(keys))
                    if len(token.strip().split("-")[0]) == 14
                })
                if not blocks:
                    key = connectivity(_clean(smiles).split("|")[0])
                    blocks = [key] if key else []
                for block in blocks:
                    record = LibraryRecord(str(chemical_id), str(name), block)
                    self.by_connectivity.setdefault(block, []).append(record)
                    self.by_name.setdefault(str(name).casefold().strip(), []).append(record)
                    self.by_chemical_id.setdefault(str(chemical_id), []).append(record)
            # A supplied library_name may be a synonym rather than the
            # canonical label. Resolve it against the same global-library
            # subset before allowing it to name a candidate structure.
            if supplied_names and connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='synonym_index'"
            ).fetchone():
                for supplied in sorted({_clean(name) for name in supplied_names if _clean(name)}):
                    for chemical_id, value in connection.execute(
                        "SELECT chemical_id, entry_value FROM synonym_index "
                        "WHERE entry_type='synonym' AND entry_value=?", (supplied,)
                    ):
                        bucket = self.by_name.setdefault(_clean(value).casefold(), [])
                        for record in self.by_chemical_id.get(str(chemical_id), []):
                            if record not in bucket:
                                bucket.append(record)

    def match_smiles(self, smiles: str) -> list[LibraryRecord]:
        return self.by_connectivity.get(connectivity(smiles), [])

    def name_conflict(self, supplied_name: str, smiles: str) -> list[LibraryRecord]:
        name = _clean(supplied_name).casefold()
        named = self.by_name.get(name, [])
        key = connectivity(smiles)
        return named if named and key and all(item.connectivity != key for item in named) else []


def curated_hypothesis(known_x: Mapping[str, Any]) -> str:
    """Use only a specific, explicitly tentative annotation from this X row."""
    tentative = _clean(known_x.get("Tentative information"))
    cleaned = _clean(known_x.get("Cleaned External Info"))
    source = tentative if tentative and _HEDGE.match(tentative) else cleaned
    source = _HEDGE.sub("", source).split(";")[0].strip().rstrip("?.")
    source = re.split(r"\s*,\s*(?:correlated|associated|see)\b", source, maxsplit=1, flags=re.I)[0]
    if not source or len(source) < 6 or _GENERIC.search(source):
        return ""
    if re.fullmatch(r"(?:sulfate|glucuronide|phosphate|glycine|taurine)\s+of\s+[A-Z][A-Za-z0-9+.-]+", source, re.I):
        return ""
    if source.lower() in {"sulfated indole", "sulfate", "glucuronide", "lipid", "amino acid"}:
        return ""
    return source


def _candidate_dict(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "model_rank": int(_num(row.get("rank")) or 0),
        "smiles": _clean(row.get("predicted_smiles") or row.get("smiles")),
        "formula": _clean(row.get("predicted_formula") or row.get("candidate_formula")),
        "name": _clean(row.get("candidate_name") or row.get("base_name") or row.get("library_name")),
        "library_name": _clean(row.get("library_name")),
        "mass_match_score": _num(row.get("mass_match_score")),
        "abs_mass_error_da": _num(row.get("abs_mass_error_da")),
        "candidate_report_score": _num(row.get("candidate_report_score")),
    }


def classify_candidates(
    rows: Sequence[Mapping[str, Any]],
    index: CollisionIndex,
    *,
    formula: str = "",
    standard_confirmed: bool = False,
    confirmed_structure: str = "",
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for row in rows:
        item = _candidate_dict(row)
        matches = index.match_smiles(item["smiles"])
        conflicts = index.name_conflict(item["library_name"], item["smiles"])
        observed_formula = _clean(formula)
        formula_ok = not observed_formula or not item["formula"] or item["formula"] == observed_formula
        mass_score = item["mass_match_score"]
        error = item["abs_mass_error_da"]
        mass_ok = (mass_score is not None and mass_score >= 0.95) or (
            error is not None and error <= 0.01
        )
        if conflicts:
            reason = "library_name_structure_conflict"
        elif matches and not (
            standard_confirmed and (
                not confirmed_structure or connectivity(item["smiles"]) == connectivity(confirmed_structure)
            )
        ):
            reason = "existing_global_library_structure"
        elif not item["smiles"]:
            reason = "missing_structure"
        elif not formula_ok:
            reason = "formula_mismatch"
        elif not mass_ok:
            reason = "mass_gate_failed"
        else:
            reason = ""
        item.update({
            "library_collision": bool(matches),
            "matched_library": [m.as_dict() for m in matches],
            "name_structure_conflicts": [m.as_dict() for m in conflicts],
            "identity_eligible": not bool(reason),
            "identity_exclusion_reason": reason,
        })
        out.append(item)
    return out


def decide(
    rows: Sequence[Mapping[str, Any]],
    index: CollisionIndex,
    known_x: Mapping[str, Any],
    *,
    assurance: str = "none",
    confirmed_structure: str = "",
) -> dict[str, Any]:
    annotated = classify_candidates(
        rows, index, formula=_clean(known_x.get("Formula")),
        standard_confirmed=assurance == "standard_confirmed",
        confirmed_structure=confirmed_structure,
    )
    raw = annotated[0] if annotated else None
    eligible = next((item for item in annotated if item["identity_eligible"]), None)
    confirmed = next((item for item in annotated if item["identity_eligible"] and
                      (not confirmed_structure or connectivity(item["smiles"]) == connectivity(confirmed_structure))), None)
    curated = "" if assurance == "standard_confirmed" else curated_hypothesis(known_x)
    if curated and index.by_name.get(curated.casefold()):
        curated = ""
    if assurance == "standard_confirmed" and confirmed:
        source, selected = "authentic_standard", confirmed
    elif curated:
        source, selected = "known_x_tentative", None
    elif eligible:
        source, selected = "ranked_candidate", eligible
    else:
        source, selected = "formula_or_mass", None
    excluded = [
        {"model_rank": item["model_rank"], "reason": item["identity_exclusion_reason"],
         "matches": item["matched_library"], "name_structure_conflicts": item["name_structure_conflicts"]}
        for item in annotated if not item["identity_eligible"]
    ]
    return {
        "naming_policy_version": POLICY_VERSION,
        "name_source": source,
        "curated_hypothesis": curated,
        "raw_rank1_candidate": raw,
        "identity_candidate": selected,
        "identity_candidate_rank": selected["model_rank"] if selected else None,
        "library_collision": any(item["library_collision"] for item in annotated),
        "identity_eligible": bool(selected) or bool(curated),
        "identity_exclusion_reason": raw["identity_exclusion_reason"] if raw else "no_candidates",
        "excluded_library_matches": excluded,
        "candidate_annotations": annotated,
        "exception_basis": assurance if source == "authentic_standard" else "",
        "review_verdict": assurance if assurance not in {"none", "standard_confirmed"} else "",
    }
