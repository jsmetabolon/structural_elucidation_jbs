"""Display-only structural descriptors for unnamed report compounds.

The resolver deliberately does not inspect or alter calibrated identification
confidence.  It synthesizes a tentative display name from rank-1 candidate names,
peptide, lipid, motif, conjugate-deconjugation, and extended lookup helpers, in
that priority order.
"""
from __future__ import annotations

import csv
import math
import re
import sqlite3
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

try:
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem, rdMolDescriptors

    HAS_RDKIT = True
except ImportError:  # pragma: no cover - optional dependency
    Chem = None
    rdBase = None
    AllChem = None
    rdMolDescriptors = None
    HAS_RDKIT = False


@dataclass(frozen=True)
class TentativeName:
    """A non-calibrated structural descriptor for report presentation."""

    kind: str
    label: str
    parent_name: str | None = None
    parent_inchikey: str | None = None
    parent_formula: str | None = None
    conjugate_count: int = 0
    motifs: tuple[str, ...] = ()
    confidence: str = "low"
    source: str = "none"
    analog_neighbor_name: str | None = None
    analog_neighbor_inchikey: str | None = None
    analog_neighbor_tanimoto: float | None = None
    class_superclass: str | None = None
    class_class: str | None = None
    class_subclass: str | None = None
    iupac_name: str | None = None
    iupac_round_trip_ok: bool | None = None

    def as_dict(self) -> dict[str, Any]:
        """Serialize the descriptor for evidence-summary sidecars."""

        payload = asdict(self)
        payload["motifs"] = list(self.motifs)
        return payload


_LIPID_CLASS_NAMES = {
    "PC": "Phosphatidylcholine",
    "LPC": "Lysophosphatidylcholine",
    "PE": "Phosphatidylethanolamine",
    "LPE": "Lysophosphatidylethanolamine",
    "PS": "Phosphatidylserine",
    "LPS": "Lysophosphatidylserine",
    "PI": "Phosphatidylinositol",
    "LPI": "Lysophosphatidylinositol",
    "PG": "Phosphatidylglycerol",
    "LPG": "Lysophosphatidylglycerol",
    "SM": "Sphingomyelin",
    "Cer": "Ceramide",
    "HexCer": "Hexosylceramide",
    "DAG": "Diacylglycerol",
    "TAG": "Triacylglycerol",
    "MAG": "Monoacylglycerol",
    "sterol": "Sterol",
}

_CONJUGATE_TERMS = {
    "glucuronide": "glucuronide",
    "sulfate": "sulfate",
    "phosphate": "phosphate",
    "glycine": "glycine",
    "taurine": "taurine",
    "carnitine": "carnitine",
    "lactoyl": "lactoyl",
    "acetyl": "acetyl",
    "methyl": "methyl",
}
_CARBON_FREE_CONJUGATES = frozenset({"sulfate", "phosphate"})
DECORATING_CONJUGATES = frozenset(
    {"glucuronide", "sulfate", "glycine", "taurine", "carnitine", "lactoyl", "acetyl", "methyl"}
)
TEXT_CONJUGATE_TRIGGERS = {
    "glucuron": "glucuronide",
    "sulfat": "sulfate",
    "sulph": "sulfate",
    "taurin": "taurine",
    "glycin": "glycine",
    "carnitin": "carnitine",
    "lactoyl": "lactoyl",
}
_TEXT_CONJUGATE_TRIGGER_ORDER = (
    "glucuron",
    "sulfat",
    "sulph",
    "taurin",
    "glycin",
    "carnitin",
    "lactoyl",
)
_STEROID_CLASS_TERMS = {
    "glucuronide": "glucuronide",
    "sulfate": "sulfate",
    "taurine": "taurine conjugate",
    "glycine": "glycine conjugate",
    "carnitine": "carnitine conjugate",
    "lactoyl": "lactoyl conjugate",
}
_STEROID_CLASS_TOKENS = frozenset(
    {
        "steroid",
        "sterol",
        "bile",
        "acid",
        "hormone",
        "organic",
        "molecule",
        "compound",
        "glucuronidated",
        "sulfated",
        "possible",
    }
)
_HEDGE_PREFIX_RE = re.compile(
    r"^(possible|possibly|likely|probably|isomer of|something similar to|a |an )\s+",
    re.I,
)
_SULFATE_SMILES_RE = re.compile(r"S\(=O\)\(=O\)")
_GLUCURONIDE_SMILES_RE = re.compile(
    r"C\(C\(=O\)O\)C\([@#=+\-\\/0-9A-Za-z]*O\)C\([@#=+\-\\/0-9A-Za-z]*O\)C",
    re.I,
)
_MULTIPLICITY_PREFIXES = {
    2: "di",
    3: "tri",
    4: "tetra",
}

_MOTIF_PHRASES = {
    "glucuronide": "glucuronidated",
    "sulfate": "sulfated",
    "sulfonate": "sulfonated",
    "sulfur_oxygen_rich": "sulfur-rich oxygenated",
    "phosphate": "phosphorylated",
    "phosphate_like": "phosphate-like",
    "glycine": "glycine-conjugated",
    "taurine": "taurine-conjugated",
    "carnitine": "carnitine-conjugated",
    "isoalloxazine": "isoalloxazine-containing",
}

_MOTIF_PRIORITY = (
    "glucuronide",
    "sulfate",
    "sulfonate",
    "phosphate",
    "glycine",
    "taurine",
    "carnitine",
    "isoalloxazine",
    "phosphate_like",
    "sulfur_oxygen_rich",
)
_SPECIFIC_MOTIFS = frozenset(_MOTIF_PHRASES)
_PEPTIDE_LENGTHS = {
    1: "Amino acid",
    2: "Dipeptide",
    3: "Tripeptide",
    4: "Tetrapeptide",
}

GENERIC_CLASS_LABELS = frozenset(
    {
        "Sterol-like compound",
        "Phosphatidylcholine",
        "Phosphatidylethanolamine",
        "Phosphatidylserine",
        "Phosphatidylinositol",
        "Phosphatidylglycerol",
        "Sphingomyelin",
        "Ceramide",
        "Hexosylceramide",
        "Diacylglycerol",
        "Triacylglycerol",
        "Monoacylglycerol",
        "Sterol",
        "Lysophosphatidylcholine",
        "Lysophosphatidylethanolamine",
        "Lysophosphatidylserine",
        "Lysophosphatidylinositol",
        "Lysophosphatidylglycerol",
    }
    | set(_LIPID_CLASS_NAMES.values())
)
_GENERIC_LABEL_LIPID_CLASS = {
    label: code
    for code, label in _LIPID_CLASS_NAMES.items()
}
_GENERIC_LABEL_LIPID_CLASS["Sterol-like compound"] = "sterol"
_GENERIC_LOWER = frozenset(label.lower() for label in GENERIC_CLASS_LABELS)

_TRAILING_PAREN_RE = re.compile(r"\s*\([^)]*\)\s*$")
_CURATOR_NOTE_RES = (
    re.compile(r"\s*\([^()]*\?\s*\)\s*$"),
    re.compile(r"\s*\(\s*\d+\s*\)\s*$"),
)
_DESCRIPTIVE_RE = re.compile(r"\b(small organic|unknown|unidentified|putative|-?like)\b", re.I)
_NAME_TAG_SUFFIX_RE = re.compile(r"\s*[-–]\s*RELEASED\s*$", re.IGNORECASE)
_RELEASED_STRUCTURE_SUFFIX_RE = re.compile(r"\s+released\s+structure\s*$", re.IGNORECASE)
_SOURCE_PRIORITY = {
    "rank1_candidate_name": 0,
    "known_x_row": 1,
    "html_alt": 2,
}
_SMILES_LIKE_RE = re.compile(r"^[\[\]@#=+\-\\/.()0-9A-Za-z]+$")
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}(-[A-Z]{10}(-[A-Z])?)?$")
_FORMULA_ONLY_RE = re.compile(r"^C\d+H\d+[A-Za-z0-9]*$")


def _clean(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "nan", "none", "null"} else text


def _num(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _boolish(value: Any) -> bool:
    text = _clean(value).lower()
    if text in {"true", "t", "yes", "y", "1", "pass", "passed", "supported"}:
        return True
    if text in {"false", "f", "no", "n", "0", ""}:
        return False
    number = _num(value)
    return number is not None and number > 0


def _format_label(base: str) -> str:
    return " ".join(base.split())


# --- MSI-consistent naming ---------------------------------------------------
# What a report is allowed to call a compound follows from its MSI level, not
# from the confidence score. The score does not order correctness well enough to
# carry a naming decision (see docs/paper/data/thresholds/), whereas the level
# already encodes what kind of identification the evidence supports: a name at
# Levels 1-2, a hedged name at Level 3, formula at Level 4, mass at Level 5.

HEDGE_SUFFIX = "(proposed structure)"
FORMULA_ONLY_NOTE = "structure unresolved"
FORMULA_ONLY_LABEL = "Molecular formula only; structure unresolved"
MASS_ONLY_LABEL = "Exact mass of interest; unidentified"


def msi_int(value: Any) -> int | None:
    """Coerce an MSI level that may arrive as an int, a float, or ``"3"``."""

    try:
        level = int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None
    return level if 1 <= level <= 5 else None


def hedged_label(label: str) -> str:
    """Mark a name as a proposal, unless it already reads as one."""

    text = _format_label(label)
    if not text or HEDGE_SUFFIX in text:
        return text
    return f"{text} {HEDGE_SUFFIX}"


def display_name_for_level(
    name: str,
    *,
    msi_level: Any,
    formula: str = "",
    observed_mz: float | None = None,
    name_is_reference: bool = False,
) -> str:
    """Return the name a report of this MSI level is entitled to display.

    ``name_is_reference`` distinguishes a name copied from the reference database
    from a descriptor the pipeline synthesized. Only the former needs hedging at
    Level 3: a synthesized descriptor such as "Steroid glucuronide" already
    states its own uncertainty, and appending a hedge to it would be noise.

    Levels 4 and 5 withhold the structure entirely, so they discard the name
    rather than qualify it. A name that survives to the reader implies a
    structural claim, and at those levels there is none to make.
    """

    text = _format_label(name)
    level = msi_int(msi_level)
    if level is None or level <= 2:
        return text
    if level == 3:
        return hedged_label(text) if name_is_reference else text
    if level == 4:
        formula_text = _format_label(formula)
        return f"{formula_text}; {FORMULA_ONLY_NOTE}" if formula_text else FORMULA_ONLY_LABEL
    if observed_mz is not None:
        try:
            return f"Exact mass {float(observed_mz):.4f}; unidentified"
        except (TypeError, ValueError):
            pass
    return MASS_ONLY_LABEL


def _sentence_case_label(base: str) -> str:
    label = _format_label(base)
    if label and label[0].isalpha():
        return f"{label[0].upper()}{label[1:]}"
    return label


def descriptor_fields(name: TentativeName | None) -> dict[str, Any]:
    """Return the stable descriptor subset stored in confidence sidecars."""

    source_row = None
    if name and name.kind == "candidate_name":
        if name.source == "rank1_candidate_name":
            source_row = "rank1"
        elif name.source == "known_x_row":
            source_row = "known_x_row"
        elif name.source == "html_alt":
            source_row = "html_alt"
    return {
        "descriptor_kind": name.kind if name else None,
        "descriptor_parent_name": name.parent_name if name else None,
        "descriptor_parent_inchikey": name.parent_inchikey if name else None,
        "descriptor_parent_formula": name.parent_formula if name else None,
        "descriptor_conjugate_count": name.conjugate_count if name else 0,
        "descriptor_motifs": list(name.motifs) if name else [],
        "descriptor_confidence": name.confidence if name else None,
        "descriptor_source": name.source if name else "none",
        "descriptor_source_row": source_row,
        "analog_neighbor_name": name.analog_neighbor_name if name else None,
        "analog_neighbor_inchikey": name.analog_neighbor_inchikey if name else None,
        "analog_neighbor_tanimoto": name.analog_neighbor_tanimoto if name else None,
        "class_superclass": name.class_superclass if name else None,
        "class_class": name.class_class if name else None,
        "class_subclass": name.class_subclass if name else None,
        "iupac_name": name.iupac_name if name else None,
        "iupac_round_trip_ok": name.iupac_round_trip_ok if name else None,
    }


def _strip_trailing_parenthetical(name: str) -> str:
    text = _clean(name)
    while text:
        match = _TRAILING_PAREN_RE.search(text)
        if not match:
            break
        text = text[: match.start()].strip()
    return text


def _strip_curator_notes(name: str) -> str:
    text = _clean(name)
    changed = True
    while text and changed:
        changed = False
        for pattern in _CURATOR_NOTE_RES:
            match = pattern.search(text)
            if match:
                text = text[: match.start()].strip()
                changed = True
    return text


def _strip_stray_trailing_paren(name: str) -> str:
    text = _clean(name)
    opens = text.count("(")
    closes = text.count(")")
    while closes > opens and text.endswith(")"):
        text = text[:-1].rstrip()
        closes -= 1
    return text


def _clean_candidate_display_name(name: str) -> str:
    """Strip upstream workflow tags from adopted candidate names."""

    text = _strip_curator_notes(name)
    text = _NAME_TAG_SUFFIX_RE.sub("", text)
    text = _RELEASED_STRUCTURE_SUFFIX_RE.sub("", text)
    text = _strip_stray_trailing_paren(text)
    return _format_label(text)


def _clean_known_x_name(raw: str) -> str:
    text = _clean(raw)
    if text.upper().startswith("KNOWN "):
        text = text[6:].strip()
    return _clean_candidate_display_name(text)


def _is_descriptive_phrase(name: str) -> bool:
    text = _clean(name).lower()
    if _DESCRIPTIVE_RE.search(text):
        return True
    return " of " in text and any(
        token in text for token in ("compound", "molecule", "alcohol", "acid derivative")
    )


def _candidate_specificity(name: str) -> int:
    text = _format_label(name)
    if not text or _is_descriptive_phrase(text) or text.lower() in _GENERIC_LOWER:
        return 0
    return 1


def _looks_like_chemical_name(name: str, *, control_id: str = "") -> bool:
    text = _clean(name)
    if len(text) < 3:
        return False
    if control_id and text.upper() == control_id.upper():
        return False
    if text.startswith("De novo structure"):
        return False
    if text.lower() in _GENERIC_LOWER:
        return False
    if not re.search(r"[A-Za-z]", text):
        return False
    if _INCHIKEY_RE.fullmatch(text):
        return False
    if _FORMULA_ONLY_RE.fullmatch(text.replace(" ", "")):
        return False
    if _SMILES_LIKE_RE.fullmatch(text) and any(ch in text for ch in "[]@=#\\/"):
        return False
    return True


def infer_lipid_class_from_generic_label(label: str) -> str:
    """Map a bare generic lipid display label back to a lipid_class code."""

    return _GENERIC_LABEL_LIPID_CLASS.get(_format_label(label), "")


def specificity_rank(descriptor_kind: str | None, label: str) -> int:
    """Return a monotonic specificity score for display-name comparison."""

    kind = _clean(descriptor_kind)
    text = _format_label(label)
    lower = text.lower()
    if kind == "candidate_name":
        return 100
    if kind == "external_db_exact":
        return 100
    if kind == "iupac_generated":
        return 80
    if kind == "external_db_skeleton":
        return 70
    if kind == "conjugate_resolved":
        return 60
    if kind == "known_x_name":
        return 55
    if kind == "peptide" and "composition" not in lower and ("–" in text or "-" in text):
        return 55
    if kind == "lipid" and "(" in text:
        return 50
    if kind == "conjugate_unresolved" and " of " in lower:
        return 50
    if text in GENERIC_CLASS_LABELS or lower.startswith("de novo structure"):
        return 10
    if kind == "peptide":
        return 40
    if kind == "lipid":
        return 40
    if kind == "analog":
        return 30
    if kind == "chemical_class":
        return 25
    if kind == "motif_only":
        return 20
    if text:
        return 85
    return 0


_TENTATIVE_SUFFIX_RE = re.compile(r"\s*[\u2014-]\s*tentative\s*$", re.IGNORECASE)
_GENERATED_TENTATIVE_RE = re.compile(r"\(generated,\s*tentative\)", re.IGNORECASE)


def normalize_descriptor_label(label: str) -> str:
    """Strip legacy inline tentative suffixes from a display label."""

    text = _format_label(label)
    if not text:
        return text
    text = _GENERATED_TENTATIVE_RE.sub("(generated)", text)
    text = _TENTATIVE_SUFFIX_RE.sub("", text).strip()
    return _format_label(text)


def _read_known_x_row(bundle_dir: Path | str | None) -> dict[str, Any]:
    if bundle_dir is None:
        return {}
    path = Path(bundle_dir)
    slug = path.name.lower()
    known_x_path = path / "data" / f"{slug}_known_x_row.csv"
    if not known_x_path.is_file():
        return {}
    try:
        with known_x_path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
    except (OSError, csv.Error):
        return {}
    if not rows:
        return {}
    return dict(rows[0])


def _read_known_x_name(bundle_dir: Path | str | None) -> str:
    row = _read_known_x_row(bundle_dir)
    if not row:
        return ""
    level = _clean(row.get("Level") or row.get("level"))
    if level != "1":
        return ""
    for key in (
        "External Information (Class or modification)",
        "Tentative information",
        "Cleaned External Info",
        "KNOWN",
    ):
        raw = _clean(row.get(key))
        if raw.upper().startswith("KNOWN "):
            cleaned = _clean_known_x_name(raw)
            if cleaned:
                return cleaned
    return ""


def _smiles_conjugate_moieties(smiles: str) -> set[str]:
    text = _clean(smiles)
    if not text:
        return set()
    found: set[str] = set()
    if _SULFATE_SMILES_RE.search(text):
        found.add("sulfate")
    if _GLUCURONIDE_SMILES_RE.search(text):
        found.add("glucuronide")
    return found


def _text_conjugate_term(*texts: str) -> str:
    combined = " ".join(_clean(text).lower() for text in texts if _clean(text))
    if not combined:
        return ""
    for stem in _TEXT_CONJUGATE_TRIGGER_ORDER:
        if stem in combined:
            return TEXT_CONJUGATE_TRIGGERS[stem]
    return ""


def _clean_curator_name(text: str) -> str:
    clause = _clean(text)
    clause = re.split(r"\s*;\s*", clause, maxsplit=1)[0]
    clause = re.split(r"\s+or\s+", clause, maxsplit=1, flags=re.I)[0]
    clause = re.split(
        r",\s*(?:Possible|possible|likely|probably|isomer)\b",
        clause,
        maxsplit=1,
        flags=re.I,
    )[0]
    clause = _HEDGE_PREFIX_RE.sub("", clause).rstrip("?").strip()
    return _clean_candidate_display_name(clause)


def _is_steroid_class_name(name: str) -> bool:
    text = _clean(name).lower()
    conjugation_patterns = (
        r"glucuron\w*",
        r"sulfat\w*",
        r"sulph\w*",
        r"taurin\w*",
        r"glycin\w*",
        r"carnitin\w*",
        r"lactoyl\w*",
    )
    for pattern in conjugation_patterns:
        text = re.sub(pattern, " ", text)
    text = re.sub(r"[-–]+", " ", text)
    tokens = [token for token in re.split(r"\s+", text) if token]
    if not tokens:
        return True
    return all(token in _STEROID_CLASS_TOKENS for token in tokens)


def _try_known_x_conjugate(
    bundle_dir: Path | str | None,
    smiles: str,
    motifs: tuple[str, ...],
    formula: str,
) -> TentativeName | None:
    row = _read_known_x_row(bundle_dir)
    if not row:
        return None

    external = _clean(row.get("External Information (Class or modification)"))
    cleaned_external = _clean(row.get("Cleaned External Info"))
    tentative = _clean(row.get("Tentative information"))
    text_term = _text_conjugate_term(external, cleaned_external, tentative)
    if not text_term:
        return None

    smiles_terms = set(motifs) & set(TEXT_CONJUGATE_TRIGGERS.values())
    if not smiles_terms:
        smiles_terms = _smiles_conjugate_moieties(smiles)
    if smiles_terms and text_term not in smiles_terms:
        return None

    name_source = cleaned_external or tentative
    cleaned_name = _clean_curator_name(name_source)
    level = _clean(row.get("Level") or row.get("level"))
    class_term = _STEROID_CLASS_TERMS.get(text_term, text_term)
    is_class = (
        not cleaned_name
        or _is_steroid_class_name(cleaned_name)
        or not _looks_like_chemical_name(cleaned_name)
        or _candidate_specificity(cleaned_name) == 0
        or _is_descriptive_phrase(cleaned_name)
    )
    if is_class:
        return TentativeName(
            kind="chemical_class",
            label=_format_label(f"Steroid {class_term}"),
            motifs=(text_term,),
            confidence="low",
            source="known_x_row",
        )

    label = _format_label(cleaned_name)
    confidence = "medium" if level in {"1", "2"} else "low"
    return TentativeName(
        kind="known_x_name",
        label=label,
        parent_name=cleaned_name,
        motifs=(text_term,),
        confidence=confidence,
        source="known_x_row",
    )


def _format_candidate_label(name: str, *, msi_level: str, truth_match: bool) -> str:
    """Hedge a borrowed reference name the evidence does not support asserting.

    A rank-1 ``candidate_name`` or a known_X row supplies a real compound name,
    which reads as an assertion. At Level 3 the structure is a hypothesis, so the
    name has to say so; a control whose top candidate matches its confirmed
    structure is exempt, because there the name is not a hypothesis.
    """

    label = _format_label(name)
    if truth_match:
        return label
    level = msi_int(msi_level)
    return hedged_label(label) if level is not None and level >= 3 else label


def _try_candidate_name(
    rank1_row: Mapping[str, Any] | Any,
    *,
    bundle_dir: Path | str | None,
    control_id: str,
    msi_level: str = "",
    truth_match: bool = False,
    html_alt_name: str = "",
) -> TentativeName | None:
    row = _row_dict(rank1_row)
    raw_candidates: list[tuple[str, str]] = []
    rank1_name = _clean_candidate_display_name(_clean(row.get("candidate_name")))
    if rank1_name:
        raw_candidates.append(("rank1_candidate_name", rank1_name))
    known_x_name = _read_known_x_name(bundle_dir)
    if known_x_name:
        raw_candidates.append(("known_x_row", known_x_name))
    alt_name = _clean_candidate_display_name(html_alt_name)
    if alt_name:
        raw_candidates.append(("html_alt", alt_name))

    ranked = [
        (source, name, _candidate_specificity(name))
        for source, name in raw_candidates
        if _looks_like_chemical_name(name, control_id=control_id) and _candidate_specificity(name) > 0
    ]
    if not ranked:
        return None
    source, name, _ = max(
        ranked,
        key=lambda item: (item[2], -_SOURCE_PRIORITY[item[0]]),
    )
    label = _format_candidate_label(name, msi_level=msi_level, truth_match=truth_match)
    return TentativeName(
        kind="candidate_name",
        label=label,
        parent_name=name,
        confidence="high",
        source=source,
    )


def _row_dict(rank1_row: Any) -> dict[str, Any]:
    if isinstance(rank1_row, Mapping):
        return dict(rank1_row)
    to_dict = getattr(rank1_row, "to_dict", None)
    if callable(to_dict):
        value = to_dict()
        return dict(value) if isinstance(value, Mapping) else {}
    return {}


def _peptide_tokens(value: Any, *, ordered: bool) -> list[str]:
    text = _clean(value)
    if not text or text.lower() == "order uncertain":
        return []
    try:
        from structural_elucidation.pipeline.enrichment import peptide_composition as peptide

        normalized = (
            peptide.normalize_ordered_sequence(text)
            if ordered
            else peptide.normalize_peptide_composition(text)
        )
        return [token for token in normalized.split("-") if token]
    except Exception:
        residue_pattern = (
            r"Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|"
            r"Pro|Pyr|Ser|Thr|Trp|Tyr|Val|Xle"
        )
        return [
            ("Xle" if token in {"Ile", "Leu"} else token)
            for token in re.findall(residue_pattern, text, flags=re.I)
        ]


def _try_peptide(smiles: str, formula: str, row: dict[str, Any]) -> TentativeName | None:
    parsed: dict[str, Any] = {}
    if HAS_RDKIT and smiles:
        try:
            from structural_elucidation.pipeline.enrichment.peptide_composition import (
                infer_sequence_from_smiles,
            )

            candidate = infer_sequence_from_smiles(smiles)
            parsed = candidate if isinstance(candidate, dict) else {}
        except Exception:
            parsed = {}

    parsed_ok = _clean(parsed.get("sequence_parse_status")).lower() == "parsed"
    row_parsed = _clean(row.get("sequence_parse_status")).lower() == "parsed"
    row_candidate = _boolish(row.get("is_peptide_candidate"))
    if not (parsed_ok or row_parsed or row_candidate):
        return None

    ordered: list[str] = []
    for value in (
        parsed.get("ordered_peptide_sequence"),
        row.get("peptide_ordering_label"),
        row.get("ordered_peptide_sequence"),
        row.get("candidate_ordered_sequence"),
        row.get("actual_ordered_sequence"),
    ):
        ordered = _peptide_tokens(value, ordered=True)
        if ordered:
            break

    composition: list[str] = []
    for value in (
        row.get("normalized_peptide_composition"),
        row.get("peptide_composition"),
        row.get("candidate_sequence_composition"),
        row.get("peptide_composition_label"),
        parsed.get("ordered_peptide_sequence"),
    ):
        composition = _peptide_tokens(value, ordered=False)
        if composition:
            break
    if not composition and ordered:
        composition = sorted(ordered)
    declared_match = _clean(row.get("sequence_matches_declared_composition"))
    if ordered and composition and (
        sorted(ordered) != sorted(composition)
        or (declared_match and not _boolish(declared_match))
    ):
        ordered = []

    length = int(_num(row.get("peptide_length")) or len(ordered) or len(composition) or 0)
    if length <= 0:
        return None
    peptide_type = _PEPTIDE_LENGTHS.get(length, f"{length}-residue peptide")
    penalty = _num(row.get("peptide_branch_or_cycle_penalty")) or 0.0
    cyclic = penalty >= 0.5 or _boolish(row.get("is_cyclic_peptide"))

    if ordered and length <= 4:
        sequence = "\N{EN DASH}".join(ordered)
        base = (
            f"Cyclic {peptide_type.lower()} ({sequence})"
            if cyclic
            else f"{peptide_type} ({sequence})"
        )
    elif composition or ordered:
        tokens = composition or ordered
        if length <= 5:
            detail = "/".join(tokens)
        else:
            unique_tokens = list(dict.fromkeys(tokens))
            detail = "/".join(unique_tokens) if len(unique_tokens) <= 5 else "mixed amino-acid"
        if length == 1:
            base = f"Amino acid ({detail})"
        else:
            prefix = f"Cyclic {peptide_type.lower()}" if cyclic else peptide_type
            base = f"{prefix} ({detail} composition)"
    else:
        base = f"Cyclic {peptide_type.lower()}" if cyclic else peptide_type

    connectivity = max(
        [
            value
            for value in (
                _num(row.get("peptide_connectivity_confidence")),
                _num(row.get("peptide_ordering_confidence")),
                _num(row.get("peptide_composition_confidence")),
                _num(row.get("composition_confidence")),
            )
            if value is not None
        ]
        or [0.0]
    )
    if connectivity >= 0.8 and penalty < 0.25:
        confidence = "high"
    elif connectivity >= 0.5 or ordered or composition:
        confidence = "medium"
    else:
        confidence = "low"
    if penalty >= 0.5:
        confidence = "low" if confidence == "medium" else "medium"
    source = "smiles_rdkit" if parsed_ok else "rank1_csv"
    return TentativeName(
        kind="peptide",
        label=_format_label(base),
        motifs=(),
        confidence=confidence,
        source=source,
    )


def _try_lipid(smiles: str, formula: str, row: dict[str, Any]) -> TentativeName | None:
    lipid_class = _clean(row.get("lipid_class"))
    shorthand = _clean(row.get("lipid_shorthand"))
    source = "rank1_csv" if lipid_class else ""
    if HAS_RDKIT and smiles:
        try:
            from structural_elucidation.pipeline.enrichment.lipid_class import (
                annotate_from_smiles,
            )

            annotation = annotate_from_smiles(smiles)
            if not annotation.is_empty:
                lipid_class = _clean(annotation.lipid_class)
                shorthand = _clean(annotation.lipid_shorthand)
                source = "smiles_rdkit"
        except Exception:
            pass
    if not lipid_class:
        return None

    class_name = _LIPID_CLASS_NAMES.get(lipid_class, lipid_class)
    detail = ""
    if shorthand:
        match = re.match(
            rf"^(?:{re.escape(lipid_class)}|{re.escape(class_name)})\s*(?:\(([^)]+)\)|\s+(.+))$",
            shorthand,
            flags=re.I,
        )
        if match:
            detail = _clean(match.group(1) or match.group(2))
    if lipid_class.lower() == "sterol":
        base = "Sterol-like compound"
    elif detail:
        base = f"{class_name} ({detail})"
    else:
        base = class_name
    return TentativeName(
        kind="lipid",
        label=_format_label(base),
        motifs=(),
        confidence="high" if bool(detail) else "medium",
        source=source or "rank1_csv",
    )


def _motif_labels(smiles: str) -> tuple[str, ...]:
    if not HAS_RDKIT or not smiles:
        return ()
    try:
        from structural_elucidation.nce.evidence import motif_labels_for_smiles

        return tuple(motif_labels_for_smiles(smiles))
    except Exception:
        return ()


def _reaction_resources() -> tuple[dict[str, tuple[str, ...]], dict[str, str]]:
    try:
        from structural_elucidation.chemistry.conjugates import CONJUGATE_REACTIONS
        from structural_elucidation.nce.evidence import TARGETED_REACTION_SMARTS

        return dict(CONJUGATE_REACTIONS), dict(TARGETED_REACTION_SMARTS)
    except Exception:
        return {}, {}


def _mol_formula(smiles: str) -> str:
    if not HAS_RDKIT or Chem is None or rdMolDescriptors is None:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        return str(rdMolDescriptors.CalcMolFormula(mol))
    except Exception:
        return ""


def _mol_inchikey(smiles: str) -> str:
    if not HAS_RDKIT or Chem is None:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        return str(Chem.MolToInchiKey(mol))
    except Exception:
        return ""


def _reverse_products(smiles: str, motif: str) -> list[str]:
    """Reverse existing targeted-conjugation reactions for one motif."""

    if not HAS_RDKIT or Chem is None or AllChem is None:
        return []
    reaction_map, smarts_map = _reaction_resources()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    input_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    input_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    products: list[str] = []
    seen: set[str] = set()
    for reaction_name in reaction_map.get(motif, ()):
        smarts = smarts_map.get(reaction_name, "")
        if ">>" not in smarts:
            continue
        reactant, product = smarts.split(">>", 1)
        try:
            if rdBase is not None and hasattr(rdBase, "BlockLogs"):
                with rdBase.BlockLogs():
                    reaction = AllChem.ReactionFromSmarts(f"{product}>>{reactant}")
                    product_sets = reaction.RunReactants((mol,))
            else:  # pragma: no cover - older RDKit
                reaction = AllChem.ReactionFromSmarts(f"{product}>>{reactant}")
                product_sets = reaction.RunReactants((mol,))
        except Exception:
            continue
        for product_set in product_sets:
            for candidate in product_set:
                try:
                    Chem.SanitizeMol(candidate)
                    parent_smiles = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=True)
                except Exception:
                    continue
                if (
                    not parent_smiles
                    or "." in parent_smiles
                    or parent_smiles == input_smiles
                    or parent_smiles in seen
                ):
                    continue
                parent_carbons = sum(
                    1 for atom in candidate.GetAtoms() if atom.GetAtomicNum() == 6
                )
                if parent_carbons <= 0 or parent_carbons > input_carbons:
                    continue
                if (
                    motif not in _CARBON_FREE_CONJUGATES
                    and parent_carbons >= input_carbons
                ):
                    continue
                seen.add(parent_smiles)
                products.append(parent_smiles)
    return products


def _conjugate_group_count(smiles: str, motif: str) -> int:
    """Count sulfur/phosphorus centers used to rank deconjugated products."""

    if (
        motif not in _CARBON_FREE_CONJUGATES
        or not HAS_RDKIT
        or Chem is None
    ):
        return 0
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    center_atomic_number = 16 if motif == "sulfate" else 15
    return sum(
        atom.GetAtomicNum() == center_atomic_number
        for atom in mol.GetAtoms()
    )


def _fully_reverse_carbon_free_conjugate(
    smiles: str,
    motif: str,
) -> tuple[list[str], int]:
    """Return the least-conjugated sulfate/phosphate parents and group count."""

    initial_count = _conjugate_group_count(smiles, motif)
    seen = {smiles}
    frontier = {smiles}
    for _ in range(max(initial_count, 1) + 1):
        next_frontier = {
            product
            for candidate in frontier
            for product in _reverse_products(candidate, motif)
        }
        next_frontier.difference_update(seen)
        if not next_frontier:
            break
        seen.update(next_frontier)
        frontier = next_frontier
    products = seen - {smiles}
    if not products:
        return [], 0
    remaining = {
        product: _conjugate_group_count(product, motif)
        for product in products
    }
    minimum = min(remaining.values())
    least_conjugated = [
        product for product, count in remaining.items() if count == minimum
    ]
    carbon_counts = {
        product: sum(
            atom.GetAtomicNum() == 6
            for atom in Chem.MolFromSmiles(product).GetAtoms()
        )
        for product in least_conjugated
    }
    maximum_carbons = max(carbon_counts.values())
    parents = sorted(
        product
        for product in least_conjugated
        if carbon_counts[product] == maximum_carbons
    )
    removed_count = max(1, initial_count - minimum)
    return parents, removed_count


def _db_name(db_names: Mapping[str, Any] | Any, inchikey: str) -> str | None:
    if not inchikey or db_names is None:
        return None
    getter = getattr(db_names, "get", None)
    if not callable(getter):
        return None
    record = getter(inchikey.split("-")[0])
    if isinstance(record, Mapping):
        return _clean(record.get("name") or record.get("canonical_name")) or None
    return _clean(record) or None


def _multiplicity_term(term: str, count: int) -> str:
    if count <= 1:
        return term
    prefix = _MULTIPLICITY_PREFIXES.get(count)
    return f"{prefix}{term}" if prefix else f"{count}-{term}"


def _resolved_conjugate_base(
    motifs: tuple[str, ...],
    parent: str,
    conjugate_count: int = 1,
) -> str:
    terms = [_CONJUGATE_TERMS.get(motif, motif.replace("_", " ")) for motif in motifs]
    if len(terms) == 1:
        term = terms[0]
        if term in {"glucuronide", "sulfate", "phosphate"}:
            if conjugate_count == 1 and parent.casefold().endswith(f" {term}"):
                stem = parent[: -(len(term) + 1)].rstrip()
                return f"{stem} di{term}"
            return f"{parent} {_multiplicity_term(term, conjugate_count)}"
        return f"{parent}\N{EN DASH}{term} conjugate"
    return f"{parent} {'/'.join(terms)} conjugate"


def _unresolved_conjugate_base(
    motifs: tuple[str, ...],
    parent_formula: str,
    conjugate_count: int = 1,
) -> str:
    terms = [_CONJUGATE_TERMS.get(motif, motif.replace("_", " ")) for motif in motifs]
    if len(terms) == 1:
        terms[0] = _multiplicity_term(terms[0], conjugate_count)
    descriptor = "/".join(terms)
    return (
        f"{descriptor[:1].upper()}{descriptor[1:]} conjugate "
        f"of {parent_formula}"
    )


def _conjugate_hint_term(row: Mapping[str, Any]) -> str:
    for key in ("conjugate_name", "transform_name"):
        raw = _clean(row.get(key)).lower()
        if not raw:
            continue
        if "glucuron" in raw:
            return "glucuronide"
        for term in DECORATING_CONJUGATES:
            if term in raw:
                return term
    return ""


def _try_conjugate_hint(row: Mapping[str, Any], formula: str) -> TentativeName | None:
    term = _conjugate_hint_term(row)
    if not term:
        return None
    parent_formula = _clean(row.get("core_formula"))
    if parent_formula:
        base = _unresolved_conjugate_base((term,), parent_formula, 1)
        return TentativeName(
            kind="conjugate_unresolved",
            label=_sentence_case_label(base),
            parent_formula=parent_formula,
            conjugate_count=1,
            motifs=(term,),
            confidence="medium",
            source="rank1_transform_hint",
        )
    phrase = _MOTIF_PHRASES.get(term)
    if phrase:
        base = f"{phrase[:1].upper()}{phrase[1:]} compound"
        return TentativeName(
            kind="motif_only",
            label=_format_label(base),
            motifs=(term,),
            confidence="low",
            source="rank1_transform_hint",
        )
    return None


def _conjugate_is_decorating(conj: TentativeName | None) -> bool:
    return (
        conj is not None
        and conj.kind in {"conjugate_resolved", "conjugate_unresolved"}
        and bool(conj.motifs)
        and set(conj.motifs) <= DECORATING_CONJUGATES
    )


def _try_conjugate(
    smiles: str,
    formula: str,
    motifs: tuple[str, ...],
    db_names: Mapping[str, Any] | Any,
) -> TentativeName | None:
    if not HAS_RDKIT or not smiles:
        return None
    reaction_map, _ = _reaction_resources()
    conjugate_motifs = tuple(motif for motif in motifs if motif in reaction_map)
    if not conjugate_motifs:
        return None
    carbon_bearing_motifs = tuple(
        motif
        for motif in conjugate_motifs
        if motif not in _CARBON_FREE_CONJUGATES
    )
    carbon_free_motifs = tuple(
        motif
        for motif in conjugate_motifs
        if motif in _CARBON_FREE_CONJUGATES
    )
    conjugate_motifs = carbon_bearing_motifs + carbon_free_motifs

    candidates = [smiles]
    removed: list[str] = []
    conjugate_count = 0
    for motif in conjugate_motifs:
        if (
            motif in _CARBON_FREE_CONJUGATES
            and set(removed) - _CARBON_FREE_CONJUGATES
        ):
            continue
        next_candidates: list[str] = []
        removed_counts: list[int] = []
        for candidate in candidates:
            if motif in _CARBON_FREE_CONJUGATES:
                products, count = _fully_reverse_carbon_free_conjugate(
                    candidate,
                    motif,
                )
                next_candidates.extend(products)
                if products:
                    removed_counts.append(count)
            else:
                products = _reverse_products(candidate, motif)
                next_candidates.extend(products)
                if products:
                    removed_counts.append(1)
        next_candidates = list(dict.fromkeys(next_candidates))[:64]
        if next_candidates:
            candidates = next_candidates
            removed.append(motif)
            conjugate_count += max(removed_counts, default=1)
    if not removed:
        return None

    resolved_options: list[tuple[str, str, str]] = []
    unresolved_options: list[tuple[str, str]] = []
    for parent_smiles in candidates:
        parent_ik = _mol_inchikey(parent_smiles)
        parent_formula = _mol_formula(parent_smiles)
        parent_name = _db_name(db_names, parent_ik)
        if parent_name:
            resolved_options.append((parent_name, parent_ik, parent_formula))
        elif parent_formula:
            unresolved_options.append((parent_formula, parent_ik))

    removed_tuple = tuple(removed)
    carbon_free_resolution = bool(
        set(removed_tuple) & _CARBON_FREE_CONJUGATES
    )
    resolved = (
        sorted(
            resolved_options,
            key=lambda item: (item[0].casefold(), item[1]),
        )
        if carbon_free_resolution
        else resolved_options
    )
    if resolved:
        parent_name, parent_ik, parent_formula = resolved[0]
        return TentativeName(
            kind="conjugate_resolved",
            label=_sentence_case_label(
                _resolved_conjugate_base(
                    removed_tuple,
                    parent_name,
                    conjugate_count,
                )
            ),
            parent_name=parent_name,
            parent_inchikey=parent_ik or None,
            parent_formula=parent_formula or None,
            conjugate_count=conjugate_count,
            motifs=removed_tuple,
            confidence="high",
            source="smiles_rdkit+db",
        )
    parent_formulas = (
        {
            parent_formula
            for parent_formula, _parent_ik in unresolved_options
            if parent_formula
        }
        if carbon_free_resolution
        else {
            unresolved_options[0][0]
        }
        if unresolved_options
        else set()
    )
    if len(parent_formulas) == 1:
        parent_formula = next(iter(parent_formulas))
        if carbon_free_resolution:
            parent_inchikeys = {
                parent_ik
                for candidate_formula, parent_ik in unresolved_options
                if candidate_formula == parent_formula and parent_ik
            }
            parent_ik = (
                next(iter(parent_inchikeys))
                if len(parent_inchikeys) == 1
                else ""
            )
        else:
            parent_ik = unresolved_options[0][1]
        return TentativeName(
            kind="conjugate_unresolved",
            label=_format_label(
                _unresolved_conjugate_base(
                    removed_tuple,
                    parent_formula,
                    conjugate_count,
                )
            ),
            parent_inchikey=parent_ik or None,
            parent_formula=parent_formula or None,
            conjugate_count=conjugate_count,
            motifs=removed_tuple,
            confidence="medium",
            source="smiles_rdkit",
        )
    return None


def _try_motif(formula: str, motifs: tuple[str, ...]) -> TentativeName | None:
    specific_set = {motif for motif in motifs if motif in _SPECIFIC_MOTIFS}
    if {"sulfate", "sulfonate"} & specific_set:
        specific_set.discard("sulfur_oxygen_rich")
    if "phosphate" in specific_set:
        specific_set.discard("phosphate_like")
    specific = tuple(motif for motif in _MOTIF_PRIORITY if motif in specific_set)
    if not specific:
        return None
    phrases = list(dict.fromkeys(_MOTIF_PHRASES[motif] for motif in specific if motif in _MOTIF_PHRASES))
    if not phrases:
        return None
    visible = phrases[:2]
    phrase = " and ".join(visible)
    base = f"{phrase[:1].upper()}{phrase[1:]} compound"
    return TentativeName(
        kind="motif_only",
        label=_format_label(base),
        motifs=specific,
        confidence="low",
        source="smiles_rdkit",
    )


def resolve(
    smiles: str,
    formula: str,
    rank1_row: Mapping[str, Any] | Any,
    db_names: Mapping[str, Any] | Any,
    *,
    name_resolved: bool,
    bundle_dir: Path | str | None = None,
    extended_ctx: Any | None = None,
    inchikey: str = "",
    control_id: str = "",
    msi_level: str = "",
    truth_match: bool = False,
    html_alt_name: str = "",
) -> TentativeName | None:
    """Resolve the first matching tentative descriptor branch.

    Named compounds are an explicit no-op.  ``None`` for an unnamed compound
    means the caller should retain its existing ``De novo structure`` fallback.
    Layer 0 prefers rank-1 ``candidate_name``, known_X Level-1 KNOWN rows, and
    HTML structure alt text before peptide/lipid/conjugate/motif branches.
    """

    if name_resolved:
        return None
    row = _row_dict(rank1_row)
    smiles_text = _clean(smiles)
    formula_text = _clean(formula) or _clean(row.get("predicted_formula"))

    candidate = _try_candidate_name(
        row,
        bundle_dir=bundle_dir,
        control_id=control_id,
        msi_level=msi_level,
        truth_match=truth_match,
        html_alt_name=html_alt_name,
    )
    if candidate:
        return candidate

    peptide = _try_peptide(smiles_text, formula_text, row)
    if peptide:
        return peptide
    lipid = _try_lipid(smiles_text, formula_text, row)
    if lipid is not None and lipid.label not in GENERIC_CLASS_LABELS:
        return lipid
    if lipid is not None:
        motifs = _motif_labels(smiles_text)
        conjugate = _try_conjugate(smiles_text, formula_text, motifs, db_names)
        hint = _try_conjugate_hint(row, formula_text)
        known_x_conj = _try_known_x_conjugate(
            bundle_dir,
            smiles_text,
            motifs,
            formula_text,
        )
        contenders = [lipid]
        if _conjugate_is_decorating(conjugate):
            contenders.append(conjugate)
        if hint is not None:
            contenders.append(hint)
        if known_x_conj is not None:
            contenders.append(known_x_conj)
        upgrades = [item for item in contenders if item.kind != "lipid"]
        if upgrades:
            return max(upgrades, key=lambda item: specificity_rank(item.kind, item.label))
        return lipid
    motifs = _motif_labels(smiles_text)
    conjugate = _try_conjugate(smiles_text, formula_text, motifs, db_names)
    if conjugate:
        return conjugate
    motif = _try_motif(formula_text, motifs)
    if motif:
        return motif
    known_x_conj = _try_known_x_conjugate(
        bundle_dir,
        smiles_text,
        motifs,
        formula_text,
    )
    if known_x_conj is not None:
        return known_x_conj
    if bundle_dir is not None and extended_ctx is not None:
        try:
            from structural_elucidation.reports import tentative_naming_ext as ext
        except ImportError:  # pragma: no cover
            import tentative_naming_ext as ext  # type: ignore

        inchikey_text = _clean(inchikey) or ext.inchikey_from_smiles(smiles_text)
        return ext.apply_extended_layers(
            smiles_text,
            formula_text,
            inchikey_text,
            Path(bundle_dir),
            extended_ctx,
        )
    return None


def load_db_names(db_path: Path | str | None) -> dict[str, dict[str, Any]]:
    """Load canonical names keyed by the first InChIKey block."""

    path = Path(db_path) if db_path else None
    if path is None or not path.exists():
        return {}
    out: dict[str, dict[str, Any]] = {}
    with sqlite3.connect(str(path)) as connection:
        rows = connection.execute(
            "SELECT inchikeys, canonical_name, short_names, chemical_names FROM metabolites"
        )
        for inchikeys, canonical_name, short_names, chemical_names in rows:
            name = _clean(canonical_name)
            if not name:
                continue
            synonyms: list[str] = []
            for source in (short_names, chemical_names):
                synonyms.extend(
                    token.strip()
                    for token in str(source or "").replace(";", ",").split(",")
                    if token.strip()
                )
            for inchikey in str(inchikeys or "").replace(";", ",").split(","):
                skeleton = inchikey.strip().split("-")[0]
                if skeleton and skeleton not in out:
                    out[skeleton] = {"name": name, "synonyms": synonyms[:8]}
    return out


class SqliteNameLookup:
    """Small mapping-like SQLite lookup used by single-report fresh builds."""

    def __init__(self, db_path: Path | str | None):
        self.db_path = Path(db_path) if db_path else Path()
        self.cache: dict[str, dict[str, str] | None] = {}

    def get(self, skeleton: str, default: Any = None) -> Any:
        key = _clean(skeleton).split("-")[0]
        if not key or not self.db_path.is_file():
            return default
        if key not in self.cache:
            with sqlite3.connect(str(self.db_path)) as connection:
                row = connection.execute(
                    "SELECT canonical_name FROM metabolites "
                    "WHERE inchikeys LIKE ? AND canonical_name IS NOT NULL "
                    "AND trim(canonical_name) <> '' LIMIT 1",
                    (f"%{key}%",),
                ).fetchone()
            self.cache[key] = {"name": str(row[0]).strip()} if row else None
        return self.cache[key] if self.cache[key] is not None else default
