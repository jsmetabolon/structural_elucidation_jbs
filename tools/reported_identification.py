"""Publication name and level for one metabolite, from the white-paper rules.

The result is one string, ``reported_identification``, written on every candidate
row of the executive-summary sheet. It is a reporting call, separate from the
curated Metabolon level and from the pipeline MSI level.

Rules, first match wins:

1. Authentic standard (Metabolon 1, MSI 1, or a standard-confirmed exception)
   → ``Level 1: {name}``.
2. A named-library entry that is not one structure, such as a glucuronide of a
   formula → ``Level 2*: {name}*``.
3. An ``isomer of`` label (exact match to a library entry besides retention
   time) → ``Level 3: isomer of {library name}``. The library name is enough.
4. A specific structure name → ``Level 2: [X-######], Probable {name}``.
   A name is enough; the candidates do not have to agree.
5. A class phrase such as ``Glucoside`` → ``Level 3``.
6. No supported formula → ``Level 5: exact mass of interest``.
7. Among the formula-supporting candidates, a structurally related family
   that clears every unrelated candidate by a distinctive Match Score gap,
   with no name found another way → ``Level 3: Probable isomer of {top
   candidate}``.
8. A conjugate, sugar, or peptide shared by a strict majority of the
   formula-supporting candidates → ``Level 3: {class}``.
9. A supported formula and no name or shared class → ``Level 4: {formula}``.

``candidate_identity_confidence`` is a separate, 0-100 number answering a different
question than Match Score. Match Score (``evidence_confidence`` /
``predicted_spectrum_score_v2``) says how well one candidate's predicted
spectrum fits the experimental data, and the Match Score margin between
structural families is frequently uninformative on its own: many unrelated
candidates score near the overlap-score fallback tie, and many others all
saturate near 1.0, so the raw margin alone does not separate a confident
"Probable [name]" call from an unsupported molecular formula.

So the scale is banded by what the publication call itself asserts, and the
candidate margin only places a compound within its band. The margin uses the
pure NCE score (``pure_nce_score``) when a compound has one, since it does not
saturate. When ``identity_confidence_calibration.json`` (written by
``tools/calibrate_identity_confidence.py``) sits beside this module, its
anchors and bands, fitted on the truth controls, replace the hand-set
``_CONFIDENCE_BANDS``. Level 1 is always 100. Level 2* is floored at
``LEVEL_2_STAR_FLOOR`` (99) and a higher margin is kept, since the confidence
there is the manual review, not the spectrum.
"""

from __future__ import annotations

import contextlib
import io
import json
import math
import re
import warnings
from pathlib import Path
from typing import Any

try:
    # RDKit prints a NumPy traceback on some environments before recovering.
    with contextlib.redirect_stderr(io.StringIO()):
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")
    HAS_RDKIT = True
except Exception:  # pragma: no cover - optional dependency
    Chem = None
    HAS_RDKIT = False
try:
    with contextlib.redirect_stderr(io.StringIO()):
        from rdkit import DataStructs
        from rdkit.Chem import RDKFingerprint
    HAS_FINGERPRINT = True
except Exception:  # pragma: no cover - optional dependency
    DataStructs = None
    RDKFingerprint = None
    HAS_FINGERPRINT = False


# A family must clear unrelated candidates by this much evidence to count as
# distinctively higher. Tied, saturated scores do not.
DISTINCTIVE_GAP = 0.05
# Curator Level 2* names are confident by review. 100 stays reserved for an
# authentic standard (Level 1). A higher cluster score is kept.
LEVEL_2_STAR_FLOOR = 99.0
# Mass-match floor for the candidates allowed to vote on the formula.
MASS_CONSISTENT = 0.95
RELATED_TANIMOTO = 0.60

# candidate_identity_confidence bands, by the publication call that was made
# (see _confidence_tier). The 0-100 Match Score cluster margin only places a
# compound within its band, so Level 2/3 calls stay well above Level 4/5 even
# when the margin itself is uninformative. Level 1 (100) and Level 2* (floor
# 99) are handled separately, not through this table.
_CONFIDENCE_BANDS: dict[str, tuple[float, float]] = {
    "2_probable": (90.0, 99.0),
    "3_isomer_exact": (75.0, 92.0),
    "3_class": (55.0, 80.0),
    "3_probable_isomer_group": (35.0, 65.0),
    "4": (15.0, 45.0),
    "5": (0.0, 20.0),
}

MSI_LABELS = {
    "1": "Confirmed (authentic standard)",
    "2": "Probable structure",
    "3": "Tentative candidate structure",
    "4": "Molecular formula",
    "5": "Exact mass of interest",
}

_HEDGE_PREFIX = re.compile(
    r"^(?:known|released|possible|possibly|probable|probably|putative|likely)\s+",
    re.I,
)
_ISOMER_PREFIX = re.compile(r"^isomer of\s+", re.I)
_CLASS_OF = re.compile(
    r"^(?:sulfates?|sulphates?|sulfated|glucuronides?|peptides?|ribonucleosides?)\s+of\b",
    re.I,
)
_GENERIC_CLASS = re.compile(
    r"^(?:[\w][\w\-/]*\s+){0,4}"
    r"(?:amino acid derivatives?|bile acids?|ribonucleosides?|derivatives?|"
    r"glucosides?|indole compounds?)$",
    re.I,
)
_PEPTIDE_COMPOSITION = re.compile(
    r"\b(?:di|tri|tetra|penta)?peptide composed of\b",
    re.I,
)
_FORMULA_PARENT = re.compile(r"\bof\s+C\d", re.I)
_BARE_FORMULA = re.compile(r"^(?:[A-Z][a-z]?\d*)+$")
_NOTE = re.compile(
    r"\b(?:similar to|elutes|may be|version on|overlapping|indicates|spectrum|"
    r"noisy|containing|proprietary|linked|fragmentation|somehow|approx)\b",
    re.I,
)
_VAGUE_CLASS_WORD = re.compile(
    r"^(?:sulfates?|sulphates?|sulfated|glucuronides?|peptides?|derivatives?|unknown)$",
    re.I,
)
_GENERIC_WORDS = frozenset(
    {
        "a",
        "an",
        "the",
        "feature",
        "structure",
        "unknown",
        "tentative",
        "compound",
        "metabolite",
        "name",
        "unresolved",
    }
)
_NOTE_PREFIX = re.compile(r"^(?:connection to|not in\b|see )", re.I)


_LEVEL_CALL = re.compile(r"^Level (\d(?:\*)?):")


def reported_level(call: str) -> str:
    """The level token from a publication call, such as ``2*`` or ``4``."""

    match = _LEVEL_CALL.match(_clean(call))
    return match.group(1) if match else ""


def reported_identification(
    compound_id: str,
    candidates: list[dict[str, Any]],
    fields: dict[str, Any] | None = None,
) -> str:
    """Return the publication level and name for one compound."""

    fields = fields or {}
    if _is_level1(fields):
        formula, _competitive = _supported_formula(candidates)
        name = _identity_name(fields) or formula
        return f"Level 1: {name}" if name else "Level 1: confirmed"
    library_name = _nondefinitive_library_name(fields)
    if library_name:
        star = "" if library_name.endswith("*") else "*"
        return f"Level 2*: {library_name}{star}"
    parent = _isomer_parent(fields)
    if parent:
        return f"Level 3: isomer of {parent}"
    specific = _specific_name(fields)
    if specific:
        return f"Level 2: [{_clean(compound_id)}], Probable {specific}"
    klass = _class_phrase(fields)
    if klass:
        return f"Level 3: {klass}"
    formula, competitive = _supported_formula(candidates)
    if not formula or not competitive:
        return "Level 5: exact mass of interest"
    group_parent = _group_isomer_name(competitive)
    if group_parent:
        return f"Level 3: Probable isomer of {group_parent}"
    motif = _shared_motif(competitive)
    if motif:
        return f"Level 3: {motif}"
    return f"Level 4: {formula}"


def _confidence_tier(identification: str) -> str:
    """Which ``_CONFIDENCE_BANDS`` entry a publication call falls in.

    Mirrors the branches of :func:`reported_identification`, read back off
    its own output text, so the two never disagree about which rule fired.
    """

    text = _clean(identification)
    if text.startswith("Level 1:"):
        return "1"
    if text.startswith("Level 2*:"):
        return "2_star"
    if text.startswith("Level 2:"):
        return "2_probable"
    lowered = text.lower()
    if lowered.startswith("level 3: probable isomer of"):
        return "3_probable_isomer_group"
    if lowered.startswith("level 3: isomer of"):
        return "3_isomer_exact"
    if text.startswith("Level 3:"):
        return "3_class"
    if text.startswith("Level 4:"):
        return "4"
    if text.startswith("Level 5:"):
        return "5"
    return ""


def candidate_identity_confidence(
    candidates: list[dict[str, Any]], fields: dict[str, Any] | None = None
) -> float | None:
    """0-100: confidence that the reported identity, not a rival, is right.

    This is not Match Score. Match Score (``evidence_confidence`` /
    ``predicted_spectrum_score_v2``) is how well one candidate's predicted
    spectrum fits the data, and the Match Score margin between structural
    families is frequently uninformative on its own (see the module
    docstring). So the scale is banded by the publication call itself
    (:func:`reported_identification`), and the margin only places a compound
    within its band:

    - Level 1 (authentic standard): 100, always.
    - Level 2* (curator-named modification): at least
      ``LEVEL_2_STAR_FLOOR`` (99); a higher margin is kept.
    - Level 2 "Probable [name]": 90-99.
    - Level 3 "isomer of" (exact match besides retention time): 75-92.
    - Level 3 class/motif (a determined conjugate or chemical class): 55-80.
    - Level 3 "Probable isomer of" (a related family only, by Match Score):
      35-65.
    - Level 4 (bare molecular formula): 15-45.
    - Level 5 (exact mass only): 0-20.

    Those are the hand-set bands; a valid ``CALIBRATION_PATH`` file replaces
    them (see :func:`load_calibration`).

    Blank when there is neither a Level 2* name nor a candidate to be
    confident about.

    This is unrelated to ``structural_elucidation.reports.confidence``'s own
    ``identity_confidence`` field, which is a curator/auto-rule flag (0-1,
    often ``None``) for "named library entry pending an authentic standard,"
    not a candidate-cluster score. Do not conflate the two.
    """

    fields = fields or {}
    if _is_level1(fields):
        return 100.0
    margin = _cluster_identity_confidence(candidates)
    tier = _confidence_tier(reported_identification("", candidates, fields))
    if tier == "2_star":
        # The confidence here is the manual review, not the spectrum, so a
        # missing Match Score still floors at 99 rather than going blank.
        if margin is None:
            return LEVEL_2_STAR_FLOOR
        return max(margin, LEVEL_2_STAR_FLOOR)
    if margin is None:
        return None
    calibration = active_calibration()
    if calibration is not None:
        return _calibrated_score(calibration, tier, margin)
    band = _CONFIDENCE_BANDS.get(tier)
    if band is None:
        return margin
    floor, ceiling = band
    return round(floor + (ceiling - floor) * (margin / 100.0), 1)


def _calibrated_score(calibration: dict[str, Any], tier: str, margin: float) -> float:
    """Place a margin inside a calibrated tier.

    The tier's ``center`` margin maps to its ``anchor`` (the observed
    correct-call rate). With a fitted ``slope`` the offset is logistic,
    otherwise linear across the half-width; either way it stays in
    ``[floor, ceiling]`` so tiers keep their order.
    """

    spec = calibration["tiers"].get(tier)
    if spec is None:
        return margin
    floor, ceiling, anchor = spec["floor"], spec["ceiling"], spec["anchor"]
    offset = (margin - spec.get("center", 50.0)) / 100.0
    slope = calibration.get("slope")
    if slope:
        p = min(max(anchor / 100.0, 1e-4), 1 - 1e-4)
        logit = math.log(p / (1 - p)) + float(slope) * offset
        value = 100.0 / (1.0 + math.exp(-logit))
    else:
        half = max(anchor - floor, ceiling - anchor)
        value = anchor + 2.0 * half * offset
    return round(min(max(value, floor), ceiling), 1)


CALIBRATION_PATH = Path(__file__).with_name("identity_confidence_calibration.json")
# Highest to lowest; calibrated anchors must not increase along this order.
CALIBRATED_TIER_ORDER = (
    "2_probable",
    "3_isomer_exact",
    "3_class",
    "3_probable_isomer_group",
    "4",
    "5",
)
_CALIBRATION_CACHE: dict[str, Any] = {}


def load_calibration(path: Path | str) -> dict[str, Any]:
    """Read and validate a calibration file; raise ``ValueError`` if it is unusable.

    Every tier in ``CALIBRATED_TIER_ORDER`` needs ``floor <= anchor <=
    ceiling`` inside ``[0, LEVEL_2_STAR_FLOOR]``, and anchors must not
    increase down the order, so no lower call can outscore a higher one at
    its anchor and nothing reaches the Level 2* floor or Level 1.
    """

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    tiers = payload.get("tiers") or {}
    previous = None
    for tier in CALIBRATED_TIER_ORDER:
        spec = tiers.get(tier)
        if not isinstance(spec, dict):
            raise ValueError(f"calibration is missing tier {tier}")
        floor, anchor, ceiling = (float(spec[key]) for key in ("floor", "anchor", "ceiling"))
        if not 0.0 <= floor <= anchor <= ceiling <= LEVEL_2_STAR_FLOOR:
            raise ValueError(f"tier {tier}: need 0 <= floor <= anchor <= ceiling <= {LEVEL_2_STAR_FLOOR}")
        if previous is not None and anchor > previous:
            raise ValueError(f"tier {tier}: anchor {anchor} exceeds the tier above ({previous})")
        previous = anchor
        spec.update(floor=floor, anchor=anchor, ceiling=ceiling, center=float(spec.get("center", 50.0)))
    return payload


def active_calibration() -> dict[str, Any] | None:
    """The calibration in ``CALIBRATION_PATH``, or ``None`` to use ``_CONFIDENCE_BANDS``."""

    key = str(CALIBRATION_PATH)
    if key not in _CALIBRATION_CACHE:
        calibration = None
        if CALIBRATION_PATH.exists():
            try:
                calibration = load_calibration(CALIBRATION_PATH)
            except (OSError, ValueError, KeyError, TypeError) as exc:
                warnings.warn(f"ignoring {CALIBRATION_PATH.name}: {exc}", stacklevel=2)
        _CALIBRATION_CACHE[key] = calibration
    return _CALIBRATION_CACHE[key]


def _margin_scores(ranked: list[dict[str, Any]]) -> dict[int, float | None]:
    """Per-candidate score behind the identity margin, keyed by ``id(candidate)``.

    The pure NCE score is used when rank 1 has a successful finite score,
    including zero, divided by the compound's best NCE score so the margin reads
    rank 1 against its rivals rather than the matcher's absolute scale
    (a calibration with ``"nce_normalization": "raw"`` skips the division).
    An all-zero successful result stays zero. Otherwise ``evidence_confidence``
    is used. One compound never mixes the two.
    """

    nce = [
        None if candidate.get("nce_status") in {"error", "no_observed_bins"}
        else _num(candidate.get("pure_nce_score"))
        for candidate in ranked
    ]
    present = [score for score in nce if score is not None]
    if nce and nce[0] is not None and present:
        calibration = active_calibration() or {}
        peak = 1.0 if calibration.get("nce_normalization") == "raw" else (max(present) or 1.0)
        return {id(candidate): None if score is None else score / peak for candidate, score in zip(ranked, nce)}
    return {id(candidate): _num(candidate.get("evidence_confidence")) for candidate in ranked}


def _best_margin_score(candidates: list[dict[str, Any]], scores: dict[int, float | None]) -> float | None:
    present = [scores.get(id(candidate)) for candidate in candidates]
    present = [score for score in present if score is not None]
    return max(present) if present else None


def _cluster_identity_confidence(candidates: list[dict[str, Any]]) -> float | None:
    """Cluster-margin score, before the Level 1 and Level 2* floors."""

    if not candidates:
        return None
    ranked = sorted(candidates, key=_rank)
    scores = _margin_scores(ranked)
    top_score = scores.get(id(ranked[0]))
    if top_score is None:
        return None
    if not HAS_RDKIT:
        return round(100.0 * top_score, 1)
    families = _cluster(ranked)
    top_family = _family_containing(families, ranked[0])
    unrelated = [member for family in families if family is not top_family for member in family]
    rival_score = _best_margin_score(unrelated, scores) if unrelated else None
    if rival_score is None:
        return round(100.0 * top_score, 1)
    denominator = top_score + rival_score
    if denominator <= 0:
        return 0.0
    return round(100.0 * (top_score * top_score) / denominator, 1)


_MOTIF_SMARTS = (
    ("glucuronide", "OC1OC(C(=O)O)C(O)C(O)C1O"),
    ("sulfate", "[SX4](=O)(=O)(O)O"),
    ("phosphate", "[PX4](=O)(O)(O)O"),
    ("taurine conjugate", "NCCS(=O)(=O)O"),
    ("carnitine conjugate", "C[N+](C)(C)CC(O)CC(=O)[O;H,-]"),
    ("hexose", "OC1OC(CO)C(O)C(O)C1O"),
)
_GLYCINE_SMARTS = "[CX3](=O)NCC(=O)[O;H,-]"
_AMIDE_SMARTS = "C(=O)N"


def _shared_motif(candidates: list[dict[str, Any]]) -> str:
    """A conjugate, sugar, or peptide on a strict majority of the candidates."""

    if not HAS_RDKIT or not candidates:
        return ""
    patterns = [(label, Chem.MolFromSmarts(smarts)) for label, smarts in _MOTIF_SMARTS]
    glycine = Chem.MolFromSmarts(_GLYCINE_SMARTS)
    amide = Chem.MolFromSmarts(_AMIDE_SMARTS)
    counts = {label: 0 for label, _pattern in patterns}
    glycine_single = 0
    peptide = 0
    for candidate in candidates:
        molecule = Chem.MolFromSmiles(_smiles(candidate))
        if molecule is None:
            continue
        for label, pattern in patterns:
            if pattern is not None and molecule.HasSubstructMatch(pattern):
                counts[label] += 1
        if glycine is not None and molecule.HasSubstructMatch(glycine):
            amides = len(molecule.GetSubstructMatches(amide)) if amide is not None else 0
            if amides <= 1:
                glycine_single += 1
            else:
                peptide += 1
    needed = len(candidates)
    for label in ("glucuronide", "sulfate", "phosphate", "taurine conjugate", "carnitine conjugate"):
        if counts[label] * 2 > needed:
            return label
    if glycine_single * 2 > needed:
        return "glycine conjugate"
    if counts["hexose"] * 2 > needed:
        return "hexose"
    if peptide * 2 > needed:
        return "peptide"
    return ""


def _supported_formula(
    candidates: list[dict[str, Any]],
) -> tuple[str, list[dict[str, Any]]]:
    """Formula shared by a strict majority of the mass-consistent candidates."""

    eligible: list[tuple[dict[str, Any], str]] = []
    for candidate in candidates:
        formula = _clean(candidate.get("formula") or candidate.get("candidate_formula"))
        if not formula:
            continue
        mass = _num(candidate.get("mass_match_score"))
        if mass is not None and mass < MASS_CONSISTENT:
            continue
        eligible.append((candidate, formula))
    if not eligible:
        return "", []
    counts: dict[str, int] = {}
    for _, formula in eligible:
        counts[formula] = counts.get(formula, 0) + 1
    ranked = sorted(eligible, key=lambda item: _rank(item[0]))
    rank1_formula = ranked[0][1]
    best_count = max(counts.values())
    chosen = rank1_formula if counts.get(rank1_formula) == best_count else max(counts, key=counts.get)
    if len(counts) > 1 and counts[chosen] * 2 <= len(eligible):
        return "", []
    return chosen, [candidate for candidate, formula in eligible if formula == chosen]


def _cluster(candidates: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    families: list[list[dict[str, Any]]] = []
    for candidate in sorted(candidates, key=_rank):
        placed = False
        for family in families:
            if _related(candidate, family[0]):
                family.append(candidate)
                placed = True
                break
        if not placed:
            families.append([candidate])
    return families


def _related(left: dict[str, Any], right: dict[str, Any]) -> bool:
    left_smiles = _smiles(left)
    right_smiles = _smiles(right)
    if not left_smiles or not right_smiles:
        return False
    if left_smiles == right_smiles:
        return True
    left_key = _connectivity(left_smiles)
    right_key = _connectivity(right_smiles)
    if left_key and left_key == right_key:
        return True
    if not HAS_FINGERPRINT:
        return False
    left_fp = _fingerprint(left_smiles)
    right_fp = _fingerprint(right_smiles)
    if left_fp is None or right_fp is None:
        return False
    return DataStructs.TanimotoSimilarity(left_fp, right_fp) >= RELATED_TANIMOTO


def _fingerprint(smiles: str):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    return RDKFingerprint(molecule)


def _connectivity(smiles: str) -> str:
    if not smiles or not HAS_RDKIT:
        return smiles
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return smiles
    try:
        return Chem.MolToInchiKey(molecule).split("-", 1)[0]
    except Exception:
        return smiles


def _one_connectivity(family: list[dict[str, Any]]) -> bool:
    keys = {_connectivity(_smiles(candidate)) for candidate in family}
    return len(keys) == 1 and "" not in keys


def _family_containing(
    families: list[list[dict[str, Any]]], candidate: dict[str, Any]
) -> list[dict[str, Any]]:
    rank = _rank(candidate)
    for family in families:
        if any(_rank(member) == rank for member in family):
            return family
    return [candidate]


def _outranks(family: list[dict[str, Any]], unrelated: list[dict[str, Any]]) -> bool:
    family_score = _best_evidence(family)
    other_score = _best_evidence(unrelated)
    if family_score is None or other_score is None:
        return False
    return family_score - other_score >= DISTINCTIVE_GAP


def _best_evidence(candidates: list[dict[str, Any]]) -> float | None:
    scores = [_num(candidate.get("evidence_confidence")) for candidate in candidates]
    present = [score for score in scores if score is not None]
    return max(present) if present else None


def _group_isomer_name(candidates: list[dict[str, Any]]) -> str:
    """Name of the dominant structural family, when nothing else named it.

    Clusters the candidates by structural relatedness. If the family holding
    rank 1 clears every other, unrelated family by ``DISTINCTIVE_GAP`` or
    more Match Score, that family is the metabolite's probable isomer group,
    named for its own top-ranked member.
    """

    if not HAS_RDKIT or not candidates:
        return ""
    ranked = sorted(candidates, key=_rank)
    families = _cluster(ranked)
    if len(families) < 2:
        return ""
    top = ranked[0]
    top_family = _family_containing(families, top)
    unrelated = [member for family in families if family is not top_family for member in family]
    if not unrelated or not _outranks(top_family, unrelated):
        return ""
    leader = min(top_family, key=_rank)
    name = _as_specific(_clean(leader.get("name")))
    if name:
        return name
    return _clean(leader.get("formula") or leader.get("candidate_formula"))


def msi_from_library_level(level: str) -> tuple[str, str]:
    """MSI level and label for a known-library Metabolon level.

    ``4?`` is reported as MSI 4. Labels that are not a level, such as
    ``Exclude``, stay blank.
    """

    text = _clean(level)
    if text == "4?":
        text = "4"
    label = MSI_LABELS.get(text, "")
    return (text, label) if label else ("", "")


def known_library_msi_rows(path) -> list[dict[str, str]]:
    """One row per known-library record, with MSI level and label added."""

    import csv

    columns = (
        "Unknown",
        "LIB_ID",
        "CompID",
        "KNOWN",
        "Formula",
        "Cleaned External Info",
        "Tentative information",
        "Level",
        "SubLevel",
        "msi_level",
        "msi_label",
        "SuperPathway",
        "Draft Subpathway",
    )
    rows: list[dict[str, str]] = []
    with open(path, newline="", encoding="utf-8-sig") as handle:
        for source in csv.DictReader(handle):
            msi_level, msi_label = msi_from_library_level(source.get("Level", ""))
            rows.append(
                {
                    "Unknown": _clean(source.get("Unknown")),
                    "LIB_ID": _clean(source.get("LIB_ID")),
                    "CompID": _clean(source.get("CompID")),
                    "KNOWN": _clean(source.get("KNOWN")),
                    "Formula": _clean(source.get("Formula")),
                    "Cleaned External Info": _clean(source.get("Cleaned External Info")),
                    "Tentative information": _clean(source.get("Tentative information")),
                    "Level": _clean(source.get("Level")),
                    "SubLevel": _clean(source.get("SubLevel")),
                    "msi_level": msi_level,
                    "msi_label": msi_label,
                    "SuperPathway": _clean(source.get("SuperPathway")),
                    "Draft Subpathway": _clean(source.get("Draft Subpathway")),
                }
            )
    rows.sort(key=lambda row: (row["Unknown"], row["LIB_ID"], row["CompID"]))
    return [{column: row[column] for column in columns} for row in rows]


def _nondefinitive_library_name(fields: dict[str, Any]) -> str:
    """A known-library name that names a modification, not one structure.

    ``glucuronide of C20H26O4`` is the case: the library named it, and the
    parent is a formula rather than a definitive structure.
    """

    for key in ("known_x_cleaned_external_info", "known_x_tentative_information"):
        for clause in _clauses(fields.get(key)):
            phrase = _as_class(clause)
            if phrase and _FORMULA_PARENT.search(phrase):
                return phrase
    return ""


def _is_level1(fields: dict[str, Any]) -> bool:
    if _clean(fields.get("exception_basis")) == "standard_confirmed":
        return True
    if _clean(fields.get("metabolon_level")) == "1":
        return True
    return _clean(fields.get("msi_level")) in {"1", "1.0"}


def _isomer_parent(fields: dict[str, Any]) -> str:
    feature = _clean(fields.get("feature_name"))
    match = _ISOMER_PREFIX.match(feature)
    if match:
        return feature[match.end() :].strip()
    excluded = _clean(fields.get("excluded_library_reference"))
    reason = _clean(fields.get("identity_exclusion_reason"))
    collision = _clean(fields.get("library_collision")).lower()
    if excluded and (reason == "existing_global_library_structure" or collision in {"true", "1", "yes"}):
        return excluded
    return ""


def _identity_name(fields: dict[str, Any]) -> str:
    return _report_specific(fields) or _specific_name(fields) or _class_phrase(fields)


def _report_specific(fields: dict[str, Any]) -> str:
    """A name the report itself is willing to state without a hedge."""

    raw = _clean(fields.get("report_display_name"))
    if _HEDGE_PREFIX.match(raw) or _ISOMER_PREFIX.match(raw):
        return ""
    return _as_specific(raw)


def _specific_name(fields: dict[str, Any]) -> str:
    for key in ("report_display_name", "feature_name", "known_x_tentative_information"):
        for clause in _clauses(fields.get(key)):
            specific = _as_specific(clause)
            if specific:
                return specific
    return ""


def _class_phrase(fields: dict[str, Any]) -> str:
    for key in ("known_x_cleaned_external_info", "known_x_tentative_information", "feature_name"):
        for clause in _clauses(fields.get(key)):
            phrase = _as_class(clause)
            if phrase:
                return phrase
    return ""


def _as_specific(text: str) -> str:
    cleaned = _strip_hedges(text)
    cleaned = re.sub(r"\s+MARS PROPRIETARY$", "", cleaned, flags=re.I).strip()
    match = _ISOMER_PREFIX.match(cleaned)
    if match:
        cleaned = _strip_hedges(cleaned[match.end() :])
    if (
        not cleaned
        or _is_placeholder(cleaned)
        or _is_class(cleaned)
        or _is_formulaish(cleaned)
        or _NOTE_PREFIX.match(cleaned)
        or _NOTE.search(cleaned)
        or _VAGUE_CLASS_WORD.match(cleaned)
        or re.match(r"(?i)^(?:a|an)\s+", cleaned)
    ):
        return ""
    words = re.findall(r"[A-Za-z]+", cleaned.lower())
    if not words or all(word in _GENERIC_WORDS for word in words):
        return ""
    if not re.search(r"[A-Za-z]", cleaned):
        return ""
    return cleaned


def _as_class(text: str) -> str:
    cleaned = _strip_hedges(text)
    if (
        _is_placeholder(cleaned)
        or _NOTE.search(cleaned)
        or _VAGUE_CLASS_WORD.match(cleaned)
        or cleaned.lower() == "unknown"
        or not _is_class(cleaned)
    ):
        return ""
    return cleaned


def _is_formulaish(text: str) -> bool:
    compact = text.replace(" ", "")
    if _BARE_FORMULA.match(compact):
        return True
    if not compact.startswith("C"):
        return False
    remainder = re.sub(r"[A-Z][a-z]?\d*", "", compact)
    remainder = re.sub(r"[\+\-\(\)\[\],;/]+", "", remainder)
    return remainder == ""


def _is_class(text: str) -> bool:
    if not text:
        return False
    if _CLASS_OF.match(text):
        return True
    if _GENERIC_CLASS.match(text):
        return True
    if _PEPTIDE_COMPOSITION.search(text):
        return True
    return bool(_FORMULA_PARENT.search(text) and re.search(r"(?i)sulfate|glucuronide|peptide|conjugate", text))


def _strip_hedges(text: str) -> str:
    cleaned = _clean(text)
    while True:
        shortened = _HEDGE_PREFIX.sub("", cleaned).strip()
        if shortened == cleaned:
            return cleaned
        cleaned = shortened


def _clauses(value: Any) -> list[str]:
    return [part.strip() for part in _clean(value).split(";") if part.strip()]


def _is_placeholder(text: str) -> bool:
    lower = text.lower()
    return (
        not text
        or "structure unresolved" in lower
        or lower.startswith("exact mass ")
        or lower.startswith("tentative structure")
        or lower.startswith("de novo structure")
    )


def _smiles(candidate: dict[str, Any]) -> str:
    return _clean(candidate.get("smiles") or candidate.get("candidate_smiles"))


def _rank(candidate: dict[str, Any]) -> int:
    number = _num(candidate.get("rank"))
    return int(number) if number is not None else 10**9


def _num(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _clean(value: Any) -> str:
    if value is None:
        return ""
    return " ".join(str(value).split())
