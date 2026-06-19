#!/usr/bin/env python3
"""Parse Metabolon HTML feature reports into structured formula/isotope signals.

Pipeline context:
    Truth-blind enrichment and isotope audits read per-feature HTML exports
    to extract formula constraints, buddy labels, mass-library ion names, and
    isotope ratio tables without re-opening reports in downstream steps.

Key I/O:
    ``parse_html_signals``: HTML path → feature signal dict for scoring joins.
    ``audit_candidate_formula_against_html_signals``: formula vs isotope audit.
"""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

from bs4 import BeautifulSoup

# ── Formula parsing patterns ───────────────────────────────────────────
INLINE_FORMULA_RE = re.compile(r"(?<![A-Za-z0-9])((?:[A-Z][a-z]?\d*){2,}(?:[+-]\d*)?)(?![A-Za-z0-9])")
FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*)")
VALID_ELEMENTS = {
    "Ac",
    "Ag",
    "Al",
    "Am",
    "Ar",
    "As",
    "At",
    "Au",
    "B",
    "Ba",
    "Be",
    "Bh",
    "Bi",
    "Bk",
    "Br",
    "C",
    "Ca",
    "Cd",
    "Ce",
    "Cf",
    "Cl",
    "Cm",
    "Cn",
    "Co",
    "Cr",
    "Cs",
    "Cu",
    "Db",
    "Ds",
    "Dy",
    "Er",
    "Es",
    "Eu",
    "F",
    "Fe",
    "Fl",
    "Fm",
    "Fr",
    "Ga",
    "Gd",
    "Ge",
    "H",
    "He",
    "Hf",
    "Hg",
    "Ho",
    "Hs",
    "I",
    "In",
    "Ir",
    "K",
    "Kr",
    "La",
    "Li",
    "Lr",
    "Lu",
    "Lv",
    "Mc",
    "Md",
    "Mg",
    "Mn",
    "Mo",
    "Mt",
    "N",
    "Na",
    "Nb",
    "Nd",
    "Ne",
    "Nh",
    "Ni",
    "No",
    "Np",
    "O",
    "Og",
    "Os",
    "P",
    "Pa",
    "Pb",
    "Pd",
    "Pm",
    "Po",
    "Pr",
    "Pt",
    "Pu",
    "Ra",
    "Rb",
    "Re",
    "Rf",
    "Rg",
    "Rh",
    "Rn",
    "Ru",
    "S",
    "Sb",
    "Sc",
    "Se",
    "Sg",
    "Si",
    "Sm",
    "Sn",
    "Sr",
    "Ta",
    "Tb",
    "Tc",
    "Te",
    "Th",
    "Ti",
    "Tl",
    "Tm",
    "Ts",
    "U",
    "V",
    "W",
    "Xe",
    "Y",
    "Yb",
    "Zn",
    "Zr",
}
# ── HTML table typing ────────────────────────────────────────────────
FORMULA_COLUMNS = {"mol_formula", "neutral", "formula"}
INT_COLUMNS = {"comp_id", "chemical_id", "lib_id", "library_matches", "tally"}
FLOAT_COLUMNS = {
    "c_predict",
    "dbe",
    "delta",
    "harmony",
    "inferred_m",
    "ion_mass",
    "mass",
    "match_score",
    "mz_buddy",
    "parent_mass",
    "ppm",
    "ratio",
    "ri",
}
TABLE_SIGNATURES: dict[tuple[str, ...], str] = {
    ("mz", "norm_int"): "ms2",
    ("mz_buddy", "tally", "ratio", "inferred_m", "buddy", "feature"): "buddy",
    ("comp_id", "ri", "chemical_id", "mol_formula", "name"): "inferred_mass_matches",
    (
        "chemical_id",
        "comp_id",
        "ion_mass",
        "ri",
        "mol_formula",
        "smiles",
        "ion_name",
        "name",
    ): "mass_library",
    (
        "mass",
        "ppm",
        "neutral",
        "dbe",
        "harmony",
        "ri_pred",
        "source",
        "id",
        "smiles",
        "library_matches",
    ): "candidate_formula",
    ("match_score", "lib_id", "comp_id", "parent_mass", "ri", "formula", "smiles", "name"): "parent_mass_matches",
}
# ── Isotopologue detection ─────────────────────────────────────────────
COMMENT_ISOTOPOLOGUE_WARNING_RE = re.compile(
    r"(\bisotope of\b|\bpotential overlap\b.*\bisotope\b|\bappears to be\b.*\bisotope\b)",
    flags=re.IGNORECASE,
)
ISOTOPE_TOKEN_RE = re.compile(r"\[(13C|18O|34S|15N|37Cl|81Br)(\d+)\]")
FEATURE_ISOTOPE_TOKEN_RE = re.compile(r"\[(C13|O18|S34|N15|Cl37|Br81)-(\d+)\]")
ISOTOPE_TO_ELEMENT = {
    "13C": "C",
    "18O": "O",
    "34S": "S",
    "15N": "N",
    "37Cl": "Cl",
    "81Br": "Br",
}
FEATURE_ISOTOPE_TO_ELEMENT = {
    "C13": "C",
    "O18": "O",
    "S34": "S",
    "N15": "N",
    "Cl37": "Cl",
    "Br81": "Br",
}
ISOTOPE_ABUNDANCE = {
    "13C": 0.0107,
    "18O": 0.00205,
    "34S": 0.0421,
    "15N": 0.00364,
    "37Cl": 0.2422,
    "81Br": 0.4931,
}
ELEMENT_PREDICT_HEADERS = {
    "C": "c_predict",
    "O": "o_predict",
    "S": "s_predict",
    "N": "n_predict",
    "Cl": "cl_predict",
    "Br": "br_predict",
}
# ── Default empty signal record ────────────────────────────────────────
EMPTY_SIGNAL_RECORD = {
    "comment_count": 0,
    "comment_texts": [],
    "comment_formulas": [],
    "inferred_mass_match_row_count": 0,
    "inferred_mass_x_row_count": 0,
    "inferred_mass_formulas": [],
    "inferred_mass_x_formulas": [],
    "mass_library_row_count": 0,
    "mass_library_x_row_count": 0,
    "mass_library_formulas": [],
    "mass_library_x_formulas": [],
    "mass_library_ion_names": [],
    "mass_library_x_ion_names": [],
    "mass_library_x_ion_masses": [],
    "mass_library_ion_name_unique_count": 0,
    "mass_library_x_ion_name_unique_count": 0,
    "candidate_formula_row_count": 0,
    "candidate_formula_unique_count": 0,
    "candidate_formulas": [],
    "candidate_formula_top_weight": 0,
    "candidate_formula_top_formulas": [],
    "parent_mass_row_count": 0,
    "parent_mass_x_row_count": 0,
    "parent_mass_formulas": [],
    "parent_mass_x_formulas": [],
    "buddy_row_count": 0,
    "buddy_labels": [],
    "non_isotope_buddy_labels": [],
    "mass_library_isotope_labels": [],
    "buddy_isotope_labels": [],
    "isotope_element_hints": {},
    "feature_is_likely_isotopologue": False,
    "feature_isotopologue_reason": None,
    "inferred_neutral_masses": [],
    "feature_level_formula_set": [],
    "candidate_level_formula_set": [],
    "usable_formula_source": "none",
    "usable_formula_set": [],
    "mass_library_top_ion_name": None,
    "mass_library_top_ion_count": 0,
    "mass_library_top_ion_share": 0.0,
    "mass_library_x_top_ion_name": None,
    "mass_library_x_top_ion_count": 0,
    "mass_library_x_top_ion_share": 0.0,
    "top_buddy_label": None,
    "top_buddy_count": 0,
    "top_buddy_share": 0.0,
    "confident_dominant_ion": False,
    "mass_library_ion_name_counter": {},
    "mass_library_x_ion_name_counter": {},
    "buddy_counter": {},
    "observed_isotope_rows": [],
    "observed_isotope_rows_json": [],
    "isotope_signal_count": 0,
    "c_predict": None,
    "parse_status": "ok",
}


# ── Text and formula normalization ───────────────────────────────────
def normalize_header(text: object) -> str:
    """Normalize HTML table header text to snake_case identifiers."""
    value = str(text).replace("\xa0", " ").strip().lower()
    return re.sub(r"[^a-z0-9]+", "_", value).strip("_")


def normalize_text(value: object) -> str | None:
    """Collapse whitespace and map empty/nan/null sentinels to ``None``."""
    if value is None:
        return None
    if isinstance(value, float) and not math.isfinite(value):
        return None
    text = re.sub(r"\s+", " ", str(value)).strip()
    if not text or text.lower() in {"none", "nan", "null"}:
        return None
    return text


def normalize_formula(value: object) -> str | None:
    """Validate and compact a molecular formula string."""
    text = normalize_text(value)
    if text is None:
        return None
    compact = text.replace(" ", "")
    body = re.sub(r"(?:[+-]\d*)$", "", compact)
    tokens = FORMULA_TOKEN_RE.findall(body)
    if not tokens:
        return None
    reconstructed = "".join(element + count for element, count in tokens)
    if reconstructed != body:
        return None
    if any(element not in VALID_ELEMENTS for element, _ in tokens):
        return None
    if not any(count for _, count in tokens):
        return None
    return compact


def unique_preserve_order(values: Iterable[object]) -> list[str]:
    """Deduplicate normalized strings while preserving first-seen order."""
    out: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = normalize_text(value)
        if text is None or text in seen:
            continue
        out.append(text)
        seen.add(text)
    return out


def unique_formulas(values: Iterable[object]) -> list[str]:
    """Deduplicate validated molecular formulas while preserving order."""
    out: list[str] = []
    seen: set[str] = set()
    for value in values:
        formula = normalize_formula(value)
        if formula is None or formula in seen:
            continue
        out.append(formula)
        seen.add(formula)
    return out


def formula_element_counts(value: object) -> dict[str, int]:
    """Parse a formula into per-element integer atom counts."""
    formula = normalize_formula(value)
    if formula is None:
        return {}
    body = re.sub(r"(?:[+-]\d*)$", "", formula.replace(" ", ""))
    counts: dict[str, int] = {}
    for element, count_text in FORMULA_TOKEN_RE.findall(body):
        counts[element] = counts.get(element, 0) + (int(count_text) if count_text else 1)
    return counts


def safe_int(value: object) -> int | None:
    """Parse an integer from messy table text, returning ``None`` on failure."""
    try:
        parsed = int(float(str(value)))
    except Exception:
        return None
    return parsed


def safe_float(value: object) -> float | None:
    """Parse a finite float from messy table text, returning ``None`` on failure."""
    try:
        parsed = float(str(value))
    except Exception:
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def extract_inline_formulas(text: object) -> list[str]:
    """Pull validated molecular formulas embedded in comment text."""
    normalized = normalize_text(text)
    if normalized is None:
        return []
    hits = [normalize_formula(match.group(1)) for match in INLINE_FORMULA_RE.finditer(normalized)]
    return [value for value in unique_preserve_order(hits) if value is not None]


def isotope_element_hints_from_label(label: object) -> dict[str, int]:
    """Parse isotope bracket tokens into minimum heavy-atom counts."""
    text = normalize_text(label)
    if text is None:
        return {}

    hints: dict[str, int] = {}
    for isotope, count_text in ISOTOPE_TOKEN_RE.findall(text):
        element = ISOTOPE_TO_ELEMENT.get(isotope)
        if element is None:
            continue
        hints[element] = max(hints.get(element, 0), int(count_text))

    for isotope, count_text in FEATURE_ISOTOPE_TOKEN_RE.findall(text):
        element = FEATURE_ISOTOPE_TO_ELEMENT.get(isotope)
        if element is None:
            continue
        hints[element] = max(hints.get(element, 0), int(count_text))
    return hints


def collect_isotope_labels(values: Iterable[object]) -> list[str]:
    """Collect unique labels that encode isotope enrichment hints."""
    out: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = normalize_text(value)
        if text is None:
            continue
        if not isotope_element_hints_from_label(text):
            continue
        if text in seen:
            continue
        out.append(text)
        seen.add(text)
    return out


def merge_isotope_element_hints(*label_groups: Iterable[object]) -> dict[str, int]:
    """Merge isotope hints from multiple label groups by max count."""
    hints: dict[str, int] = {}
    for labels in label_groups:
        for label in labels:
            for element, count in isotope_element_hints_from_label(label).items():
                hints[element] = max(hints.get(element, 0), int(count))
    return dict(sorted(hints.items()))


def formula_satisfies_isotope_hints(
    formula: object,
    isotope_element_hints: dict[str, int] | None,
) -> bool | None:
    """Return whether a formula meets minimum isotope element counts."""
    if not isotope_element_hints:
        return None
    counts = formula_element_counts(formula)
    if not counts:
        return None
    return all(int(counts.get(element, 0)) >= int(required) for element, required in isotope_element_hints.items())


def annotate_formula_against_html_signals(
    formula: object,
    feature_signal: dict[str, object] | None,
) -> dict[str, object]:
    """Attach isotope-consistency flags for one formula vs HTML signals."""
    feature_signal = feature_signal or {}
    isotope_element_hints = dict(feature_signal.get("isotope_element_hints") or {})
    consistent = formula_satisfies_isotope_hints(formula, isotope_element_hints)
    return {
        "feature_is_likely_isotopologue": bool(feature_signal.get("feature_is_likely_isotopologue", False)),
        "feature_isotopologue_reason": feature_signal.get("feature_isotopologue_reason"),
        "isotope_element_hints": isotope_element_hints,
        "isotope_formula_consistent": consistent,
        "isotope_formula_score": None if consistent is None else float(1.0 if consistent else 0.0),
    }


def detect_feature_isotopologue_reason(
    buddy_labels: Iterable[object],
    comment_texts: Iterable[object],
    *,
    correct_isotopologue_detection: bool = False,
) -> str | None:
    """Infer why a feature is flagged as a likely isotopologue."""
    for label in buddy_labels:
        text = normalize_text(label)
        if text is None:
            continue
        lowered = text.lower()
        if correct_isotopologue_detection:
            if lowered.startswith("feature[") or lowered.startswith("self["):
                return f"buddy_label:{text}"
            continue
        if lowered.startswith("x["):
            return f"buddy_label:{text}"
    for comment in comment_texts:
        text = normalize_text(comment)
        if text is None:
            continue
        match = COMMENT_ISOTOPOLOGUE_WARNING_RE.search(text)
        if match:
            return f"comment_warning:{match.group(1)}"
    return None


def parse_isotope_label(label: object) -> dict[str, object] | None:
    """Parse buddy or isotope-table labels into element, key, and multiplicity."""
    text = normalize_text(label)
    if text is None:
        return None

    feature_match = FEATURE_ISOTOPE_TOKEN_RE.search(text)
    if feature_match is not None:
        isotope, count_text = feature_match.groups()
        element = FEATURE_ISOTOPE_TO_ELEMENT.get(isotope)
        if element is None:
            return None
        isotope_key = {
            "C13": "13C",
            "O18": "18O",
            "S34": "34S",
            "N15": "15N",
            "Cl37": "37Cl",
            "Br81": "81Br",
        }.get(isotope)
        return {
            "label": text,
            "element": element,
            "isotope_key": isotope_key,
            "k": int(count_text),
        }

    mass_match = ISOTOPE_TOKEN_RE.search(text)
    if mass_match is not None:
        isotope_key, count_text = mass_match.groups()
        element = ISOTOPE_TO_ELEMENT.get(isotope_key)
        if element is None:
            return None
        return {
            "label": text,
            "element": element,
            "isotope_key": isotope_key,
            "k": int(count_text),
        }

    return None


def expected_isotope_ratio(element_count: int, isotope_key: str, k: int) -> float | None:
    """Compute binomial expected isotopologue intensity ratio for ``k`` labels."""
    if element_count < 0 or k <= 0 or element_count < k:
        return None
    abundance = ISOTOPE_ABUNDANCE.get(str(isotope_key))
    if abundance is None or abundance <= 0.0 or abundance >= 1.0:
        return None
    return float(math.comb(int(element_count), int(k)) * ((abundance / (1.0 - abundance)) ** int(k)))


def isotope_ratio_score(observed_ratio: float, expected_ratio: float) -> float | None:
    """Score observed vs expected isotope ratio on a log-distance scale."""
    if observed_ratio <= 0.0 or expected_ratio <= 0.0:
        return None
    distance = abs(math.log(float(observed_ratio) / float(expected_ratio)))
    return float(max(0.0, 1.0 - (distance / math.log(3.0))))


def audit_candidate_formula_against_html_signals(
    formula: object,
    feature_signal: dict[str, object] | None,
    ratio_supported_min_score: float = 0.5,
) -> dict[str, object]:
    """Run full isotope audit including ratio support for one formula."""
    feature_signal = feature_signal or {}
    base = annotate_formula_against_html_signals(formula, feature_signal)
    signal_rows = list(feature_signal.get("observed_isotope_rows") or [])
    signal_count = int(feature_signal.get("isotope_signal_count", len(signal_rows)) or 0)
    parse_status = str(feature_signal.get("parse_status") or "unavailable")
    formula_counts = formula_element_counts(formula)
    ratio_score_value: float | None = None
    audit_status = "unavailable"

    if parse_status != "ok":
        audit_status = "unavailable"
    elif not signal_count and not bool(base["feature_is_likely_isotopologue"]):
        audit_status = "no_signal"
    elif not formula_counts and (signal_count or base["isotope_element_hints"]):
        audit_status = "unavailable"
    elif base["isotope_formula_consistent"] is False:
        audit_status = "formula_impossible"
    elif bool(base["feature_is_likely_isotopologue"]):
        audit_status = "feature_isotopologue"
    else:
        ratio_scores: list[float] = []
        for row in signal_rows:
            observed_ratio = safe_float(row.get("ratio"))
            if observed_ratio is None or observed_ratio <= 0.0:
                continue
            element = str(row.get("element") or "").strip()
            isotope_key = str(row.get("isotope_key") or "").strip()
            k = safe_int(row.get("k"))
            if not element or not isotope_key or k is None:
                label_details = parse_isotope_label(row.get("label"))
                if label_details is None:
                    continue
                element = str(label_details["element"])
                isotope_key = str(label_details["isotope_key"])
                k = int(label_details["k"])
            element_count = int(formula_counts.get(element, 0))
            expected_ratio = expected_isotope_ratio(element_count, isotope_key, int(k))
            if expected_ratio is None:
                continue
            score_value = isotope_ratio_score(observed_ratio, expected_ratio)
            if score_value is not None:
                ratio_scores.append(float(score_value))
        if ratio_scores:
            ratio_score_value = float(sum(ratio_scores) / len(ratio_scores))
            audit_status = (
                "ratio_supported" if ratio_score_value >= float(ratio_supported_min_score) else "ratio_mismatch"
            )
        else:
            audit_status = "no_signal"

    return {
        **base,
        "isotope_ratio_score": ratio_score_value,
        "isotope_audit_status": audit_status,
        "isotope_signal_count": signal_count,
    }


def table_headers(table) -> tuple[str, ...]:
    """Normalize ``<thead>`` cell text into a tuple of snake_case headers."""
    thead = table.find("thead")
    if thead is None:
        return ()
    header_cells = thead.find_all("th")
    if not header_cells:
        return ()
    headers = tuple(normalize_header(cell.get_text(" ", strip=True)) for cell in header_cells)
    return tuple(header for header in headers if header)


def is_candidate_isotope_header_set(headers: tuple[str, ...]) -> bool:
    """Return whether headers match the candidate-isotope table signature."""
    if len(headers) < 6:
        return False
    if headers[:5] != ("tally", "delta", "ratio", "isotope", "has_is"):
        return False
    return all(header.endswith("_predict") for header in headers[5:])


def identify_table_kind(headers: tuple[str, ...]) -> str | None:
    """Map normalized headers to a known Metabolon HTML table kind."""
    if is_candidate_isotope_header_set(headers):
        return "candidate_isotopes"
    return TABLE_SIGNATURES.get(headers)


def table_rows(table, headers: tuple[str, ...]) -> list[dict[str, object]]:
    """Parse HTML table body rows into typed dict records."""
    tbody = table.find("tbody")
    if tbody is None:
        return []

    out: list[dict[str, object]] = []
    for tr in tbody.find_all("tr", recursive=False):
        cells = tr.find_all(["td", "th"], recursive=False)
        if len(cells) < len(headers):
            continue
        row: dict[str, object] = {}
        for idx, header in enumerate(headers):
            cell = cells[idx]
            if header == "smiles":
                row[header] = None
                continue

            raw_text = cell.get_text("" if header in FORMULA_COLUMNS else " ", strip=True)
            if header in FORMULA_COLUMNS:
                row[header] = normalize_formula(raw_text)
            elif header in INT_COLUMNS:
                row[header] = safe_int(raw_text)
            elif header in FLOAT_COLUMNS or header.endswith("_predict"):
                row[header] = safe_float(raw_text)
            else:
                row[header] = normalize_text(raw_text)
        out.append(row)
    return out


def extract_comments_block(html: str) -> list[str]:
    """Extract normalized comment paragraphs from the Comments section."""
    match = re.search(
        r"<h2>\s*Comments\s*</h2>(.*?)(?:<button class=\"accordion\">|</body>)",
        html,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if not match:
        return []
    block = match.group(1)
    soup = BeautifulSoup(block, "html.parser")
    texts: list[str] = []
    for tag in soup.find_all(["h3", "p"]):
        text = normalize_text(tag.get_text(" ", strip=True))
        if text is None or text.lower() == "comments":
            continue
        texts.append(text)
    return unique_preserve_order(texts)


def strip_embedded_image_payloads(html: str) -> str:
    """Remove inline base64 ``<img>`` payloads to shrink HTML before parsing."""
    return re.sub(r"<img\b[^>]*\bsrc=[\"']data:image/[^\"']+[\"'][^>]*>", "", html, flags=re.IGNORECASE)


def is_x_row(row: dict[str, object], xid: str) -> bool:
    """Return whether a table row name column matches the feature X identifier."""
    return normalize_text(row.get("name")) == xid


def build_formula_weight_counter(rows: Iterable[dict[str, object]], formula_key: str) -> dict[str, int]:
    """Tally candidate formulas weighted by library match counts."""
    weights: dict[str, int] = defaultdict(int)
    for row in rows:
        formula = normalize_formula(row.get(formula_key))
        if formula is None:
            continue
        weight = safe_int(row.get("library_matches"))
        weights[formula] += weight if weight is not None else 1
    return dict(sorted(weights.items()))


def counter_stats(counter: Counter[str]) -> tuple[str | None, int, float]:
    """Return top counter value, count, and share for dominance heuristics."""
    if not counter:
        return None, 0, 0.0
    top_value, top_count = counter.most_common(1)[0]
    total = sum(counter.values())
    share = float(top_count / max(total, 1))
    return top_value, top_count, share


def empty_signal_record(
    *,
    feature_key: str,
    feature_label: str,
    source_html: object,
    parse_status: str,
) -> dict[str, object]:
    """Build a parse-status placeholder signal dict for missing HTML input."""
    source_text = normalize_text(source_html) or ""
    return {
        "feature_key": feature_key,
        "feature_label": feature_label,
        "source_html": source_text,
        **EMPTY_SIGNAL_RECORD,
        "usable_formula_source": parse_status if parse_status != "ok" else "none",
        "parse_status": parse_status,
    }


def normalize_observed_isotope_row(
    row: dict[str, object],
    *,
    source_table: str,
) -> dict[str, object] | None:
    """Normalize one buddy or candidate-isotope row for scoring joins."""
    label_text = normalize_text(row.get("isotope") if source_table == "candidate_isotopes" else row.get("buddy"))
    parsed = parse_isotope_label(label_text)
    if parsed is None:
        return None
    element = str(parsed["element"])
    predictor_header = ELEMENT_PREDICT_HEADERS.get(element)
    predictor_value = safe_float(row.get(predictor_header)) if predictor_header else None
    return {
        "label": str(parsed["label"]),
        "ratio": safe_float(row.get("ratio")),
        "tally": safe_int(row.get("tally")),
        "has_is": normalize_text(row.get("has_is")),
        "source_table": source_table,
        "element": element,
        "isotope_key": str(parsed["isotope_key"]),
        "k": int(parsed["k"]),
        "predictor_header": predictor_header,
        "predictor_value": predictor_value,
    }


def artifact_observed_isotope_rows(rows: Iterable[dict[str, object]]) -> list[dict[str, object]]:
    """Trim isotope rows to stable JSON-serializable fields for artifacts."""
    out: list[dict[str, object]] = []
    for row in rows:
        out.append(
            {
                "label": row.get("label"),
                "ratio": row.get("ratio"),
                "tally": row.get("tally"),
                "has_is": row.get("has_is"),
                "source_table": row.get("source_table"),
            }
        )
    return out


def parse_html_signals(
    source_html: Path | str | None,
    feature_key: str,
    feature_label: str,
    candidate_formula_max: int,
    correct_isotopologue_detection: bool = False,
) -> dict[str, object]:
    """Parse one feature HTML report into a structured signal record.

    Extracts comments, mass-library rows, buddy labels, candidate formulas,
    and isotope tables; derives usable formula sets for downstream scoring.

    Args:
        source_html: Path to HTML file (or missing/unavailable sentinel).
        feature_key: Internal feature identifier.
        feature_label: Display X ID used to detect X-specific table rows.
        candidate_formula_max: Max unique candidate formulas for narrow usable set.
        correct_isotopologue_detection: Use feature/self buddy labels for isotopologue flag.

    Returns:
        Signal dict matching ``EMPTY_SIGNAL_RECORD`` keys plus metadata.
    """
    source_text = normalize_text(source_html)
    if source_text is None:
        return empty_signal_record(
            feature_key=feature_key,
            feature_label=feature_label,
            source_html=source_html,
            parse_status="unavailable",
        )

    source_path = Path(source_text).expanduser()
    if not source_path.exists():
        return empty_signal_record(
            feature_key=feature_key,
            feature_label=feature_label,
            source_html=source_text,
            parse_status="unavailable",
        )

    html = strip_embedded_image_payloads(source_path.read_text(encoding="utf-8", errors="ignore"))
    soup = BeautifulSoup(html, "html.parser")

    tables: dict[str, list[dict[str, object]]] = defaultdict(list)
    for table in soup.find_all("table"):
        headers = table_headers(table)
        if not headers:
            continue
        kind = identify_table_kind(headers)
        if kind is None:
            continue
        rows = table_rows(table, headers)
        if kind == "buddy":
            rows = [
                row
                for row in rows
                if any(normalize_text(row.get(key)) for key in ("buddy", "feature"))
                or row.get("inferred_m") is not None
            ]
        elif kind == "candidate_formula":
            rows = [row for row in rows if normalize_formula(row.get("neutral")) is not None]
        elif kind == "mass_library":
            rows = [
                row
                for row in rows
                if any(normalize_text(row.get(key)) for key in ("name", "ion_name")) or row.get("ion_mass") is not None
            ]
        elif kind == "inferred_mass_matches":
            rows = [
                row
                for row in rows
                if any(normalize_text(row.get(key)) for key in ("name", "mol_formula"))
                or row.get("comp_id") is not None
            ]
        elif kind == "parent_mass_matches":
            rows = [
                row
                for row in rows
                if any(normalize_text(row.get(key)) for key in ("name", "formula"))
                or row.get("parent_mass") is not None
            ]
        elif kind == "candidate_isotopes":
            rows = [
                row
                for row in rows
                if normalize_text(row.get("isotope")) is not None and safe_float(row.get("ratio")) is not None
            ]
        tables[kind].extend(rows)

    comment_texts = extract_comments_block(html)
    comment_formulas = unique_formulas(
        formula for comment in comment_texts for formula in extract_inline_formulas(comment)
    )

    inferred_rows = tables.get("inferred_mass_matches", [])
    inferred_x_rows = [row for row in inferred_rows if is_x_row(row, feature_label)]
    inferred_formulas = unique_formulas(row.get("mol_formula") for row in inferred_rows)
    inferred_x_formulas = unique_formulas(row.get("mol_formula") for row in inferred_x_rows)

    mass_library_rows = tables.get("mass_library", [])
    mass_library_x_rows = [row for row in mass_library_rows if is_x_row(row, feature_label)]
    mass_library_formulas = unique_formulas(row.get("mol_formula") for row in mass_library_rows)
    mass_library_x_formulas = unique_formulas(row.get("mol_formula") for row in mass_library_x_rows)

    mass_library_ion_names = unique_preserve_order(row.get("ion_name") for row in mass_library_rows)
    mass_library_x_ion_names = unique_preserve_order(row.get("ion_name") for row in mass_library_x_rows)
    mass_library_x_ion_masses = [
        float(value)
        for value in unique_preserve_order(row.get("ion_mass") for row in mass_library_x_rows)
        if safe_float(value) is not None
    ]

    candidate_formula_rows = tables.get("candidate_formula", [])
    candidate_formula_weights = build_formula_weight_counter(candidate_formula_rows, "neutral")
    candidate_formulas = list(candidate_formula_weights.keys())
    candidate_top_weight = max(candidate_formula_weights.values(), default=0)
    candidate_top_formulas = [
        formula
        for formula, weight in candidate_formula_weights.items()
        if weight == candidate_top_weight and candidate_top_weight > 0
    ]

    parent_mass_rows = tables.get("parent_mass_matches", [])
    parent_mass_x_rows = [row for row in parent_mass_rows if is_x_row(row, feature_label)]
    parent_mass_formulas = unique_formulas(row.get("formula") for row in parent_mass_rows)
    parent_mass_x_formulas = unique_formulas(row.get("formula") for row in parent_mass_x_rows)

    buddy_rows = tables.get("buddy", [])
    buddy_labels = unique_preserve_order(row.get("buddy") for row in buddy_rows)
    non_isotope_buddy_labels = [
        label for label in buddy_labels if normalize_text(label) is not None and not str(label).startswith("x[")
    ]
    mass_library_isotope_labels = collect_isotope_labels(mass_library_ion_names)
    buddy_isotope_labels = collect_isotope_labels(buddy_labels)
    candidate_isotope_rows = tables.get("candidate_isotopes", [])
    candidate_isotope_labels = collect_isotope_labels(row.get("isotope") for row in candidate_isotope_rows)

    observed_isotope_rows = [
        normalized
        for normalized in (
            normalize_observed_isotope_row(row, source_table="candidate_isotopes") for row in candidate_isotope_rows
        )
        if normalized is not None
    ]
    if not observed_isotope_rows:
        observed_isotope_rows = [
            normalized
            for normalized in (normalize_observed_isotope_row(row, source_table="buddy") for row in buddy_rows)
            if normalized is not None
        ]

    if candidate_isotope_labels:
        isotope_element_hints = merge_isotope_element_hints(candidate_isotope_labels)
    else:
        isotope_element_hints = merge_isotope_element_hints(mass_library_isotope_labels, buddy_isotope_labels)

    feature_isotopologue_reason = detect_feature_isotopologue_reason(
        buddy_labels=buddy_labels,
        comment_texts=comment_texts,
        correct_isotopologue_detection=bool(correct_isotopologue_detection),
    )
    feature_is_likely_isotopologue = feature_isotopologue_reason is not None
    inferred_masses = [
        value for value in (safe_float(row.get("inferred_m")) for row in buddy_rows) if value is not None
    ]

    feature_level_formula_set = unique_formulas(
        list(comment_formulas) + inferred_x_formulas + mass_library_x_formulas + parent_mass_x_formulas
    )
    candidate_level_formula_set = unique_formulas(
        candidate_formulas + inferred_formulas + mass_library_formulas + parent_mass_formulas
    )
    scoring_feature_level_formula_set = [] if feature_is_likely_isotopologue else feature_level_formula_set

    if scoring_feature_level_formula_set:
        usable_formula_source = "feature_level"
        usable_formula_set = scoring_feature_level_formula_set
    elif 0 < len(candidate_level_formula_set) <= int(candidate_formula_max):
        usable_formula_source = "candidate_level_narrow"
        usable_formula_set = candidate_level_formula_set
    elif candidate_level_formula_set:
        usable_formula_source = "candidate_level_broad"
        usable_formula_set = []
    else:
        usable_formula_source = "none"
        usable_formula_set = []

    ion_counter = Counter(
        ion_name
        for ion_name in (normalize_text(row.get("ion_name")) for row in mass_library_rows)
        if ion_name is not None
    )
    x_ion_counter = Counter(
        ion_name
        for ion_name in (normalize_text(row.get("ion_name")) for row in mass_library_x_rows)
        if ion_name is not None
    )
    buddy_counter = Counter(
        buddy for buddy in (normalize_text(row.get("buddy")) for row in buddy_rows) if buddy is not None
    )

    top_ion_name, top_ion_count, top_ion_share = counter_stats(ion_counter)
    x_top_ion_name, x_top_ion_count, x_top_ion_share = counter_stats(x_ion_counter)
    top_buddy, top_buddy_count, top_buddy_share = counter_stats(buddy_counter)

    x_row_unique_ion_name_count = len(set(mass_library_x_ion_names))
    all_unique_ion_name_count = len(set(mass_library_ion_names))
    confident_dominant_ion = bool(
        mass_library_x_ion_names
        and x_row_unique_ion_name_count == 1
        and top_ion_name is not None
        and mass_library_x_ion_names[0] == top_ion_name
        and top_ion_share >= 0.5
    )

    c_predict = None
    for row in observed_isotope_rows:
        if str(row.get("element")) == "C" and safe_float(row.get("predictor_value")) is not None:
            c_predict = float(row["predictor_value"])
            break

    return {
        "feature_key": feature_key,
        "feature_label": feature_label,
        "source_html": str(source_path.resolve()),
        "comment_count": len(comment_texts),
        "comment_texts": comment_texts,
        "comment_formulas": comment_formulas,
        "inferred_mass_match_row_count": len(inferred_rows),
        "inferred_mass_x_row_count": len(inferred_x_rows),
        "inferred_mass_formulas": inferred_formulas,
        "inferred_mass_x_formulas": inferred_x_formulas,
        "mass_library_row_count": len(mass_library_rows),
        "mass_library_x_row_count": len(mass_library_x_rows),
        "mass_library_formulas": mass_library_formulas,
        "mass_library_x_formulas": mass_library_x_formulas,
        "mass_library_ion_names": mass_library_ion_names,
        "mass_library_x_ion_names": mass_library_x_ion_names,
        "mass_library_x_ion_masses": mass_library_x_ion_masses,
        "mass_library_ion_name_unique_count": all_unique_ion_name_count,
        "mass_library_x_ion_name_unique_count": x_row_unique_ion_name_count,
        "candidate_formula_row_count": len(candidate_formula_rows),
        "candidate_formula_unique_count": len(candidate_formulas),
        "candidate_formulas": candidate_formulas,
        "candidate_formula_top_weight": candidate_top_weight,
        "candidate_formula_top_formulas": candidate_top_formulas,
        "parent_mass_row_count": len(parent_mass_rows),
        "parent_mass_x_row_count": len(parent_mass_x_rows),
        "parent_mass_formulas": parent_mass_formulas,
        "parent_mass_x_formulas": parent_mass_x_formulas,
        "buddy_row_count": len(buddy_rows),
        "buddy_labels": buddy_labels,
        "non_isotope_buddy_labels": non_isotope_buddy_labels,
        "mass_library_isotope_labels": mass_library_isotope_labels,
        "buddy_isotope_labels": buddy_isotope_labels,
        "isotope_element_hints": isotope_element_hints,
        "feature_is_likely_isotopologue": feature_is_likely_isotopologue,
        "feature_isotopologue_reason": feature_isotopologue_reason,
        "inferred_neutral_masses": sorted(set(inferred_masses)),
        "feature_level_formula_set": feature_level_formula_set,
        "candidate_level_formula_set": candidate_level_formula_set,
        "usable_formula_source": usable_formula_source,
        "usable_formula_set": usable_formula_set,
        "mass_library_top_ion_name": top_ion_name,
        "mass_library_top_ion_count": top_ion_count,
        "mass_library_top_ion_share": top_ion_share,
        "mass_library_x_top_ion_name": x_top_ion_name,
        "mass_library_x_top_ion_count": x_top_ion_count,
        "mass_library_x_top_ion_share": x_top_ion_share,
        "top_buddy_label": top_buddy,
        "top_buddy_count": top_buddy_count,
        "top_buddy_share": top_buddy_share,
        "confident_dominant_ion": confident_dominant_ion,
        "mass_library_ion_name_counter": dict(sorted(ion_counter.items())),
        "mass_library_x_ion_name_counter": dict(sorted(x_ion_counter.items())),
        "buddy_counter": dict(sorted(buddy_counter.items())),
        "observed_isotope_rows": observed_isotope_rows,
        "observed_isotope_rows_json": artifact_observed_isotope_rows(observed_isotope_rows),
        "isotope_signal_count": int(len(observed_isotope_rows)),
        "c_predict": c_predict,
        "parse_status": "ok",
    }
