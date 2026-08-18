#!/usr/bin/env python3
"""Check which of the documented identification rules actually bind.

Read-only, pure standard library, and unlike ``tools/audit_mass_policy.py`` it
runs against inputs committed to this repository rather than a delivered report
bundle, so it can join the CI regenerate-and-diff set.

The rules themselves are written down in ``docs/identification_rules.md``; this
replays the ones that are checkable from ``data/feature_rows/*_neutral_mass_hypotheses.csv``
and reports, per rule, whether the shipped data is consistent with it:

1. Adduct arithmetic. Every hypothesis row names an ion (``2m-H``, ``(m-2H)/2``,
   ``m-H[13C1]``) and stores a neutral mass. Re-deriving the neutral mass from the
   label says whether the label was applied or merely recorded. A ladder whose
   rows all collapse onto ``[M+/-H]`` is not a ladder, and no downstream scorer
   can tell the difference from the score alone.
2. Charge and multimer handling, split out because those are the two label
   families whose error is a clean multiplicative factor rather than a tolerance
   miss, which is what makes them invisible to a mass window.
3. Isotope and halogen evidence carried per feature. The isotope tokens name
   elements the molecular formula must contain, which is the only halogen
   evidence in the pipeline. This counts how much of it exists.

Usage:
    python3 tools/audit_identification_rules.py
    python3 tools/audit_identification_rules.py --feature-rows-dir data/feature_rows
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

PROJECT_ROOT = Path(__file__).resolve().parent.parent

csv.field_size_limit(10**7)

# --- physical constants ------------------------------------------------------

ELECTRON_MASS = 0.000548579909065

# Monoisotopic masses for the elements that appear in the shipped adduct labels.
ELEMENT_MASS = {
    "H": 1.00782503223,
    "C": 12.0,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "Na": 22.9897692820,
    "S": 31.9720711744,
    "Cl": 34.968852682,
    "K": 38.9637064864,
    "Br": 78.9183376,
}

# Mass a heavy isotope adds to the ion relative to the monoisotopic element it
# replaces. An ion carrying k of these is k deltas heavier than the molecule it
# came from, so the underlying neutral mass is that much lighter.
ISOTOPE_DELTA = {
    "2H": 2.01410177812 - ELEMENT_MASS["H"],
    "13C": 13.00335483507 - ELEMENT_MASS["C"],
    "15N": 15.00010889888 - ELEMENT_MASS["N"],
    "18O": 17.99915961286 - ELEMENT_MASS["O"],
    "34S": 33.967867004 - ELEMENT_MASS["S"],
    "37Cl": 36.965902602 - ELEMENT_MASS["Cl"],
    "81Br": 80.9162906 - ELEMENT_MASS["Br"],
}

ISOTOPE_ELEMENT = {
    "2H": "H",
    "13C": "C",
    "15N": "N",
    "18O": "O",
    "34S": "S",
    "37Cl": "Cl",
    "81Br": "Br",
}

HALOGENS = ("F", "Cl", "Br", "I")

# Agreement bands for "did the label get applied". 1 mDa is well inside any
# instrument window; 10 mDa is already twice the pipeline's own 0.005 Da
# tolerance, so anything past it is a different molecule, not a rounding
# difference.
AGREE_DA = 0.001
NEAR_DA = 0.010

ISOTOPE_TOKEN_RE = re.compile(r"\[(\d+)([A-Z][a-z]?)(\d*)\]")
LABEL_PREFIX_RE = re.compile(r"^[A-Za-z_]+:")
CHARGE_SUFFIX_RE = re.compile(r"^\((?P<body>.+)\)/(?P<z>\d+)$")
MULTIMER_RE = re.compile(r"^(?P<n>\d*)m")
# The numeric branch has to come first: a greedy count would eat the "62" of
# "+62.9959" and leave ".9959" as the formula.
TERM_RE = re.compile(
    r"(?P<sign>[+-])(?:(?P<number>\d+\.\d+)|(?P<count>\d*)(?P<token>[A-Za-z][A-Za-z0-9]*))"
)


# --- label parsing -----------------------------------------------------------


@dataclass(frozen=True)
class Adduct:
    """One parsed ion label, in the form needed to invert m/z to a neutral mass.

    ``multiplier`` is the ``n`` of ``nM``; ``charge`` the ``z`` of ``/z``. Both
    default to 1, which is what the great majority of labels mean. ``gain_mass``
    and ``loss_mass`` are neutral atom masses, with the electron accounted for
    separately by polarity, so ``m-H`` in negative mode and ``m+H`` in positive
    mode both come out at the proton mass without either needing a special case.
    """

    label: str
    multiplier: int = 1
    charge: int = 1
    gain_mass: float = 0.0
    loss_mass: float = 0.0
    isotopes: tuple[tuple[str, int], ...] = ()
    named_terms: tuple[str, ...] = ()
    numeric_terms: tuple[float, ...] = ()

    @property
    def isotope_shift(self) -> float:
        return sum(ISOTOPE_DELTA[iso] * k for iso, k in self.isotopes)

    @property
    def is_plain_protonation(self) -> bool:
        """``m-H`` / ``m+H``: one proton, one molecule, one charge, no isotopes."""

        return (
            self.multiplier == 1
            and self.charge == 1
            and not self.isotopes
            and not self.numeric_terms
            and self.named_terms == ("H",)
        )

    @property
    def element_hints(self) -> dict[str, int]:
        """Minimum atom counts the neutral formula must have to explain the label."""

        hints: dict[str, int] = {}
        for iso, k in self.isotopes:
            element = ISOTOPE_ELEMENT[iso]
            hints[element] = max(hints.get(element, 0), int(k))
        return hints


def is_derivationless_label(label: object) -> bool:
    """True for the rows that state a neutral mass without naming an ion.

    ``observed_consensus`` and ``feature_neutral_mass_fallback`` are not failures
    to parse: they carry no derivation to check, so they are excluded from the
    arithmetic rather than counted against it.
    """

    text = str(label or "").strip()
    return not text or text.endswith(("observed_consensus", "fallback"))


def token_mass(token: str) -> float | None:
    """Mass of one adduct term: a bare number, or a small neutral formula."""

    text = token.strip()
    if not text:
        return None
    if re.fullmatch(r"\d+\.\d+", text):
        return float(text)
    total = 0.0
    consumed = 0
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", text):
        if element not in ELEMENT_MASS:
            return None
        total += ELEMENT_MASS[element] * int(count or 1)
        consumed += len(element) + len(count)
    if consumed != len(text):
        return None
    return total


def parse_adduct_label(label: str) -> Adduct | None:
    """Parse a CHRO_LIB ion label into the terms of the neutral-mass inversion.

    Returns ``None`` for labels that name no ion at all -- the consensus and
    fallback rows, which carry a neutral mass with no stated derivation and so
    cannot be checked against one.
    """

    text = str(label or "").strip()
    if not text:
        return None
    text = LABEL_PREFIX_RE.sub("", text)
    if "m" not in text:
        return None

    isotopes: Counter[str] = Counter()

    def strip_isotopes(source: str) -> str:
        for mass_number, element, count in ISOTOPE_TOKEN_RE.findall(source):
            key = f"{mass_number}{element}"
            if key not in ISOTOPE_DELTA:
                return ""
            isotopes[key] += int(count or 1)
        return ISOTOPE_TOKEN_RE.sub("", source)

    stripped = strip_isotopes(text)
    if not stripped:
        return None

    charge = 1
    charged = CHARGE_SUFFIX_RE.match(stripped)
    if charged:
        charge = int(charged.group("z"))
        stripped = charged.group("body")

    head = MULTIMER_RE.match(stripped)
    if head is None:
        return None
    multiplier = int(head.group("n") or 1)
    rest = stripped[head.end() :]

    gain = 0.0
    loss = 0.0
    position = 0
    named: list[str] = []
    numeric: list[float] = []
    for term in TERM_RE.finditer(rest):
        if term.start() != position:
            return None
        position = term.end()
        if term.group("number") is not None:
            mass = float(term.group("number"))
            numeric.append(mass)
        else:
            token = term.group("token")
            unit = token_mass(token)
            if unit is None:
                return None
            mass = unit * int(term.group("count") or 1)
            named.append(token)
        if term.group("sign") == "+":
            gain += mass
        else:
            loss += mass
    if position != len(rest):
        return None

    return Adduct(
        label=text,
        multiplier=multiplier,
        charge=charge,
        gain_mass=gain,
        loss_mass=loss,
        isotopes=tuple(sorted(isotopes.items())),
        named_terms=tuple(named),
        numeric_terms=tuple(numeric),
    )


def neutral_mass_from_adduct(parent_mz: float, adduct: Adduct, ion_mode: str) -> float:
    """Invert an observed m/z to the neutral monoisotopic mass the label implies.

    Forward: ``mz = (n*M + gains - losses + isotope_shift - s*z*m_e) / z``, with
    ``s`` +1 for cations (an electron short) and -1 for anions. Solved for M.
    """

    sign = 1.0 if str(ion_mode).strip().lower().startswith("pos") else -1.0
    ion_mass = adduct.charge * float(parent_mz) + sign * adduct.charge * ELECTRON_MASS
    return (ion_mass - adduct.gain_mass + adduct.loss_mass - adduct.isotope_shift) / adduct.multiplier


def naive_neutral_mass(parent_mz: float, ion_mode: str) -> float:
    """The [M+/-H] neutral mass, which is what the pipeline assumes everywhere."""

    proton = ELEMENT_MASS["H"] - ELECTRON_MASS
    if str(ion_mode).strip().lower().startswith("pos"):
        return float(parent_mz) - proton
    return float(parent_mz) + proton


def label_family(adduct: Adduct | None, label: str) -> str:
    """Group labels by the property whose mishandling would matter.

    Ordered so that the most consequential property wins: a doubly charged
    isotopologue is reported under charge, because that is the term whose error
    is largest.
    """

    if adduct is None:
        return "consensus_or_fallback" if is_derivationless_label(label) else "unparsed"
    if adduct.charge > 1:
        return "multiply_charged"
    if adduct.multiplier > 1:
        return "multimer"
    if adduct.isotopes:
        return "isotopologue"
    if adduct.numeric_terms:
        return "numeric_offset"
    if adduct.is_plain_protonation:
        return "protonated"
    return "adduct_or_cluster"


def merge_element_hints(*sources: str) -> dict[str, int]:
    """Union isotope tokens found anywhere in a row into per-element minimums."""

    hints: dict[str, int] = {}
    for source in sources:
        for mass_number, element, count in ISOTOPE_TOKEN_RE.findall(str(source or "")):
            key = f"{mass_number}{element}"
            if key not in ISOTOPE_ELEMENT:
                continue
            target = ISOTOPE_ELEMENT[key]
            hints[target] = max(hints.get(target, 0), int(count or 1))
    return hints


# --- audit -------------------------------------------------------------------


@dataclass
class RowCheck:
    """One hypothesis row, re-derived."""

    control_id: str
    feature_key: str
    lib_id: str
    ion_mode: str
    label: str
    family: str
    parent_mz: float
    stored_neutral: float
    expected_neutral: float | None
    naive_neutral: float
    declared_charge: int | None
    declared_multiplier: int | None
    parsed_charge: int | None
    parsed_multiplier: int | None

    @property
    def error_da(self) -> float | None:
        if self.expected_neutral is None:
            return None
        return abs(self.stored_neutral - self.expected_neutral)

    @property
    def naive_error_da(self) -> float:
        return abs(self.stored_neutral - self.naive_neutral)


@dataclass
class FeatureEvidence:
    """Isotope evidence attached to one feature across all of its rows."""

    control_id: str
    lib_id: str
    hints: dict[str, int] = field(default_factory=dict)

    def absorb(self, hints: dict[str, int]) -> None:
        for element, count in hints.items():
            self.hints[element] = max(self.hints.get(element, 0), count)

    @property
    def halogen_hints(self) -> dict[str, int]:
        return {e: n for e, n in self.hints.items() if e in HALOGENS}


def read_hypotheses(paths: Iterable[Path]) -> Iterator[dict[str, str]]:
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                row["__source__"] = path.name
                yield row


def as_float(value: object) -> float | None:
    text = str(value if value is not None else "").strip()
    if not text:
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    return parsed


def as_int(value: object) -> int | None:
    parsed = as_float(value)
    if parsed is None:
        return None
    return int(round(parsed))


def audit(paths: list[Path]) -> dict[str, object]:
    checks: list[RowCheck] = []
    features: dict[tuple[str, str], FeatureEvidence] = {}
    round_trip_echo = 0
    round_trip_blank = 0
    round_trip_disagree = 0
    unparsed_labels: Counter[str] = Counter()

    for row in read_hypotheses(paths):
        label = str(row.get("adduct_label") or "")
        ion_mode = str(row.get("ion_mode") or "negative")
        lib_id = str(row.get("lib_id") or "")
        control_id = str(row.get("control_id") or "")
        feature_key = str(row.get("feature_key") or "")
        parent_mz = as_float(row.get("parent_mz"))
        stored = as_float(row.get("neutral_mass"))

        key = (lib_id, feature_key)
        evidence = features.setdefault(key, FeatureEvidence(control_id=control_id, lib_id=lib_id))
        evidence.absorb(merge_element_hints(label, row.get("note", "")))

        recomputed = as_float(row.get("recomputed_mz"))
        if recomputed is None:
            round_trip_blank += 1
        elif parent_mz is not None and abs(recomputed - parent_mz) <= 1e-9:
            round_trip_echo += 1
        else:
            round_trip_disagree += 1

        if parent_mz is None or stored is None:
            continue

        adduct = parse_adduct_label(label)
        family = label_family(adduct, label)
        if family == "unparsed":
            unparsed_labels[label] += 1
        checks.append(
            RowCheck(
                control_id=control_id,
                feature_key=feature_key,
                lib_id=lib_id,
                ion_mode=ion_mode,
                label=label,
                family=family,
                parent_mz=parent_mz,
                stored_neutral=stored,
                expected_neutral=(
                    neutral_mass_from_adduct(parent_mz, adduct, ion_mode) if adduct else None
                ),
                naive_neutral=naive_neutral_mass(parent_mz, ion_mode),
                declared_charge=as_int(row.get("charge")),
                declared_multiplier=as_int(row.get("multiplier")),
                parsed_charge=adduct.charge if adduct else None,
                parsed_multiplier=adduct.multiplier if adduct else None,
            )
        )

    return {
        "inputs": sorted(path.name for path in paths),
        "totals": _totals(checks, features, round_trip_echo, round_trip_blank, round_trip_disagree),
        "adduct_arithmetic": _adduct_arithmetic(checks),
        "charge_and_multimer": _charge_and_multimer(checks),
        "isotope_evidence": _isotope_evidence(checks, features),
        "unparsed_labels": dict(unparsed_labels.most_common(20)),
    }


def _totals(
    checks: list[RowCheck],
    features: dict[tuple[str, str], FeatureEvidence],
    echo: int,
    blank: int,
    disagree: int,
) -> dict[str, object]:
    return {
        "hypothesis_rows": len(checks),
        "features": len(features),
        "controls": len({c.control_id for c in checks if c.control_id}),
        "rows_with_parseable_ion_label": sum(1 for c in checks if c.expected_neutral is not None),
        "recomputed_mz_echoes_parent_mz": echo,
        "recomputed_mz_blank": blank,
        "recomputed_mz_disagrees_with_parent_mz": disagree,
    }


def _band(error: float) -> str:
    if error <= AGREE_DA:
        return "agrees_within_1mda"
    if error <= NEAR_DA:
        return "within_10mda"
    return "over_10mda"


def _adduct_arithmetic(checks: list[RowCheck]) -> dict[str, object]:
    """Does the stored neutral mass follow from the label, or from [M+/-H]?"""

    per_family: dict[str, dict[str, object]] = {}
    grouped: dict[str, list[RowCheck]] = defaultdict(list)
    for check in checks:
        grouped[check.family].append(check)

    for family, rows in sorted(grouped.items()):
        derivable = [r for r in rows if r.error_da is not None]
        bands = Counter(_band(r.error_da) for r in derivable if r.error_da is not None)
        errors = sorted(r.error_da for r in derivable if r.error_da is not None)
        naive_errors = sorted(r.naive_error_da for r in rows)
        per_family[family] = {
            "rows": len(rows),
            "derivable_rows": len(derivable),
            "vs_label": dict(bands),
            "median_error_vs_label_da": round(statistics.median(errors), 6) if errors else None,
            "max_error_vs_label_da": round(errors[-1], 6) if errors else None,
            "matches_plain_protonation_within_10mda": sum(1 for e in naive_errors if e <= NEAR_DA),
        }

    derivable = [r for r in checks if r.error_da is not None]
    follows_label = sum(1 for r in derivable if (r.error_da or 0.0) <= NEAR_DA)
    follows_naive = sum(1 for r in checks if r.naive_error_da <= NEAR_DA)
    return {
        "by_family": per_family,
        "rows_following_their_own_label": follows_label,
        "rows_following_plain_protonation": follows_naive,
        "rows_total": len(checks),
        "share_following_plain_protonation": (
            round(follows_naive / len(checks), 4) if checks else None
        ),
    }


def _charge_and_multimer(checks: list[RowCheck]) -> dict[str, object]:
    """The two label families whose error is a factor, not a tolerance miss.

    A doubly charged precursor read as singly charged is out by roughly the whole
    molecular mass, and a dimer by half of it. No mass window catches either.
    """

    multiply_charged = [c for c in checks if (c.parsed_charge or 1) > 1]
    multimers = [c for c in checks if (c.parsed_multiplier or 1) > 1]

    def summarize(rows: list[RowCheck]) -> dict[str, object]:
        errors = sorted(r.error_da for r in rows if r.error_da is not None)
        return {
            "rows": len(rows),
            "controls": len({r.control_id for r in rows if r.control_id}),
            "stored_mass_follows_label_within_10mda": sum(1 for e in errors if e <= NEAR_DA),
            "median_error_vs_label_da": round(statistics.median(errors), 4) if errors else None,
            "max_error_vs_label_da": round(errors[-1], 4) if errors else None,
        }

    declared_z = Counter(str(c.declared_charge) for c in checks)
    agree_z = sum(
        1
        for c in checks
        if c.declared_charge is not None
        and c.parsed_charge is not None
        and c.declared_charge == c.parsed_charge
    )
    agree_n = sum(
        1
        for c in checks
        if c.declared_multiplier is not None
        and c.parsed_multiplier is not None
        and c.declared_multiplier == c.parsed_multiplier
    )
    return {
        "multiply_charged": summarize(multiply_charged),
        "multimer": summarize(multimers),
        "declared_charge_values": dict(sorted(declared_z.items())),
        "declared_charge_matches_label": agree_z,
        "declared_multiplier_matches_label": agree_n,
    }


def _isotope_evidence(
    checks: list[RowCheck], features: dict[tuple[str, str], FeatureEvidence]
) -> dict[str, object]:
    """How much isotope and halogen evidence the committed data carries."""

    with_hints = {k: v for k, v in features.items() if v.hints}
    by_element: Counter[str] = Counter()
    for evidence in with_hints.values():
        for element in evidence.hints:
            by_element[element] += 1

    halogen_features = {k: v for k, v in with_hints.items() if v.halogen_hints}
    halogen_controls = sorted(
        {v.control_id for v in halogen_features.values() if v.control_id}
    )

    isotopologue_rows = [c for c in checks if c.family == "isotopologue"]
    shift_applied = sum(
        1 for c in isotopologue_rows if c.error_da is not None and c.error_da <= NEAR_DA
    )
    return {
        "features_total": len(features),
        "features_with_isotope_evidence": len(with_hints),
        "features_by_hinted_element": dict(by_element.most_common()),
        "features_with_halogen_evidence": len(halogen_features),
        "halogen_controls": halogen_controls,
        "isotopologue_hypothesis_rows": len(isotopologue_rows),
        "isotopologue_rows_with_shift_applied": shift_applied,
    }


# --- rendering ---------------------------------------------------------------


def render_markdown(result: dict[str, object]) -> str:
    totals = result["totals"]
    arithmetic = result["adduct_arithmetic"]
    charge = result["charge_and_multimer"]
    isotopes = result["isotope_evidence"]

    lines: list[str] = []
    add = lines.append
    add("# Identification-rule audit")
    add("")
    add(
        "Replays the adduct, charge, and isotope arithmetic in "
        "`data/feature_rows/*_neutral_mass_hypotheses.csv` against the rules written down in "
        "[docs/identification_rules.md](../../../identification_rules.md). Generated by "
        "`tools/audit_identification_rules.py`; read-only, and reproducible from committed inputs."
    )
    add("")
    add(
        f"{totals['hypothesis_rows']:,} hypothesis rows over {totals['features']:,} features "
        f"and {totals['controls']:,} controls."
    )
    add("")

    add("## Does a hypothesis follow the ion it names")
    add("")
    add(
        "Each row names an ion and stores a neutral mass. Re-deriving the neutral mass from the "
        "label distinguishes a ladder that was applied from one that was only recorded."
    )
    add("")
    add(
        "| Ion family | Rows | Follows its own label (<=10 mDa) | Median error vs label "
        "| Equals plain [M+/-H] (<=10 mDa) |"
    )
    add("|---|---:|---:|---:|---:|")
    for family, stats in arithmetic["by_family"].items():
        bands = stats["vs_label"]
        follows = bands.get("agrees_within_1mda", 0) + bands.get("within_10mda", 0)
        median = stats["median_error_vs_label_da"]
        median_text = "n/a" if median is None else f"{median:.4f} Da"
        add(
            f"| `{family}` | {stats['rows']:,} | {follows:,} | {median_text} | "
            f"{stats['matches_plain_protonation_within_10mda']:,} |"
        )
    add("")
    share = arithmetic["share_following_plain_protonation"]
    add(
        f"{arithmetic['rows_following_plain_protonation']:,} of {arithmetic['rows_total']:,} rows "
        f"({share:.1%}) carry the neutral mass of a singly charged, monomeric, monoisotopic "
        "protonation, whatever ion their label names."
    )
    add("")

    add("## Charge state and multimers")
    add("")
    add(
        "These two are separated out because their error is a multiplicative factor rather than a "
        "tolerance miss: a doubly charged precursor read as singly charged is out by roughly a whole "
        "molecular mass, a dimer by half of one. No mass window catches either."
    )
    add("")
    add("| Family | Rows | Controls | Stored mass follows the label | Median error vs label |")
    add("|---|---:|---:|---:|---:|")
    for name in ("multiply_charged", "multimer"):
        stats = charge[name]
        median = stats["median_error_vs_label_da"]
        median_text = "n/a" if median is None else f"{median:.2f} Da"
        add(
            f"| `{name}` | {stats['rows']:,} | {stats['controls']:,} | "
            f"{stats['stored_mass_follows_label_within_10mda']:,} | {median_text} |"
        )
    add("")
    add(
        f"The `charge` column agrees with the label on {charge['declared_charge_matches_label']:,} rows "
        f"and `multiplier` on {charge['declared_multiplier_matches_label']:,}, so the metadata is right "
        "and only the mass derived from it is not."
    )
    add("")
    add(
        f"`recomputed_mz` equals `parent_mz` exactly on "
        f"{totals['recomputed_mz_echoes_parent_mz']:,} rows, is blank on "
        f"{totals['recomputed_mz_blank']:,}, and disagrees on "
        f"{totals['recomputed_mz_disagrees_with_parent_mz']:,}. It echoes its input rather than "
        "round-tripping the neutral mass back through the adduct, so it cannot fail."
    )
    add("")

    add("## Isotope and halogen evidence carried per feature")
    add("")
    add(
        "Isotope tokens in the ion labels and in the consensus `note` name elements the molecular "
        "formula has to contain. This is the pipeline's only halogen evidence: chlorine and bromine "
        "are inferred from their isotope peaks and from nothing else."
    )
    add("")
    add(
        f"{isotopes['features_with_isotope_evidence']:,} of {isotopes['features_total']:,} features "
        "carry at least one isotope hint."
    )
    add("")
    add("| Element hinted | Features |")
    add("|---|---:|")
    for element, count in isotopes["features_by_hinted_element"].items():
        add(f"| {element} | {count:,} |")
    add("")
    add(
        f"{isotopes['features_with_halogen_evidence']:,} features across "
        f"{len(isotopes['halogen_controls']):,} controls carry direct halogen evidence: "
        + ", ".join(f"`{c}`" for c in isotopes["halogen_controls"])
        + "."
    )
    add("")
    add(
        f"Of {isotopes['isotopologue_hypothesis_rows']:,} isotopologue hypothesis rows, "
        f"{isotopes['isotopologue_rows_with_shift_applied']:,} have the isotope shift applied to the "
        "stored neutral mass."
    )
    add("")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--feature-rows-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "feature_rows",
        help="Directory holding *_neutral_mass_hypotheses.csv (default: data/feature_rows).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=PROJECT_ROOT / "docs" / "paper" / "data" / "identification",
        help="Where to write audit_identification_rules.json and .md.",
    )
    args = parser.parse_args(argv)

    paths = sorted(args.feature_rows_dir.glob("*_neutral_mass_hypotheses.csv"))
    if not paths:
        print(f"No *_neutral_mass_hypotheses.csv under {args.feature_rows_dir}", file=sys.stderr)
        return 1

    result = audit(paths)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.out_dir / "audit_identification_rules.json"
    md_path = args.out_dir / "audit_identification_rules.md"
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    md_path.write_text(render_markdown(result), encoding="utf-8")
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
