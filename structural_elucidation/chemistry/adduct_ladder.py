#!/usr/bin/env python3
"""Canonical adduct ladder for phase-0 and analog candidate mass matching.

Loads ``data/adducts/adduct_ladder.csv`` once and inverts precursor m/z to a
neutral mass under each tabulated hypothesis. Carbon-13 overlays stay in
``phase0.best_adduct_neutral_mass_match``; they are not extra CSV rows.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from structural_elucidation.paths import REPO_ROOT

LADDER_CSV: Path = REPO_ROOT / "data" / "adducts" / "adduct_ladder.csv"
VALID_ION_MODES = frozenset({"negative", "positive"})
FORMULA_KINDS = frozenset({"offset", "multiply_charged", "dimer"})


@dataclass(frozen=True)
class AdductRow:
    """One row of the committed adduct ladder."""

    ion_mode: str
    adduct_label: str
    chro_lib_alias: str
    charge: int
    n_molecules: int
    formula_kind: str
    mz_minus_M: float
    is_default: bool
    source: str
    platform_hint: str
    notes: str
    csv_order: int


def _normalize_ion_mode(ion_mode: str) -> str:
    mode = str(ion_mode).strip().lower()
    if mode not in VALID_ION_MODES:
        raise ValueError(
            f"Unsupported ion mode {ion_mode!r}. "
            f"Expected one of {sorted(VALID_ION_MODES)}."
        )
    return mode


def _parse_bool(value: str) -> bool:
    text = str(value).strip().upper()
    if text in {"TRUE", "1", "YES"}:
        return True
    if text in {"FALSE", "0", "NO", ""}:
        return False
    raise ValueError(f"Invalid boolean {value!r}")


@lru_cache(maxsize=1)
def load_ladder() -> tuple[AdductRow, ...]:
    """Load and cache the canonical adduct table."""
    if not LADDER_CSV.is_file():
        raise FileNotFoundError(f"Adduct ladder CSV not found: {LADDER_CSV}")
    rows: list[AdductRow] = []
    with LADDER_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "ion_mode",
            "adduct_label",
            "chro_lib_alias",
            "charge",
            "n_molecules",
            "formula_kind",
            "mz_minus_M",
            "is_default",
            "source",
            "platform_hint",
            "notes",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{LADDER_CSV} missing columns: {sorted(missing)}"
            )
        for csv_order, raw in enumerate(reader):
            kind = str(raw["formula_kind"]).strip()
            if kind not in FORMULA_KINDS:
                raise ValueError(
                    f"Unknown formula_kind {kind!r} in {LADDER_CSV}"
                )
            rows.append(
                AdductRow(
                    ion_mode=_normalize_ion_mode(raw["ion_mode"]),
                    adduct_label=str(raw["adduct_label"]).strip(),
                    chro_lib_alias=str(
                        raw.get("chro_lib_alias") or ""
                    ).strip(),
                    charge=int(raw["charge"]),
                    n_molecules=int(raw["n_molecules"]),
                    formula_kind=kind,
                    mz_minus_M=float(raw["mz_minus_M"]),
                    is_default=_parse_bool(raw["is_default"]),
                    source=str(raw["source"]).strip(),
                    platform_hint=str(raw.get("platform_hint") or "").strip(),
                    notes=str(raw.get("notes") or "").strip(),
                    csv_order=csv_order,
                )
            )
    return tuple(rows)


def implied_neutral_mass(row: AdductRow, parent_mz: float) -> float:
    """Invert precursor m/z to a neutral mass for one ladder row."""
    mz = float(parent_mz)
    if row.formula_kind == "offset":
        return mz - row.mz_minus_M
    if row.formula_kind == "multiply_charged":
        return float(row.charge) * mz - row.mz_minus_M
    if row.formula_kind == "dimer":
        return (mz - row.mz_minus_M) / float(row.n_molecules)
    raise ValueError(f"Unknown formula_kind {row.formula_kind!r}")


NEG_ONLY_LIBS = frozenset({209, 305})
SKIP_DERIVE_LABELS = frozenset(
    {"", "neutral_mass_fallback", "nce_snapshot_control", "consensus"}
)
LABEL_ALIASES = {
    "m+h+nh3": "[M+NH4]+",
    "m+nh4": "[M+NH4]+",
    "m+h+ch3oh": "[M+CH3OH+H]+",
    "m+chooh-h": "[M+HCOO]-",
    "m+fa-h": "[M+HCOO]-",
    "m+hcoo": "[M+HCOO]-",
    "m+cl35": "[M+Cl]-",
    "m+cl": "[M+Cl]-",
    "m+cl37": "[M+37Cl]-",
    "m+37cl": "[M+37Cl]-",
    "m-e": "[M]+",
    "m+na": "[M+Na]+",
    "m+k": "[M+K]+",
}


def include_polar_for_lib(lib_id: int) -> bool:
    """Polar-hinted NEG adducts (Cl/formate) are unusual on NEG-only 209/305."""
    return int(lib_id) not in NEG_ONLY_LIBS


def _fold_label(label: str) -> str:
    text = str(label or "").strip()
    if ":" in text:
        text = text.split(":", 1)[-1].strip()
    text = text.replace("[", "").replace("]", "").replace(" ", "").lower()
    return text.rstrip("+-")


def selected_rows(ion_mode: str, *, include_polar: bool = True) -> list[AdductRow]:
    """Ladder rows for one polarity, default ``[M±H]`` first."""
    mode = _normalize_ion_mode(ion_mode)
    selected = [row for row in load_ladder() if row.ion_mode == mode]
    if not include_polar:
        selected = [row for row in selected if row.platform_hint != "polar"]
    selected.sort(key=lambda row: (not row.is_default, row.csv_order))
    return selected


def hypotheses(
    parent_mz: float,
    ion_mode: str,
    *,
    include_polar: bool = True,
) -> list[tuple[str, float]]:
    """Return ``(adduct_label, implied_neutral_mass)`` for one polarity.

    Default ``[M±H]`` is first; remaining rows keep CSV order. Negative mode
    includes NEG and POLAR unless ``include_polar=False``.
    """
    return [
        (row.adduct_label, implied_neutral_mass(row, parent_mz))
        for row in selected_rows(ion_mode, include_polar=include_polar)
    ]


def row_for_observed_label(label: str, ion_mode: str) -> AdductRow | None:
    """Map a CHRO_LIB / Excel-style ion label onto one ladder row, if any."""
    raw = str(label or "").strip()
    if _fold_label(raw) in SKIP_DERIVE_LABELS or raw.lower() in SKIP_DERIVE_LABELS:
        return None
    mode = _normalize_ion_mode(ion_mode)
    folded = _fold_label(raw)
    alias = LABEL_ALIASES.get(folded) or LABEL_ALIASES.get(raw.lower())
    want = {_fold_label(alias)} if alias else set()
    want.add(folded)
    for row in load_ladder():
        if row.ion_mode != mode:
            continue
        keys = {_fold_label(row.adduct_label)}
        if row.chro_lib_alias:
            keys.add(_fold_label(row.chro_lib_alias))
        if keys & want:
            return row
    return None


def keeps_search_mass(adduct_label: str, ion_mode: str) -> bool:
    """Keep a hypothesis as a 10 ppm search target.

    Non-default ladder adducts move the neutral mass off the protonated ion.
    Those rows are not search targets. Unmapped labels and ``[M±H]`` stay.
    """

    try:
        row = row_for_observed_label(adduct_label, ion_mode)
    except ValueError:
        return True
    if row is None:
        return True
    return bool(row.is_default)


def derive_target_neutral_mass(
    parent_mz: float,
    ion_mode: str,
    adduct_label: str,
    stored_mass: float | None,
) -> tuple[float | None, str, float | None]:
    """Invert ``parent_mz`` from a mapped ladder row; else keep ``stored_mass``.

    Default ``[M±H]`` rows keep the stored mass so CHG-2 protonated values do
    not move. Returns ``(mass, source, shift_da)`` with source ``label`` or
    ``stored``.
    """
    stored = None if stored_mass is None else float(stored_mass)
    row = row_for_observed_label(adduct_label, ion_mode)
    if row is None or row.is_default:
        return stored, "stored", 0.0 if stored is not None else None
    derived = implied_neutral_mass(row, parent_mz)
    shift = None if stored is None else derived - stored
    return derived, "label", shift
