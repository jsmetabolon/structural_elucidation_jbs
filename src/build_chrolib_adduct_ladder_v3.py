#!/usr/bin/env python3
"""Build observed CHRO adduct ladder evidence from parent rows."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

from build_recovered_ms2_from_chrolib_v1 import (
    adduct_token_mass,
    clean_text,
    ion_mode_for_lib,
    is_isotope_parent_name,
    is_parent_row,
    math_isfinite,
    normalize_id,
)
from x_domain_analog_utils import write_csv, write_json


DEFAULT_CHROLIB = Path("ions_with_ri.csv")
DEFAULT_OUT = Path("merged_encoder/data/chrolib_adduct_ladder.csv")
DEFAULT_SUMMARY = Path("merged_encoder/data/chrolib_adduct_ladder_summary.json")
DEFAULT_KNOWN_X = Path("known_X.csv")
CHUNK_SIZE = 500_000

ISOTOPE_SHIFT = {
    "13C": 1.00335483507,
    "18O": 2.004245,
    "15N": 0.9970349,
    "34S": 1.995796,
    "37CL": 1.99705,
    "81BR": 1.9979535,
}
ISOTOPE_TOKEN = re.compile(r"\[(13C|18O|15N|34S|37Cl|81Br)(\d*)\]", re.IGNORECASE)
UNKNOWN_LABEL_TO_XNUM = re.compile(r"Unknown\s*-\s*(\d+)", re.I)
ADDUCT_COLUMNS = [
    "lib_id",
    "ion_mode",
    "comp_id",
    "x_id",
    "adduct_name",
    "parent_mz",
    "implied_neutral_mass",
    "multiplier",
    "charge",
    "weighting",
    "quant_mass",
    "is_isotope",
    "is_neutral_loss",
    "neutral_loss_tokens",
    "adduct_class",
]


def xid_from_text(value: Any) -> str:
    match = re.search(r"(X-\d+)", str(value or ""))
    return match.group(1) if match else ""


def resolve_canonical_control_id(compound_name: Any, report_name: Any = "") -> tuple[str, str, str]:
    for raw in (compound_name, report_name):
        match = UNKNOWN_LABEL_TO_XNUM.search(str(raw or "").strip())
        if match:
            xid = f"X-{match.group(1)}"
            return xid, xid, "unknown_to_x"
    for raw in (compound_name, report_name):
        xid = xid_from_text(raw)
        if xid:
            return xid, xid, "already_x"
    original = str(compound_name or "").strip()
    return "", original, "named" if original else "empty"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrolib-csv", type=Path, default=DEFAULT_CHROLIB)
    parser.add_argument("--known-x", type=Path, default=DEFAULT_KNOWN_X)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--summary-json", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--chunksize", type=int, default=CHUNK_SIZE)
    parser.add_argument("--force", action="store_true", help="Allow overwriting existing output files.")
    return parser.parse_args(argv)


def load_known_x_map(path: Path) -> dict[tuple[str, str], str]:
    if not path.exists():
        return {}
    known = pd.read_csv(path, low_memory=False)
    known.columns = [str(col).strip().strip('"').lstrip("\ufeff") for col in known.columns]
    out: dict[tuple[str, str], str] = {}
    for _, row in known.iterrows():
        lib_id = normalize_id(row.get("LIB_ID"))
        comp_id = normalize_id(row.get("CompID"))
        x_id = clean_text(row.get("Unknown"))
        if lib_id and comp_id and x_id.startswith("X-"):
            out[(lib_id, comp_id)] = x_id
    return out


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(col).strip().strip('"').lstrip("\ufeff") for col in out.columns]
    if "ION_NAME" in out.columns and "NAME" not in out.columns:
        out = out.rename(columns={"ION_NAME": "NAME"})
    for col in ("LIB_ID", "COMP_ID", "PARENT_ID", "NAME", "MASS", "QUANT_MASS", "WEIGHTING"):
        if col not in out.columns:
            out[col] = ""
    return out


def numeric_float(value: Any, default: float = float("nan")) -> float:
    try:
        parsed = float(pd.to_numeric(value, errors="coerce"))
    except Exception:
        return default
    return parsed if math_isfinite(parsed) else default


def strip_isotopes(text: str) -> tuple[str, float]:
    shift = 0.0

    def repl(match: re.Match[str]) -> str:
        nonlocal shift
        token = match.group(1).upper()
        count = int(match.group(2) or "1")
        shift += ISOTOPE_SHIFT.get(token, 0.0) * count
        return ""

    return ISOTOPE_TOKEN.sub(repl, text), float(shift)


def split_body(text: str) -> tuple[int, int, str, float] | None:
    body, isotope_shift = strip_isotopes(clean_text(text).replace(" ", ""))
    if not body:
        return None
    charge = 1
    wrapped = re.fullmatch(r"\((.+)\)/(\d+)", body, flags=re.IGNORECASE)
    if wrapped:
        body = wrapped.group(1)
        charge = int(wrapped.group(2))

    match = re.match(r"(?:(\d+)[mM]|[mM])", body)
    if not match:
        return None
    multiplier = int(match.group(1) or "1")
    rest = body[match.end() :]
    return multiplier, charge, rest, isotope_shift


def parse_adduct_expression_v3(name: Any) -> tuple[int, int, float, float] | None:
    split = split_body(clean_text(name))
    if split is None:
        return None
    multiplier, charge, rest, isotope_shift = split
    offset = 0.0
    tokens = re.findall(r"([+-])([^+-]+)", rest)
    reconstructed = "".join(sign + token for sign, token in tokens)
    if reconstructed != rest:
        return None
    for sign, token in tokens:
        mass = adduct_token_mass(token)
        if not math_isfinite(mass):
            return None
        offset += mass if sign == "+" else -mass
    return multiplier, charge, float(offset), float(isotope_shift)


def source_neutral_mass_v3(parent_mz: float, adduct_name: Any) -> float:
    parsed = parse_adduct_expression_v3(adduct_name)
    if parsed is None or not math_isfinite(parent_mz):
        return float("nan")
    multiplier, charge, offset, isotope_shift = parsed
    neutral = (float(parent_mz) * float(charge) - float(offset) - float(isotope_shift)) / float(multiplier)
    return neutral if neutral > 0 and math_isfinite(neutral) else float("nan")


def neutral_loss_tokens(adduct_name: Any) -> list[str]:
    split = split_body(clean_text(adduct_name))
    if split is None:
        return []
    _, _, rest, _ = split
    losses: list[str] = []
    for sign, token in re.findall(r"([+-])([^+-]+)", rest):
        token_clean = clean_text(token)
        if sign != "-":
            continue
        if token_clean.upper() in {"H", "2H"}:
            continue
        losses.append(token_clean)
    return losses


def adduct_class(adduct_name: str, multiplier: int, is_isotope: bool, losses: list[str]) -> str:
    upper = adduct_name.upper()
    if is_isotope:
        return "isotopologue"
    if losses:
        return "loss"
    if int(multiplier) > 1:
        return "dimer"
    if "+NA" in upper or "+K" in upper:
        return "alkali"
    if "+NH4" in upper:
        return "ammonium"
    if "CH3OH" in upper or "CH3CN" in upper:
        return "solvent"
    if upper in {"M+H", "M-H"}:
        return "primary"
    return "other"


def row_x_id(row: pd.Series, known_map: dict[tuple[str, str], str]) -> str:
    lib_id = normalize_id(row.get("LIB_ID"))
    comp_id = normalize_id(row.get("COMP_ID"))
    known = known_map.get((lib_id, comp_id), "")
    if known:
        return known
    canonical = clean_text(row.get("CANONICAL_CONTROL_ID"))
    if canonical.startswith("X-"):
        return canonical
    inferred, _, _ = resolve_canonical_control_id(row.get("COMPOUND_NAME"), row.get("REPORT_NAME"))
    return inferred if inferred.startswith("X-") else ""


def build_ladder(chrolib_csv: Path, known_map: dict[tuple[str, str], str], chunksize: int) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for chunk in pd.read_csv(chrolib_csv, chunksize=int(chunksize), dtype=str, low_memory=False):
        chunk = normalize_columns(chunk)
        parents = chunk.loc[chunk["PARENT_ID"].map(is_parent_row)].copy()
        if parents.empty:
            continue
        parents = parents[parents["NAME"].map(clean_text).ne("")].copy()
        if parents.empty:
            continue
        for _, row in parents.iterrows():
            lib_raw = normalize_id(row.get("LIB_ID"))
            comp_id = normalize_id(row.get("COMP_ID"))
            if not lib_raw or not comp_id:
                continue
            lib_id = int(lib_raw)
            mode = ion_mode_for_lib(lib_id)
            if mode is None:
                continue
            _, ion_mode_label = mode
            adduct_name = clean_text(row.get("NAME"))
            parent_mz = numeric_float(row.get("MASS"))
            parsed = parse_adduct_expression_v3(adduct_name)
            if parsed is None or not math_isfinite(parent_mz):
                multiplier, charge = 0, 0
                neutral = float("nan")
            else:
                multiplier, charge, _, _ = parsed
                neutral = source_neutral_mass_v3(parent_mz, adduct_name)
            losses = neutral_loss_tokens(adduct_name)
            isotope = bool(is_isotope_parent_name(adduct_name))
            rows.append(
                {
                    "lib_id": int(lib_id),
                    "ion_mode": ion_mode_label,
                    "comp_id": comp_id,
                    "x_id": row_x_id(row, known_map),
                    "adduct_name": adduct_name,
                    "parent_mz": parent_mz,
                    "implied_neutral_mass": neutral,
                    "multiplier": int(multiplier),
                    "charge": int(charge),
                    "weighting": clean_text(row.get("WEIGHTING")),
                    "quant_mass": clean_text(row.get("QUANT_MASS")),
                    "is_isotope": bool(isotope),
                    "is_neutral_loss": bool(losses),
                    "neutral_loss_tokens": "|".join(losses),
                    "adduct_class": adduct_class(adduct_name, int(multiplier), bool(isotope), losses),
                }
            )
    out = pd.DataFrame(rows, columns=ADDUCT_COLUMNS)
    if out.empty:
        return out
    out = out.drop_duplicates(["lib_id", "comp_id", "adduct_name", "parent_mz"], keep="first")
    return out.sort_values(["lib_id", "comp_id", "adduct_class", "adduct_name"], kind="mergesort").reset_index(drop=True)


def run(args: argparse.Namespace) -> dict[str, Any]:
    if not args.force:
        existing = [path for path in (args.out, args.summary_json) if path.exists()]
        if existing:
            joined = "\n".join(f"  {path}" for path in existing)
            raise FileExistsError(f"Refusing to overwrite existing output files without --force:\n{joined}")
    known_map = load_known_x_map(args.known_x)
    ladder = build_ladder(args.chrolib_csv, known_map, int(args.chunksize))
    write_csv(args.out, ladder)
    payload = {
        "chrolib_csv": str(args.chrolib_csv),
        "out_csv": str(args.out),
        "n_rows": int(ladder.shape[0]),
        "n_comp_lib_pairs": int(ladder[["lib_id", "comp_id"]].drop_duplicates().shape[0]) if not ladder.empty else 0,
        "adduct_class_counts": ladder["adduct_class"].value_counts().to_dict() if not ladder.empty else {},
        "n_neutral_loss_rows": int(ladder["is_neutral_loss"].astype(bool).sum()) if not ladder.empty else 0,
        "n_isotope_rows": int(ladder["is_isotope"].astype(bool).sum()) if not ladder.empty else 0,
    }
    write_json(args.summary_json, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
