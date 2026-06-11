#!/usr/bin/env python3
"""Convert CHRO_LIB-style recovered MS2 CSV into ions_unknown_merged peak schema."""

from __future__ import annotations

import argparse
import json
import re
import shutil
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from x_domain_analog_utils import write_csv, write_json

DEFAULT_CHROLIB = Path("previously_no_ms2_metabolites.csv")
DEFAULT_KNOWN_X = Path("known_X.csv")
DEFAULT_OUT = Path("merged_encoder/data/ions_recovered_chrolib_v1.csv")
DEFAULT_AUDIT = Path("merged_encoder/data/ions_recovered_chrolib_v1_audit.csv")
DEFAULT_MERGED = Path("merged_encoder/data/ions_unknown_merged.csv")

NEGATIVE_LIBS = {209, 305}
POSITIVE_LIBS = {400, 402}
SUPPORTED_LIBS = NEGATIVE_LIBS | POSITIVE_LIBS
PROTON_MASS = 1.007276466812
ADDUCT_MASS = {
    "H": PROTON_MASS,
    "NA": 22.989218,
    "K": 38.9632,
    "NH4": 18.03382,
}
ATOM_MASS = {
    "C": 12.0,
    "H": 1.00782503223,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "P": 30.97376199842,
    "S": 31.9720711744,
    "F": 18.99840316273,
    "Cl": 34.968852682,
    "Br": 78.9183376,
    "I": 126.9044719,
}
ISOTOPE_PATTERN = re.compile(r"(13C|C13|18O|O18|15N|34S|37Cl|Cl37|35Cl|Cl35|81Br|Br81|D-1|\[D|\[2H)", re.IGNORECASE)

PEAK_COLUMNS = [
    "comp_id",
    "mass",
    "intensity",
    "ri",
    "rt",
    "compound_name",
    "parent_mz",
    "feature_key",
    "lib_id",
    "ion_mode",
    "ion_mode_label",
    "dataset_split",
    "chrolib_parent_id",
    "source_parent_adduct",
    "source_neutral_mass",
    "recovered_ms2_version",
]

AUDIT_COLUMNS = [
    "lib_id",
    "comp_id",
    "x_id",
    "chrolib_parent_id",
    "feature_key",
    "status",
    "parent_adduct",
    "parent_mz",
    "source_neutral_mass",
    "n_fragments",
    "n_fragments_trimmed",
    "quant_mass",
    "weighting",
    "parent_scope",
    "recovered_ms2_version",
    "skip_reason",
]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrolib-csv", type=Path, default=DEFAULT_CHROLIB)
    parser.add_argument("--known-x", type=Path, default=DEFAULT_KNOWN_X)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--audit-out", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--merged-out", type=Path, default=DEFAULT_MERGED)
    parser.add_argument("--merge-base", type=Path, default=None, help="Optional CSV to use as merge base instead of current --merged-out.")
    parser.add_argument("--raw-out", type=Path, default=None, help="Optional path for untrimmed recovered peaks.")
    parser.add_argument("--parent-scope", choices=("primary", "non_isotope", "all_mapped"), default="primary")
    parser.add_argument("--recovered-ms2-version", default="", help="Metadata version label; defaults to v1 for primary and v2 otherwise.")
    parser.add_argument("--trim-output", action="store_true", help="Apply NCE-style consensus trimming before writing --out and merging.")
    parser.add_argument("--min-relative-intensity", type=float, default=0.01)
    parser.add_argument("--max-consensus-peaks", type=int, default=80)
    parser.add_argument(
        "--merge-strategy",
        choices=("append", "sidecar"),
        default="append",
        help="append: backup merged-out and concat recovered peaks; sidecar: only write --out",
    )
    parser.add_argument("--summary-json", type=Path, default=None)
    return parser.parse_args(argv)


def clean_text(value: Any) -> str:
    text = "" if value is None else str(value)
    return "" if text.lower() == "nan" else text.strip()


def normalize_id(value: Any) -> str:
    text = clean_text(value)
    if not text:
        return ""
    try:
        numeric = float(text)
        if numeric.is_integer():
            return str(int(numeric))
    except ValueError:
        pass
    return text


def is_parent_row(parent_id: Any) -> bool:
    return clean_text(parent_id) == ""


def load_known_x_map(path: Path) -> dict[tuple[str, str], dict[str, Any]]:
    known = pd.read_csv(path, low_memory=False)
    lib_col = "LIB_ID"
    if lib_col not in known.columns:
        for col in known.columns:
            if col.replace("\ufeff", "") == "LIB_ID":
                lib_col = col
                break
    out: dict[tuple[str, str], dict[str, Any]] = {}
    for _, row in known.iterrows():
        lib_id = normalize_id(row.get(lib_col))
        comp_id = normalize_id(row.get("CompID"))
        x_id = clean_text(row.get("Unknown"))
        if not lib_id or not comp_id or not x_id.startswith("X-"):
            continue
        out[(lib_id, comp_id)] = {
            "x_id": x_id,
            "level": clean_text(row.get("Level")),
            "sublevel": clean_text(row.get("SubLevel")),
        }
    return out


def ion_mode_for_lib(lib_id: int) -> tuple[float, str] | None:
    if lib_id in NEGATIVE_LIBS:
        return 0.0, "negative"
    if lib_id in POSITIVE_LIBS:
        return 1.0, "positive"
    return None


def primary_adduct_name(lib_id: int) -> str:
    return "m-H" if lib_id in NEGATIVE_LIBS else "m+H"


def select_parent(parent_candidates: pd.DataFrame) -> pd.Series | None:
    if parent_candidates.empty:
        return None
    ranked = parent_candidates.copy()
    ranked["_quant"] = ranked["QUANT_MASS"].map(lambda v: 1 if clean_text(v).upper() == "Y" else 0)
    ranked["_weight"] = pd.to_numeric(ranked["WEIGHTING"], errors="coerce").fillna(-1.0)
    ranked["_mass_id"] = ranked["CHRO_LIB_MASSES_ID"].map(clean_text)
    ranked = ranked.sort_values(
        by=["_quant", "_weight", "_mass_id"],
        ascending=[False, False, True],
    )
    return ranked.iloc[0]


def build_feature_key(lib_id: int, parent_mz: float, x_id: str) -> str:
    return f"feature_X_{lib_id}_mz{parent_mz:.4f}_ri0_{x_id}"


def build_v2_feature_key(lib_id: int, parent_id: str, parent_mz: float, x_id: str) -> str:
    return f"feature_X_{lib_id}_chrolib{parent_id}_mz{parent_mz:.4f}_ri0_{x_id}"


def is_isotope_parent_name(name: Any) -> bool:
    return bool(ISOTOPE_PATTERN.search(clean_text(name)))


def formula_mass(token: str) -> float:
    text = clean_text(token)
    if not text:
        return float("nan")
    total = 0.0
    pos = 0
    for match in re.finditer(r"([A-Z][a-z]?)(\d*)", text):
        if match.start() != pos:
            return float("nan")
        elem = match.group(1)
        if elem not in ATOM_MASS:
            return float("nan")
        count = int(match.group(2) or "1")
        total += ATOM_MASS[elem] * count
        pos = match.end()
    return total if pos == len(text) else float("nan")


def adduct_token_mass(token: str) -> float:
    text = clean_text(token)
    if not text:
        return float("nan")
    counted = re.fullmatch(r"(\d+)([A-Za-z].*)", text)
    if counted:
        mass = adduct_token_mass(counted.group(2))
        return float(int(counted.group(1)) * mass) if math_isfinite(mass) else float("nan")
    upper = text.upper()
    if upper in ADDUCT_MASS:
        return float(ADDUCT_MASS[upper])
    try:
        return float(text)
    except ValueError:
        pass
    return formula_mass(text)


def parse_offset(rest: str) -> float:
    if not rest:
        return 0.0
    offset = 0.0
    for sign, token in re.findall(r"([+-])([^+-]+)", rest):
        mass = adduct_token_mass(token)
        if not math_isfinite(mass):
            return float("nan")
        offset += mass if sign == "+" else -mass
    reconstructed = "".join(sign + token for sign, token in re.findall(r"([+-])([^+-]+)", rest))
    if reconstructed != rest:
        return float("nan")
    return float(offset)


def parse_adduct_expression(name: Any) -> tuple[int, int, float] | None:
    text = clean_text(name).replace(" ", "")
    if not text:
        return None
    body = text
    charge = 1
    wrapped = re.fullmatch(r"\((.+)\)/(\d+)", text, flags=re.IGNORECASE)
    if wrapped:
        body = wrapped.group(1)
        charge = int(wrapped.group(2))
    multiplier = 1
    lower = body.lower()
    if lower.startswith("2m"):
        multiplier = 2
        rest = body[2:]
    elif lower.startswith("m"):
        rest = body[1:]
    elif body.startswith("M"):
        rest = body[1:]
    else:
        return None
    offset = parse_offset(rest)
    if not math_isfinite(offset) or multiplier <= 0 or charge <= 0:
        return None
    return multiplier, charge, offset


def source_neutral_mass(parent_mz: float, adduct_name: Any) -> float:
    parsed = parse_adduct_expression(adduct_name)
    if parsed is None or not math_isfinite(parent_mz):
        return float("nan")
    multiplier, charge, offset = parsed
    neutral = (float(parent_mz) * float(charge) - float(offset)) / float(multiplier)
    return neutral if neutral > 0 and math_isfinite(neutral) else float("nan")


def base_audit(
    *,
    lib_id: int,
    comp_id: str,
    x_id: str = "",
    parent: pd.Series | None = None,
    parent_scope: str,
    version: str,
    skip_reason: str = "",
) -> dict[str, Any]:
    parent_id = clean_text(parent.get("CHRO_LIB_MASSES_ID")) if parent is not None else ""
    parent_adduct = clean_text(parent.get("NAME")) if parent is not None else ""
    parent_mz = clean_text(parent.get("MASS")) if parent is not None else ""
    neutral = ""
    if parent is not None:
        mz = float(pd.to_numeric(parent.get("MASS"), errors="coerce"))
        parsed = source_neutral_mass(mz, parent_adduct)
        neutral = parsed if math_isfinite(parsed) else ""
    return {
        "lib_id": lib_id,
        "comp_id": comp_id,
        "x_id": x_id,
        "chrolib_parent_id": parent_id,
        "feature_key": "",
        "status": "skipped",
        "parent_adduct": parent_adduct,
        "parent_mz": parent_mz,
        "source_neutral_mass": neutral,
        "n_fragments": 0,
        "n_fragments_trimmed": 0,
        "quant_mass": clean_text(parent.get("QUANT_MASS")) if parent is not None else "",
        "weighting": clean_text(parent.get("WEIGHTING")) if parent is not None else "",
        "parent_scope": parent_scope,
        "recovered_ms2_version": version,
        "skip_reason": skip_reason,
    }


def convert_group(
    group: pd.DataFrame,
    lib_id: int,
    comp_id: str,
    known_map: dict[tuple[str, str], dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    peaks, audits = convert_group_parents(
        group,
        lib_id,
        comp_id,
        known_map,
        parent_scope="primary",
        version="v1",
    )
    return peaks, audits[0] if audits else base_audit(lib_id=lib_id, comp_id=comp_id, parent_scope="primary", version="v1")


def candidate_parents(group: pd.DataFrame, lib_id: int, parent_scope: str) -> tuple[pd.DataFrame, str]:
    parents = group.loc[group["PARENT_ID"].map(is_parent_row)].copy()
    if parent_scope == "primary":
        adduct = primary_adduct_name(lib_id)
        candidates = parents.loc[parents["NAME"].map(clean_text) == adduct]
        parent = select_parent(candidates)
        if parent is None:
            return pd.DataFrame(), "no_primary_adduct"
        return parent.to_frame().T.reset_index(drop=True), ""
    if parent_scope == "non_isotope":
        return parents.reset_index(drop=True), ""
    if parent_scope == "all_mapped":
        return parents.reset_index(drop=True), ""
    raise ValueError(f"Unsupported parent_scope={parent_scope}")


def convert_parent(
    parent: pd.Series,
    group: pd.DataFrame,
    lib_id: int,
    comp_id: str,
    x_id: str,
    parent_scope: str,
    version: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    audit = base_audit(
        lib_id=lib_id,
        comp_id=comp_id,
        x_id=x_id,
        parent=parent,
        parent_scope=parent_scope,
        version=version,
    )

    mode = ion_mode_for_lib(lib_id)
    if mode is None:
        audit["skip_reason"] = "unsupported_lib"
        return [], audit

    ion_mode, ion_mode_label = mode
    parent_mass_id = clean_text(parent["CHRO_LIB_MASSES_ID"])
    parent_mz = float(pd.to_numeric(parent["MASS"], errors="coerce"))
    if not math_isfinite(parent_mz):
        audit["skip_reason"] = "invalid_parent_mass"
        return [], audit

    fragments = group.loc[group["PARENT_ID"].map(clean_text) == parent_mass_id].copy()
    if fragments.empty:
        audit["skip_reason"] = "no_fragments"
        return [], audit

    neutral = source_neutral_mass(parent_mz, parent.get("NAME"))
    if not math_isfinite(neutral):
        audit["skip_reason"] = "neutral_mass_parse_failed"
        return [], audit
    audit["source_neutral_mass"] = neutral

    feature_key = build_feature_key(lib_id, parent_mz, x_id) if version == "v1" else build_v2_feature_key(lib_id, parent_mass_id, parent_mz, x_id)
    audit["feature_key"] = feature_key
    peaks: list[dict[str, Any]] = []
    for _, frag in fragments.iterrows():
        mass = float(pd.to_numeric(frag["MASS"], errors="coerce"))
        intensity = float(pd.to_numeric(frag["WEIGHTING"], errors="coerce"))
        if not math_isfinite(mass) or not math_isfinite(intensity):
            continue
        peaks.append(
            {
                "comp_id": str(comp_id),
                "mass": mass,
                "intensity": intensity,
                "ri": 0.0,
                "rt": 0.0,
                "compound_name": x_id,
                "parent_mz": parent_mz,
                "feature_key": feature_key,
                "lib_id": lib_id,
                "ion_mode": ion_mode,
                "ion_mode_label": ion_mode_label,
                "dataset_split": "unknown",
                "chrolib_parent_id": parent_mass_id,
                "source_parent_adduct": clean_text(parent.get("NAME")),
                "source_neutral_mass": neutral,
                "recovered_ms2_version": version,
            }
        )

    if not peaks:
        audit["skip_reason"] = "no_valid_fragment_rows"
        return [], audit

    audit.update(
        {
            "status": "emitted",
            "parent_adduct": clean_text(parent.get("NAME")),
            "parent_mz": parent_mz,
            "n_fragments": len(peaks),
            "n_fragments_trimmed": len(peaks),
            "quant_mass": clean_text(parent["QUANT_MASS"]),
            "weighting": clean_text(parent["WEIGHTING"]),
            "skip_reason": "",
        }
    )
    return peaks, audit


def convert_group_parents(
    group: pd.DataFrame,
    lib_id: int,
    comp_id: str,
    known_map: dict[tuple[str, str], dict[str, Any]],
    *,
    parent_scope: str,
    version: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    mode = ion_mode_for_lib(lib_id)
    if mode is None:
        return [], [base_audit(lib_id=lib_id, comp_id=comp_id, parent_scope=parent_scope, version=version, skip_reason="unsupported_lib")]

    parents, missing_reason = candidate_parents(group, lib_id, parent_scope)
    if parents.empty:
        return [], [base_audit(lib_id=lib_id, comp_id=comp_id, parent_scope=parent_scope, version=version, skip_reason=missing_reason or "no_parent_rows")]

    peaks: list[dict[str, Any]] = []
    audits: list[dict[str, Any]] = []
    meta = known_map.get((normalize_id(lib_id), normalize_id(comp_id)))
    for _, parent in parents.iterrows():
        audit = base_audit(lib_id=lib_id, comp_id=comp_id, parent=parent, parent_scope=parent_scope, version=version)
        if parent_scope == "non_isotope" and is_isotope_parent_name(parent.get("NAME")):
            audit["skip_reason"] = "isotope_parent_excluded"
            audits.append(audit)
            continue
        if meta is None:
            audit["skip_reason"] = "not_in_known_X"
            audits.append(audit)
            continue
        parent_peaks, parent_audit = convert_parent(
            parent,
            group,
            lib_id,
            comp_id,
            str(meta["x_id"]),
            parent_scope,
            version,
        )
        peaks.extend(parent_peaks)
        audits.append(parent_audit)
    return peaks, audits


def math_isfinite(value: float) -> bool:
    return value == value and abs(value) != float("inf")


def load_chrolib(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False, dtype=str)


def convert_chrolib(
    chrolib: pd.DataFrame,
    known_map: dict[tuple[str, str], dict[str, Any]],
    *,
    parent_scope: str = "primary",
    version: str = "v1",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    peaks: list[dict[str, Any]] = []
    audits: list[dict[str, Any]] = []

    grouped = chrolib.groupby(["LIB_ID", "COMP_ID"], sort=False)
    for (lib_raw, comp_raw), group in grouped:
        lib_id = int(normalize_id(lib_raw))
        comp_id = normalize_id(comp_raw)
        group_peaks, group_audits = convert_group_parents(
            group,
            lib_id,
            comp_id,
            known_map,
            parent_scope=parent_scope,
            version=version,
        )
        peaks.extend(group_peaks)
        audits.extend(group_audits)

    peak_df = pd.DataFrame(peaks, columns=PEAK_COLUMNS) if peaks else pd.DataFrame(columns=PEAK_COLUMNS)
    audit_df = pd.DataFrame(audits, columns=AUDIT_COLUMNS)
    return peak_df, audit_df


def existing_x_ids(merged_path: Path) -> set[str]:
    if not merged_path.exists():
        return set()
    merged = pd.read_csv(merged_path, usecols=["compound_name"], low_memory=False)
    return {clean_text(v) for v in merged["compound_name"] if clean_text(v).startswith("X-")}


def merge_peaks(existing: pd.DataFrame, recovered: pd.DataFrame) -> pd.DataFrame:
    if existing.empty:
        return recovered.copy()
    if recovered.empty:
        return existing.copy()
    combined = pd.concat([existing, recovered], ignore_index=True)
    dedup_keys = ["feature_key", "mass"] if "feature_key" in combined.columns else ["comp_id", "lib_id", "parent_mz", "mass"]
    for key in dedup_keys:
        if key not in combined.columns:
            combined[key] = ""
    before = len(combined)
    combined = combined.drop_duplicates(subset=dedup_keys, keep="first")
    combined.attrs["dedup_dropped"] = before - len(combined)
    return combined


def trim_recovered_peaks(peak_df: pd.DataFrame, min_relative_intensity: float, max_peaks: int) -> pd.DataFrame:
    if peak_df.empty:
        return peak_df.copy()
    work = peak_df.copy()
    work["intensity"] = pd.to_numeric(work.get("intensity"), errors="coerce").fillna(0.0)
    max_by_feature = work.groupby("feature_key", dropna=False)["intensity"].transform("max")
    work["_relative_intensity"] = work["intensity"] / max_by_feature.replace(0, pd.NA)
    if float(min_relative_intensity) > 0:
        work = work[work["_relative_intensity"].ge(float(min_relative_intensity))].copy()
    work = work.sort_values(["feature_key", "_relative_intensity", "mass"], ascending=[True, False, True], kind="mergesort")
    if int(max_peaks) > 0:
        work = work.groupby("feature_key", dropna=False).head(int(max_peaks)).copy()
    return work.drop(columns=["_relative_intensity"]).reset_index(drop=True)


def annotate_trimmed_counts(audit_df: pd.DataFrame, trimmed_df: pd.DataFrame) -> pd.DataFrame:
    if audit_df.empty or "feature_key" not in audit_df.columns:
        return audit_df
    out = audit_df.copy()
    counts = trimmed_df.groupby("feature_key", dropna=False).size().to_dict() if not trimmed_df.empty else {}
    emitted = out["status"].astype(str).eq("emitted")
    out.loc[emitted, "n_fragments_trimmed"] = out.loc[emitted, "feature_key"].map(lambda key: int(counts.get(key, 0)))
    return out


def summarize(audit_df: pd.DataFrame, peak_df: pd.DataFrame, merged_path: Path) -> dict[str, Any]:
    emitted = audit_df[audit_df["status"] == "emitted"] if not audit_df.empty else pd.DataFrame()
    skipped = audit_df[audit_df["status"] != "emitted"] if not audit_df.empty else pd.DataFrame()
    skip_counts: dict[str, int] = defaultdict(int)
    for reason in skipped.get("skip_reason", pd.Series(dtype=str)).fillna("").astype(str):
        skip_counts[reason or "unknown"] += 1

    per_lib: dict[str, int] = defaultdict(int)
    for lib in emitted.get("lib_id", pd.Series(dtype=int)):
        per_lib[str(int(lib))] += 1

    new_x_ids = sorted({clean_text(v) for v in emitted.get("x_id", []) if clean_text(v).startswith("X-")})
    prior = existing_x_ids(merged_path)
    gained = sorted(set(new_x_ids) - prior)

    return {
        "n_audit_rows": int(len(audit_df)),
        "n_emitted_spectra": int(len(emitted)),
        "n_skipped_spectra": int(len(skipped)),
        "n_peak_rows": int(len(peak_df)),
        "skip_reason_counts": dict(sorted(skip_counts.items())),
        "emitted_per_lib": dict(sorted(per_lib.items())),
        "n_unique_x_ids_emitted": len(new_x_ids),
        "n_x_ids_new_vs_merged": len(gained),
        "sample_new_x_ids": gained[:20],
    }


def backup_merged(path: Path) -> Path:
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    backup = path.with_name(f"{path.name}.bak.{ts}")
    shutil.copy2(path, backup)
    return backup


def run(args: argparse.Namespace) -> dict[str, Any]:
    chrolib = load_chrolib(args.chrolib_csv)
    known_map = load_known_x_map(args.known_x)
    version = clean_text(args.recovered_ms2_version) or ("v1" if args.parent_scope == "primary" else "v2")
    raw_peak_df, audit_df = convert_chrolib(
        chrolib,
        known_map,
        parent_scope=str(args.parent_scope),
        version=version,
    )
    peak_df = (
        trim_recovered_peaks(raw_peak_df, float(args.min_relative_intensity), int(args.max_consensus_peaks))
        if bool(args.trim_output)
        else raw_peak_df.copy()
    )
    audit_df = annotate_trimmed_counts(audit_df, peak_df)

    if args.raw_out is not None:
        write_csv(args.raw_out, raw_peak_df)
    write_csv(args.out, peak_df)
    write_csv(args.audit_out, audit_df)

    summary = summarize(audit_df, peak_df, args.merged_out)
    summary["out_csv"] = str(args.out)
    summary["raw_out_csv"] = str(args.raw_out) if args.raw_out is not None else ""
    summary["audit_csv"] = str(args.audit_out)
    summary["merge_strategy"] = args.merge_strategy
    summary["parent_scope"] = str(args.parent_scope)
    summary["recovered_ms2_version"] = version
    summary["trim_output"] = bool(args.trim_output)
    summary["min_relative_intensity"] = float(args.min_relative_intensity)
    summary["max_consensus_peaks"] = int(args.max_consensus_peaks)
    summary["n_raw_peak_rows"] = int(len(raw_peak_df))
    summary["n_trimmed_peak_rows"] = int(len(peak_df))

    if args.merge_strategy == "append" and not peak_df.empty:
        if args.merged_out.exists():
            backup = backup_merged(args.merged_out)
            summary["merged_backup"] = str(backup)
        merge_base = args.merge_base if args.merge_base is not None else args.merged_out
        if merge_base.exists():
            existing = pd.read_csv(merge_base, low_memory=False)
            summary["merged_base"] = str(merge_base)
        else:
            existing = pd.DataFrame(columns=PEAK_COLUMNS)
            args.merged_out.parent.mkdir(parents=True, exist_ok=True)
        merged = merge_peaks(existing, peak_df)
        dedup_dropped = merged.attrs.get("dedup_dropped", 0)
        summary["merged_rows_before"] = int(len(existing))
        summary["merged_rows_recovered"] = int(len(peak_df))
        summary["merged_rows_after"] = int(len(merged))
        summary["merged_dedup_dropped"] = int(dedup_dropped)
        write_csv(args.merged_out, merged)
        summary["merged_out"] = str(args.merged_out)
    else:
        summary["merged_out"] = ""

    print(json.dumps(summary, indent=2))

    summary_path = args.summary_json or args.audit_out.with_name("ions_recovered_chrolib_v1_summary.json")
    write_json(summary_path, summary)
    summary["summary_json"] = str(summary_path)
    return summary


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
