#!/usr/bin/env python3
"""Add peptide composition summary columns to equalization rerank CSVs.

Derives human-readable peptide composition labels, confidence scores, and
ordering labels from existing peptide-track columns and optional SMILES inference.

Pipeline context:
    Third enrichment stage; reads lipid-class-enriched top-10 CSVs from
    ``pipeline_equalization/rerank_lib*/``.

Key I/O:
    Inputs: ``*_lipid_class.csv`` per library.
    Outputs: ``*_peptide_composition.csv``, ``peptide_composition_enrichment_summary.json``.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

import pandas as pd

from structural_elucidation.paths import REPO_ROOT, resolve

try:
    from peptide_permutation_v1_utils import infer_sequence_from_smiles as _infer_sequence_from_smiles
except Exception:  # pragma: no cover - only available in the full EC2 workspace
    _infer_sequence_from_smiles = None


# ── Defaults and output schema ───────────────────────────────────────
DEFAULT_EQUALIZATION_ROOT: Path = Path("models/pipeline_equalization")
DEFAULT_LIBS: str = "209,305,400,402"
DEFAULT_INPUT_NAME: str = "unified_rerank_all_top10_with_pathways_lipid_class.csv"
DEFAULT_OUTPUT_NAME: str = (
    "unified_rerank_all_top10_with_pathways_lipid_class_peptide_composition.csv"
)
PEPTIDE_COLUMNS: list[str] = [
    "peptide_composition_label",
    "peptide_composition_confidence",
    "peptide_ordering_confidence",
    "peptide_ordering_label",
    "peptide_length",
]
RESIDUES: dict[str, str] = {  # three-letter code → display token
    "ALA": "Ala",
    "ARG": "Arg",
    "ASN": "Asn",
    "ASP": "Asp",
    "CYS": "Cys",
    "GLN": "Gln",
    "GLU": "Glu",
    "GLY": "Gly",
    "HIS": "His",
    "ILE": "Xle",
    "LEU": "Xle",
    "LYS": "Lys",
    "MET": "Met",
    "PHE": "Phe",
    "PRO": "Pro",
    "PYR": "Pyr",
    "SER": "Ser",
    "THR": "Thr",
    "TRP": "Trp",
    "TYR": "Tyr",
    "VAL": "Val",
    "XLE": "Xle",
}
FULL_RESIDUES: dict[str, str] = {  # full residue name → display token
    "ALANINE": "Ala",
    "ARGININE": "Arg",
    "ASPARAGINE": "Asn",
    "ASPARTIC": "Asp",
    "CYSTEINE": "Cys",
    "GLUTAMINE": "Gln",
    "GLUTAMIC": "Glu",
    "GLYCINE": "Gly",
    "HISTIDINE": "His",
    "ISOLEUCINE": "Xle",
    "LEUCINE": "Xle",
    "LYSINE": "Lys",
    "METHIONINE": "Met",
    "PHENYLALANINE": "Phe",
    "PROLINE": "Pro",
    "PYROGLUTAMATE": "Pyr",
    "PYROGLUTAMIC": "Pyr",
    "SERINE": "Ser",
    "THREONINE": "Thr",
    "TRYPTOPHAN": "Trp",
    "TYROSINE": "Tyr",
    "VALINE": "Val",
}
RESIDUE_ORDER: dict[str, int] = {
    residue: idx
    for idx, residue in enumerate(
        [
            "Ala",
            "Arg",
            "Asn",
            "Asp",
            "Cys",
            "Gln",
            "Glu",
            "Gly",
            "His",
            "Lys",
            "Met",
            "Phe",
            "Pro",
            "Pyr",
            "Ser",
            "Thr",
            "Trp",
            "Tyr",
            "Val",
            "Xle",
        ]
    )
}
LENGTH_LABELS: dict[int, str] = {1: "amino acid", 2: "dipeptide", 3: "tripeptide", 4: "tetrapeptide"}


# ── CLI and parsing helpers ──────────────────────────────────────────
def parse_args() -> argparse.Namespace:
    """Build the peptide-composition enrichment CLI argument parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--equalization-root", type=Path, default=DEFAULT_EQUALIZATION_ROOT)
    parser.add_argument("--libs", nargs="+", default=[DEFAULT_LIBS])
    parser.add_argument("--input-name", default=DEFAULT_INPUT_NAME)
    parser.add_argument("--output-name", default=DEFAULT_OUTPUT_NAME)
    parser.add_argument("--summary-json", type=Path, default=None)
    parser.add_argument("--ordering-threshold", type=float, default=0.65)
    parser.add_argument("--allow-missing", action="store_true")
    return parser.parse_args()


def parse_libs(values: list[str]) -> list[int]:
    libs: list[int] = []
    for value in values:
        for token in str(value).split(","):
            token = token.strip()
            if token:
                libs.append(int(token))
    return libs


def clean_text(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return str(value).strip()


def numeric_value(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
        return out if math.isfinite(out) else default
    except Exception:
        return default


def clip01(value: float) -> float:
    return round(float(min(1.0, max(0.0, value))), 6)


def boolish(value: object) -> bool:
    text = clean_text(value).lower()
    if text in {"true", "t", "yes", "y", "1", "pass", "passed", "tight", "supported"}:
        return True
    if text in {"false", "f", "no", "n", "0", "", "nan", "none"}:
        return False
    return numeric_value(value, 0.0) > 0.0


def parse_residue_tokens(value: Any) -> list[str]:
    text = clean_text(value)
    if not text:
        return []
    upper = text.upper()
    upper = re.sub(r"\bKNOWN\b|\bRELEASED\b|\bPUTATIVE\b|\bTENTATIVE\b", " ", upper)
    upper = re.sub(r"\b[DL]-", "", upper)
    pattern = "|".join(sorted([*FULL_RESIDUES.keys(), *RESIDUES.keys()], key=len, reverse=True))
    tokens: list[str] = []
    for match in re.finditer(rf"(?<![A-Z])({pattern})(?![A-Z])", upper):
        raw = match.group(1)
        token = FULL_RESIDUES.get(raw, RESIDUES.get(raw, ""))
        if token:
            tokens.append(token)
    return tokens


def normalize_peptide_composition(value: Any) -> str:
    tokens = parse_residue_tokens(value)
    if not tokens:
        return ""
    return "-".join(sorted(tokens, key=lambda token: (RESIDUE_ORDER.get(token, 999), token)))


def normalize_ordered_sequence(value: Any) -> str:
    return "-".join(parse_residue_tokens(value))


def composition_tokens(value: Any) -> list[str]:
    return [token for token in normalize_peptide_composition(value).split("-") if token]


def peptide_length_label(length: int) -> str:
    return LENGTH_LABELS.get(int(length), f"{int(length)}-residue peptide")


def peptide_composition_label(composition: Any) -> str:
    tokens = composition_tokens(composition)
    if not tokens:
        return ""
    return f"{peptide_length_label(len(tokens))}: {', '.join(tokens)}"


def infer_sequence_from_smiles(smiles: Any) -> dict[str, Any]:
    if _infer_sequence_from_smiles is None:
        return {}
    try:
        result = _infer_sequence_from_smiles(smiles)
        return result if isinstance(result, dict) else {}
    except Exception:
        return {}


def row_is_peptide_candidate(row: pd.Series) -> bool:
    if clean_text(row.get("normalized_peptide_composition")):
        return True
    if boolish(row.get("is_peptide_candidate")):
        return True
    return numeric_value(row.get("peptide_bond_count"), 0.0) > 0.0


def extract_composition_and_order(row: pd.Series) -> tuple[str, str, dict[str, Any]]:
    composition = normalize_peptide_composition(row.get("normalized_peptide_composition"))
    ordered = normalize_ordered_sequence(row.get("candidate_ordered_sequence"))
    parsed = infer_sequence_from_smiles(row.get("predicted_smiles"))
    parsed_ordered = normalize_ordered_sequence(parsed.get("ordered_peptide_sequence", ""))
    if not composition and parsed_ordered:
        composition = normalize_peptide_composition(parsed_ordered)
    if not ordered and parsed_ordered:
        ordered = parsed_ordered
    if not composition and numeric_value(row.get("peptide_bond_count"), numeric_value(parsed.get("peptide_bond_count"), 0.0)) > 0.0:
        composition = normalize_peptide_composition(row.get("candidate_ordered_sequence"))
    if not ordered and composition:
        tokens = composition_tokens(composition)
        if len(set(tokens)) == 1:
            ordered = "-".join(tokens)
    return composition, ordered, parsed


def immonium_evidence_score(row: pd.Series) -> float:
    if boolish(row.get("immonium_evidence")) or boolish(row.get("high_confidence_residue_composition_support")):
        return 1.0
    score = max(
        numeric_value(row.get("spectrum_immonium_score")),
        numeric_value(row.get("immonium_score")),
        numeric_value(row.get("neutral_loss_evidence_score")),
    )
    count = max(
        numeric_value(row.get("spectrum_immonium_count")),
        numeric_value(row.get("immonium_count")),
        numeric_value(row.get("spectrum_neutral_loss_count")),
    )
    return clip01(max(score, min(1.0, count / 3.0)))


def mass_fit_score(row: pd.Series) -> float:
    raw = row.get("peptide_neutral_mass_fit")
    if boolish(raw):
        numeric = numeric_value(raw, math.nan)
        return clip01(numeric) if math.isfinite(numeric) and numeric <= 1.0 else 1.0
    scores = [
        numeric_value(row.get("mass_match_score"), math.nan),
        numeric_value(row.get("score_mass_base"), math.nan),
    ]
    valid_scores = [value for value in scores if math.isfinite(value)]
    errors = [
        abs(numeric_value(row.get("peptide_mass_error_da"), math.nan)),
        abs(numeric_value(row.get("composition_mass_error_da"), math.nan)),
        abs(numeric_value(row.get("abs_mass_error_da"), math.nan)),
        abs(numeric_value(row.get("primary_abs_mass_error_da"), math.nan)),
    ]
    valid_errors = [value for value in errors if math.isfinite(value)]
    if valid_errors:
        valid_scores.append(max(0.0, 1.0 - min(valid_errors) / 0.05))
    return clip01(max(valid_scores or [0.0]))


def peptide_composition_confidence(row: pd.Series) -> float:
    """Blend hypothesis, support, immonium, and mass-fit evidence into [0, 1]."""
    hypothesis = numeric_value(row.get("composition_confidence"), 0.0)
    support = max(
        numeric_value(row.get("composition_support_score"), 0.0),
        numeric_value(row.get("composition_support_score_norm"), 0.0),
    )
    score = (
        0.40 * clip01(hypothesis)
        + 0.30 * clip01(support)
        + 0.15 * immonium_evidence_score(row)
        + 0.15 * mass_fit_score(row)
    )  # mass fit weighted lightly because peptide rows often lack tight DA errors
    return clip01(score)


def b_y_ion_score(row: pd.Series) -> float:
    b_score = numeric_value(row.get("spectrum_b_ion_score"), math.nan)
    y_score = numeric_value(row.get("spectrum_y_ion_score"), math.nan)
    scores = [value for value in (b_score, y_score) if math.isfinite(value)]
    if scores:
        return clip01(sum(scores) / max(len(scores), 1))
    b_count = numeric_value(row.get("spectrum_b_ion_count"), 0.0)
    y_count = numeric_value(row.get("spectrum_y_ion_count"), 0.0)
    return clip01((b_count + y_count) / 4.0)


def only_one_ordering(tokens: list[str]) -> bool:
    return bool(tokens) and len(set(tokens)) == 1


def peptide_ordering_confidence(row: pd.Series, composition: str, parsed: dict[str, Any]) -> float:
    tokens = composition_tokens(composition)
    if only_one_ordering(tokens):
        return 1.0
    rank = numeric_value(row.get("ordered_sequence_within_composition_rank"), math.nan)
    rank_score = 1.0 if math.isfinite(rank) and rank <= 1.0 else 0.0
    connectivity = max(
        numeric_value(row.get("peptide_connectivity_confidence"), 0.0),
        numeric_value(parsed.get("peptide_connectivity_confidence"), 0.0),
    )
    support = max(
        numeric_value(row.get("ordered_sequence_support_score"), 0.0),
        numeric_value(row.get("sequence_support_score"), 0.0),
    )
    score = 0.30 * rank_score + 0.30 * b_y_ion_score(row) + 0.20 * clip01(connectivity) + 0.20 * clip01(support)
    return clip01(score)


def peptide_ordering_label(ordered: str, composition: str, confidence: float, threshold: float) -> str:
    """Return ordered sequence when confidence clears threshold, else 'order uncertain'."""
    if not composition:
        return ""
    tokens = composition_tokens(composition)
    if not ordered and only_one_ordering(tokens):
        ordered = "-".join(tokens)
    if ordered and confidence >= float(threshold):
        return ordered
    return "order uncertain"


def annotate_row(row: pd.Series, ordering_threshold: float = 0.65) -> dict[str, Any]:
    """Build peptide summary columns for one rerank row."""
    if not row_is_peptide_candidate(row):
        return {
            "peptide_composition_label": "",
            "peptide_composition_confidence": "",
            "peptide_ordering_confidence": "",
            "peptide_ordering_label": "",
            "peptide_length": "",
        }
    composition, ordered, parsed = extract_composition_and_order(row)
    if not composition:
        return {
            "peptide_composition_label": "",
            "peptide_composition_confidence": "",
            "peptide_ordering_confidence": "",
            "peptide_ordering_label": "",
            "peptide_length": "",
        }
    comp_conf = peptide_composition_confidence(row)
    order_conf = peptide_ordering_confidence(row, composition, parsed)
    return {
        "peptide_composition_label": peptide_composition_label(composition),
        "peptide_composition_confidence": comp_conf,
        "peptide_ordering_confidence": order_conf,
        "peptide_ordering_label": peptide_ordering_label(ordered, composition, order_conf, ordering_threshold),
        "peptide_length": int(len(composition_tokens(composition))),
    }


# ── I/O and enrichment driver ────────────────────────────────────────
def rerank_path(root: Path, lib_id: int, name: str) -> Path:
    """Path to a named CSV under ``rerank_lib{lib_id}/``."""
    return root / f"rerank_lib{int(lib_id)}" / name


def load_inputs(root: Path, libs: list[int], input_name: str, allow_missing: bool) -> dict[int, pd.DataFrame]:
    frames: dict[int, pd.DataFrame] = {}
    missing: list[Path] = []
    for lib_id in libs:
        path = rerank_path(root, lib_id, input_name)
        if not path.exists():
            missing.append(path)
            continue
        frames[int(lib_id)] = pd.read_csv(path, low_memory=False)
    if missing and not allow_missing:
        joined = "\n".join(f"  {path}" for path in missing)
        raise FileNotFoundError(f"Missing rerank input CSVs:\n{joined}")
    return frames


def ordered_columns(original_columns: list[str], enriched: pd.DataFrame) -> list[str]:
    columns = list(original_columns)
    for column in PEPTIDE_COLUMNS:
        if column not in columns and column in enriched.columns:
            columns.append(column)
    return columns + [column for column in enriched.columns if column not in columns]


def enrich_frame(frame: pd.DataFrame, ordering_threshold: float = 0.65) -> pd.DataFrame:
    """Add peptide composition columns to all rows in a rerank top-10 frame."""
    original_columns = list(frame.columns)
    out = frame.copy()
    annotations = [annotate_row(row, ordering_threshold=ordering_threshold) for _, row in out.iterrows()]
    annotation_frame = pd.DataFrame(annotations, index=out.index)
    for column in PEPTIDE_COLUMNS:
        out[column] = annotation_frame[column]
    return out.loc[:, ordered_columns(original_columns, out)]


def summarize_lib(lib_id: int, input_path: Path, output_path: Path, frame: pd.DataFrame) -> dict[str, Any]:
    labels = frame["peptide_composition_label"].fillna("").astype(str)
    ordering = frame["peptide_ordering_label"].fillna("").astype(str)
    length = pd.to_numeric(frame["peptide_length"], errors="coerce")
    return {
        "lib_id": int(lib_id),
        "input_csv": str(input_path),
        "output_csv": str(output_path),
        "n_rows": int(len(frame)),
        "n_rows_with_peptide_composition": int(labels.ne("").sum()),
        "n_rows_with_ordered_peptide": int(ordering.ne("").sum() - ordering.eq("order uncertain").sum()),
        "peptide_length_counts": {str(int(k)): int(v) for k, v in length.dropna().astype(int).value_counts().sort_index().items()},
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    """Run peptide-composition enrichment for all requested libraries."""
    repo_root = REPO_ROOT
    root = resolve(args.equalization_root, repo_root)
    libs = parse_libs(args.libs)
    frames = load_inputs(root, libs, args.input_name, args.allow_missing)
    summaries: list[dict[str, Any]] = []
    for lib_id, frame in frames.items():
        input_path = rerank_path(root, int(lib_id), args.input_name)
        output_path = rerank_path(root, int(lib_id), args.output_name)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        enriched = enrich_frame(frame, ordering_threshold=float(args.ordering_threshold))
        enriched.to_csv(output_path, index=False)
        summaries.append(summarize_lib(int(lib_id), input_path, output_path, enriched))
    summary = {
        "equalization_root": str(root),
        "input_name": str(args.input_name),
        "output_name": str(args.output_name),
        "n_libs": int(len(frames)),
        "n_rows": int(sum(item["n_rows"] for item in summaries)),
        "n_rows_with_peptide_composition": int(sum(item["n_rows_with_peptide_composition"] for item in summaries)),
        "n_rows_with_ordered_peptide": int(sum(item["n_rows_with_ordered_peptide"] for item in summaries)),
        "per_lib": summaries,
    }
    summary_path = resolve(args.summary_json, repo_root) if args.summary_json else root / "peptide_composition_enrichment_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def main() -> int:
    """CLI entry point."""
    run(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
