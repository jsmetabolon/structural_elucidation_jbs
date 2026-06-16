#!/usr/bin/env python3
"""Add pathway/subpathway predictions to equalization rerank top-10 CSVs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from build_v14_all_library_prediction_review_v1 import (
    clean_text,
    enrich_pathway_annotations,
    generate_pathway_assignment_map,
    normalize_assignment_map,
)
from build_v19b_updated_summary_v1 import apply_amino_acid_pathway_overrides
from run_retrieval_anchored_denovo_assembly_v4 import resolve
from x209_phase0_utils import canonicalize_smiles


DEFAULT_EQUALIZATION_ROOT = Path("models/pipeline_equalization_v21_clean")
DEFAULT_REFERENCE_DB = Path("Metabolite Synonym Database/metabolite_synonyms_merged_with_mtx_giga.sqlite")
PATHWAY_COLUMNS = [
    "predicted_pathway",
    "predicted_subpathway",
    "pathway_confidence_tier",
    "pathway_assignment_status",
    "pathway_vote_share",
    "pathway_best_reference_similarity",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--equalization-root", type=Path, default=DEFAULT_EQUALIZATION_ROOT)
    parser.add_argument("--reference-db", type=Path, default=DEFAULT_REFERENCE_DB)
    parser.add_argument("--libs", nargs="+", default=["209,305,400,402"])
    parser.add_argument("--input-name", default="unified_rerank_all_top10.csv")
    parser.add_argument("--output-name", default="unified_rerank_all_top10_with_pathways.csv")
    parser.add_argument("--summary-json", type=Path, default=None)
    parser.add_argument("--assignment-cache", type=Path, default=None)
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


def canonicalize_for_pathway(value: object) -> str:
    text = clean_text(value)
    if not text:
        return ""
    return clean_text(canonicalize_smiles(text))


def rerank_path(root: Path, lib_id: int, name: str) -> Path:
    return root / f"rerank_lib{int(lib_id)}" / name


def load_inputs(root: Path, libs: list[int], input_name: str, allow_missing: bool) -> dict[int, pd.DataFrame]:
    frames: dict[int, pd.DataFrame] = {}
    missing: list[Path] = []
    for lib_id in libs:
        path = rerank_path(root, lib_id, input_name)
        if not path.exists():
            missing.append(path)
            continue
        frame = pd.read_csv(path, low_memory=False)
        if "predicted_smiles" not in frame.columns:
            raise ValueError(f"{path} missing required column: predicted_smiles")
        frames[int(lib_id)] = frame
    if missing and not allow_missing:
        joined = "\n".join(f"  {path}" for path in missing)
        raise FileNotFoundError(f"Missing rerank top10 input CSVs:\n{joined}")
    return frames


def load_assignment_cache(path: Path | None) -> pd.DataFrame:
    if path is None or not path.exists():
        return normalize_assignment_map(None)
    return normalize_assignment_map(pd.read_csv(path, low_memory=False))


def build_assignment_map(
    *,
    frames: dict[int, pd.DataFrame],
    reference_db: Path,
    assignment_cache: Path | None,
) -> tuple[pd.DataFrame, dict[str, int]]:
    unique_smiles = sorted(
        {
            canonical
            for frame in frames.values()
            for canonical in frame["predicted_smiles"].map(canonicalize_for_pathway)
            if canonical
        }
    )
    cached = load_assignment_cache(assignment_cache)
    cached_smiles = set(cached["canonical_smiles"].map(clean_text)) if not cached.empty else set()
    missing = [smiles for smiles in unique_smiles if smiles not in cached_smiles]
    generated = normalize_assignment_map(None)
    if missing:
        if not reference_db.exists():
            raise FileNotFoundError(reference_db)
        generated = generate_pathway_assignment_map(missing, reference_db)
    assignment_map = (
        pd.concat([cached, generated], ignore_index=True, sort=False)
        .pipe(normalize_assignment_map)
        .drop_duplicates(subset=["canonical_smiles"], keep="first")
        .reset_index(drop=True)
    )
    if assignment_cache is not None:
        assignment_cache.parent.mkdir(parents=True, exist_ok=True)
        assignment_map.to_csv(assignment_cache, index=False)
    stats = {
        "n_unique_input_smiles": int(len(unique_smiles)),
        "n_cached_assignment_rows": int(len(cached)),
        "n_missing_assignment_rows_requested": int(len(missing)),
        "n_generated_assignment_rows": int(len(generated)),
        "n_assignment_rows": int(len(assignment_map)),
    }
    return assignment_map, stats


def ordered_columns(original_columns: list[str], enriched: pd.DataFrame) -> list[str]:
    columns = list(original_columns)
    for column in PATHWAY_COLUMNS:
        if column not in columns and column in enriched.columns:
            columns.append(column)
    return columns + [column for column in enriched.columns if column not in columns and column not in {"pathway_match", "subpathway_match"}]


def enrich_one(
    *,
    frame: pd.DataFrame,
    assignment_map: pd.DataFrame,
) -> pd.DataFrame:
    original_columns = list(frame.columns)
    out = enrich_pathway_annotations(frame, assignment_map=assignment_map, reference_db=None)
    out = apply_amino_acid_pathway_overrides(out)
    if "pathway_match" not in original_columns:
        out = out.drop(columns=["pathway_match"], errors="ignore")
    if "subpathway_match" not in original_columns:
        out = out.drop(columns=["subpathway_match"], errors="ignore")
    return out.loc[:, ordered_columns(original_columns, out)]


def summarize_lib(lib_id: int, input_path: Path, output_path: Path, frame: pd.DataFrame) -> dict[str, Any]:
    status = frame["pathway_assignment_status"].map(clean_text) if "pathway_assignment_status" in frame.columns else pd.Series([], dtype=object)
    tier = frame["pathway_confidence_tier"].map(clean_text) if "pathway_confidence_tier" in frame.columns else pd.Series([], dtype=object)
    path = frame["predicted_pathway"].map(clean_text) if "predicted_pathway" in frame.columns else pd.Series([], dtype=object)
    sub = frame["predicted_subpathway"].map(clean_text) if "predicted_subpathway" in frame.columns else pd.Series([], dtype=object)
    smiles = frame["predicted_smiles"].map(clean_text) if "predicted_smiles" in frame.columns else pd.Series([], dtype=object)
    return {
        "lib_id": int(lib_id),
        "input_csv": str(input_path),
        "output_csv": str(output_path),
        "n_rows": int(len(frame)),
        "n_rows_with_smiles": int(smiles.ne("").sum()),
        "n_unique_smiles": int(len({canonicalize_for_pathway(value) for value in smiles if clean_text(value)} - {""})),
        "n_rows_with_predicted_pathway": int(path.ne("").sum()),
        "n_rows_with_predicted_subpathway": int(sub.ne("").sum()),
        "pathway_assignment_status_counts": {str(k): int(v) for k, v in status.value_counts(dropna=False).sort_index().items()},
        "pathway_confidence_tier_counts": {str(k): int(v) for k, v in tier.value_counts(dropna=False).sort_index().items()},
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parent
    root = resolve(args.equalization_root, repo_root)
    reference_db = resolve(args.reference_db, repo_root)
    libs = parse_libs(args.libs)
    if not libs:
        raise ValueError("No libraries requested")
    frames = load_inputs(root, libs, args.input_name, bool(args.allow_missing))
    if not frames:
        raise FileNotFoundError(f"No rerank inputs found under {root}")
    assignment_cache = args.assignment_cache
    if assignment_cache is None:
        assignment_cache = root / "pathway_assignment_map_v1.csv"
    else:
        assignment_cache = resolve(assignment_cache, repo_root)
    assignment_map, assignment_stats = build_assignment_map(
        frames=frames,
        reference_db=reference_db,
        assignment_cache=assignment_cache,
    )

    lib_summaries: list[dict[str, Any]] = []
    for lib_id in libs:
        if lib_id not in frames:
            continue
        input_path = rerank_path(root, lib_id, args.input_name)
        output_path = rerank_path(root, lib_id, args.output_name)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        enriched = enrich_one(frame=frames[lib_id], assignment_map=assignment_map)
        enriched.to_csv(output_path, index=False)
        lib_summaries.append(summarize_lib(lib_id, input_path, output_path, enriched))

    summary_path = args.summary_json
    if summary_path is None:
        summary_path = root / "pathway_enrichment_summary.json"
    else:
        summary_path = resolve(summary_path, repo_root)
    summary = {
        "equalization_root": str(root),
        "reference_db": str(reference_db),
        "assignment_cache_csv": str(assignment_cache),
        "input_name": args.input_name,
        "output_name": args.output_name,
        **assignment_stats,
        "libs": lib_summaries,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def main() -> int:
    run(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
