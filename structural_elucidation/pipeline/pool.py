#!/usr/bin/env python3
"""Build a V16 unified pool while leaving the V1 pool builder unchanged.

Thin overlay on ``pool_v1``: reuses V1 merge/dedup logic, then attaches NCE
draw-count features from scored candidates and refreshes S0/non-S0 splits.

Pipeline context:
    Called when a prior unified pool CSV exists (lib-specific defaults in
    ``default_base_pool_for_scope``). Falls back to ``pool_v1.run`` when no
    base pool is available.

Key I/O:
    Inputs: base unified pool, scored candidates, extra enumerator CSVs,
    integrated per-control best for winner-preservation check.
    Outputs: ``unified_candidates_per_control.csv`` (+ S0 splits),
    ``unified_candidate_pool_v16_summary.json``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

import structural_elucidation.pipeline.pool_v1 as base
from structural_elucidation.utils.denovo import canonicalize_smiles
from structural_elucidation.paths import REPO_ROOT, resolve
from structural_elucidation.utils.domain import write_csv, write_json
from structural_elucidation.utils.timing import stage_timer

# ── NCE feature columns overlaid onto the unified pool ───────────────
NCE_COLS: list[str] = ["nce_draw_count", "nce_max_log_prob", "nce_mean_log_prob"]


# ── Scope-specific base pool defaults ────────────────────────────────
def default_base_pool_for_scope(lib_id: int, ion_mode: str) -> Path | None:
    """Return the canonical prior unified-pool CSV for a lib/ion scope."""
    lookup = {
        (209, "negative"): Path("models/rerank_lib209_step7_targeted_v2_auto/unified_candidates_scored.csv"),
        (305, "negative"): Path("models/rerank_v9_no_wgcna_lib305_negative/unified_candidates_scored.csv"),
        (400, "positive"): Path("models/rerank_v9_no_wgcna_lib400_positive/unified_candidates_scored.csv"),
        (402, "positive"): Path("models/rerank_v9_no_wgcna_lib402_positive/unified_candidates_scored.csv"),
    }
    return lookup.get((int(lib_id), str(ion_mode)))


def nce_feature_table(path: Path, repo_root: Path) -> pd.DataFrame:
    """Extract per-(control, SMILES) NCE draw statistics from scored candidates."""
    resolved = resolve(path, repo_root)
    if not resolved.exists() or resolved.stat().st_size == 0:
        return pd.DataFrame(columns=["control_id", "canonical_2d_smiles", *NCE_COLS])
    try:
        raw = pd.read_csv(resolved, low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["control_id", "canonical_2d_smiles", *NCE_COLS])
    if raw.empty:
        return pd.DataFrame(columns=["control_id", "canonical_2d_smiles", *NCE_COLS])
    out = raw.copy()
    if "control_id" not in out.columns:
        out["control_id"] = (
            out.get("feature_key", pd.Series("", index=out.index))
            .astype(str)
            .str.extract(r"(X-\d+)", expand=False)
            .fillna("")
        )
    if "canonical_2d_smiles" not in out.columns:
        out["canonical_2d_smiles"] = out.get("predicted_smiles", pd.Series("", index=out.index)).map(
            lambda value: canonicalize_smiles(value) or ""
        )
    for col in NCE_COLS:
        if col not in out.columns:
            out[col] = 0.0
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0.0)
    out = out[out["control_id"].astype(str).str.startswith("X-") & out["canonical_2d_smiles"].astype(str).ne("")].copy()
    if out.empty:
        return pd.DataFrame(columns=["control_id", "canonical_2d_smiles", *NCE_COLS])
    return out.groupby(["control_id", "canonical_2d_smiles"], dropna=False)[NCE_COLS].max().reset_index()


def overlay_nce_features(pool: pd.DataFrame, nce_features: pd.DataFrame) -> pd.DataFrame:
    """Left-join NCE draw features; flag rows sourced from nce_direct."""
    if pool.empty or nce_features.empty:
        return pool
    out = pool.drop(columns=[col for col in NCE_COLS if col in pool.columns], errors="ignore")
    out = out.merge(nce_features, on=["control_id", "canonical_2d_smiles"], how="left")
    for col in NCE_COLS:
        out[col] = pd.to_numeric(out.get(col), errors="coerce").fillna(0.0)
    if "source_labels" in out.columns:
        out["is_nce_direct"] = (
            out["source_labels"].fillna("").astype(str).str.contains("nce_direct", regex=False).astype(int)
        )
    return out


def safe_read_csv(path: Path) -> pd.DataFrame:
    """Read CSV or return empty frame on missing/empty files."""
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def numeric_col(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    """Coerce column to float with fillna."""
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def normalize_existing_pool(pool: pd.DataFrame) -> pd.DataFrame:
    """Align legacy pool column names before delegating to V1 normalize."""
    if pool.empty:
        return pool
    out = pool.copy()
    if "source_label" not in out.columns and "source_labels" in out.columns:
        out["source_label"] = out["source_labels"]
    return base.normalize(out)


def load_extra_only(extra_text: str, repo_root: Path) -> pd.DataFrame:
    """Load comma-separated extra candidate CSVs without the scored-candidates frame."""
    frames: list[pd.DataFrame] = []
    for part in str(extra_text or "").split(","):
        text = part.strip()
        if not text:
            continue
        path = base.prefer_verifier_csv(Path(text), repo_root)
        df = safe_read_csv(path)
        if df.empty:
            continue
        if "source_label" not in df.columns:
            df = df.copy()
            df["source_label"] = "mlims_anchored_v1" if "mlims_anchored" in str(path) else path.parent.name
        frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def refresh_splits(output_root: Path, pool: pd.DataFrame, args: argparse.Namespace, repo_root: Path) -> None:
    """Write full, S0-only, and non-S0 pool CSV splits."""
    s0_controls = base.s0_controls_for_scope(resolve(args.all_x_catalog, repo_root), int(args.lib_id))
    pool_s0 = (
        pool[pool.get("control_id", pd.Series(dtype=str)).astype(str).isin(s0_controls)].copy()
        if not pool.empty and s0_controls
        else pd.DataFrame(columns=pool.columns)
    )
    pool_non_s0 = (
        pool[~pool.get("control_id", pd.Series(dtype=str)).astype(str).isin(s0_controls)].copy()
        if not pool.empty and s0_controls
        else pd.DataFrame(columns=pool.columns)
    )
    write_csv(output_root / "unified_candidates_per_control.csv", pool)
    write_csv(output_root / "unified_candidates_per_control_s0.csv", pool_s0)
    write_csv(output_root / "unified_candidates_per_control_non_s0.csv", pool_non_s0)


def run_from_base_pool(
    args: argparse.Namespace, repo_root: Path, output_root: Path, base_pool_path: Path
) -> dict[str, Any]:
    """Merge base pool + scored + extras, overlay NCE features, check winner preservation."""
    current = safe_read_csv(resolve(args.current_integrated_best, repo_root))
    if "lib_id" in current.columns:
        current = current[pd.to_numeric(current["lib_id"], errors="coerce").eq(int(args.lib_id))].copy()
    if "ion_mode" in current.columns:
        current = current[current["ion_mode"].astype(str).eq(str(args.ion_mode))].copy()
    allowed_controls = set(current["control_id"].astype(str)) if "control_id" in current.columns else set()
    if bool(args.include_all_x):
        catalog_path = resolve(args.all_x_catalog, repo_root)
        all_x = base.all_x_controls_for_scope(catalog_path, int(args.lib_id))
        if not all_x:
            # Without the catalog the pool silently collapses to the
            # integrated-best controls, which is not what --include-all-x asks for.
            raise SystemExit(
                f"--include-all-x given but no X controls for lib {args.lib_id} in {catalog_path}"
            )
        allowed_controls |= all_x
    base_pool = normalize_existing_pool(safe_read_csv(resolve(base_pool_path, repo_root)))
    scored = base.normalize(safe_read_csv(resolve(args.scored_candidates, repo_root)))
    extra_text = str(args.extra_inputs or "")
    if extra_text == ",".join(base.DEFAULT_EXTRA_INPUTS) and (int(args.lib_id), str(args.ion_mode)) != (
        209,
        "negative",
    ):
        extra_text = base.default_extra_inputs_for_scope(int(args.lib_id), str(args.ion_mode))
    extras = base.normalize(load_extra_only(extra_text, repo_root))
    frames = [frame for frame in (base_pool, scored, extras) if not frame.empty]
    raw = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()
    norm = base.restrict_controls(base.normalize(raw), allowed_controls)
    pool = base.dedup(norm)
    nce_features = nce_feature_table(args.scored_candidates, repo_root)
    pool = overlay_nce_features(pool, nce_features)
    if not pool.empty:
        sort_cols = [
            col
            for col in ["control_id", "predicted_spectrum_score_v2", "mass_match_score", "min_source_rank"]
            if col in pool.columns
        ]
        ascending_lookup = {
            "control_id": True,
            "predicted_spectrum_score_v2": False,
            "mass_match_score": False,
            "min_source_rank": True,
        }
        pool = pool.sort_values(sort_cols, ascending=[ascending_lookup[col] for col in sort_cols], kind="mergesort")
        pool["unified_pool_rank"] = pool.groupby("control_id", dropna=False).cumcount() + 1
    refresh_splits(output_root, pool, args, repo_root)
    winners = (
        current[pd.to_numeric(current.get("top10_best_tanimoto"), errors="coerce").fillna(0.0).ge(0.7)].copy()
        if not current.empty
        else pd.DataFrame()
    )
    pool_keys = (
        set(
            zip(
                pool.get("control_id", pd.Series(dtype=str)).astype(str),
                pool.get("canonical_2d_smiles", pd.Series(dtype=str)).astype(str),
            )
        )
        if not pool.empty
        else set()
    )
    preserve_rows = []
    for _, row in winners.iterrows():
        key = (str(row.get("control_id", "")), base.canonicalize_smiles(row.get("top10_smiles", "")) or "")
        preserve_rows.append(
            {"control_id": key[0], "winner_present": key in pool_keys, "winner_smiles": row.get("top10_smiles", "")}
        )
    preserve = pd.DataFrame(preserve_rows)
    write_csv(output_root / "preserve_25_pool_check.csv", preserve)
    n_expected = int(preserve.shape[0])
    n_present = int(preserve["winner_present"].sum()) if not preserve.empty else 0
    # Require at least 25 of current high-tanimoto winners to survive merge.
    summary = {
        "base_unified_pool_csv": str(resolve(base_pool_path, repo_root)),
        "n_base_rows": int(base_pool.shape[0]),
        "n_extra_rows": int(extras.shape[0]),
        "n_raw_rows": int(raw.shape[0]),
        "n_normalized_rows": int(norm.shape[0]),
        "n_unified_rows": int(pool.shape[0]),
        "n_controls": int(pool["control_id"].nunique()) if not pool.empty else 0,
        "n_s0_rows": int(safe_read_csv(output_root / "unified_candidates_per_control_s0.csv").shape[0]),
        "n_nce_feature_rows": int(nce_features.shape[0]),
        "n_pool_rows_with_nce_draw_count": (
            int(numeric_col(pool, "nce_draw_count").gt(0.0).sum()) if not pool.empty else 0
        ),
        "n_current_winners_present": n_present,
        "n_current_winners_expected": n_expected,
        "acceptance_pass": bool(n_expected == 0 or n_present >= min(25, n_expected)),
        "v16_overlay": True,
    }
    write_json(output_root / "unified_candidate_pool_summary.json", summary)
    write_json(output_root / "unified_candidate_pool_v16_summary.json", summary)
    return summary


def run(args: argparse.Namespace) -> dict[str, Any]:
    """Run V16 overlay from base pool or delegate to V1 then overlay NCE features."""
    repo_root = REPO_ROOT
    output_root = resolve(args.output_root, repo_root)
    current_path = resolve(args.current_integrated_best, repo_root)
    if not current_path.exists():
        nested_current = (
            Path("models")
            / "integrated_denovo_audit_v15"
            / f"lib_{int(args.lib_id)}_{str(args.ion_mode)}"
            / "integrated_per_control_best.csv"
        )
        if resolve(nested_current, repo_root).exists():
            args.current_integrated_best = nested_current
    base_pool_path = args.base_unified_pool_csv or default_base_pool_for_scope(int(args.lib_id), str(args.ion_mode))
    if base_pool_path is not None and resolve(base_pool_path, repo_root).exists():
        summary = run_from_base_pool(args, repo_root, output_root, base_pool_path)
        print(json.dumps(summary, indent=2, sort_keys=True))
        return summary
    args.base_unified_pool_csv = None
    summary = base.run(args)
    pool_path = output_root / "unified_candidates_per_control.csv"
    pool = pd.read_csv(pool_path, low_memory=False) if pool_path.exists() else pd.DataFrame()
    nce_features = nce_feature_table(args.scored_candidates, repo_root)
    pool = overlay_nce_features(pool, nce_features)
    refresh_splits(output_root, pool, args, repo_root)
    summary = {
        **summary,
        "v16_overlay": True,
        "n_nce_feature_rows": int(nce_features.shape[0]),
        "n_pool_rows_with_nce_draw_count": (
            int(pd.to_numeric(pool.get("nce_draw_count"), errors="coerce").fillna(0.0).gt(0.0).sum())
            if not pool.empty
            else 0
        ),
    }
    write_json(output_root / "unified_candidate_pool_summary.json", summary)
    write_json(output_root / "unified_candidate_pool_v16_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def main() -> None:
    """Parse CLI arguments and invoke :func:`run`."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=Path("models/unified_candidate_pool_nce_sweep_v16"))
    parser.add_argument("--scored-candidates", type=Path, required=True)
    parser.add_argument("--extra-inputs", default=",".join(base.DEFAULT_EXTRA_INPUTS))
    parser.add_argument(
        "--current-integrated-best",
        type=Path,
        default=Path("models/integrated_denovo_audit_v15/integrated_per_control_best.csv"),
    )
    parser.add_argument("--include-all-x", action="store_true")
    parser.add_argument(
        "--all-x-catalog", type=Path, default=Path("models/x_known_curation_v1/known_x_evidence_catalog.csv")
    )
    parser.add_argument("--isotope-signals-csv", type=Path, default=None)
    parser.add_argument("--base-unified-pool-csv", type=Path, default=None)
    parser.add_argument("--isotope-filter-impossible", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--lib-id", type=int, required=True)
    parser.add_argument("--ion-mode", choices=["negative", "positive"], required=True)
    args = parser.parse_args()
    with stage_timer("pool", lib_id=int(args.lib_id), ion_mode=str(args.ion_mode)):
        run(args)
    raise SystemExit(0)


if __name__ == "__main__":
    main()
