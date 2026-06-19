#!/usr/bin/env python3
"""Reusable control benchmark harness for de novo prediction evaluation.

Compares prediction outputs against known controls from known_X.csv and
known_X_all_available_structures_simple.csv, reporting exact-match and
top-k Tanimoto metrics per control X-ID.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFingerprintGenerator
    from rdkit.DataStructs import TanimotoSimilarity

    HAS_RDKIT = True
    MORGAN_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
except ImportError:
    HAS_RDKIT = False
    MORGAN_GENERATOR = None


def canonicalize(smiles: str) -> Optional[str]:
    """Canonicalize a SMILES string with RDKit sanitization."""
    if not HAS_RDKIT or not smiles or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return None


def tanimoto(smi_a: str, smi_b: str, radius: int = 2) -> float:
    """Compute Morgan fingerprint Tanimoto similarity between two SMILES."""
    if not HAS_RDKIT:
        return 0.0
    try:
        mol_a = Chem.MolFromSmiles(smi_a)
        mol_b = Chem.MolFromSmiles(smi_b)
        if mol_a is None or mol_b is None:
            return 0.0
        if MORGAN_GENERATOR is not None and int(radius) == 2:
            fp_a = MORGAN_GENERATOR.GetFingerprint(mol_a)
            fp_b = MORGAN_GENERATOR.GetFingerprint(mol_b)
        else:
            fp_a = AllChem.GetMorganFingerprintAsBitVect(mol_a, int(radius), nBits=2048)
            fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, int(radius), nBits=2048)
        return float(TanimotoSimilarity(fp_a, fp_b))
    except Exception:
        return 0.0


def extract_xid(feature_key: str) -> Optional[str]:
    """Extract X-compound ID (e.g. 'X-12345') from a feature_key string."""
    m = re.search(r"(X-\d+)", feature_key)
    return m.group(1) if m else None


def load_controls(
    known_x_path: Path,
    known_x_all_path: Optional[Path] = None,
    control_lib_id: Optional[int] = None,
    fail_on_mixed_lib_id: bool = False,
) -> Tuple[Dict[str, List[str]], dict]:
    """Load control X-IDs mapped to their canonical SMILES (may have multiple).

    Returns {x_compound_id: [canonical_smiles, ...]}.
    """
    known_x_df = pd.DataFrame()
    control_scope: set[str] | None = None
    scope_summary = {
        "known_x_path": str(known_x_path),
        "control_lib_id": int(control_lib_id) if control_lib_id is not None else None,
        "known_x_rows": 0,
        "known_x_unique_lib_ids": [],
        "mixed_lib_ids_detected": False,
        "control_scope_xids": 0,
    }
    if known_x_path.exists():
        known_x_df = pd.read_csv(known_x_path)
        known_x_df.columns = [c.strip().strip('"').lstrip("\ufeff") for c in known_x_df.columns]
        scope_summary["known_x_rows"] = int(known_x_df.shape[0])
        if {"LIB_ID", "Unknown"}.issubset(set(known_x_df.columns)):
            lib_ids = sorted(
                {int(v) for v in pd.to_numeric(known_x_df["LIB_ID"], errors="coerce").dropna().astype(int).tolist()}
            )
            scope_summary["known_x_unique_lib_ids"] = lib_ids
            scope_summary["mixed_lib_ids_detected"] = len(lib_ids) > 1
            if control_lib_id is not None:
                scoped_df = known_x_df[
                    pd.to_numeric(known_x_df["LIB_ID"], errors="coerce").eq(int(control_lib_id))
                ].copy()
                control_scope = {
                    str(value).strip() for value in scoped_df["Unknown"].astype(str).tolist() if str(value).strip()
                }
                scope_summary["control_scope_xids"] = int(len(control_scope))
            else:
                if scope_summary["mixed_lib_ids_detected"]:
                    message = (
                        f"Warning: mixed LIB_ID controls detected in {known_x_path}. "
                        "Set --control_lib_id for mode-specific evaluation."
                    )
                    print(message)
                    if fail_on_mixed_lib_id:
                        raise RuntimeError(message)
                control_scope = {
                    str(value).strip() for value in known_x_df["Unknown"].astype(str).tolist() if str(value).strip()
                }
                scope_summary["control_scope_xids"] = int(len(control_scope))

    controls: Dict[str, List[str]] = {}

    for path in [known_x_path, known_x_all_path]:
        if path is None or not path.exists():
            continue
        df = pd.read_csv(path)
        df.columns = [c.strip().strip('"').lstrip("\ufeff") for c in df.columns]
        for _, row in df.iterrows():
            xid = str(row.get("x_compound_id", "")).strip()
            raw_smi = str(row.get("smiles", "")).strip()
            if not xid or not raw_smi or raw_smi.lower() == "nan":
                continue
            if control_scope is not None and xid not in control_scope:
                continue
            can = canonicalize(raw_smi)
            if can is None:
                continue
            controls.setdefault(xid, [])
            if can not in controls[xid]:
                controls[xid].append(can)

    scope_summary["controls_loaded"] = int(len(controls))
    return controls, scope_summary


def load_predictions(
    pred_path: Path,
    smiles_col: str = "predicted_smiles",
    feature_key_col: str = "feature_key",
    rank_col: str = "rank",
) -> pd.DataFrame:
    """Load and normalize a prediction CSV for control evaluation."""
    if not pred_path.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(pred_path)
    except pd.errors.EmptyDataError:
        print(f"  Warning: {pred_path} is empty; treating as no predictions.")
        return pd.DataFrame()
    df.columns = [c.strip() for c in df.columns]
    if smiles_col == "predicted_smiles" and smiles_col not in df.columns:
        if "smiles" in df.columns:
            smiles_col = "smiles"
    needed = {feature_key_col, smiles_col}
    if not needed.issubset(set(df.columns)):
        print(f"  Warning: {pred_path} missing columns {needed - set(df.columns)}")
        return pd.DataFrame()
    df = df.rename(columns={smiles_col: "_smiles", feature_key_col: "_feature_key"})
    if rank_col in df.columns:
        df = df.rename(columns={rank_col: "_rank"})
    elif rank_col == "merged_rank" and "rank" in df.columns:
        df = df.rename(columns={"rank": "_rank"})
    else:
        df["_rank"] = range(1, len(df) + 1)
    df["_xid"] = df["_feature_key"].apply(extract_xid)
    return df


def evaluate_controls(
    controls: Dict[str, List[str]],
    preds_df: pd.DataFrame,
    track_name: str,
    top_k: int = 10,
) -> Tuple[List[dict], dict]:
    """Evaluate predictions against controls.

    Returns (per_control_rows, summary_dict).
    """
    per_control: List[dict] = []
    n_exact = 0
    n_topk_tanimoto_above = {0.5: 0, 0.7: 0, 0.9: 0}
    n_controls_covered = 0
    all_best_tanimotos: List[float] = []

    for xid, control_smiles_list in sorted(controls.items()):
        xid_preds = preds_df[preds_df["_xid"] == xid].copy()
        if xid_preds.empty:
            per_control.append(
                {
                    "x_compound_id": xid,
                    "track": track_name,
                    "n_predictions": 0,
                    "exact_match": False,
                    "exact_match_rank": None,
                    "best_tanimoto": 0.0,
                    "best_tanimoto_rank": None,
                    "top10_mean_tanimoto": 0.0,
                    "failure_bucket": "no_predictions",
                }
            )
            continue

        n_controls_covered += 1
        xid_preds = xid_preds.sort_values("_rank").head(top_k)

        exact_found = False
        exact_rank = None
        best_tan = 0.0
        best_tan_rank = None
        tan_scores: List[float] = []

        for _, row in xid_preds.iterrows():
            pred_smi = str(row["_smiles"]).strip()
            pred_can = canonicalize(pred_smi)
            rank = int(row["_rank"])

            if pred_can is None:
                tan_scores.append(0.0)
                continue

            max_tan_for_this = 0.0
            is_exact = False
            for ctrl_smi in control_smiles_list:
                if pred_can == ctrl_smi:
                    is_exact = True
                t = tanimoto(pred_can, ctrl_smi)
                max_tan_for_this = max(max_tan_for_this, t)

            tan_scores.append(max_tan_for_this)

            if is_exact and not exact_found:
                exact_found = True
                exact_rank = rank
            if max_tan_for_this > best_tan:
                best_tan = max_tan_for_this
                best_tan_rank = rank

        if exact_found:
            n_exact += 1

        for thresh in list(n_topk_tanimoto_above.keys()):
            if best_tan >= thresh:
                n_topk_tanimoto_above[thresh] += 1

        all_best_tanimotos.append(best_tan)
        mean_tan = float(np.mean(tan_scores)) if tan_scores else 0.0

        if exact_found:
            bucket = "exact_match"
        elif best_tan >= 0.7:
            bucket = "structural_hit"
        elif best_tan >= 0.3:
            bucket = "low_similarity"
        else:
            bucket = "structural_miss"

        per_control.append(
            {
                "x_compound_id": xid,
                "track": track_name,
                "n_predictions": len(xid_preds),
                "exact_match": exact_found,
                "exact_match_rank": exact_rank,
                "best_tanimoto": round(best_tan, 4),
                "best_tanimoto_rank": best_tan_rank,
                "top10_mean_tanimoto": round(mean_tan, 4),
                "failure_bucket": bucket,
            }
        )

    n_total = len(controls)
    summary = {
        "track": track_name,
        "n_controls": n_total,
        "n_controls_covered": n_controls_covered,
        "n_exact_match": n_exact,
        "exact_match_rate": round(n_exact / max(n_total, 1), 4),
        "n_topk_tanimoto_ge_0.5": n_topk_tanimoto_above[0.5],
        "n_topk_tanimoto_ge_0.7": n_topk_tanimoto_above[0.7],
        "n_topk_tanimoto_ge_0.9": n_topk_tanimoto_above[0.9],
        "topk_tanimoto_ge_0.7_rate": round(n_topk_tanimoto_above[0.7] / max(n_total, 1), 4),
        "mean_best_tanimoto": round(float(np.mean(all_best_tanimotos)) if all_best_tanimotos else 0.0, 4),
        "median_best_tanimoto": round(float(np.median(all_best_tanimotos)) if all_best_tanimotos else 0.0, 4),
    }
    return per_control, summary


def main() -> None:
    """CLI entry point for control-benchmark evaluation."""
    parser = argparse.ArgumentParser(description="Evaluate de novo and retrieval predictions against known controls.")
    parser.add_argument(
        "--prediction_dir",
        type=str,
        default="models/x209_smiles_predictions",
        help="Directory containing prediction CSV outputs.",
    )
    parser.add_argument("--known_x", type=str, default="known_X.csv")
    parser.add_argument("--known_x_all", type=str, default="known_X_all_available_structures_simple.csv")
    parser.add_argument(
        "--denovo_csv",
        type=str,
        default="",
        help="Override denovo report path (default: auto-detect in prediction_dir).",
    )
    parser.add_argument("--retrieval_csv", type=str, default="", help="Override retrieval report path.")
    parser.add_argument("--merged_csv", type=str, default="", help="Override merged report path.")
    parser.add_argument(
        "--output_dir", type=str, default="", help="Output directory (default: same as prediction_dir)."
    )
    parser.add_argument("--top_k", type=int, default=10)
    parser.add_argument(
        "--control_lib_id",
        type=int,
        default=None,
        help="Optional LIB_ID filter applied through known_X.csv scope selection.",
    )
    parser.add_argument(
        "--fail_on_mixed_lib_id",
        action="store_true",
        help="Fail fast when known_X.csv contains multiple LIB_ID values and --control_lib_id is not set.",
    )
    args = parser.parse_args()

    pred_dir = Path(args.prediction_dir).expanduser().resolve()
    out_dir = Path(args.output_dir or args.prediction_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading controls...")
    controls, control_scope_summary = load_controls(
        Path(args.known_x).expanduser().resolve(),
        Path(args.known_x_all).expanduser().resolve(),
        control_lib_id=args.control_lib_id,
        fail_on_mixed_lib_id=bool(args.fail_on_mixed_lib_id),
    )
    print(f"  {len(controls)} control X-IDs with canonical SMILES")

    denovo_path = (
        Path(args.denovo_csv).expanduser().resolve() if args.denovo_csv else pred_dir / "x209_denovo_report.csv"
    )
    retrieval_path = (
        Path(args.retrieval_csv).expanduser().resolve()
        if args.retrieval_csv
        else pred_dir / "x209_retrieval_report.csv"
    )
    merged_path = (
        Path(args.merged_csv).expanduser().resolve() if args.merged_csv else pred_dir / "x209_merged_predictions.csv"
    )

    all_per_control: List[dict] = []
    all_summaries: List[dict] = []

    for path, track, smi_col, rank_col in [
        (denovo_path, "de_novo", "predicted_smiles", "rank"),
        (retrieval_path, "retrieval", "smiles", "rank"),
        (merged_path, "merged", "smiles", "merged_rank"),
    ]:
        if not path.exists():
            print(f"  Skipping {track}: {path} not found")
            continue
        print(f"\nEvaluating {track} from {path.name}...")
        preds = load_predictions(path, smiles_col=smi_col, rank_col=rank_col)
        if preds.empty:
            print(f"  No predictions loaded for {track}")
            continue
        per_ctrl, summary = evaluate_controls(
            controls,
            preds,
            track_name=track,
            top_k=args.top_k,
        )
        all_per_control.extend(per_ctrl)
        all_summaries.append(summary)
        print(
            f"  {track}: exact={summary['n_exact_match']}, "
            f"tan>=0.7={summary['n_topk_tanimoto_ge_0.7']}, "
            f"mean_best_tan={summary['mean_best_tanimoto']:.4f}"
        )

    per_ctrl_df = pd.DataFrame(all_per_control)
    per_ctrl_path = out_dir / "control_eval_per_xid.csv"
    per_ctrl_df.to_csv(per_ctrl_path, index=False)

    combined_summary = {
        "controls_loaded": len(controls),
        "control_scope": control_scope_summary,
        "tracks": all_summaries,
    }
    summary_path = out_dir / "control_eval_summary.json"
    summary_path.write_text(json.dumps(combined_summary, indent=2))

    print(f"\nWrote: {per_ctrl_path}")
    print(f"Wrote: {summary_path}")

    if per_ctrl_df.empty or "track" not in per_ctrl_df.columns or "failure_bucket" not in per_ctrl_df.columns:
        print("\nFailure buckets:\n  no evaluated control rows")
        return

    bucket_counts = per_ctrl_df.groupby(["track", "failure_bucket"]).size()
    print(f"\nFailure buckets:\n{bucket_counts.to_string()}")


if __name__ == "__main__":
    main()
