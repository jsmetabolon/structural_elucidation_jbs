#!/usr/bin/env python3
"""Rescore every exec-summary top-10 candidate with the pure NCE matcher.

The delivered ``evidence_confidence`` saturates near 1.0, so rival candidate
families tie and the identity margin carries little signal. This scores each
(compound, SMILES) pair directly with the NCE matcher against the compound's
observed bins. Successful zero scores and unavailable scores have distinct
statuses; the portable overlap fallback is never substituted.

Output columns: x_compound_id, candidate_smiles, pure_nce_score, nce_status,
where nce_status is one of scored, zero, no_observed_bins, error.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing
import math
import os
import sys
from functools import lru_cache
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

def load_observed(lookup_dir: Path, ids: set[str]) -> pd.DataFrame:
    """Load observed bins, retaining the first available library per compound."""

    columns = ["control_id", "ion_json", "loss_json", "lib_id"]
    frames = []
    for lib, ion in ((209, "negative"), (305, "negative"), (400, "positive"), (402, "positive")):
        path = lookup_dir / f"lib{lib}_{ion}_observed_ms2_lookup.csv"
        if not path.exists():
            continue
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
        if "control_id" not in frame:
            continue
        frame = frame.loc[frame["control_id"].isin(ids)].copy()
        frame["lib_id"] = lib
        for original, renamed in (("observed_ion_bins_json", "ion_json"), ("observed_loss_bins_json", "loss_json")):
            frame[renamed] = frame[original] if original in frame else "[]"
        frames.append(frame[columns])
    if not frames:
        return pd.DataFrame(columns=columns)
    return pd.concat(frames, ignore_index=True).drop_duplicates("control_id", keep="first")


OUT_COLUMNS = ["x_compound_id", "candidate_smiles", "pure_nce_score", "nce_status"]
SCORED = "scored"
ZERO = "zero"
NO_BINS = "no_observed_bins"
ERROR = "error"


def _has_bins(ion_json: str, loss_json: str) -> bool:
    for text in (ion_json, loss_json):
        try:
            values = json.loads(text or "[]")
        except Exception:
            continue
        if isinstance(values, list) and any(str(v).strip() not in {"", "0"} for v in values):
            return True
    return False


def _valid_smiles(smiles: str) -> bool:
    try:
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")
        return Chem.MolFromSmiles(smiles) is not None
    except ImportError:
        return bool(smiles)


@lru_cache(maxsize=8)
def matcher_available(lib_id: int) -> bool:
    """predicted_spectrum_score_v2 turns a matcher load failure into 0, so check it first."""
    from structural_elucidation.nce.runtime import DEFAULT_REFERENCE_DIR, get_predicted_spectrum_matcher

    reference_dir = os.environ.get("NCE_REFERENCE_DIR", "").strip() or str(DEFAULT_REFERENCE_DIR)
    try:
        get_predicted_spectrum_matcher(int(lib_id), reference_dir=reference_dir)
    except Exception as exc:
        print(f"[rescore] matcher unavailable for lib {lib_id}: {type(exc).__name__}: {exc}", flush=True)
        return False
    return True


def score_key(key: tuple[str, str, str, int]) -> tuple[float, str]:
    """Return (score, status) for one (smiles, ion_json, loss_json, lib_id) key."""
    from structural_elucidation.nce.runtime import (
        DEFAULT_REFERENCE_DIR,
        get_predicted_spectrum_matcher,
        observed_bin_spectrum,
    )
    from structural_elucidation.nce.evidence import parse_json_list

    smiles, ion_json, loss_json, lib_id = key
    spectrum = observed_bin_spectrum(parse_json_list(ion_json), parse_json_list(loss_json))
    if not spectrum:
        return 0.0, NO_BINS
    if not _valid_smiles(smiles):
        return 0.0, ERROR
    reference_dir = os.environ.get("NCE_REFERENCE_DIR", "").strip() or str(DEFAULT_REFERENCE_DIR)
    try:
        # Both verifier_v2 and its runtime helper turn matcher exceptions into
        # zero. Calling the matcher here preserves that distinction for scoring.
        matcher = get_predicted_spectrum_matcher(int(lib_id), reference_dir=reference_dir)
        model_bins = set(spectrum) & set(matcher.model_dict)
        score = float(matcher.score_a_smiles_for_spectrum(smiles, spectrum, model_bins))
    except Exception:
        return 0.0, ERROR
    if not math.isfinite(score) or not 0.0 <= score <= 1.0:
        return 0.0, ERROR
    return score, SCORED if score > 0 else ZERO


def score_keys(keys: list[tuple[str, str, str, int]], workers: int) -> dict:
    """Score keys serially or on a fork pool; the parent warms the matcher first."""
    results: dict = {}
    if not keys:
        return results
    head, rest = keys[0], keys[1:]
    results[head] = score_key(head)
    if workers <= 1 or not rest:
        for key in rest:
            results[key] = score_key(key)
        return results
    ctx = multiprocessing.get_context("fork")
    n = min(workers, len(rest))
    with ctx.Pool(processes=n) as pool:
        for key, value in zip(rest, pool.map(score_key, rest, chunksize=max(1, len(rest) // (n * 4)))):
            results[key] = value
    return results


def build_pairs(top10: Path) -> pd.DataFrame:
    frame = pd.read_csv(top10, dtype=str, keep_default_na=False, usecols=["x_compound_id", "candidate_smiles"])
    frame = frame[frame["candidate_smiles"].str.strip().ne("")]
    return frame.drop_duplicates(["x_compound_id", "candidate_smiles"]).reset_index(drop=True)


def run(args: argparse.Namespace) -> pd.DataFrame:
    pairs = build_pairs(args.top10)
    done = pd.DataFrame(columns=OUT_COLUMNS)
    if args.resume and args.out.exists():
        done = pd.read_csv(args.out, dtype=str, keep_default_na=False)
        seen = set(zip(done["x_compound_id"], done["candidate_smiles"]))
        pairs = pairs[[pair not in seen for pair in zip(pairs["x_compound_id"], pairs["candidate_smiles"])]]
    if args.limit:
        pairs = pairs.head(args.limit)
    observed = load_observed(args.lookup_dir, set(pairs["x_compound_id"]))
    pairs = pairs.merge(observed[["control_id", "ion_json", "loss_json", "lib_id"]],
                        left_on="x_compound_id", right_on="control_id", how="left")
    pairs[["ion_json", "loss_json"]] = pairs[["ion_json", "loss_json"]].fillna("[]")
    pairs["lib_id"] = pd.to_numeric(pairs["lib_id"], errors="coerce").fillna(209).astype(int)
    keys = list(zip(pairs["candidate_smiles"], pairs["ion_json"], pairs["loss_json"], pairs["lib_id"]))
    unique = list(dict.fromkeys(keys))
    print(f"[rescore] pairs={len(pairs)} unique_keys={len(unique)} workers={args.workers}", flush=True)
    scored = score_keys(unique, args.workers)
    pairs["pure_nce_score"] = [scored[key][0] for key in keys]
    pairs["nce_status"] = [scored[key][1] for key in keys]
    out = pd.concat([done, pairs[OUT_COLUMNS]], ignore_index=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print("[rescore] status", json.dumps(pairs["nce_status"].value_counts().to_dict()), flush=True)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--top10", type=Path, required=True, help="exec_summary_top10.csv")
    parser.add_argument("--lookup-dir", type=Path, required=True, help="dir with lib*_observed_ms2_lookup.csv")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=0, help="score only the first N pairs (smoke test)")
    parser.add_argument("--resume", action="store_true", help="skip pairs already present in --out")
    parser.add_argument("--workers", type=int, default=1)
    run(parser.parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
