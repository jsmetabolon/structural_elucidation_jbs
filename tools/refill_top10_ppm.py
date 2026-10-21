#!/usr/bin/env python3
"""Rebuild each library's displayed top-10 from the full scored pools.

Drops candidates whose mass error is missing or exceeds a reporting ceiling
(default 10 ppm), then keeps the next-best rerank survivors up to ``--top-k``.
Does not overwrite ``unified_rerank_all_top10.csv``.

Ppm uses the same conversion as ``tools/build_exec_summary.py``:
``abs(abs_mass_error_da) / mass_reference * 1e6``. Ranking uses the same
tie-breakers as ``pipeline.rerank.rank_frame``.

Usage (on the host that has the scored pools):

    python3 tools/refill_top10_ppm.py \\
        --equalization-root models/pipeline_equalization_v21
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd

REPORTING_PPM = 10.0
DEFAULT_LIBS = ("209", "305", "400", "402")
SCORE_COLUMNS = ("unified_final_score", "rerank_score", "rerank_score_calibrated")
SMILES_COLUMNS = ("predicted_smiles", "candidate_smiles", "smiles", "canonical_2d_smiles")
RANK_SORT_COLUMNS = (
    "predicted_spectrum_score_v2",
    "mass_match_score",
    "formula_prior_score",
    "min_source_rank",
    "predicted_smiles",
)
MASS_REFERENCE_KEYS = (
    "predicted_exact_mass",
    "candidate_exact_mass",
    "theoretical_neutral_mass",
    "mass_hypothesis_neutral_mass",
    "target_neutral_mass",
    "parent_mz",
)


def _num_series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(float("nan"), index=frame.index, dtype="float64")
    return pd.to_numeric(frame[column], errors="coerce")


def mass_reference(frame: pd.DataFrame) -> pd.Series:
    """First positive mass among the exec-summary reference keys."""

    out = pd.Series(float("nan"), index=frame.index, dtype="float64")
    for key in MASS_REFERENCE_KEYS:
        values = _num_series(frame, key)
        take = out.isna() & values.notna() & (values > 0)
        out = out.where(~take, values)
    return out


def abs_mass_error_ppm(error_da: pd.Series, reference: pd.Series) -> pd.Series:
    error = pd.to_numeric(error_da, errors="coerce").abs()
    mass = pd.to_numeric(reference, errors="coerce")
    ok = error.notna() & mass.notna() & (mass > 0)
    ppm = pd.Series(float("nan"), index=error.index, dtype="float64")
    return ppm.where(~ok, error / mass * 1e6)


def _proton_adduct(label: object) -> bool:
    text = str(label or "").replace("[", "").replace("]", "").replace(" ", "").lower().rstrip("+-")
    return text in {"", "nan", "none", "m+h", "m-h", "mh"}


def fill_missing_mass_error(frame: pd.DataFrame) -> pd.DataFrame:
    """Fill a blank mass error from the proton hypothesis.

    A candidate can match ``[M±H]`` and still have an empty
    ``abs_mass_error_da``. The ppm gate treats that blank as a miss. Non-default
    adduct labels are left blank so they do not become a 10 ppm hit.
    """

    out = frame.copy()
    if "abs_mass_error_da" not in out.columns:
        out["abs_mass_error_da"] = pd.NA
    error = pd.to_numeric(out["abs_mass_error_da"], errors="coerce")
    exact = _num_series(out, "predicted_exact_mass")
    if exact.isna().all():
        exact = _num_series(out, "candidate_exact_mass")
    hypothesis = _num_series(out, "mass_hypothesis_neutral_mass")
    target = _num_series(out, "target_neutral_mass")
    reference = hypothesis.where(hypothesis.notna() & (hypothesis > 0), target)
    missing = error.isna() & exact.notna() & reference.notna()
    if "mass_hypothesis_adduct" in out.columns:
        proton = out["mass_hypothesis_adduct"].map(_proton_adduct)
        missing = missing & proton
    out.loc[missing, "abs_mass_error_da"] = (exact - reference).abs()[missing]
    return out


def ensure_rank_columns(frame: pd.DataFrame, score_col: str) -> pd.DataFrame:
    """Fill columns ``rank_frame`` sorts on so heterogeneous lib CSVs still rank."""

    out = frame.copy()
    if "control_id" not in out.columns:
        if "x_compound_id" in out.columns:
            out["control_id"] = out["x_compound_id"]
        else:
            raise SystemExit("scored CSV needs control_id or x_compound_id")
    smiles_col = next((c for c in SMILES_COLUMNS if c in out.columns), None)
    if smiles_col is None:
        raise SystemExit("scored CSV needs a SMILES column")
    if "predicted_smiles" not in out.columns:
        out["predicted_smiles"] = out[smiles_col]
    if score_col not in out.columns:
        raise SystemExit(f"scored CSV needs {score_col}")
    for column in RANK_SORT_COLUMNS:
        if column not in out.columns:
            out[column] = 0.0 if column != "predicted_smiles" else ""
        elif column != "predicted_smiles":
            out[column] = pd.to_numeric(out[column], errors="coerce")
            if column == "min_source_rank":
                out[column] = out[column].fillna(10**9)
            else:
                out[column] = out[column].fillna(0.0)
    if "feature_key" not in out.columns:
        out["feature_key"] = ""
    out["feature_key"] = out["feature_key"].fillna("").astype(str)
    return out


def assign_ranks(frame: pd.DataFrame, score_col: str, policy: str) -> pd.DataFrame:
    """Sort like ``pipeline.rerank.rank_frame`` and assign 1-based per-control ranks."""

    out = ensure_rank_columns(frame, score_col)
    out["rerank_policy"] = policy
    out["rerank_score"] = pd.to_numeric(out[score_col], errors="coerce").fillna(0.0)
    missing_key = out["feature_key"].eq("")
    out.loc[missing_key, "feature_key"] = out.loc[missing_key, "control_id"].map(
        lambda value: f"unified_rerank_{value}"
    )
    out = out.sort_values(
        [
            "control_id",
            "rerank_score",
            "predicted_spectrum_score_v2",
            "mass_match_score",
            "formula_prior_score",
            "min_source_rank",
            "predicted_smiles",
        ],
        ascending=[True, False, False, False, False, True, True],
        kind="mergesort",
    )
    out["rank"] = out.groupby("control_id", dropna=False).cumcount() + 1
    return out.reset_index(drop=True)


def read_scored_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size <= 1:
        return pd.DataFrame()
    return pd.read_csv(path, low_memory=False)


def concat_scored(lib_dir: Path) -> pd.DataFrame:
    s0 = read_scored_csv(lib_dir / "unified_candidates_scored.csv")
    non_s0 = read_scored_csv(lib_dir / "non_s0_candidates_scored.csv")
    frames = [frame for frame in (s0, non_s0) if not frame.empty]
    if not frames:
        raise SystemExit(f"no scored CSVs under {lib_dir}")
    return pd.concat(frames, ignore_index=True, sort=False)


def score_column(frame: pd.DataFrame) -> str:
    for column in SCORE_COLUMNS:
        if column in frame.columns:
            return column
    raise SystemExit(f"none of {SCORE_COLUMNS} present")


def refill_lib(
    lib_dir: Path,
    *,
    ppm_ceiling: float,
    top_k: int,
    out_name: str,
    min_spectrum_score: float | None = None,
    out_path: Path | None = None,
) -> dict[str, Any]:
    """Filter one lib's scored pool and write a ppm-gated top-k CSV.

    When ``min_spectrum_score`` is set, a row must also have
    ``predicted_spectrum_score_v2`` strictly above that value. A mass match
    with no ion-intensity support is dropped.
    """

    scored = fill_missing_mass_error(concat_scored(lib_dir))
    score_col = score_column(scored)
    ranked = assign_ranks(scored, score_col, policy="unified_rerank_v1_ppm10")
    ranked["prefilter_rerank_rank"] = ranked["rank"]
    ranked["abs_mass_error_ppm"] = abs_mass_error_ppm(
        ranked.get("abs_mass_error_da", pd.Series(dtype="float64")),
        mass_reference(ranked),
    )
    ppm = pd.to_numeric(ranked["abs_mass_error_ppm"], errors="coerce")
    passing = ranked[ppm.notna() & (ppm <= float(ppm_ceiling))].copy()
    if min_spectrum_score is not None:
        spectrum = pd.to_numeric(passing.get("predicted_spectrum_score_v2"), errors="coerce")
        passing = passing[spectrum.notna() & (spectrum > float(min_spectrum_score))].copy()
    passing["rank"] = passing.groupby("control_id", dropna=False).cumcount() + 1
    gated = passing[passing["rank"].le(int(top_k))].reset_index(drop=True)

    old_path = lib_dir / "unified_rerank_all_top10.csv"
    old = read_scored_csv(old_path)
    old_rank1: dict[str, str] = {}
    old_counts: dict[str, int] = {}
    if not old.empty:
        id_col = "control_id" if "control_id" in old.columns else "x_compound_id"
        old[id_col] = old[id_col].astype(str)
        old["rank"] = pd.to_numeric(old.get("rank"), errors="coerce")
        smiles_col = next((c for c in SMILES_COLUMNS if c in old.columns), None)
        old_counts = old.groupby(id_col).size().to_dict()
        if smiles_col:
            for control_id, group in old.groupby(id_col, dropna=False):
                row = group.sort_values("rank", kind="mergesort").iloc[0]
                old_rank1[str(control_id)] = str(row[smiles_col])

    new_counts = gated.groupby("control_id").size()
    n_compounds = int(ranked["control_id"].nunique())
    n_new = int(gated["control_id"].nunique()) if not gated.empty else 0
    n_short = int((new_counts < top_k).sum()) if n_new else n_compounds
    n_empty = n_compounds - n_new
    n_from_deeper = int((gated["prefilter_rerank_rank"] > top_k).sum()) if not gated.empty else 0
    n_rank1_changed = 0
    if old_rank1 and not gated.empty:
        new_rank1 = (
            gated[gated["rank"].eq(1)]
            .set_index("control_id")["predicted_smiles"]
            .astype(str)
            .to_dict()
        )
        for control_id, smiles in old_rank1.items():
            if new_rank1.get(control_id, "") != smiles:
                n_rank1_changed += 1

    dest = out_path if out_path is not None else lib_dir / out_name
    dest.parent.mkdir(parents=True, exist_ok=True)
    gated.to_csv(dest, index=False)
    out_path = dest
    summary = {
        "lib_dir": str(lib_dir),
        "out_path": str(out_path),
        "n_scored_rows": int(len(ranked)),
        "n_passing_rows": int(len(passing)),
        "n_written_rows": int(len(gated)),
        "n_compounds_scored": n_compounds,
        "n_compounds_written": n_new,
        "n_short_of_top_k": n_short,
        "n_empty": n_empty,
        "n_rows_from_rank_gt_top_k": n_from_deeper,
        "n_rank1_smiles_changed": n_rank1_changed,
        "n_old_top10_rows": int(len(old)),
        "n_old_compounds": len(old_counts),
        "ppm_ceiling": float(ppm_ceiling),
        "top_k": int(top_k),
    }
    return summary


def verify_full_passers(
    lib_dir: Path,
    gated: pd.DataFrame,
    verify_path: Path,
    top_k: int,
) -> list[str]:
    """Compounds that already had k passers in the old top-10 should not change."""

    import csv

    problems: list[str] = []
    with verify_path.open(newline="") as handle:
        by_compound = list(csv.DictReader(handle))
    full = {
        row["x_compound_id"]
        for row in by_compound
        if row.get("n_le10_ppm") == str(top_k)
    }
    if not full:
        return problems
    old = read_scored_csv(lib_dir / "unified_rerank_all_top10.csv")
    if old.empty:
        return ["missing unified_rerank_all_top10.csv for verify"]
    id_col = "control_id" if "control_id" in old.columns else "x_compound_id"
    smiles_col = next((c for c in SMILES_COLUMNS if c in old.columns), "predicted_smiles")
    old[id_col] = old[id_col].astype(str)
    gated = gated.copy()
    gated["control_id"] = gated["control_id"].astype(str)
    present = set(old[id_col].unique())
    for xid in sorted(full & present):
        old_smiles = tuple(
            old.loc[old[id_col].eq(xid)]
            .sort_values("rank", kind="mergesort")[smiles_col]
            .astype(str)
            .tolist()
        )
        new_smiles = tuple(
            gated.loc[gated["control_id"].eq(xid)]
            .sort_values("rank", kind="mergesort")[smiles_col]
            .astype(str)
            .tolist()
        )
        if old_smiles != new_smiles:
            problems.append(f"{xid}: expected {old_smiles!r} got {new_smiles!r}")
    return problems


def unique_coverage(root: Path, libs: list[str], out_name: str, top_k: int) -> dict[str, int]:
    """First-library-wins counts of compounds with 10, fewer than 10, or no rows.

    A compound present in an earlier library is not counted again. Empty means
    the scored pool had no candidate at or under the ppm ceiling.
    """

    seen: set[str] = set()
    full = short = empty = 0
    for lib in libs:
        path = root / f"rerank_lib{lib}" / out_name
        frame = read_scored_csv(path)
        if frame.empty or "control_id" not in frame.columns:
            continue
        frame["control_id"] = frame["control_id"].astype(str)
        counts = frame.groupby("control_id").size()
        for control_id, count in counts.items():
            if control_id in seen:
                continue
            seen.add(str(control_id))
            if int(count) >= top_k:
                full += 1
            elif int(count) > 0:
                short += 1
        # Compounds scored in this lib but absent from the gated file are empty
        # only when no earlier lib already claimed them. The scored pool is the
        # source of that set.
        scored = concat_scored(root / f"rerank_lib{lib}")
        id_col = "control_id" if "control_id" in scored.columns else "x_compound_id"
        if id_col not in scored.columns:
            continue
        for control_id in scored[id_col].astype(str).unique():
            if control_id in seen:
                continue
            seen.add(control_id)
            empty += 1
    return {"full": full, "short": short, "empty": empty, "compounds": full + short + empty}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--equalization-root",
        type=Path,
        required=True,
        help="Directory holding rerank_lib*/ scored CSVs",
    )
    parser.add_argument("--libs", default=",".join(DEFAULT_LIBS))
    parser.add_argument("--ppm-ceiling", type=float, default=REPORTING_PPM)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument(
        "--min-spectrum-score",
        type=float,
        default=None,
        help="Require predicted_spectrum_score_v2 strictly above this value",
    )
    parser.add_argument(
        "--out-name",
        default="unified_rerank_all_top10_ppm10.csv",
        help="Filename written next to the ungated top-10 (not overwritten)",
    )
    parser.add_argument(
        "--verify-against",
        type=Path,
        default=None,
        help="Optional audit_ppm_qc_filtered_by_compound.csv",
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=None,
        help="Write the per-lib summaries as JSON",
    )
    args = parser.parse_args(argv)

    root = args.equalization_root.expanduser().resolve()
    if not root.is_dir():
        raise SystemExit(f"equalization root not found: {root}")

    summaries: list[dict[str, Any]] = []
    for lib in [item.strip() for item in str(args.libs).split(",") if item.strip()]:
        lib_dir = root / f"rerank_lib{lib}"
        if not lib_dir.is_dir():
            raise SystemExit(f"missing {lib_dir}")
        print(f"[refill] lib{lib} {lib_dir}", flush=True)
        summary = refill_lib(
            lib_dir,
            ppm_ceiling=args.ppm_ceiling,
            top_k=args.top_k,
            out_name=args.out_name,
            min_spectrum_score=args.min_spectrum_score,
        )
        summaries.append(summary)
        print(
            f"  compounds {summary['n_compounds_written']}/{summary['n_compounds_scored']} "
            f"short={summary['n_short_of_top_k']} empty={summary['n_empty']} "
            f"deeper_rows={summary['n_rows_from_rank_gt_top_k']} "
            f"rank1_changed={summary['n_rank1_smiles_changed']}",
            flush=True,
        )
        if args.verify_against:
            gated = pd.read_csv(lib_dir / args.out_name, low_memory=False)
            problems = verify_full_passers(
                lib_dir, gated, args.verify_against, args.top_k
            )
            if problems:
                preview = "\n".join(problems[:8])
                raise SystemExit(
                    f"verify failed for lib{lib} ({len(problems)} compounds):\n{preview}"
                )
            print(f"  verify-against ok ({args.verify_against.name})", flush=True)

    totals = {
        "n_written_rows": sum(s["n_written_rows"] for s in summaries),
        "n_compounds_written": sum(s["n_compounds_written"] for s in summaries),
        "n_short_of_top_k": sum(s["n_short_of_top_k"] for s in summaries),
        "n_empty": sum(s["n_empty"] for s in summaries),
        "n_rows_from_rank_gt_top_k": sum(s["n_rows_from_rank_gt_top_k"] for s in summaries),
        "n_rank1_smiles_changed": sum(s["n_rank1_smiles_changed"] for s in summaries),
        "libs": summaries,
    }
    libs = [item.strip() for item in str(args.libs).split(",") if item.strip()]
    coverage = unique_coverage(root, libs, args.out_name, args.top_k)
    totals["coverage"] = coverage
    print(
        "[refill] coverage "
        f"full={coverage['full']} short={coverage['short']} empty={coverage['empty']} "
        f"compounds={coverage['compounds']}",
        flush=True,
    )
    print("[refill] TOTAL", json.dumps({k: v for k, v in totals.items() if k != "libs"}))
    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(totals, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
