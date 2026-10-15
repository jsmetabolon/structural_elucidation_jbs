#!/usr/bin/env python3
"""Cohort the unpublished compounds and pilot same-formula isomers.

A compound is unpublished when the merged delivery marks it ``no_prediction``.
An on-mass seed is the best stored-mass candidate within 10 ppm, or the best
adduct-fallback candidate within 10 ppm when stored mass has none. Isomers are
scored with ``predicted_spectrum_score_v2`` against the compound's observed bins.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import pandas as pd

from refill_top10_ppm import MASS_REFERENCE_KEYS, abs_mass_error_ppm, mass_reference

LIBS = (209, 305, 400, 402)
ION = {209: "negative", 305: "negative", 400: "positive", 402: "positive"}
SCORED_NAMES = ("unified_candidates_scored.csv", "non_s0_candidates_scored.csv")
POOL_COLUMNS = {
    "control_id",
    "predicted_smiles",
    "candidate_smiles",
    "smiles",
    "predicted_formula",
    "candidate_formula",
    "abs_mass_error_da",
    "abs_mass_error_ppm",
    "predicted_spectrum_score_v2",
    "lib_id",
    *MASS_REFERENCE_KEYS,
}
NO_OBSERVED = "no_observed_ions"
ON_MASS = "on_mass_seed"
OFF_MASS = "off_mass_only"
NO_CANDIDATE = "no_candidate"
PPM_CEILING = 10.0

# Formula-preserving transforms. Products that fail sanitization or change the
# formula are dropped, so a bad match cannot enter the pilot.
FORMULA_PRESERVING_SMIRKS = (
    "[C:1]=[C:2]-[C:3]-[C:4]>>[C:1]-[C:2]=[C:3]-[C:4]",
    "[C:1]-[C:2]=[C:3]-[C:4]>>[C:1]=[C:2]-[C:3]-[C:4]",
    "[CH2:1][CH2:2][OH:3]>>[CH3:1][O:3][CH3:2]",
    "[CH3:1][O:2][CH3:3]>>[CH3:1][CH2:3][OH:2]",
    "[C:1][CH2:2][OH:3]>>[C:1]([OH:3])[CH3:2]",
    "[C:1]([OH:2])[CH2:3][C:4]>>[C:1][CH:3]([OH:2])[C:4]",
    "[C:1][CH2:2][C:3](=[O:4])>>[C:1][C:2](=[O:4])[CH3:3]",
    "[C:1][C:2](=[O:3])[CH3:4]>>[C:1][CH2:2][CH:4]=[O:3]",
    "[C:1][NH:2][C:3]=[O:4]>>[C:1][N:2]=[C:3][OH:4]",
    "[C:1](=[O:2])[O:3][S:4](=[O:5])(=[O:6])[OH:7].[C:8][OH:9]>>[C:1](=[O:2])[OH:9].[C:8][O:3][S:4](=[O:5])(=[O:6])[OH:7]",
    "[C:1][OH:2].[C:3][OH:4]>>[C:1][O:2][C:3].[OH:4]",
)


def assign_cohort(
    *,
    no_inputs: bool,
    has_ion_bins: bool,
    on_mass: bool,
    has_candidate: bool,
) -> str:
    """One label per compound. Compounds with no feature row stay ``no_candidate``."""

    if no_inputs:
        return NO_CANDIDATE
    if not has_ion_bins:
        return NO_OBSERVED
    if on_mass:
        return ON_MASS
    if has_candidate:
        return OFF_MASS
    return NO_CANDIDATE


def seed_pool(stored_ppm: float | None, fallback_ppm: float | None) -> str:
    """Stored mass wins when it already has a candidate within 10 ppm."""

    if stored_ppm is not None and stored_ppm <= PPM_CEILING:
        return "stored"
    if fallback_ppm is not None and fallback_ppm <= PPM_CEILING:
        return "adduct_fallback"
    return ""


def _finite(value: float | None) -> float | None:
    if value is None or value != value:
        return None
    return float(value)


def _smiles_of(row: pd.Series) -> str:
    for column in ("predicted_smiles", "candidate_smiles", "smiles"):
        text = str(row.get(column, "") or "").strip()
        if text and text.lower() not in {"nan", "none"}:
            return text
    return ""


def _formula_of(row: pd.Series) -> str:
    for column in ("predicted_formula", "candidate_formula"):
        text = str(row.get(column, "") or "").strip()
        if text and text.lower() not in {"nan", "none"}:
            return text
    return ""


def _ppm_series(frame: pd.DataFrame) -> pd.Series:
    if "abs_mass_error_ppm" in frame.columns:
        parsed = pd.to_numeric(frame["abs_mass_error_ppm"], errors="coerce").abs()
        if parsed.notna().any():
            return parsed
    error = pd.to_numeric(frame.get("abs_mass_error_da"), errors="coerce")
    return abs_mass_error_ppm(error, mass_reference(frame))


def best_rows(root: Path, ids: set[str]) -> pd.DataFrame:
    """Best (lowest ppm) scored row per compound under one equalization root."""

    frames: list[pd.DataFrame] = []
    for lib in LIBS:
        for name in SCORED_NAMES:
            path = root / f"rerank_lib{lib}" / name
            if not path.exists() or path.stat().st_size < 50:
                continue
            frame = pd.read_csv(
                path,
                dtype=str,
                keep_default_na=False,
                low_memory=False,
                usecols=lambda column: column in POOL_COLUMNS,
            )
            if "control_id" not in frame.columns or frame.empty:
                continue
            frame = frame[frame["control_id"].isin(ids)]
            if frame.empty:
                continue
            frame["lib_id"] = str(lib)
            frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["control_id", "ppm", "smiles", "formula", "lib_id", "spectrum"])
    pool = pd.concat(frames, ignore_index=True)
    pool["ppm"] = _ppm_series(pool)
    pool["smiles"] = pool.apply(_smiles_of, axis=1)
    pool["formula"] = pool.apply(_formula_of, axis=1)
    pool["spectrum"] = pd.to_numeric(pool.get("predicted_spectrum_score_v2"), errors="coerce")
    pool = pool.dropna(subset=["ppm"]).sort_values(["control_id", "ppm"], kind="mergesort")
    best = pool.drop_duplicates("control_id", keep="first")
    return best[["control_id", "ppm", "smiles", "formula", "lib_id", "spectrum"]].reset_index(drop=True)


def load_observed(lookup_dir: Path, ids: set[str]) -> pd.DataFrame:
    """First library wins. Peak count is the number of observed MS2 pairs."""

    frames: list[pd.DataFrame] = []
    columns = [
        "control_id",
        "lib_id",
        "observed_ion_bins_json",
        "observed_loss_bins_json",
        "n_observed_ion_bins",
        "n_ms2_pairs",
        "ms2_pairs_json",
    ]
    for lib in LIBS:
        path = lookup_dir / f"lib{lib}_{ION[lib]}_observed_ms2_lookup.csv"
        if not path.exists():
            continue
        frame = pd.read_csv(path, dtype=str, keep_default_na=False, usecols=lambda column: column in columns)
        if "control_id" not in frame.columns:
            continue
        frame = frame[frame["control_id"].isin(ids)]
        if frame.empty:
            continue
        frame["lib_id"] = str(lib)
        frame["ord"] = LIBS.index(lib)
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["control_id", "n_ion_bins", "peak_count", "ion_json", "loss_json", "lib_id"])
    observed = pd.concat(frames, ignore_index=True).sort_values(["ord", "control_id"], kind="mergesort")
    observed = observed.drop_duplicates("control_id", keep="first")
    observed["n_ion_bins"] = pd.to_numeric(observed.get("n_observed_ion_bins"), errors="coerce").fillna(0)
    peaks = pd.to_numeric(observed.get("n_ms2_pairs"), errors="coerce")
    if "ms2_pairs_json" in observed.columns:
        counted = observed["ms2_pairs_json"].map(_pair_count)
        peaks = peaks.where(peaks.notna() & peaks.gt(0), counted)
    observed["peak_count"] = peaks.fillna(0).astype(int)
    return observed.rename(
        columns={"observed_ion_bins_json": "ion_json", "observed_loss_bins_json": "loss_json"}
    )[["control_id", "n_ion_bins", "peak_count", "ion_json", "loss_json", "lib_id"]]


def _pair_count(text: str) -> int:
    try:
        parsed = json.loads(text or "[]")
    except json.JSONDecodeError:
        return 0
    return len(parsed) if isinstance(parsed, list) else 0


def _bin_count(text: str, reported: float) -> int:
    if reported and reported == reported and reported > 0:
        return int(reported)
    try:
        parsed = json.loads(text or "[]")
    except json.JSONDecodeError:
        return 0
    return len(parsed) if isinstance(parsed, list) else 0


def build_cohort(
    ids: list[str],
    *,
    no_input_ids: set[str],
    observed: pd.DataFrame,
    stored: pd.DataFrame,
    fallback: pd.DataFrame,
) -> pd.DataFrame:
    observed = observed.set_index("control_id") if len(observed) else pd.DataFrame()
    stored = stored.set_index("control_id") if len(stored) else pd.DataFrame()
    fallback = fallback.set_index("control_id") if len(fallback) else pd.DataFrame()
    rows: list[dict[str, object]] = []
    for compound_id in ids:
        obs = observed.loc[compound_id] if compound_id in getattr(observed, "index", []) else None
        ion_json = "" if obs is None else str(obs.get("ion_json", "") or "")
        loss_json = "" if obs is None else str(obs.get("loss_json", "") or "")
        reported_bins = 0.0 if obs is None else float(obs.get("n_ion_bins", 0) or 0)
        n_bins = _bin_count(ion_json, reported_bins)
        peak_count = 0 if obs is None else int(obs.get("peak_count", 0) or 0)
        stored_ppm = _finite(float(stored.loc[compound_id, "ppm"])) if compound_id in getattr(stored, "index", []) else None
        fallback_ppm = (
            _finite(float(fallback.loc[compound_id, "ppm"])) if compound_id in getattr(fallback, "index", []) else None
        )
        source = seed_pool(stored_ppm, fallback_ppm)
        seed_frame = stored if source == "stored" else fallback if source == "adduct_fallback" else None
        seed = seed_frame.loc[compound_id] if seed_frame is not None and compound_id in seed_frame.index else None
        has_candidate = stored_ppm is not None or fallback_ppm is not None
        rows.append(
            {
                "x_compound_id": compound_id,
                "has_ion_bins": bool(n_bins > 0),
                "peak_count": peak_count,
                "best_stored_ppm": "" if stored_ppm is None else round(stored_ppm, 4),
                "best_fallback_ppm": "" if fallback_ppm is None else round(fallback_ppm, 4),
                "seed_smiles": "" if seed is None else str(seed.get("smiles", "") or ""),
                "seed_formula": "" if seed is None else str(seed.get("formula", "") or ""),
                "seed_lib": "" if seed is None else str(seed.get("lib_id", "") or ""),
                "seed_source": source,
                "seed_ppm": "" if seed is None else round(float(seed["ppm"]), 4),
                "ion_json": ion_json,
                "loss_json": loss_json,
                "cohort": assign_cohort(
                    no_inputs=compound_id in no_input_ids,
                    has_ion_bins=n_bins > 0,
                    on_mass=source != "",
                    has_candidate=has_candidate,
                ),
            }
        )
    return pd.DataFrame(rows)


def _relocate_single_bonds(mol: object, accept) -> None:
    """Move a heteroatom attachment onto a different carbon of the same formula."""

    from rdkit import Chem

    carbons = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() > 0
    ]
    moves: list[tuple[int, int]] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in {7, 8, 15, 16}:
            continue
        carbon_neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue
        moves.append((atom.GetIdx(), carbon_neighbors[0]))
    for hetero, origin in moves:
        for dest in carbons:
            if dest == origin:
                continue
            edited = Chem.RWMol(mol)
            edited.RemoveBond(hetero, origin)
            edited.AddBond(hetero, dest, Chem.BondType.SINGLE)
            accept(edited.GetMol())


def connectivity_key(mol: object) -> str:
    from rdkit import Chem

    return Chem.MolToInchiKey(mol)[:14]


def formula_of(mol: object) -> str:
    from rdkit.Chem import rdMolDescriptors

    return rdMolDescriptors.CalcMolFormula(mol)


def enumerate_constitutional_isomers(smiles: str, limit: int = 50) -> list[str]:
    """Return up to ``limit`` same-formula constitutions other than ``smiles``."""

    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    seed = Chem.MolFromSmiles(smiles)
    if seed is None:
        return []
    seed_formula = formula_of(seed)
    seen = {connectivity_key(seed)}
    found: list[str] = []
    queue = [seed]

    def accept(mol: object) -> None:
        if len(found) >= limit:
            return
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            return
        if formula_of(mol) != seed_formula:
            return
        key = connectivity_key(mol)
        if key in seen:
            return
        seen.add(key)
        text = Chem.MolToSmiles(mol)
        if text:
            found.append(text)
            queue.append(mol)

    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetMaxTautomers(limit + 1)
    try:
        for tautomer in enumerator.Enumerate(seed):
            accept(tautomer)
    except Exception:
        pass
    _relocate_single_bonds(seed, accept)
    reactions = []
    for smarts in FORMULA_PRESERVING_SMIRKS:
        try:
            reactions.append(Chem.rdChemReactions.ReactionFromSmarts(smarts))
        except Exception:
            continue
    while queue and len(found) < limit:
        current = queue.pop(0)
        _relocate_single_bonds(current, accept)
        for reaction in reactions:
            if len(found) >= limit:
                break
            try:
                products = reaction.RunReactants((current,))
            except Exception:
                continue
            for product in products:
                if len(found) >= limit:
                    break
                accept(product[0])
    return found[:limit]


def select_pilot(cohort: pd.DataFrame, n: int = 20, min_peaks: int = 3) -> pd.DataFrame:
    work = cohort[cohort["cohort"].eq(ON_MASS)].copy()
    work["peak_count"] = pd.to_numeric(work["peak_count"], errors="coerce").fillna(0)
    work = work[work["peak_count"].ge(min_peaks) & work["seed_smiles"].astype(str).str.len().gt(0)]
    return work.sort_values("x_compound_id", kind="mergesort").head(n).reset_index(drop=True)


def score_isomers(row: pd.Series, isomers: list[str]) -> dict[str, object]:
    from structural_elucidation.nce.verifier_v2 import predicted_spectrum_score_v2

    lib_id = int(float(row["seed_lib"])) if str(row.get("seed_lib", "")).strip() else 209
    reference_dir = os.environ.get("NCE_REFERENCE_DIR", "").strip() or None
    best_score = 0.0
    best_smiles = ""
    for smiles in isomers:
        score = float(
            predicted_spectrum_score_v2(
                smiles,
                row.get("ion_json", "[]"),
                row.get("loss_json", "[]"),
                lib_id=lib_id,
                reference_dir=reference_dir,
            )
        )
        if score > best_score:
            best_score = score
            best_smiles = smiles
    return {
        "x_compound_id": row["x_compound_id"],
        "seed_smiles": row["seed_smiles"],
        "seed_formula": row["seed_formula"],
        "seed_ppm": row["seed_ppm"],
        "seed_lib": row["seed_lib"],
        "n_isomers": len(isomers),
        "best_smiles": best_smiles,
        "best_score": best_score,
        "above_zero": best_score > 0,
    }


def _no_input_ids(fallback_ids: Path, selection: Path) -> set[str]:
    fallback = set(fallback_ids.read_text(encoding="utf-8").split())
    payload = json.loads(selection.read_text(encoding="utf-8"))
    searched = {str(item["control_id"]) for item in payload["assignments"]}
    return fallback - searched


def _unpublished_ids(by_compound: Path) -> list[str]:
    frame = pd.read_csv(by_compound, dtype=str, usecols=["x_compound_id", "prediction_source"])
    ids = frame.loc[frame["prediction_source"].eq("no_prediction"), "x_compound_id"].astype(str)
    return sorted(set(ids))


def write_cohort(args: argparse.Namespace) -> pd.DataFrame:
    ids = _unpublished_ids(args.by_compound)
    no_inputs = _no_input_ids(args.fallback_ids, args.fallback_selection)
    cohort = build_cohort(
        ids,
        no_input_ids=no_inputs,
        observed=load_observed(args.lookup_dir, set(ids)),
        stored=best_rows(args.stored_root, set(ids)),
        fallback=best_rows(args.fallback_root, set(ids)),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    cohort.to_csv(args.out, index=False)
    counts = cohort["cohort"].value_counts().to_dict()
    print("[cohort]", json.dumps({"compounds": len(cohort), **counts}))
    return cohort


def write_pilot(args: argparse.Namespace) -> pd.DataFrame:
    cohort = pd.read_csv(args.cohort, dtype=str, keep_default_na=False)
    chosen = select_pilot(cohort, n=args.n, min_peaks=args.min_peaks)
    if len(chosen) == 0:
        raise SystemExit("no on-mass seeds with enough observed peaks")
    rows = []
    for _, row in chosen.iterrows():
        isomers = enumerate_constitutional_isomers(str(row["seed_smiles"]), limit=args.isomer_limit)
        scored = score_isomers(row, isomers)
        rows.append(scored)
        print(
            f"[pilot] {scored['x_compound_id']} isomers={scored['n_isomers']} best={scored['best_score']}",
            flush=True,
        )
    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    above = int(out["above_zero"].astype(bool).sum()) if len(out) else 0
    print("[pilot] TOTAL", json.dumps({"scored": len(out), "above_zero": above}))
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    cohort = sub.add_parser("cohort")
    cohort.add_argument("--by-compound", type=Path, required=True)
    cohort.add_argument("--lookup-dir", type=Path, required=True)
    cohort.add_argument("--stored-root", type=Path, required=True)
    cohort.add_argument("--fallback-root", type=Path, required=True)
    cohort.add_argument("--fallback-ids", type=Path, required=True)
    cohort.add_argument("--fallback-selection", type=Path, required=True)
    cohort.add_argument("--out", type=Path, required=True)
    cohort.set_defaults(func=write_cohort)

    pilot = sub.add_parser("pilot")
    pilot.add_argument("--cohort", type=Path, required=True)
    pilot.add_argument("--out", type=Path, required=True)
    pilot.add_argument("--n", type=int, default=20)
    pilot.add_argument("--min-peaks", type=int, default=3)
    pilot.add_argument("--isomer-limit", type=int, default=50)
    pilot.set_defaults(func=write_pilot)

    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
