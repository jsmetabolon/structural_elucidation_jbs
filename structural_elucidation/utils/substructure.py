#!/usr/bin/env python3
"""Substructure and similarity utilities for expanded catalog retrieval."""

from __future__ import annotations

import argparse
import base64
import json
import math
import pickle
from pathlib import Path
from typing import Any, Dict, Iterable, Sequence

import pandas as pd

try:
    from rdkit import Chem
    from rdkit import DataStructs
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdFMCS

    HAS_RDKIT = True
except Exception:  # pragma: no cover - environment dependent
    HAS_RDKIT = False


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".pkl", ".pickle"}:
        df = pd.read_pickle(path)
    else:
        try:
            df = pd.read_csv(path, sep=None, engine="python")
        except Exception:
            df = pd.read_csv(path)
    df.columns = [str(c).strip().strip('"').lstrip("\ufeff") for c in df.columns]
    return df


def pick_col(df: pd.DataFrame, options: Sequence[str]) -> str | None:
    lookup = {str(col).lower(): str(col) for col in df.columns}
    for option in options:
        col = lookup.get(str(option).lower())
        if col is not None:
            return col
    return None


def decode_fp_binary(value: object) -> bytes | None:
    if value is None:
        return None
    if isinstance(value, bytes):
        return value
    if isinstance(value, memoryview):
        return bytes(value)
    text = str(value).strip()
    if not text:
        return None
    try:
        return base64.b64decode(text, validate=True)
    except Exception:
        return None


def canonicalize_smiles(smiles: str) -> str | None:
    if not HAS_RDKIT:
        text = str(smiles or "").strip()
        return text if text else None
    text = str(smiles or "").strip()
    if not text:
        return None
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return None


def morgan_fp_from_smiles(smiles: str):
    if not HAS_RDKIT:
        return None
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def pattern_fp_from_smiles(smiles: str):
    if not HAS_RDKIT:
        return None
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return Chem.PatternFingerprint(mol)


def smiles_matches_substructure(smiles: str, query_smarts: str) -> bool:
    if not HAS_RDKIT:
        return False
    mol = Chem.MolFromSmiles(str(smiles or "").strip())
    query = Chem.MolFromSmarts(str(query_smarts or "").strip())
    if mol is None or query is None:
        return False
    return bool(mol.HasSubstructMatch(query))


def tanimoto_similarity_smiles(smiles_a: str, smiles_b: str) -> float:
    if not HAS_RDKIT:
        return 0.0
    fp_a = morgan_fp_from_smiles(smiles_a)
    fp_b = morgan_fp_from_smiles(smiles_b)
    if fp_a is None or fp_b is None:
        return 0.0
    return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))


def build_consensus_smarts(
    smiles_values: Sequence[str],
    *,
    min_mcs_atoms: int = 6,
    timeout_sec: int = 3,
) -> str | None:
    if not HAS_RDKIT:
        return None
    mols = []
    for smiles in smiles_values:
        mol = Chem.MolFromSmiles(str(smiles or "").strip())
        if mol is not None:
            mols.append(mol)
    if len(mols) < 2:
        return None

    try:
        mcs = rdFMCS.FindMCS(
            mols,
            timeout=max(int(timeout_sec), 1),
            completeRingsOnly=True,
            ringMatchesRingOnly=True,
        )
    except Exception:
        return None

    smarts = str(getattr(mcs, "smartsString", "") or "").strip()
    num_atoms = int(getattr(mcs, "numAtoms", 0) or 0)
    if not smarts or num_atoms < int(min_mcs_atoms):
        return None
    return smarts


def build_feature_consensus_smarts_from_retrieval(
    retrieval_df: pd.DataFrame,
    *,
    top_n: int = 5,
    min_mcs_atoms: int = 6,
    timeout_sec: int = 3,
) -> tuple[dict[str, str], dict[str, str]]:
    if retrieval_df.empty:
        return {}, {}
    required = {"feature_key", "smiles"}
    if not required.issubset(set(retrieval_df.columns)):
        return {}, {}

    work = retrieval_df.copy()
    work["feature_key"] = work["feature_key"].astype(str).str.strip()
    work["smiles"] = work["smiles"].astype(str).str.strip()
    work = work[work["smiles"].ne("") & work["smiles"].str.lower().ne("nan")]
    if "status" in work.columns:
        work = work[work["status"].astype(str).eq("ok")]
    if "rank" in work.columns:
        work["rank"] = pd.to_numeric(work["rank"], errors="coerce").fillna(9999)
        work = work.sort_values(["feature_key", "rank"], kind="mergesort")
    else:
        work = work.sort_values(["feature_key"], kind="mergesort")

    consensus_by_feature: dict[str, str] = {}
    top_smiles_by_feature: dict[str, str] = {}
    for feature_key, group in work.groupby("feature_key", sort=True):
        smiles_values = []
        for smi in group["smiles"].tolist():
            can = canonicalize_smiles(smi)
            if can is None:
                continue
            if can in smiles_values:
                continue
            smiles_values.append(can)
            if len(smiles_values) >= max(int(top_n), 1):
                break
        if not smiles_values:
            continue
        top_smiles_by_feature[str(feature_key)] = str(smiles_values[0])
        smarts = build_consensus_smarts(
            smiles_values,
            min_mcs_atoms=int(min_mcs_atoms),
            timeout_sec=int(timeout_sec),
        )
        if smarts:
            consensus_by_feature[str(feature_key)] = smarts
    return consensus_by_feature, top_smiles_by_feature


def _mean(values: Sequence[float]) -> float:
    finite = [float(v) for v in values if np_is_finite(v)]
    if not finite:
        return 0.0
    return float(sum(finite) / len(finite))


def np_is_finite(value: object) -> bool:
    try:
        return bool(math.isfinite(float(value)))
    except Exception:
        return False


def _has_shared_substructure(
    mol_a,
    mol_b,
    *,
    min_mcs_atoms: int,
    timeout_sec: int,
) -> bool:
    if mol_a is None or mol_b is None:
        return False
    try:
        mcs = rdFMCS.FindMCS(
            [mol_a, mol_b],
            timeout=max(int(timeout_sec), 1),
            completeRingsOnly=True,
            ringMatchesRingOnly=True,
        )
    except Exception:
        return False
    return int(getattr(mcs, "numAtoms", 0) or 0) >= int(min_mcs_atoms)


def compute_structural_neighbor_scores(
    smiles_values: Sequence[str],
    *,
    top_n: int = 5,
    min_mcs_atoms: int = 6,
    mcs_timeout_sec: int = 1,
) -> list[dict[str, float]]:
    """Compute per-candidate structural cohesion scores.

    Returns one dict per input row with:
    - substructure_consensus_score
    - tanimoto_neighbor_score
    """
    n = len(smiles_values)
    if n == 0 or not HAS_RDKIT:
        return [
            {
                "substructure_consensus_score": 0.0,
                "tanimoto_neighbor_score": 0.0,
            }
            for _ in range(n)
        ]

    mols = []
    morgan_fps = []
    for smiles in smiles_values:
        mol = Chem.MolFromSmiles(str(smiles or "").strip())
        mols.append(mol)
        if mol is None:
            morgan_fps.append(None)
        else:
            morgan_fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048))

    top_indices = [idx for idx, mol in enumerate(mols) if mol is not None][: max(int(top_n), 1)]
    if not top_indices:
        return [
            {
                "substructure_consensus_score": 0.0,
                "tanimoto_neighbor_score": 0.0,
            }
            for _ in range(n)
        ]

    top_mols = [mols[idx] for idx in top_indices]
    consensus_smarts = build_consensus_smarts(
        [smiles_values[idx] for idx in top_indices],
        min_mcs_atoms=int(min_mcs_atoms),
        timeout_sec=max(int(mcs_timeout_sec), 1),
    )
    consensus_mol = Chem.MolFromSmarts(consensus_smarts) if consensus_smarts else None
    consensus_top_match_count = (
        sum(1 for mol in top_mols if mol is not None and consensus_mol is not None and mol.HasSubstructMatch(consensus_mol))
        if consensus_mol is not None
        else 0
    )
    pair_cache: dict[tuple[int, int], bool] = {}

    out: list[dict[str, float]] = []
    for idx in range(n):
        mol_i = mols[idx]
        fp_i = morgan_fps[idx]
        if mol_i is None or fp_i is None:
            out.append(
                {
                    "substructure_consensus_score": 0.0,
                    "tanimoto_neighbor_score": 0.0,
                }
            )
            continue

        tanis = []
        for j in top_indices:
            if j == idx:
                continue
            fp_j = morgan_fps[j]
            if fp_j is None:
                continue
            tanis.append(float(DataStructs.TanimotoSimilarity(fp_i, fp_j)))
        tanimoto_neighbor_score = _mean(tanis)

        if consensus_mol is not None:
            if mol_i.HasSubstructMatch(consensus_mol):
                substructure_consensus_score = float(consensus_top_match_count / max(len(top_indices), 1))
            else:
                substructure_consensus_score = 0.0
        else:
            shared_flags: list[float] = []
            for j in top_indices:
                if j == idx:
                    continue
                key = (min(idx, j), max(idx, j))
                shared = pair_cache.get(key)
                if shared is None:
                    shared = _has_shared_substructure(
                        mol_i,
                        mols[j],
                        min_mcs_atoms=int(min_mcs_atoms),
                        timeout_sec=max(int(mcs_timeout_sec), 1),
                    )
                    pair_cache[key] = shared
                shared_flags.append(1.0 if shared else 0.0)
            substructure_consensus_score = _mean(shared_flags)

        out.append(
            {
                "substructure_consensus_score": float(substructure_consensus_score),
                "tanimoto_neighbor_score": float(tanimoto_neighbor_score),
            }
        )
    return out


def load_fingerprint_index(path: Path) -> dict[int, dict[str, Any]]:
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    if isinstance(payload, dict):
        out = {}
        for key, value in payload.items():
            try:
                out[int(key)] = dict(value)
            except Exception:
                continue
        return out
    if isinstance(payload, list):
        out = {}
        for item in payload:
            if not isinstance(item, dict):
                continue
            try:
                cid = int(item["comp_id"])
            except Exception:
                continue
            out[cid] = dict(item)
        return out
    raise TypeError(f"Unsupported fingerprint index payload: {type(payload)}")


def run_substructure_query(
    *,
    query: str,
    catalog_csv: Path,
    fingerprint_index_path: Path,
    max_results: int,
) -> pd.DataFrame:
    if not HAS_RDKIT:
        raise RuntimeError("RDKit is required for substructure search.")

    query = str(query or "").strip()
    query_substruct_mol = Chem.MolFromSmarts(query)
    query_smiles_mol = Chem.MolFromSmiles(query)
    query_mol = query_substruct_mol or query_smiles_mol
    if query_mol is None:
        raise ValueError(f"Invalid SMARTS/SMILES query: {query!r}")

    query_pattern_fp = Chem.PatternFingerprint(query_mol)
    query_morgan_fp = (
        AllChem.GetMorganFingerprintAsBitVect(query_smiles_mol, radius=2, nBits=2048)
        if query_smiles_mol is not None
        else None
    )

    catalog_df = read_table(catalog_csv)
    id_col = pick_col(catalog_df, ["comp_id", "id"])
    smiles_col = pick_col(catalog_df, ["canonical_smiles", "smiles", "SMILES"])
    if id_col is None or smiles_col is None:
        raise KeyError(f"Catalog missing comp_id/smiles columns: {catalog_df.columns.tolist()}")
    catalog_df["comp_id"] = pd.to_numeric(catalog_df[id_col], errors="coerce").astype("Int64")
    catalog_df["smiles"] = catalog_df[smiles_col].astype(str).str.strip()
    catalog_df = catalog_df.dropna(subset=["comp_id"]).copy()
    catalog_lookup = {
        int(row.comp_id): str(row.smiles)
        for row in catalog_df[["comp_id", "smiles"]].itertuples(index=False)
    }

    fp_index = load_fingerprint_index(fingerprint_index_path)
    rows: list[dict[str, Any]] = []
    for comp_id, record in fp_index.items():
        smiles = str(record.get("smiles") or catalog_lookup.get(int(comp_id), "")).strip()
        if not smiles:
            continue
        pattern_bin = decode_fp_binary(record.get("pattern_fp_binary"))
        pattern_fp = (
            DataStructs.CreateFromBinaryText(pattern_bin)
            if pattern_bin is not None
            else pattern_fp_from_smiles(smiles)
        )
        if pattern_fp is None:
            continue
        if not DataStructs.AllProbeBitsMatch(query_pattern_fp, pattern_fp):
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or not mol.HasSubstructMatch(query_mol):
            continue
        morgan_bin = decode_fp_binary(record.get("morgan_fp_binary"))
        morgan_fp = (
            DataStructs.CreateFromBinaryText(morgan_bin)
            if morgan_bin is not None
            else morgan_fp_from_smiles(smiles)
        )
        tanimoto = 0.0
        if query_morgan_fp is not None and morgan_fp is not None:
            tanimoto = float(DataStructs.TanimotoSimilarity(query_morgan_fp, morgan_fp))
        rows.append(
            {
                "comp_id": int(comp_id),
                "smiles": smiles,
                "inchikey14": str(record.get("inchikey14") or ""),
                "tanimoto_to_query": float(tanimoto),
            }
        )

    out = pd.DataFrame(rows).sort_values(
        ["tanimoto_to_query", "comp_id"],
        ascending=[False, True],
        kind="mergesort",
    )
    if int(max_results) > 0:
        out = out.head(int(max_results)).reset_index(drop=True)
    return out


def run_tanimoto_query(
    *,
    query_smiles: str,
    catalog_csv: Path,
    fingerprint_index_path: Path,
    threshold: float,
    max_results: int,
) -> pd.DataFrame:
    if not HAS_RDKIT:
        raise RuntimeError("RDKit is required for Tanimoto search.")
    query_can = canonicalize_smiles(query_smiles)
    if query_can is None:
        raise ValueError(f"Invalid query SMILES: {query_smiles!r}")
    query_fp = morgan_fp_from_smiles(query_can)
    if query_fp is None:
        raise ValueError(f"Failed to build query fingerprint for: {query_can}")

    catalog_df = read_table(catalog_csv)
    id_col = pick_col(catalog_df, ["comp_id", "id"])
    smiles_col = pick_col(catalog_df, ["canonical_smiles", "smiles", "SMILES"])
    if id_col is None or smiles_col is None:
        raise KeyError(f"Catalog missing comp_id/smiles columns: {catalog_df.columns.tolist()}")
    catalog_df["comp_id"] = pd.to_numeric(catalog_df[id_col], errors="coerce").astype("Int64")
    catalog_df["smiles"] = catalog_df[smiles_col].astype(str).str.strip()
    catalog_df = catalog_df.dropna(subset=["comp_id"]).copy()
    catalog_lookup = {
        int(row.comp_id): str(row.smiles)
        for row in catalog_df[["comp_id", "smiles"]].itertuples(index=False)
    }

    fp_index = load_fingerprint_index(fingerprint_index_path)
    rows: list[dict[str, Any]] = []
    for comp_id, record in fp_index.items():
        smiles = str(record.get("smiles") or catalog_lookup.get(int(comp_id), "")).strip()
        if not smiles:
            continue
        morgan_bin = decode_fp_binary(record.get("morgan_fp_binary"))
        fp = (
            DataStructs.CreateFromBinaryText(morgan_bin)
            if morgan_bin is not None
            else morgan_fp_from_smiles(smiles)
        )
        if fp is None:
            continue
        tanimoto = float(DataStructs.TanimotoSimilarity(query_fp, fp))
        if tanimoto < float(threshold):
            continue
        rows.append(
            {
                "comp_id": int(comp_id),
                "smiles": smiles,
                "inchikey14": str(record.get("inchikey14") or ""),
                "tanimoto_to_query": tanimoto,
            }
        )
    out = pd.DataFrame(rows).sort_values(
        ["tanimoto_to_query", "comp_id"],
        ascending=[False, True],
        kind="mergesort",
    )
    if int(max_results) > 0:
        out = out.head(int(max_results)).reset_index(drop=True)
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run substructure or Tanimoto search against precomputed catalog fingerprints."
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=("substructure", "tanimoto"),
        required=True,
    )
    parser.add_argument("--query", type=str, required=True)
    parser.add_argument("--catalog_csv", type=Path, required=True)
    parser.add_argument("--fingerprint_index", type=Path, required=True)
    parser.add_argument("--tanimoto_threshold", type=float, default=0.7)
    parser.add_argument("--max_results", type=int, default=250)
    parser.add_argument("--output_csv", type=Path, default=None)
    parser.add_argument("--output_json", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mode = str(args.mode).strip().lower()
    if mode == "substructure":
        result_df = run_substructure_query(
            query=str(args.query),
            catalog_csv=args.catalog_csv.expanduser().resolve(),
            fingerprint_index_path=args.fingerprint_index.expanduser().resolve(),
            max_results=int(args.max_results),
        )
    else:
        result_df = run_tanimoto_query(
            query_smiles=str(args.query),
            catalog_csv=args.catalog_csv.expanduser().resolve(),
            fingerprint_index_path=args.fingerprint_index.expanduser().resolve(),
            threshold=float(args.tanimoto_threshold),
            max_results=int(args.max_results),
        )

    if args.output_csv is not None:
        output_csv = args.output_csv.expanduser().resolve()
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(output_csv, index=False)
        print(f"Wrote: {output_csv}")
    if args.output_json is not None:
        output_json = args.output_json.expanduser().resolve()
        output_json.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "mode": mode,
            "query": str(args.query),
            "rows": int(result_df.shape[0]),
            "results": result_df.to_dict(orient="records"),
        }
        output_json.write_text(json.dumps(payload, indent=2))
        print(f"Wrote: {output_json}")

    print(f"Result rows: {int(result_df.shape[0])}")
    if not result_df.empty:
        print(result_df.head(min(10, len(result_df))).to_string(index=False))


if __name__ == "__main__":
    main()
