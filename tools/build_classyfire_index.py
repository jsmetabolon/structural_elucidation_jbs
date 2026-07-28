#!/usr/bin/env python3
"""Build InChIKey -> ClassyFire/ChemOnt class SQLite index from HMDB annotations."""
from __future__ import annotations

import argparse
import csv
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

try:
    from rdkit import Chem

    HAS_RDKIT = True
except ImportError:  # pragma: no cover
    HAS_RDKIT = False

WORKSPACE_ROOT = PROJECT_ROOT.parents[1]
DEFAULT_OUTPUT = PROJECT_ROOT / "reference_data" / "classyfire_hmdb.sqlite"
DEFAULT_ANNOTATIONS = (
    WORKSPACE_ROOT
    / "NewCompoundElucidation"
    / "elucidation"
    / "table_builders"
    / "HMDB_36_classyfire_21_annotations.csv"
)
DEFAULT_HMDB = WORKSPACE_ROOT / "spec2mol_merged_fullrun" / "SMILES_DB" / "hmdb_minimal.csv"
DEFAULT_CHEMONT_DEFS = WORKSPACE_ROOT / "references" / "production" / "chemont_defs.csv"
DEFAULT_METABOLITE_DB = (
    WORKSPACE_ROOT
    / "spec2mol_merged_fullrun_metabolon_only"
    / "Metabolite Synonym Database"
    / "metabolite_synonyms.sqlite"
)
DEFAULT_LIPIDMAPS = WORKSPACE_ROOT / "spec2mol_merged_fullrun" / "SMILES_DB" / "lipidmaps_name_index.csv"
DEFAULT_MASTER = WORKSPACE_ROOT / "spec2mol_merged_fullrun_metabolon_only" / "MASTER_METABOLITE_LIST.csv"

_GENERIC_CLASS_NAMES = {
    "organic compounds",
    "organoheterocyclic compounds",
    "organic oxygen compounds",
    "organic nitrogen compounds",
    "organosulfur compounds",
    "organohalogen compounds",
    "organometallic compounds",
    "hydrocarbons",
    "lipids and lipid-like molecules",
}


def _clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "nan", "none", "null"} else text


def _inchikey_from_smiles(smiles: str) -> str:
    if not HAS_RDKIT:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    return Chem.MolToInchiKey(mol)


def _normalize_hmdb_id(value: str) -> str:
    text = _clean(value).upper()
    if text.startswith("HMDB") and text[4:].isdigit():
        return f"HMDB{text[4:].zfill(7)}"
    return text


def _load_hmdb_accession_map(hmdb_csv: Path, smiles_to_inchikey: dict[str, str]) -> dict[str, tuple[str, str]]:
    mapping: dict[str, tuple[str, str]] = {}
    with hmdb_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            accession = _normalize_hmdb_id(row.get("accession", ""))
            name = _clean(row.get("name"))
            smiles = _clean(row.get("canonical_smiles") or row.get("smiles"))
            inchikey = ""
            for key in ("inchikey", "inchi_key", "inchikeys"):
                inchikey = _clean(row.get(key))
                if inchikey:
                    break
            if not inchikey and smiles:
                inchikey = smiles_to_inchikey.get(smiles, "")
                if not inchikey and HAS_RDKIT:
                    inchikey = _inchikey_from_smiles(smiles)
            if accession and inchikey:
                mapping[accession] = (inchikey, name)
    return mapping


def _load_metabolite_hmdb_map(db_path: Path) -> dict[str, tuple[str, str]]:
    if not db_path.is_file():
        return {}
    out: dict[str, tuple[str, str]] = {}
    with sqlite3.connect(str(db_path)) as connection:
        rows = connection.execute(
            "SELECT hmdb_ids, inchikeys, canonical_name FROM metabolites "
            "WHERE hmdb_ids IS NOT NULL AND trim(hmdb_ids) <> ''"
        )
        for hmdb_ids, inchikeys, canonical_name in rows:
            inchikey = _clean(str(inchikeys or "").replace(";", ",").split(",")[0])
            if not inchikey:
                continue
            for token in str(hmdb_ids or "").replace(";", ",").split(","):
                accession = _normalize_hmdb_id(token)
                if accession:
                    out[accession] = (inchikey, _clean(canonical_name))
    return out


def _load_smiles_to_inchikey(
    *,
    metabolite_db: Path,
    lipidmaps_csv: Path,
    master_csv: Path,
) -> dict[str, str]:
    out: dict[str, str] = {}
    if metabolite_db.is_file():
        with sqlite3.connect(str(metabolite_db)) as connection:
            rows = connection.execute(
                "SELECT smiles, inchikeys FROM metabolites "
                "WHERE smiles IS NOT NULL AND inchikeys IS NOT NULL"
            )
            for smiles, inchikeys in rows:
                smiles_text = _clean(smiles)
                inchikey = _clean(str(inchikeys or "").replace(";", ",").split(",")[0])
                if smiles_text and inchikey:
                    out.setdefault(smiles_text, inchikey)
    if lipidmaps_csv.is_file():
        with lipidmaps_csv.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                smiles = _clean(row.get("canonical_smiles") or row.get("smiles"))
                inchikey = _clean(row.get("inchi_key") or row.get("inchikey"))
                if smiles and inchikey:
                    out.setdefault(smiles, inchikey)
    if master_csv.is_file():
        with master_csv.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                smiles = _clean(row.get("SMILES"))
                inchikeys = _clean(row.get("INCHIKEY") or row.get("inchikey"))
                inchikey = inchikeys.replace(";", ",").split(",")[0].strip()
                if smiles and inchikey:
                    out.setdefault(smiles, inchikey)
    return out


def _load_chemont_defs(defs_csv: Path) -> dict[str, tuple[str, str]]:
    out: dict[str, tuple[str, str]] = {}
    if not defs_csv.is_file():
        return out
    with defs_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            chemont_id = _clean(row.get("chemont_id")).zfill(7)
            name = _clean(row.get("name"))
            if chemont_id and name:
                out[f"CHEMONTID:{chemont_id}"] = (name, _clean(row.get("def")))
    return out


def _pick_best_class(classes: list[tuple[str, str, str]]) -> tuple[str, str, str] | None:
    if not classes:
        return None
    filtered = [
        item
        for item in classes
        if _clean(item[2]).lower() not in _GENERIC_CLASS_NAMES
    ]
    pool = filtered or classes
    return max(pool, key=lambda item: (len(item[2]), item[1]))


def build_index(
    *,
    annotations_csv: Path,
    hmdb_csv: Path,
    chemont_defs_csv: Path,
    output_path: Path,
    metabolite_db: Path,
    lipidmaps_csv: Path,
    master_csv: Path,
) -> dict[str, int]:
    smiles_to_inchikey = _load_smiles_to_inchikey(
        metabolite_db=metabolite_db,
        lipidmaps_csv=lipidmaps_csv,
        master_csv=master_csv,
    )
    hmdb_map = _load_metabolite_hmdb_map(metabolite_db)
    hmdb_map.update(_load_hmdb_accession_map(hmdb_csv, smiles_to_inchikey))
    chemont_defs = _load_chemont_defs(chemont_defs_csv)
    by_hmdb: dict[str, list[tuple[str, str, str]]] = defaultdict(list)

    with annotations_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            compound_id = _normalize_hmdb_id(row.get("CompoundID", ""))
            chemont_token = _clean(row.get("Smiles"))
            class_name = _clean(row.get("ChemOntID"))
            parent_name = _clean(row.get("ParentName"))
            if not compound_id or not class_name:
                continue
            superclass = parent_name or chemont_defs.get(chemont_token, ("", ""))[0]
            subclass = class_name
            by_hmdb[compound_id].append((chemont_token, superclass, subclass))

    by_inchikey: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    for compound_id, classes in by_hmdb.items():
        mapped = hmdb_map.get(compound_id)
        if not mapped:
            continue
        inchikey, _ = mapped
        by_inchikey[inchikey].extend(classes)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        output_path.unlink()
    inserted = 0
    with sqlite3.connect(str(output_path)) as connection:
        connection.execute(
            """
            CREATE TABLE classyfire (
                inchikey TEXT NOT NULL PRIMARY KEY,
                inchikey_skeleton TEXT NOT NULL,
                chemont_id TEXT,
                superclass TEXT,
                class_name TEXT,
                subclass TEXT
            )
            """
        )
        for inchikey, classes in by_inchikey.items():
            best = _pick_best_class(classes)
            if not best:
                continue
            chemont_id, superclass, subclass = best
            skeleton = inchikey.split("-")[0]
            connection.execute(
                "INSERT INTO classyfire "
                "(inchikey, inchikey_skeleton, chemont_id, superclass, class_name, subclass) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (inchikey, skeleton, chemont_id, superclass, subclass, subclass),
            )
            inserted += 1
        stats = {
            "rows": inserted,
            "hmdb_accessions": len(hmdb_map),
            "annotated_hmdb": len(by_hmdb),
        }
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    parser.add_argument("--hmdb", type=Path, default=DEFAULT_HMDB)
    parser.add_argument("--chemont-defs", type=Path, default=DEFAULT_CHEMONT_DEFS)
    parser.add_argument("--metabolite-db", type=Path, default=DEFAULT_METABOLITE_DB)
    parser.add_argument("--lipidmaps", type=Path, default=DEFAULT_LIPIDMAPS)
    parser.add_argument("--master", type=Path, default=DEFAULT_MASTER)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    stats = build_index(
        annotations_csv=args.annotations.resolve(),
        hmdb_csv=args.hmdb.resolve(),
        chemont_defs_csv=args.chemont_defs.resolve(),
        output_path=args.output.resolve(),
        metabolite_db=args.metabolite_db.resolve(),
        lipidmaps_csv=args.lipidmaps.resolve(),
        master_csv=args.master.resolve(),
    )
    print(f"[build_classyfire_index] wrote {args.output}")
    for key, value in stats.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
