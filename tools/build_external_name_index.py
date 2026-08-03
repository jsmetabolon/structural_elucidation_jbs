#!/usr/bin/env python3
"""Build a merged InChIKey -> canonical-name SQLite index for second-stage naming."""
from __future__ import annotations

import argparse
import csv
import sqlite3
import sys
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
DEFAULT_OUTPUT = PROJECT_ROOT / "reference_data" / "external_name_index.sqlite"

SOURCE_PRIORITY = {
    "metabolon": 1,
    "external_names": 2,
    "hmdb": 3,
    "lipidmaps": 4,
    "master": 5,
    "catalog_hmdb": 6,
}


def _clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "nan", "none", "null"} else text


def _inchikeys_from_smiles(smiles: str) -> tuple[str, str]:
    if not HAS_RDKIT:
        return "", ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "", ""
    inchikey = Chem.MolToInchiKey(mol)
    if not inchikey:
        return "", ""
    return inchikey, inchikey.split("-")[0]


def _register(
    rows: dict[tuple[str, str, str], int],
    *,
    name: str,
    inchikey: str,
    source: str,
) -> None:
    name = _clean(name)
    inchikey = _clean(inchikey)
    if not name or not inchikey:
        return
    priority = SOURCE_PRIORITY[source]
    full = inchikey if len(inchikey) >= 27 else ""
    skeleton = inchikey.split("-")[0] if "-" in inchikey else inchikey[:14]
    for key, kind in ((full, "full"), (skeleton, "skeleton")):
        if not key:
            continue
        slot = (key, kind, name)
        rows[slot] = min(rows.get(slot, priority), priority)


def _load_metabolon(rows: dict[tuple[str, str, str], int], db_path: Path) -> int:
    if not db_path.is_file():
        return 0
    count = 0
    with sqlite3.connect(str(db_path)) as connection:
        query = connection.execute(
            "SELECT inchikeys, canonical_name FROM metabolites "
            "WHERE canonical_name IS NOT NULL AND trim(canonical_name) <> ''"
        )
        for inchikeys, canonical_name in query:
            for token in str(inchikeys or "").replace(";", ",").split(","):
                token = token.strip()
                if token:
                    _register(rows, name=canonical_name, inchikey=token, source="metabolon")
                    count += 1
    return count


def _load_external_names(rows: dict[tuple[str, str, str], int], db_path: Path) -> int:
    if not db_path.is_file():
        return 0
    count = 0
    with sqlite3.connect(str(db_path)) as connection:
        tables = {
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            )
        }
        if "metabolites" not in tables:
            return 0
        query = connection.execute(
            "SELECT inchikeys, canonical_name FROM metabolites "
            "WHERE canonical_name IS NOT NULL AND trim(canonical_name) <> ''"
        )
        for inchikeys, canonical_name in query:
            for token in str(inchikeys or "").replace(";", ",").split(","):
                token = token.strip()
                if token:
                    _register(rows, name=canonical_name, inchikey=token, source="external_names")
                    count += 1
    return count


def _load_hmdb(rows: dict[tuple[str, str, str], int], csv_path: Path, smiles_to_inchikey: dict[str, str]) -> int:
    if not csv_path.is_file():
        return 0
    count = 0
    with csv_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
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
                    inchikey, _ = _inchikeys_from_smiles(smiles)
            if name and inchikey:
                _register(rows, name=name, inchikey=inchikey, source="hmdb")
                count += 1
    return count


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


def _load_lipidmaps(rows: dict[tuple[str, str, str], int], csv_path: Path) -> int:
    if not csv_path.is_file():
        return 0
    count = 0
    with csv_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            name = _clean(row.get("display_name") or row.get("normalized_name"))
            inchikey = _clean(row.get("inchi_key") or row.get("inchikey"))
            smiles = _clean(row.get("canonical_smiles") or row.get("smiles"))
            if not inchikey and smiles:
                inchikey, _ = _inchikeys_from_smiles(smiles)
            if name and inchikey:
                _register(rows, name=name, inchikey=inchikey, source="lipidmaps")
                count += 1
    return count


def _load_master(rows: dict[tuple[str, str, str], int], csv_path: Path) -> int:
    if not csv_path.is_file():
        return 0
    count = 0
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            name = _clean(row.get("SHORT_NAME") or row.get("CHEMICAL_NAME"))
            inchikeys = _clean(row.get("INCHIKEY") or row.get("inchikey"))
            smiles = _clean(row.get("SMILES"))
            tokens = [token.strip() for token in inchikeys.replace(";", ",").split(",") if token.strip()]
            if not tokens and smiles:
                full, skeleton = _inchikeys_from_smiles(smiles)
                tokens = [token for token in (full, skeleton) if token]
            for token in tokens:
                if name:
                    _register(rows, name=name, inchikey=token, source="master")
                    count += 1
    return count


def _source_for_priority(priority: int) -> str:
    for source, pri in SOURCE_PRIORITY.items():
        if pri == priority:
            return source
    return "unknown"


def build_index(output_path: Path) -> dict[str, int]:
    rows: dict[tuple[str, str, str], int] = {}
    smiles_db = WORKSPACE_ROOT / "spec2mol_merged_fullrun" / "SMILES_DB"
    metabolite_db = WORKSPACE_ROOT / "spec2mol_merged_fullrun_metabolon_only" / "Metabolite Synonym Database"
    master_csv = WORKSPACE_ROOT / "spec2mol_merged_fullrun_metabolon_only" / "MASTER_METABOLITE_LIST.csv"
    smiles_to_inchikey = _load_smiles_to_inchikey(
        metabolite_db=metabolite_db / "metabolite_synonyms.sqlite",
        lipidmaps_csv=smiles_db / "lipidmaps_name_index.csv",
        master_csv=master_csv,
    )
    stats = {
        "metabolon": _load_metabolon(rows, metabolite_db / "metabolite_synonyms.sqlite"),
        "external_names": _load_external_names(
            rows, metabolite_db / "metabolite_synonyms_merged_with_external_names.sqlite"
        ),
        "hmdb": _load_hmdb(rows, smiles_db / "hmdb_minimal.csv", smiles_to_inchikey),
        "lipidmaps": _load_lipidmaps(rows, smiles_db / "lipidmaps_name_index.csv"),
        "master": _load_master(rows, master_csv),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        output_path.unlink()
    with sqlite3.connect(str(output_path)) as connection:
        connection.execute(
            """
            CREATE TABLE names (
                inchikey TEXT NOT NULL,
                match_kind TEXT NOT NULL,
                canonical_name TEXT NOT NULL,
                source TEXT NOT NULL,
                priority INTEGER NOT NULL,
                PRIMARY KEY (inchikey, match_kind, canonical_name)
            )
            """
        )
        connection.executemany(
            "INSERT INTO names (inchikey, match_kind, canonical_name, source, priority) "
            "VALUES (?, ?, ?, ?, ?)",
            [
                (key, kind, name, _source_for_priority(priority), priority)
                for (key, kind, name), priority in rows.items()
            ],
        )
        connection.execute("CREATE INDEX idx_names_lookup ON names(inchikey, match_kind, priority)")
        stats["rows"] = connection.execute("SELECT COUNT(*) FROM names").fetchone()[0]
        stats["full_keys"] = connection.execute(
            "SELECT COUNT(DISTINCT inchikey) FROM names WHERE match_kind='full'"
        ).fetchone()[0]
        stats["skeleton_keys"] = connection.execute(
            "SELECT COUNT(DISTINCT inchikey) FROM names WHERE match_kind='skeleton'"
        ).fetchone()[0]
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    stats = build_index(args.output.resolve())
    print(f"[build_external_name_index] wrote {args.output}")
    for key, value in stats.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
