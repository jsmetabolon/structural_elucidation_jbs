"""Tests for second-stage tentative naming layers."""
from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "structural_elucidation" / "reports"))
import tentative_naming as TN  # noqa: E402
import tentative_naming_ext  # noqa: E402
from tentative_naming_ext import (  # noqa: E402
    ExtendedNamingContext,
    apply_extended_layers,
    load_iupac_cache,
    resolve_layer6_external_db,
    resolve_layer7_analog,
    resolve_layer8_class,
    resolve_layer9_iupac,
    truncate_for_title,
)


def _make_external_db(path: Path) -> None:
    with sqlite3.connect(str(path)) as connection:
        connection.execute(
            "CREATE TABLE names (inchikey TEXT, match_kind TEXT, canonical_name TEXT, "
            "source TEXT, priority INTEGER, PRIMARY KEY (inchikey, match_kind, canonical_name))"
        )
        connection.execute(
            "INSERT INTO names VALUES (?, ?, ?, ?, ?)",
            ("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", "full", "aspirin", "hmdb", 3),
        )
        connection.execute(
            "INSERT INTO names VALUES (?, ?, ?, ?, ?)",
            ("BSYNRYMUTXBXSQ", "skeleton", "aspirin", "hmdb", 3),
        )


def _make_classyfire_db(path: Path) -> None:
    with sqlite3.connect(str(path)) as connection:
        connection.execute(
            "CREATE TABLE classyfire (inchikey TEXT PRIMARY KEY, inchikey_skeleton TEXT, "
            "chemont_id TEXT, superclass TEXT, class_name TEXT, subclass TEXT)"
        )
        connection.execute(
            "INSERT INTO classyfire VALUES (?, ?, ?, ?, ?, ?)",
            (
                "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                "BSYNRYMUTXBXSQ",
                "CHEMONTID:0000001",
                "Benzenoids",
                "Benzene and derivatives",
                "Salicylic acids and derivatives",
            ),
        )


def test_layer6_exact_match(tmp_path: Path) -> None:
    db_path = tmp_path / "external.sqlite"
    _make_external_db(db_path)
    with sqlite3.connect(str(db_path)) as connection:
        result = resolve_layer6_external_db(
            "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            "C9H8O4",
            connection,
        )
    assert result is not None
    assert result.kind == "external_db_exact"
    assert result.label == "Aspirin (C9H8O4)"
    assert result.confidence == "high"


def test_layer6_skeleton_match(tmp_path: Path) -> None:
    db_path = tmp_path / "external.sqlite"
    _make_external_db(db_path)
    with sqlite3.connect(str(db_path)) as connection:
        result = resolve_layer6_external_db("BSYNRYMUTXBXSQ-OTHER-KEY", "C9H8O4", connection)
    assert result is not None
    assert result.kind == "external_db_skeleton"
    assert "stereo unspecified" in result.label
    assert not result.label.endswith("— tentative")


def test_layer7_high_tanimoto(tmp_path: Path) -> None:
    csv_path = tmp_path / "neighbors.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["canonical_name", "inchikeys", "similarity_to_target"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "canonical_name": "caffeic acid",
                "inchikeys": "QAOWNCQODCNURD-UHFFFAOYSA-N",
                "similarity_to_target": "0.91",
            }
        )
    result = resolve_layer7_analog(csv_path)
    assert result is not None
    assert result.kind == "analog"
    assert "≈ caffeic acid" in result.label
    assert result.analog_neighbor_tanimoto == pytest.approx(0.91)


def test_layer7_low_tanimoto(tmp_path: Path) -> None:
    csv_path = tmp_path / "neighbors.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["canonical_name", "inchikeys", "similarity_to_target"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "canonical_name": "tyrosol",
                "inchikeys": "KEY",
                "similarity_to_target": "0.55",
            }
        )
    result = resolve_layer7_analog(csv_path)
    assert result is not None
    assert "Putative analog of tyrosol" in result.label
    assert result.confidence == "low"


def test_layer8_classyfire(tmp_path: Path) -> None:
    db_path = tmp_path / "classyfire.sqlite"
    _make_classyfire_db(db_path)
    with sqlite3.connect(str(db_path)) as connection:
        result = resolve_layer8_class(
            "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            "C9H8O4",
            connection,
        )
    assert result is not None
    assert result.kind == "chemical_class"
    assert "Salicylic acids and derivatives" in result.label


def test_layer9_iupac_from_cache() -> None:
    cache = {
        "CC(=O)Oc1ccccc1C(=O)O": {
            "iupac_name": "2-acetyloxybenzoic acid",
            "round_trip_ok": True,
        }
    }
    result = resolve_layer9_iupac("CC(=O)Oc1ccccc1C(=O)O", "C9H8O4", cache)
    assert result is not None
    assert result.kind == "iupac_generated"
    assert "2-acetyloxybenzoic acid" in result.label
    assert result.iupac_round_trip_ok is True


def test_layer9_rejects_failed_round_trip() -> None:
    cache = {
        "CCO": {
            "iupac_name": "ethanol",
            "round_trip_ok": False,
        }
    }
    assert resolve_layer9_iupac("CCO", "C2H6O", cache) is None


def test_priority_exact_beats_analog(tmp_path: Path) -> None:
    external_db = tmp_path / "external.sqlite"
    classyfire_db = tmp_path / "classyfire.sqlite"
    _make_external_db(external_db)
    _make_classyfire_db(classyfire_db)
    bundle_dir = tmp_path / "X-TEST"
    data_dir = bundle_dir / "data"
    data_dir.mkdir(parents=True)
    neighbors = data_dir / "similarity_neighbors.csv"
    with neighbors.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["canonical_name", "inchikeys", "similarity_to_target"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "canonical_name": "neighbor",
                "inchikeys": "KEY",
                "similarity_to_target": "0.99",
            }
        )
    ctx = ExtendedNamingContext(
        external_db_path=external_db,
        classyfire_db_path=classyfire_db,
    )
    result = apply_extended_layers(
        "CC(=O)Oc1ccccc1C(=O)O",
        "C9H8O4",
        "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        bundle_dir,
        ctx,
    )
    assert result is not None
    assert result.kind == "external_db_exact"


def test_resolve_uses_extended_layers_when_primary_misses(tmp_path: Path) -> None:
    external_db = tmp_path / "external.sqlite"
    _make_external_db(external_db)
    bundle_dir = tmp_path / "X-TEST"
    data_dir = bundle_dir / "data"
    data_dir.mkdir(parents=True)
    ctx = ExtendedNamingContext(external_db_path=external_db)
    result = TN.resolve(
        "CC(=O)Oc1ccccc1C(=O)O",
        "C9H8O4",
        {},
        {},
        name_resolved=False,
        bundle_dir=bundle_dir,
        extended_ctx=ctx,
        inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
    )
    assert result is not None
    assert result.kind == "external_db_exact"


def test_truncate_for_title() -> None:
    long_name = "a" * 120
    truncated = truncate_for_title(long_name, max_len=90)
    assert len(truncated) == 90
    assert truncated.endswith("…")


def test_load_iupac_cache_round_trip(tmp_path: Path) -> None:
    cache_path = tmp_path / "cache.jsonl"
    cache_path.write_text(
        json.dumps({"smiles": "CCO", "iupac_name": "ethanol", "round_trip_ok": True}) + "\n",
        encoding="utf-8",
    )
    loaded = load_iupac_cache(cache_path)
    assert loaded["CCO"]["iupac_name"] == "ethanol"


@pytest.mark.parametrize(
    "module_path",
    [
        "/repo/structural_elucidation/reports/tentative_naming_ext.py",
        "/a/b/c/structural_elucidation/reports/tentative_naming_ext.py",
    ],
)
def test_default_context_does_not_index_past_the_filesystem_root(monkeypatch, module_path):
    """Every resource this discovers is optional, so finding none is not an error.

    It used to index two levels above the project root unconditionally. For a
    checkout at depth one - which is where the repository sits in CI and in the
    cloud VM - that raised IndexError and took the whole report build with it.
    """

    monkeypatch.setattr(tentative_naming_ext, "__file__", module_path)
    assert tentative_naming_ext.default_extended_context() is None
