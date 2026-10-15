"""``--include-all-x`` must not silently fall back to the integrated-best controls."""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from structural_elucidation.pipeline import pool


def _args(tmp_path: Path, catalog: Path) -> argparse.Namespace:
    (tmp_path / "ib.csv").write_text(
        "control_id,lib_id,ion_mode,top10_best_tanimoto,top10_smiles\nX-1,209,negative,0.9,CCO\n"
    )
    (tmp_path / "scored.csv").write_text("control_id,smiles\nX-1,CCO\n")
    (tmp_path / "base.csv").write_text("control_id,smiles\nX-1,CCO\nX-2,CCN\n")
    return argparse.Namespace(
        current_integrated_best=tmp_path / "ib.csv",
        lib_id=209,
        ion_mode="negative",
        include_all_x=True,
        all_x_catalog=catalog,
        scored_candidates=tmp_path / "scored.csv",
        extra_inputs="",
        output_root=tmp_path / "out",
    )


def test_missing_catalog_is_an_error(tmp_path: Path) -> None:
    args = _args(tmp_path, tmp_path / "missing.csv")
    with pytest.raises(SystemExit, match="include-all-x"):
        pool.run_from_base_pool(args, tmp_path, tmp_path / "out", tmp_path / "base.csv")


def test_empty_scored_file_yields_no_nce_features(tmp_path: Path) -> None:
    empty = tmp_path / "scored.csv"
    empty.write_text("\n")
    frame = pool.nce_feature_table(empty, tmp_path)
    assert list(frame.columns)[:2] == ["control_id", "canonical_2d_smiles"]
    assert frame.empty


def test_catalog_without_this_library_is_an_error(tmp_path: Path) -> None:
    catalog = tmp_path / "catalog.csv"
    catalog.write_text("lib_id,x_compound_id\n400,X-2\n")
    args = _args(tmp_path, catalog)
    with pytest.raises(SystemExit, match="lib 209"):
        pool.run_from_base_pool(args, tmp_path, tmp_path / "out", tmp_path / "base.csv")
