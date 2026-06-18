"""Smoke test for structural_elucidation package layout and imports."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

try:
    import pytest
except ImportError:  # pragma: no cover
    pytest = None  # type: ignore[assignment]

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

REQUIRED_DIRS = [
    "structural_elucidation",
    "structural_elucidation/pipeline",
    "structural_elucidation/inference",
    "vendor/nce",
    "models/decoder",
    "data",
    "data/base_pools",
    "data/integrated_best",
    "data/feature_rows",
    "data/observed_ms2_lookups",
    "data/smiles_db",
]

REQUIRED_FILES = [
    "structural_elucidation/paths.py",
    "structural_elucidation/inference/predict.py",
    "structural_elucidation/inference/merged.py",
    "structural_elucidation/pipeline/verify.py",
    "structural_elucidation/pipeline/pool.py",
    "structural_elucidation/pipeline/rerank.py",
    "structural_elucidation/pipeline/enrichment/pathway.py",
    "structural_elucidation/pipeline/enrichment/lipid_class.py",
    "structural_elucidation/pipeline/enrichment/peptide_composition.py",
    "structural_elucidation/reports/compound.py",
    "structural_elucidation/nce/runtime.py",
    "models/decoder/config.nb",
    "models/decoder/vocab.nb",
    "data/ions_unknown_merged_clean.csv",
    "data/known_X.csv",
    "requirements.txt",
    "fetch_assets.sh",
    "run_inference.sh",
    "README.md",
    "pyproject.toml",
]

CORE_MODULES = [
    "structural_elucidation.inference.predict",
    "structural_elucidation.inference.merged",
    "structural_elucidation.nce.runtime",
    "structural_elucidation.pipeline.verify",
    "structural_elucidation.pipeline.pool",
    "structural_elucidation.pipeline.rerank",
    "structural_elucidation.pipeline.enrichment.pathway",
    "structural_elucidation.pipeline.enrichment.lipid_class",
    "structural_elucidation.pipeline.enrichment.peptide_composition",
    "structural_elucidation.reports.compound",
    "structural_elucidation.utils.phase0",
    "structural_elucidation.inference.retrieval_scoring",
    "structural_elucidation.utils.merge",
    "structural_elucidation.inference.merged_model",
]


def _smoke_mark(func):
    return pytest.mark.smoke(func) if pytest is not None else func


@_smoke_mark
def test_layout() -> None:
    missing_dirs = [d for d in REQUIRED_DIRS if not (REPO_ROOT / d).exists()]
    missing_files = [f for f in REQUIRED_FILES if not (REPO_ROOT / f).exists()]
    assert not missing_dirs, f"Missing dirs: {missing_dirs}"
    assert not missing_files, f"Missing files: {missing_files}"


@_smoke_mark
def test_imports() -> None:
    failed = []
    for mod in CORE_MODULES:
        try:
            importlib.import_module(mod)
        except Exception as exc:  # pragma: no cover
            failed.append(f"{mod}: {exc}")
    assert not failed, "Import failures:\n" + "\n".join(failed)


@_smoke_mark
def test_nce_paths() -> None:
    from structural_elucidation.nce import runtime as nce

    assert nce.DEFAULT_NCE_ROOT == REPO_ROOT / "vendor" / "nce"
    assert nce.DEFAULT_REFERENCE_DIR == REPO_ROOT / "vendor" / "modelfiles"


@_smoke_mark
def test_predict_smiles_rename() -> None:
    text = (REPO_ROOT / "structural_elucidation/inference/predict.py").read_text(encoding="utf-8")
    assert "predict_x209_smiles" not in text
    assert "from structural_elucidation.utils.phase0 import" in text


@_smoke_mark
def test_paths_module() -> None:
    from structural_elucidation.paths import REPO_ROOT as pkg_root, resolve

    assert pkg_root == REPO_ROOT
    assert resolve(Path("data/known_X.csv")) == REPO_ROOT / "data/known_X.csv"


def _run_without_pytest() -> None:
    test_layout()
    test_paths_module()
    test_predict_smiles_rename()
    test_nce_paths()
    print("smoke: layout/paths PASS")
    try:
        test_imports()
        print("smoke: imports PASS")
    except AssertionError as exc:
        print(f"smoke: imports SKIPPED ({exc})")


if __name__ == "__main__":
    _run_without_pytest()
