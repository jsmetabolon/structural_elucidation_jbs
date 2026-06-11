#!/usr/bin/env python3
"""Smoke test for structural_elucidation_jbs repo structure and imports."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"

REQUIRED_DIRS = [
    "src",
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
    "src/predict_smiles.py",
    "src/predict_merged.py",
    "src/score_candidates_with_verifier_v1.py",
    "src/build_unified_candidate_pool_v16.py",
    "src/train_unified_rerank_v1.py",
    "src/nce_local_runtime_utils.py",
    "data/ions_unknown_merged_clean.csv",
    "data/known_X.csv",
    "requirements.txt",
    "fetch_assets.sh",
    "run_inference.sh",
    "README.md",
]

CORE_MODULES = [
    "predict_smiles",
    "predict_merged",
    "nce_local_runtime_utils",
    "score_candidates_with_verifier_v1",
    "build_unified_candidate_pool_v16",
    "train_unified_rerank_v1",
    "phase0_utils",
    "retrieval_scoring",
    "merge_policy",
    "merged_model",
]


def test_layout() -> None:
    missing_dirs = [d for d in REQUIRED_DIRS if not (REPO_ROOT / d).exists()]
    missing_files = [f for f in REQUIRED_FILES if not (REPO_ROOT / f).exists()]
    assert not missing_dirs, f"Missing dirs: {missing_dirs}"
    assert not missing_files, f"Missing files: {missing_files}"


def test_imports() -> None:
    sys.path.insert(0, str(SRC_ROOT))
    failed = []
    for mod in CORE_MODULES:
        try:
            importlib.import_module(mod)
        except Exception as exc:  # pragma: no cover - diagnostic
            failed.append(f"{mod}: {exc}")
    assert not failed, "Import failures:\n" + "\n".join(failed)


def test_nce_paths() -> None:
    sys.path.insert(0, str(SRC_ROOT))
    import nce_local_runtime_utils as nce

    assert nce.DEFAULT_NCE_ROOT == REPO_ROOT / "vendor" / "nce"
    assert nce.DEFAULT_REFERENCE_DIR == REPO_ROOT / "vendor" / "modelfiles"


def test_predict_smiles_rename() -> None:
    text = (SRC_ROOT / "predict_smiles.py").read_text(encoding="utf-8")
    assert "predict_x209_smiles" not in text
    assert "from phase0_utils import" in text or "from phase0_utils" in text


def main() -> None:
    test_layout()
    test_predict_smiles_rename()
    test_nce_paths()
    test_imports()
    print("smoke_test: PASS")


if __name__ == "__main__":
    main()
