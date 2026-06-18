#!/usr/bin/env python3
"""One-shot migration: flat src/ -> structural_elucidation package layout."""

from __future__ import annotations

import re
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"
PKG = REPO / "structural_elucidation"

MODULE_MAP: dict[str, str] = {
    "score_candidates_with_verifier_v1": "pipeline/verify.py",
    "build_unified_candidate_pool_v1": "pipeline/pool_v1.py",
    "build_unified_candidate_pool_v16": "pipeline/pool.py",
    "train_unified_rerank_v1": "pipeline/rerank.py",
    "build_library_anchored_candidates_v1": "pipeline/candidates.py",
    "predict_smiles": "inference/predict.py",
    "predict_merged": "inference/merged.py",
    "model1D2Conv": "inference/model_conv.py",
    "merged_model": "inference/merged_model.py",
    "latent_calibration": "inference/calibration.py",
    "retrieval_scoring": "inference/retrieval_scoring.py",
    "ri_audit_utils": "inference/ri_audit.py",
    "build_compound_report": "reports/compound.py",
    "nce_local_runtime_utils": "nce/runtime.py",
    "predicted_spectrum_verifier_v1": "nce/verifier_v1.py",
    "predicted_spectrum_verifier_v2": "nce/verifier_v2.py",
    "nce_evidence_utils": "nce/evidence.py",
    "evaluate_denovo_failure_modes": "evaluation/failure_modes.py",
    "evaluate_denovo_controls": "evaluation/controls.py",
    "evaluate_selfies_candidate_pools": "evaluation/candidate_pools.py",
    "benchmark_selfies_merge_policies": "evaluation/merge_policies.py",
    "run_retrieval_anchored_denovo_assembly_v1": "evaluation/assembly_v1.py",
    "run_retrieval_anchored_denovo_assembly_v4": "evaluation/assembly_v4.py",
    "run_retrieval_anchored_analog_v1": "evaluation/analog.py",
    "aggregate_x_domain_analog_v1": "evaluation/aggregate.py",
    "assess_denovo_pipeline_v3_outputs": "evaluation/assess.py",
    "build_post_run_decision_package": "evaluation/decision.py",
    "rerank_denovo_assembly_v3": "evaluation/rerank_assembly.py",
    "chemistry_class_diagnostics_v1": "chemistry/diagnostics.py",
    "build_multi_conjugate_candidates_v1": "chemistry/conjugates.py",
    "formula_prior_v1": "chemistry/formula.py",
    "denovo_improvement_utils": "utils/denovo.py",
    "x_domain_analog_utils": "utils/domain.py",
    "spectrum_preprocessing": "utils/spectrum.py",
    "substructure_search": "utils/substructure.py",
    "html_signal_utils": "utils/html.py",
    "wgcna_inputs_v1": "utils/wgcna.py",
    "device_utils": "utils/device.py",
    "merge_policy": "utils/merge.py",
    "phase0_utils": "utils/phase0.py",
    "utils": "utils/rdkit_helpers.py",
}

OLD_TO_NEW: dict[str, str] = {
    name: f"structural_elucidation.{path.replace('/', '.').removesuffix('.py')}"
    for name, path in MODULE_MAP.items()
}

RESOLVE_IMPORT = "from structural_elucidation.paths import resolve"

SUBPACKAGES = [
    "pipeline",
    "inference",
    "reports",
    "nce",
    "evaluation",
    "chemistry",
    "utils",
]


def _pkg_path(old_name: str) -> str:
    rel = MODULE_MAP[old_name]
    return f"structural_elucidation.{rel.replace('/', '.').removesuffix('.py')}"


def rewrite_imports(text: str) -> str:
    text = re.sub(r"^#\s*flake8:\s*noqa\s*\n", "", text, flags=re.MULTILINE)
    text = re.sub(
        r"from run_retrieval_anchored_denovo_assembly_v4 import resolve",
        RESOLVE_IMPORT,
        text,
    )

    for old, new in sorted(OLD_TO_NEW.items(), key=lambda x: -len(x[0])):
        text = re.sub(rf"\bimport {re.escape(old)} as ", f"import {new} as ", text)
        text = re.sub(rf"\bfrom {re.escape(old)} import ", f"from {new} import ", text)
        text = re.sub(rf"\bimport {re.escape(old)}\b", f"import {new}", text)

    text = re.sub(
        r"import wgcna_inputs_v1 as wgcna",
        "from structural_elucidation.utils import wgcna",
        text,
    )
    text = re.sub(
        r"import merge_policy as merge_policy",
        "from structural_elucidation.utils import merge as merge_policy",
        text,
    )
    return text


def create_paths() -> None:
    (PKG / "paths.py").write_text(
        '''"""Repository path resolution helpers."""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def resolve(path: Path, repo_root: Path | None = None) -> Path:
    """Resolve a path relative to the repository root."""
    root = REPO_ROOT if repo_root is None else repo_root
    path = path.expanduser()
    return path if path.is_absolute() else root / path
''',
        encoding="utf-8",
    )


def create_utils_io() -> None:
    (PKG / "utils" / "io.py").write_text(
        '''"""Shared table and file IO helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception:
        df = pd.read_csv(path)
    df.columns = [str(col).strip().strip(\'"\').lstrip("\\ufeff") for col in df.columns]
    return df


def pick_col(df: pd.DataFrame, options: list[str]) -> str | None:
    lookup = {str(col).strip().lower(): str(col) for col in df.columns}
    for option in options:
        found = lookup.get(str(option).strip().lower())
        if found is not None:
            return found
    return None


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    cols = [str(col) for col in df.columns]
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in df.iterrows():
        values = []
        for col in df.columns:
            value = row.get(col, "")
            if isinstance(value, float):
                value = round(value, 4) if np.isfinite(value) else ""
            values.append(str(value).replace("|", "\\\\|"))
        lines.append("| " + " | ".join(values) + " |")
    return "\\n".join(lines)


def write_csv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\\n")
''',
        encoding="utf-8",
    )


def create_utils_chemistry() -> None:
    (PKG / "utils" / "chemistry.py").write_text(
        '''"""Shared chemistry helpers."""

from __future__ import annotations

import re
from typing import Any

import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors

    HAS_RDKIT = True
except ImportError:  # pragma: no cover
    Chem = None
    Descriptors = None
    rdMolDescriptors = None
    HAS_RDKIT = False


FORMULA_ELEMENTS = ["C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"]
FORMULA_FEATURE_COLUMNS = [f"formula_count_{element}" for element in FORMULA_ELEMENTS]
FORMULA_FEATURE_SCALE = 64.0


def canonicalize_smiles(smiles: Any) -> str | None:
    if not HAS_RDKIT:
        return None
    text = str(smiles or "").strip()
    if not text or text.lower() == "nan":
        return None
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return None


def mol_from_smiles(smiles: Any):
    if not HAS_RDKIT:
        return None
    canonical = canonicalize_smiles(smiles)
    if canonical is None:
        return None
    return Chem.MolFromSmiles(canonical)


def inchikey_from_smiles(smiles: Any) -> str:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return ""
    try:
        return str(Chem.MolToInchiKey(mol))
    except Exception:
        return ""


def formula_from_smiles(smiles: Any) -> str:
    mol = mol_from_smiles(smiles)
    if mol is None or rdMolDescriptors is None:
        return ""
    try:
        return str(rdMolDescriptors.CalcMolFormula(mol))
    except Exception:
        return ""


def exact_mass_from_smiles(smiles: Any) -> float:
    mol = mol_from_smiles(smiles)
    if mol is None or Descriptors is None:
        return float("nan")
    try:
        return float(Descriptors.ExactMolWt(mol))
    except Exception:
        return float("nan")


_FORMULA_RE = re.compile(r"([A-Z][a-z]?)(\\d*)")


def formula_counts(formula: Any) -> dict[str, int]:
    text = str(formula or "").strip()
    if not text or text.lower() == "nan":
        return {}
    counts: dict[str, int] = {}
    for element, raw_count in _FORMULA_RE.findall(text):
        if element not in FORMULA_ELEMENTS:
            continue
        count = int(raw_count) if raw_count else 1
        counts[element] = counts.get(element, 0) + count
    return counts


def formula_vector_from_counts(counts: dict[str, int]) -> list[float]:
    return [
        float(np.clip(float(counts.get(element, 0)) / FORMULA_FEATURE_SCALE, 0.0, 1.0))
        for element in FORMULA_ELEMENTS
    ]


def formula_vector_from_formula(formula: Any) -> list[float]:
    return formula_vector_from_counts(formula_counts(formula))


def formula_vector_from_smiles(smiles: Any) -> list[float]:
    return formula_vector_from_formula(formula_from_smiles(smiles))


def zero_formula_vector() -> list[float]:
    return [0.0 for _ in FORMULA_ELEMENTS]
''',
        encoding="utf-8",
    )


def patch_denovo(text: str) -> str:
    """Remove duplicated io/chemistry helpers; import from utils submodules."""
    remove_blocks = [
        r"def read_table\(path: Path\).*?(?=\n\ndef |\nclass |\Z)",
        r"def pick_col\(df: pd\.DataFrame.*?return None\n\n",
        r"def canonicalize_smiles\(smiles: Any\).*?return None\n\n",
        r"def mol_from_smiles\(smiles: Any\):.*?return Chem\.MolFromSmiles\(canonical\)\n\n",
        r"def inchikey_from_smiles\(smiles: Any\).*?return \"\"\n\n",
        r"def formula_from_smiles\(smiles: Any\).*?return \"\"\n\n",
        r"def exact_mass_from_smiles\(smiles: Any\).*?return float\(\"nan\"\)\n\n",
        r"FORMULA_ELEMENTS = .*?\nFORMULA_FEATURE_SCALE = .*?\n\n",
        r"_FORMULA_RE = re\.compile.*?\n\n",
        r"def formula_counts\(formula: Any\).*?return counts\n\n",
        r"def formula_vector_from_counts\(counts: dict\[str, int\]\).*?return \[.*?\]\n\n",
        r"def formula_vector_from_formula\(formula: Any\).*?return formula_vector_from_counts\(formula_counts\(formula\)\)\n\n",
        r"def formula_vector_from_smiles\(smiles: Any\).*?return formula_vector_from_formula\(formula_from_smiles\(smiles\)\)\n\n",
        r"def zero_formula_vector\(\).*?return \[0\.0 for _ in FORMULA_ELEMENTS\]\n\n",
    ]
    for pattern in remove_blocks:
        text = re.sub(pattern, "", text, flags=re.DOTALL)

    insert = '''from structural_elucidation.utils.chemistry import (
    FORMULA_ELEMENTS,
    FORMULA_FEATURE_COLUMNS,
    FORMULA_FEATURE_SCALE,
    canonicalize_smiles,
    exact_mass_from_smiles,
    formula_counts,
    formula_from_smiles,
    formula_vector_from_counts,
    formula_vector_from_formula,
    formula_vector_from_smiles,
    inchikey_from_smiles,
    mol_from_smiles,
    zero_formula_vector,
)
from structural_elucidation.utils.io import pick_col, read_table

'''
    text = re.sub(
        r"(import pandas as pd\n)",
        r"\1\n" + insert,
        text,
        count=1,
    )
    return text


def patch_domain(text: str) -> str:
    text = re.sub(
        r"from structural_elucidation\.utils\.denovo import \(\n.*?\)\n",
        "from structural_elucidation.utils.chemistry import canonicalize_smiles, inchikey_from_smiles\n"
        "from structural_elucidation.utils.io import read_table\n",
        text,
        flags=re.DOTALL,
    )
    remove_blocks = [
        r"def dataframe_to_markdown\(df: pd\.DataFrame\).*?return \"\\n\"\.join\(lines\)\n\n",
        r"def write_csv\(path: Path, df: pd\.DataFrame\).*?df\.to_csv\(path, index=False\)\n\n",
        r"def write_json\(path: Path, payload: dict\[str, Any\]\).*?path\.write_text\(.*?\)\n\n",
    ]
    for pattern in remove_blocks:
        text = re.sub(pattern, "", text, flags=re.DOTALL)

    insert = "from structural_elucidation.utils.io import dataframe_to_markdown, write_csv, write_json\n\n"
    text = re.sub(r"(import pandas as pd\n)", r"\1\n" + insert, text, count=1)
    return text


def patch_assembly_v4(text: str) -> str:
    text = re.sub(
        r"def resolve\(path: Path, repo_root: Path\) -> Path:.*?return path if path\.is_absolute\(\) else repo_root / path\n\n",
        "",
        text,
        flags=re.DOTALL,
    )
    if RESOLVE_IMPORT not in text:
        text = re.sub(
            r"(import pandas as pd\n)",
            r"\1\n" + RESOLVE_IMPORT + "\n",
            text,
            count=1,
        )
    return text


def patch_repo_root(text: str, depth: int = 2) -> str:
    """Replace Path(__file__).resolve().parents[N] with REPO_ROOT where N matches package depth."""
    text = re.sub(
        rf"Path\(__file__\)\.resolve\(\)\.parents\[{depth}\]",
        "REPO_ROOT",
        text,
    )
    if "REPO_ROOT" in text and "from structural_elucidation.paths import REPO_ROOT" not in text:
        if "from structural_elucidation.paths import resolve" in text:
            text = text.replace(
                "from structural_elucidation.paths import resolve",
                "from structural_elucidation.paths import REPO_ROOT, resolve",
            )
        elif "import pandas as pd" in text:
            text = re.sub(
                r"(import pandas as pd\n)",
                r"\1\nfrom structural_elucidation.paths import REPO_ROOT\n",
                text,
                count=1,
            )
        else:
            text = "from structural_elucidation.paths import REPO_ROOT\n\n" + text
    return text


def migrate_file(old_name: str, rel_path: str) -> None:
    src_path = SRC / f"{old_name}.py"
    dst_path = PKG / rel_path
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    text = src_path.read_text(encoding="utf-8")
    text = rewrite_imports(text)

    if old_name == "denovo_improvement_utils":
        text = patch_denovo(text)
    elif old_name == "x_domain_analog_utils":
        text = patch_domain(text)
    elif old_name == "run_retrieval_anchored_denovo_assembly_v4":
        text = patch_assembly_v4(text)

    depth = len(Path(rel_path).parts)
    if depth >= 2:
        text = patch_repo_root(text, depth=2)

    dst_path.write_text(text, encoding="utf-8")


def create_skeleton() -> None:
    if PKG.exists():
        shutil.rmtree(PKG)
    PKG.mkdir()
    (PKG / "__init__.py").write_text(
        '''"""Structural elucidation inference and reranking pipeline."""

__name__ = "structural_elucidation"
__description__ = "Mass-spectrometry structural elucidation inference and reranking"
''',
        encoding="utf-8",
    )
    for sub in SUBPACKAGES:
        (PKG / sub).mkdir()
        (PKG / sub / "__init__.py").write_text("", encoding="utf-8")
    create_paths()
    create_utils_io()
    create_utils_chemistry()


def main() -> None:
    create_skeleton()
    for old_name, rel_path in MODULE_MAP.items():
        migrate_file(old_name, rel_path)
    print(f"Migrated {len(MODULE_MAP)} modules to {PKG}")


if __name__ == "__main__":
    main()
