#!/usr/bin/env python3
"""Copy and rename Python modules from parent repo into structural_elucidation_jbs/src/."""

from __future__ import annotations

import re
import shutil
from pathlib import Path

PARENT = Path(__file__).resolve().parents[2]
DEST = Path(__file__).resolve().parents[1] / "src"

# source path (relative to PARENT) -> dest filename in src/
# Only modules on the inference / report critical path are packaged here.
COPY_MAP: dict[str, str] = {
    "predict_x209_smiles.py": "predict_smiles.py",
    "merged_encoder/predict_merged.py": "predict_merged.py",
    "merged_encoder/merged_model.py": "merged_model.py",
    "x209_retrieval_scoring.py": "retrieval_scoring.py",
    "x209_merge_policy.py": "merge_policy.py",
    "x209_phase0_utils.py": "phase0_utils.py",
    "x400_ri_audit_utils.py": "ri_audit_utils.py",
    "device_utils.py": "device_utils.py",
    "model1D2Conv.py": "model1D2Conv.py",
    "latent_calibration.py": "latent_calibration.py",
    "html_signal_utils.py": "html_signal_utils.py",
    "spectrum_preprocessing.py": "spectrum_preprocessing.py",
    "substructure_search.py": "substructure_search.py",
    "denovo_improvement_utils.py": "denovo_improvement_utils.py",
    "x_domain_analog_utils.py": "x_domain_analog_utils.py",
    "nce_local_runtime_utils.py": "nce_local_runtime_utils.py",
    "predicted_spectrum_verifier_v2.py": "predicted_spectrum_verifier_v2.py",
    "predicted_spectrum_verifier_v1.py": "predicted_spectrum_verifier_v1.py",
    "score_candidates_with_verifier_v1.py": "score_candidates_with_verifier_v1.py",
    "build_unified_candidate_pool_v1.py": "build_unified_candidate_pool_v1.py",
    "build_unified_candidate_pool_v16.py": "build_unified_candidate_pool_v16.py",
    "train_unified_rerank_v1.py": "train_unified_rerank_v1.py",
    "build_library_anchored_candidates_v1.py": "build_library_anchored_candidates_v1.py",
    "chemistry_class_diagnostics_v1.py": "chemistry_class_diagnostics_v1.py",
    "evaluate_denovo_failure_modes.py": "evaluate_denovo_failure_modes.py",
    "rerank_denovo_assembly_v3.py": "rerank_denovo_assembly_v3.py",
    "run_retrieval_anchored_denovo_assembly_v4.py": "run_retrieval_anchored_denovo_assembly_v4.py",
    "benchmark_selfies_merge_policies.py": "benchmark_selfies_merge_policies.py",
    "nce_evidence_utils.py": "nce_evidence_utils.py",
    "formula_prior_v1.py": "formula_prior_v1.py",
    "build_post_run_decision_package.py": "build_post_run_decision_package.py",
    "evaluate_denovo_controls.py": "evaluate_denovo_controls.py",
    "evaluate_selfies_candidate_pools.py": "evaluate_selfies_candidate_pools.py",
    "aggregate_x_domain_analog_v1.py": "aggregate_x_domain_analog_v1.py",
    "assess_denovo_pipeline_v3_outputs.py": "assess_denovo_pipeline_v3_outputs.py",
    "build_multi_conjugate_candidates_v1.py": "build_multi_conjugate_candidates_v1.py",
    "run_retrieval_anchored_denovo_assembly_v1.py": "run_retrieval_anchored_denovo_assembly_v1.py",
    "run_retrieval_anchored_analog_v1.py": "run_retrieval_anchored_analog_v1.py",
    "utils.py": "utils.py",
}

IMPORT_REPLACEMENTS = [
    (r"\bpredict_x209_smiles\b", "predict_smiles"),
    (r"\bx209_retrieval_scoring\b", "retrieval_scoring"),
    (r"\bx209_merge_policy\b", "merge_policy"),
    (r"\bx209_phase0_utils\b", "phase0_utils"),
    (r"\bx400_ri_audit_utils\b", "ri_audit_utils"),
    (r"from merged_encoder\.merged_model import", "from merged_model import"),
    (r"from merged_encoder import merged_model", "import merged_model"),
    (r"from merged_model import", "from merged_model import"),
    (r"except ImportError:\s*\n\s*from merged_model import", "except ImportError:\n    from merged_model import"),
    (r"REPO_ROOT = Path\(__file__\)\.resolve\(\)\.parents\[1\]", "REPO_ROOT = Path(__file__).resolve().parents[1]"),
]


def rewrite_content(text: str, dest_name: str) -> str:
    for pattern, repl in IMPORT_REPLACEMENTS:
        text = re.sub(pattern, repl, text)
    if dest_name == "predict_merged.py":
        text = text.replace("import predict_x209_smiles as shared_predict", "import predict_smiles as shared_predict")
        text = text.replace('"""Merged-ionization inference wrapper over predict_x209_smiles.py."""',
                            '"""Merged-ionization inference wrapper over predict_smiles.py."""')
        text = re.sub(
            r"REPO_ROOT = Path\(__file__\)\.resolve\(\)\.parents\[1\]",
            "REPO_ROOT = Path(__file__).resolve().parents[1]",
            text,
        )
    if dest_name == "predict_smiles.py":
        text = text.replace('"""Batch-predict SMILES with scores for X_209 features."""',
                            '"""Batch-predict SMILES with scores for merged-library features."""')
        text = re.sub(r"from x209_phase0_utils import", "from phase0_utils import", text)
        text = re.sub(r"from x209_retrieval_scoring import", "from retrieval_scoring.py import", text)
        text = re.sub(r"from x209_merge_policy import", "from merge_policy import", text)
        text = re.sub(r"from x400_ri_audit_utils import", "from ri_audit_utils import", text)
        text = text.replace("from retrieval_scoring.py import", "from retrieval_scoring import")
    if dest_name == "nce_local_runtime_utils.py":
        text = re.sub(
            r'DEFAULT_NCE_ROOT = Path\("~/dev/NewCompoundElucidation"\)\.expanduser\(\)\n'
            r'DEFAULT_REFERENCE_DIR = Path\("~/dev/modelfiles"\)\.expanduser\(\)\n'
            r'DEFAULT_MLIMS_TABLE_CANDIDATES = \(\n'
            r'    Path\("~/dev/modelfiles/mlims_tables\.pkl"\)\.expanduser\(\),\n'
            r'    Path\("/data/josh/spec2mol/shared/references/production/mlims_tables\.pkl"\),\n'
            r'\)',
            'REPO_ROOT = Path(__file__).resolve().parents[1]\n'
            'DEFAULT_NCE_ROOT = REPO_ROOT / "vendor" / "nce"\n'
            'DEFAULT_REFERENCE_DIR = REPO_ROOT / "vendor" / "modelfiles"\n'
            'DEFAULT_MLIMS_TABLE_CANDIDATES = (\n'
            '    REPO_ROOT / "vendor" / "modelfiles" / "mlims_tables.pkl",\n'
            '    REPO_ROOT / "vendor" / "modelfiles" / "mlims_tables.pkl",\n'
            ')',
            text,
            count=1,
        )
    return text


def main() -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    copied = []
    missing = []
    for src_rel, dest_name in COPY_MAP.items():
        src = PARENT / src_rel
        if not src.exists():
            missing.append(src_rel)
            continue
        dest = DEST / dest_name
        content = src.read_text(encoding="utf-8")
        content = rewrite_content(content, dest_name)
        dest.write_text(content, encoding="utf-8")
        copied.append((src_rel, dest_name))
    print(f"Copied {len(copied)} files to {DEST}")
    if missing:
        print("Missing sources (will need manual add):")
        for m in missing:
            print(f"  - {m}")


if __name__ == "__main__":
    main()
