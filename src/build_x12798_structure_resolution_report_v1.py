#!/usr/bin/env python3
"""Build the X-12798 structure-resolution static HTML report."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import re
import sqlite3
import sys
from pathlib import Path
from typing import Any

import pandas as pd

try:
    from rdkit import Chem
    from rdkit import rdBase
    from rdkit.Chem import Descriptors, Draw, rdDepictor, rdFingerprintGenerator, rdFMCS, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.DataStructs import TanimotoSimilarity

    HAS_RDKIT = True
    rdBase.DisableLog("rdApp.error")
    MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
except ImportError:  # pragma: no cover
    Chem = None
    Descriptors = None
    Draw = None
    MORGAN = None
    HAS_RDKIT = False


DEFAULT_OUTPUT_DIR = Path("x_compound_reports/X-12798_structure_resolution")
DEFAULT_RESOLUTION_MD = Path(
    "models/target_pure_denovo_ms2_embedding_v1/"
    "x12798_encoder_comparison_v1/x12798_real_latent_structure_resolution.md"
)
DEFAULT_ENCODER_COMPARISON = Path(
    "models/target_pure_denovo_ms2_embedding_v1/x12798_encoder_comparison_v1/encoder_branch_comparison.csv"
)
DEFAULT_DEEP_DIVE_CSV = Path("models/target_unknown_deep_dive_v1/x12798/candidate_evidence.csv")
DEFAULT_DEEP_DIVE_FINAL = Path("models/target_unknown_deep_dive_v1/x12798/final_call.json")
DEFAULT_LIB400_PURE_TOP10 = Path(
    "models/target_pure_denovo_ms2_embedding_v1/x12798_lib400_only_all_candidate_latent/pure_denovo_top10.csv"
)
DEFAULT_MERGED_PURE_TOP10 = Path(
    "models/target_pure_denovo_ms2_embedding_v1/x12798_merged_libcond_all_candidate_latent/pure_denovo_top10.csv"
)
DEFAULT_LIB400_FINAL = Path(
    "models/target_pure_denovo_ms2_embedding_v1/x12798_lib400_only_all_candidate_latent/final_denovo_call.json"
)
DEFAULT_MERGED_FINAL = Path(
    "models/target_pure_denovo_ms2_embedding_v1/x12798_merged_libcond_all_candidate_latent/final_denovo_call.json"
)
DEFAULT_V15_NAMES = Path("models/v15_predicted_names_v1/all_libraries_top10_with_predicted_names.csv")
DEFAULT_LIBRARY_MATCH_MD = Path("models/target_unknown_deep_dive_v1/x12798/x12798_metabolon_library_matches.md")
DEFAULT_CONSENSUS_GUIDED_ROOT = Path("models/x12798_consensus_guided_structure_v2")
DEFAULT_REAL_LATENT_CALIBRATION_ROOT = Path("models/x12798_real_latent_calibration_v1/gate_calibration")
DEFAULT_V17_SPECTRUM_LATENT_ROOT = Path("models/x12798_v17_spectrum_latent_alignment_v1")
DEFAULT_REFERENCE_DB = Path("models/v14_ranked_pathway_annotations_v1/reference/metabolite_synonyms_merged_with_mtx_giga.sqlite")
ANCHOR_CANDIDATE_NAME = "N-(4-acetylphenyl)benzamide"
ANCHOR_CANDIDATE_SMILES = "CC(=O)c1ccc(NC(=O)c2ccccc2)cc1"
PREDICTION_SUMMARY_SENTENCE = (
    "Likely a benzamide-like xenobiotic/benzoate-family compound with two aromatic rings connected by an amide linkage, "
    "plus an acetyl/enol-like substituent whose exact position and keto/enol form remain uncertain."
)
DEFAULT_METABOLON_LIBRARY_CATALOGS = [
    Path("merged_encoder/data/decoder_train_pos_400_402_smiles.csv"),
    Path("merged_encoder/data/decoder_train_neg_209_305_smiles.csv"),
]

RERUN_COMMANDS = """Required X-12798 source artifacts are missing. Recreate pure de novo branch outputs with:

.venv/bin/python build_target_pure_denovo_ms2_embedding_v1.py \\
  --scope lib_400_positive \\
  --control-id X-12798 \\
  --output-root models/target_pure_denovo_ms2_embedding_v1/x12798_merged_libcond_real_latent \\
  --observed-latent-pt models/target_pure_denovo_ms2_embedding_v1/x12798/observed_latents/46624_merged_libcond.pt \\
  --embedding-branch-label merged_libcond \\
  --no-attempt-ec2

.venv/bin/python build_target_pure_denovo_ms2_embedding_v1.py \\
  --scope lib_400_positive \\
  --control-id X-12798 \\
  --output-root models/target_pure_denovo_ms2_embedding_v1/x12798_lib400_only_real_latent \\
  --observed-latent-pt models/target_pure_denovo_ms2_embedding_v1/x12798/observed_latents/46624_lib400_only.pt \\
  --embedding-branch-label lib400_only \\
  --no-attempt-ec2
"""

CANDIDATE_COLUMNS = [
    "candidate_id",
    "canonical_smiles",
    "display_name",
    "candidate_family",
    "source_section",
    "source_branch",
    "rank",
    "score",
    "ms2_embedding_cosine",
    "ms2_embedding_rank",
    "ms2_embedding_status",
    "formula",
    "exact_mass",
    "adduct_consistent",
    "isotope_status",
    "predicted_pathway",
    "predicted_subpathway",
    "promotion_status",
    "gate_status",
    "candidate_consensus_asset",
    "candidate_consensus_mean_weighted_fraction",
    "candidate_consensus_atoms_ge_0p5_fraction",
    "provenance",
    "metabolon_library_filter_status",
    "metabolon_detected_library_ids",
    "metabolon_library_match_basis",
]

NODE_COLUMNS = [
    "id",
    "label",
    "canonical_smiles",
    "candidate_family",
    "source_branch",
    "rank",
    "score",
    "ms2_embedding_cosine",
    "predicted_subpathway",
    "gate_status",
    "x",
    "y",
    "radius",
    "color",
]

EDGE_COLUMNS = ["source", "target", "similarity", "similarity_bucket", "label"]

TREE_NODE_COLUMNS = [
    "id",
    "label",
    "node_type",
    "x",
    "y",
    "color",
    "motif_id",
    "motif_present_fraction",
    "same_position_fraction",
    "candidate_count",
    "top_candidate_id",
    "top_candidate_label",
    "top_candidate_smiles",
    "interpretation",
]

TREE_EDGE_COLUMNS = ["source", "target", "label", "edge_type", "support_fraction", "color", "stroke_style"]

POSSIBILITY_CLASSES = [
    {
        "id": "benzamide_like_acetyl_para",
        "label": "Benzamide-like acetyl para",
        "color": "#2E7D32",
        "reader_interpretation": "Candidates preserve the benzamide-like core, amide bridge, acetyl/enol feature, and same-ring para placement.",
    },
    {
        "id": "benzamide_like_enol_or_tautomer",
        "label": "Benzamide-like enol or tautomer",
        "color": "#4F8F61",
        "reader_interpretation": "Candidates retain the benzamide-like scaffold and acetyl/enol feature, but not the same para placement.",
    },
    {
        "id": "benzamide_like_substitution_variant",
        "label": "Benzamide-like substitution variant",
        "color": "#7AA95C",
        "reader_interpretation": "Candidates retain the benzamide-like scaffold and amide bridge while changing or losing the acetyl/enol substituent.",
    },
    {
        "id": "aryl_amide_analog",
        "label": "Aryl amide analog",
        "color": "#F28E2B",
        "reader_interpretation": "Candidates retain an aryl amide pattern, but not the full benzamide-like scaffold.",
    },
    {
        "id": "fused_or_cyclized_analog",
        "label": "Fused or cyclized analog",
        "color": "#D73A49",
        "reader_interpretation": "Pure de novo candidates form fused or cyclized aromatic analogs without retaining the amide bridge.",
    },
    {
        "id": "non_benzamide_outlier",
        "label": "Non-benzamide outlier",
        "color": "#7A8691",
        "reader_interpretation": "Candidates do not fit the benzamide, aryl-amide, or fused/cyclized analog branches.",
    },
]

POSSIBILITY_BY_ID = {item["id"]: item for item in POSSIBILITY_CLASSES}

FAMILY_PRIORITY = {
    "provisional_non_pure": 0,
    "encoder_branch_top": 1,
    "pure_denovo_lib400_only": 2,
    "pure_denovo_merged_libcond": 3,
    "non_pure_deep_dive": 4,
}

FAMILY_COLORS = {
    "provisional_non_pure": "#1F4E79",
    "non_pure_deep_dive": "#4E79A7",
    "pure_denovo_lib400_only": "#D73A49",
    "pure_denovo_merged_libcond": "#F28E2B",
    "encoder_branch_top": "#6F4EAA",
}

MOTIF_DEFINITIONS = [
    {
        "motif_id": "benzamide_core",
        "motif_label": "benzamide-like core",
        "reference_atom_selector": "left aryl + amide + right phenyl",
        "motif_smarts": "c1ccc(NC(=O)c2ccccc2)cc1",
        "reader_label": "Benzamide-like scaffold",
    },
    {
        "motif_id": "left_aryl_ring",
        "motif_label": "substituted left aryl ring",
        "reference_atom_selector": "left aryl ring",
        "motif_smarts": "c1ccccc1",
        "reader_label": "Aromatic ring on acetyl/enol side",
    },
    {
        "motif_id": "amide_linker",
        "motif_label": "amide linker",
        "reference_atom_selector": "N-C(=O)",
        "motif_smarts": "NC(=O)",
        "reader_label": "Amide bridge",
    },
    {
        "motif_id": "right_phenyl",
        "motif_label": "right phenyl ring",
        "reference_atom_selector": "right phenyl ring",
        "motif_smarts": "c1ccccc1",
        "reader_label": "Phenyl ring attached through amide",
    },
    {
        "motif_id": "acetyl_enol_substituent",
        "motif_label": "acetyl/enol substituent",
        "reference_atom_selector": "CC(=O) / C=C(O)",
        "motif_smarts": "[CX3](=O)[CH3]",
        "reader_label": "Acetyl/enol-like substituent",
    },
    {
        "motif_id": "para_placement",
        "motif_label": "para placement on left aryl",
        "reference_atom_selector": "relative left-aryl placement",
        "motif_smarts": "Nc1ccc([C,c])cc1",
        "reader_label": "Same-ring para placement",
    },
]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--resolution-md", type=Path, default=DEFAULT_RESOLUTION_MD)
    parser.add_argument("--encoder-comparison-csv", type=Path, default=DEFAULT_ENCODER_COMPARISON)
    parser.add_argument("--deep-dive-csv", type=Path, default=DEFAULT_DEEP_DIVE_CSV)
    parser.add_argument("--deep-dive-final-json", type=Path, default=DEFAULT_DEEP_DIVE_FINAL)
    parser.add_argument("--lib400-pure-top10", type=Path, default=DEFAULT_LIB400_PURE_TOP10)
    parser.add_argument("--merged-pure-top10", type=Path, default=DEFAULT_MERGED_PURE_TOP10)
    parser.add_argument("--lib400-final-json", type=Path, default=DEFAULT_LIB400_FINAL)
    parser.add_argument("--merged-final-json", type=Path, default=DEFAULT_MERGED_FINAL)
    parser.add_argument("--v15-names-csv", type=Path, default=DEFAULT_V15_NAMES)
    parser.add_argument("--library-match-md", type=Path, default=DEFAULT_LIBRARY_MATCH_MD)
    parser.add_argument("--consensus-guided-root", type=Path, default=DEFAULT_CONSENSUS_GUIDED_ROOT)
    parser.add_argument("--real-latent-calibration-root", type=Path, default=DEFAULT_REAL_LATENT_CALIBRATION_ROOT)
    parser.add_argument("--v17-spectrum-latent-root", type=Path, default=DEFAULT_V17_SPECTRUM_LATENT_ROOT)
    parser.add_argument("--reference-db", type=Path, default=DEFAULT_REFERENCE_DB)
    parser.add_argument(
        "--metabolon-library-catalogs",
        type=Path,
        nargs="*",
        default=DEFAULT_METABOLON_LIBRARY_CATALOGS,
    )
    return parser.parse_args(argv)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = columns or (list(rows[0]) if rows else [])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in fieldnames})


def safe_json_records(df: pd.DataFrame) -> list[dict[str, Any]]:
    return json.loads(df.to_json(orient="records"))


def is_missing(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    text = str(value)
    return text == "" or text.lower() in {"nan", "none", "null", "<na>"}


def clean(value: Any, default: str = "") -> str:
    if is_missing(value):
        return default
    return str(value)


def as_float(value: Any) -> float | None:
    if is_missing(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def as_bool_text(value: Any) -> str:
    if is_missing(value):
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    text = str(value).strip()
    if text.lower() in {"true", "1", "yes"}:
        return "true"
    if text.lower() in {"false", "0", "no"}:
        return "false"
    return text


def canonicalize_smiles(smiles: Any) -> str:
    if not HAS_RDKIT or is_missing(smiles):
        return ""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return ""
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def computed_formula(smiles: str) -> str:
    if not HAS_RDKIT or not smiles:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    return rdMolDescriptors.CalcMolFormula(mol)


def computed_exact_mass(smiles: str) -> float | None:
    if not HAS_RDKIT or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return float(Descriptors.ExactMolWt(mol))


def computed_inchikey(smiles: str) -> str:
    if not HAS_RDKIT or not smiles:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:
        return ""


def candidate_id_for_smiles(canonical_smiles: str) -> str:
    digest = hashlib.sha1(canonical_smiles.encode("utf-8")).hexdigest()[:10]
    return f"cand_{digest}"


def load_v15_name_map(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, low_memory=False)
    if "control_id" in df.columns:
        df = df.loc[df["control_id"].astype(str).eq("X-12798")].copy()
    out: dict[str, dict[str, Any]] = {}
    for _, row in df.iterrows():
        smiles = row.get("canonical_smiles_for_name") or row.get("predicted_smiles")
        canonical = canonicalize_smiles(smiles)
        if not canonical or canonical in out:
            continue
        out[canonical] = {
            "display_name": clean(row.get("predicted_name")) or clean(row.get("candidate_name")),
            "name_source": clean(row.get("predicted_name_source")),
            "name_tier": clean(row.get("predicted_name_tier")),
            "candidate_source": clean(row.get("candidate_source")),
            "source_labels": clean(row.get("source_labels")),
        }
    return out


def name_for_candidate(row: dict[str, Any], v15_names: dict[str, dict[str, Any]]) -> str:
    canonical = clean(row.get("canonical_smiles"))
    for key in ("display_name", "predicted_name", "candidate_name", "final_smiles_name"):
        if clean(row.get(key)):
            return clean(row.get(key))
    if canonical in v15_names and clean(v15_names[canonical].get("display_name")):
        return clean(v15_names[canonical]["display_name"])
    family = clean(row.get("candidate_family"))
    rank = clean(row.get("rank"))
    if family.startswith("pure_denovo"):
        branch = clean(row.get("source_branch"))
        return f"Pure de novo {branch} rank {rank}".strip()
    return f"Candidate {candidate_id_for_smiles(canonical)[-6:]}" if canonical else "Invalid SMILES"


def base_candidate(
    *,
    smiles: Any,
    family: str,
    source_section: str,
    source_branch: str = "",
    rank: Any = "",
    score: Any = "",
    display_name: Any = "",
    ms2_embedding_cosine: Any = "",
    ms2_embedding_rank: Any = "",
    ms2_embedding_status: Any = "",
    formula: Any = "",
    exact_mass: Any = "",
    adduct_consistent: Any = "",
    isotope_status: Any = "",
    predicted_pathway: Any = "",
    predicted_subpathway: Any = "",
    promotion_status: Any = "",
    gate_status: Any = "",
    provenance: Any = "",
) -> dict[str, Any]:
    canonical = canonicalize_smiles(smiles)
    candidate_id = candidate_id_for_smiles(canonical) if canonical else ""
    exact_mass_value = as_float(exact_mass)
    if exact_mass_value is None and canonical:
        exact_mass_value = computed_exact_mass(canonical)
    formula_value = clean(formula) or computed_formula(canonical)
    return {
        "candidate_id": candidate_id,
        "canonical_smiles": canonical,
        "display_name": clean(display_name),
        "candidate_family": family,
        "source_section": clean(source_section),
        "source_branch": clean(source_branch),
        "rank": clean(rank),
        "score": clean(score),
        "ms2_embedding_cosine": clean(ms2_embedding_cosine),
        "ms2_embedding_rank": clean(ms2_embedding_rank),
        "ms2_embedding_status": clean(ms2_embedding_status),
        "formula": formula_value,
        "exact_mass": "" if exact_mass_value is None else f"{exact_mass_value:.6f}",
        "adduct_consistent": as_bool_text(adduct_consistent),
        "isotope_status": clean(isotope_status),
        "predicted_pathway": clean(predicted_pathway),
        "predicted_subpathway": clean(predicted_subpathway),
        "promotion_status": clean(promotion_status),
        "gate_status": clean(gate_status),
        "candidate_consensus_asset": "",
        "candidate_consensus_mean_weighted_fraction": "",
        "candidate_consensus_atoms_ge_0p5_fraction": "",
        "provenance": clean(provenance),
        "metabolon_library_filter_status": "not_audited",
        "metabolon_detected_library_ids": "",
        "metabolon_library_match_basis": "",
    }


def add_deduped_candidate(
    candidates: dict[str, dict[str, Any]],
    row: dict[str, Any],
    v15_names: dict[str, dict[str, Any]],
) -> None:
    canonical = clean(row.get("canonical_smiles"))
    if not canonical:
        return
    row = {column: row.get(column, "") for column in CANDIDATE_COLUMNS}
    if not clean(row.get("display_name")):
        row["display_name"] = name_for_candidate(row, v15_names)
    existing = candidates.get(canonical)
    if existing is None:
        candidates[canonical] = row
        return

    existing_priority = FAMILY_PRIORITY.get(clean(existing.get("candidate_family")), 99)
    new_priority = FAMILY_PRIORITY.get(clean(row.get("candidate_family")), 99)
    winner, loser = (row, existing) if new_priority < existing_priority else (existing, row)
    source_sections = sorted(
        {part for value in [existing.get("source_section"), row.get("source_section")] for part in clean(value).split("|") if part}
    )
    provenance = sorted(
        {part for value in [existing.get("provenance"), row.get("provenance")] for part in clean(value).split("|") if part}
    )
    branches = sorted(
        {part for value in [existing.get("source_branch"), row.get("source_branch")] for part in clean(value).split("|") if part}
    )
    merged = dict(winner)
    merged["source_section"] = "|".join(source_sections)
    merged["source_branch"] = "|".join(branches)
    merged["provenance"] = "|".join(provenance)
    winner_name = clean(merged.get("display_name"))
    loser_name = clean(loser.get("display_name"))
    if (not winner_name or winner_name.startswith("Candidate ")) and loser_name and not loser_name.startswith("Candidate "):
        merged["display_name"] = clean(loser.get("display_name"))
    candidates[canonical] = merged


def markdown_smiles_candidates(markdown_path: Path) -> list[str]:
    if not markdown_path.exists():
        return []
    found: list[str] = []
    seen: set[str] = set()
    for line in markdown_path.read_text().splitlines():
        if "|" not in line:
            continue
        for cell in line.split("|"):
            text = cell.strip().strip("`")
            if len(text) < 5 or "/" in text or " " in text:
                continue
            if not re.search(r"[=#()\[\]@]", text) and not re.search(r"[cnosp][0-9]", text):
                continue
            canonical = canonicalize_smiles(text)
            if canonical and canonical not in seen:
                seen.add(canonical)
                found.append(canonical)
    return found


def collect_candidates(args: argparse.Namespace) -> pd.DataFrame:
    v15_names = load_v15_name_map(args.v15_names_csv)
    candidates: dict[str, dict[str, Any]] = {}
    deep_final = read_json(args.deep_dive_final_json)
    lib400_final = read_json(args.lib400_final_json)
    merged_final = read_json(args.merged_final_json)

    add_deduped_candidate(
        candidates,
        base_candidate(
            smiles=deep_final.get("final_smiles"),
            display_name=deep_final.get("final_smiles_name"),
            family="provisional_non_pure",
            source_section="deep_dive_final_call",
            rank=1,
            score=deep_final.get("final_smiles_score"),
            formula=deep_final.get("target_formula"),
            predicted_subpathway=deep_final.get("final_subpathway"),
            promotion_status="not_pure_denovo",
            gate_status="provisional_non_pure_structure",
            provenance="models/target_unknown_deep_dive_v1/x12798/final_call.json",
        ),
        v15_names,
    )

    deep = pd.read_csv(args.deep_dive_csv, low_memory=False)
    deep["deep_dive_rank_numeric"] = pd.to_numeric(deep.get("deep_dive_rank"), errors="coerce")
    for _, row in deep.loc[deep["deep_dive_rank_numeric"].le(10)].sort_values("deep_dive_rank_numeric").iterrows():
        add_deduped_candidate(
            candidates,
            base_candidate(
                smiles=row.get("canonical_smiles") or row.get("predicted_smiles"),
                display_name=row.get("predicted_name") or row.get("candidate_name"),
                family="non_pure_deep_dive",
                source_section="non_pure_deep_dive_top10",
                rank=row.get("deep_dive_rank"),
                score=row.get("deep_dive_score"),
                formula=row.get("predicted_formula_computed") or row.get("target_formula"),
                exact_mass=row.get("exact_mass_computed"),
                adduct_consistent=row.get("mass_within_5mda"),
                isotope_status=row.get("isotope_status") or row.get("isotope_audit_status"),
                predicted_pathway=row.get("predicted_pathway") or row.get("pathway"),
                predicted_subpathway=row.get("predicted_subpathway") or row.get("subpathway"),
                promotion_status="not_pure_denovo",
                gate_status=row.get("confidence_notes") or row.get("formula_status"),
                provenance=row.get("source_provenance") or row.get("source_labels"),
            ),
            v15_names,
        )

    comparison = pd.read_csv(args.encoder_comparison_csv, low_memory=False)
    for _, row in comparison.iterrows():
        add_deduped_candidate(
            candidates,
            base_candidate(
                smiles=row.get("top_smiles"),
                family="encoder_branch_top",
                source_section="encoder_branch_comparison",
                source_branch=row.get("branch_label"),
                rank=row.get("top_rank"),
                score=row.get("top_score"),
                ms2_embedding_cosine=row.get("top_ms2_embedding_cosine"),
                ms2_embedding_rank=row.get("top_ms2_embedding_rank"),
                ms2_embedding_status=row.get("top_ms2_embedding_status"),
                adduct_consistent=row.get("top_adduct_consistent"),
                isotope_status=row.get("top_isotope_status"),
                predicted_pathway=row.get("top_predicted_pathway"),
                predicted_subpathway=row.get("top_predicted_subpathway"),
                promotion_status=(
                    "isotope_adduct_promoted"
                    if str(row.get("top_isotope_adduct_promoted")).lower() == "true"
                    else "not_promoted"
                ),
                gate_status=(
                    "ms2_gate_pass"
                    if str(row.get("ms2_embedding_gate_pass")).lower() == "true"
                    else "ms2_embedding_gate_failed"
                ),
                provenance=f"encoder_branch_comparison:{clean(row.get('branch_label'))}",
            ),
            v15_names,
        )

    for path, branch, family, final in [
        (args.lib400_pure_top10, "lib400_only", "pure_denovo_lib400_only", lib400_final),
        (args.merged_pure_top10, "merged_libcond", "pure_denovo_merged_libcond", merged_final),
    ]:
        pure = pd.read_csv(path, low_memory=False)
        for _, row in pure.iterrows():
            add_deduped_candidate(
                candidates,
                base_candidate(
                    smiles=row.get("canonical_smiles"),
                    family=family,
                    source_section=f"{branch}_pure_denovo_top10",
                    source_branch=branch,
                    rank=row.get("pure_denovo_rank"),
                    score=row.get("constrained_rank_score"),
                    ms2_embedding_cosine=row.get("ms2_embedding_cosine"),
                    ms2_embedding_rank=row.get("ms2_embedding_rank"),
                    ms2_embedding_status=row.get("ms2_embedding_match_status"),
                    formula=row.get("predicted_formula") or row.get("target_formula"),
                    exact_mass=row.get("exact_mass"),
                    adduct_consistent=row.get("adduct_consistent"),
                    isotope_status=row.get("isotope_status") or row.get("isotope_audit_status"),
                    predicted_pathway=row.get("predicted_pathway") or row.get("pathway"),
                    predicted_subpathway=row.get("predicted_subpathway") or row.get("subpathway"),
                    promotion_status=(
                        "isotope_adduct_promoted"
                        if str(row.get("isotope_adduct_promoted")).lower() == "true"
                        else clean(row.get("promotion_block_reason"), "not_promoted")
                    ),
                    gate_status=(
                        "ms2_gate_pass"
                        if str(row.get("ms2_embedding_structural_gate_pass")).lower() == "true"
                        else clean(row.get("confidence_notes"), "ms2_embedding_gate_failed")
                    ),
                    provenance=(
                        f"{branch}_final:{clean(final.get('final_call_status'))}|"
                        f"operators:{clean(row.get('operator_set'))}|seeds:{clean(row.get('seed_source_provenance'))}"
                    ),
                ),
                v15_names,
            )

    for smiles in markdown_smiles_candidates(args.resolution_md):
        add_deduped_candidate(
            candidates,
            base_candidate(
                smiles=smiles,
                family="non_pure_deep_dive",
                source_section="resolution_markdown_candidate",
                promotion_status="markdown_context",
                gate_status="candidate_from_resolution_markdown",
                provenance=str(args.resolution_md),
            ),
            v15_names,
        )

    rows = list(candidates.values())
    rows.sort(
        key=lambda row: (
            FAMILY_PRIORITY.get(clean(row.get("candidate_family")), 99),
            as_float(row.get("rank")) if as_float(row.get("rank")) is not None else 9999.0,
            clean(row.get("candidate_id")),
        )
    )
    return pd.DataFrame(rows, columns=CANDIDATE_COLUMNS)


def split_library_ids(value: Any) -> list[str]:
    if is_missing(value):
        return []
    out: list[str] = []
    for part in re.split(r"[;,|]", str(value)):
        cleaned = part.strip()
        if cleaned and cleaned.lower() not in {"nan", "none"}:
            out.append(cleaned)
    return out


def parse_markdown_table(section_lines: list[str]) -> pd.DataFrame:
    table_lines = [line for line in section_lines if line.strip().startswith("|")]
    if len(table_lines) < 3:
        return pd.DataFrame()
    def split_row(line: str) -> list[str]:
        placeholder = "\u241f"
        protected = line.replace(r"\|", placeholder)
        return [cell.strip().replace(placeholder, "|") for cell in protected.strip().strip("|").split("|")]

    header = split_row(table_lines[0])
    rows: list[dict[str, Any]] = []
    for line in table_lines[2:]:
        cells = split_row(line)
        if len(cells) != len(header):
            continue
        row = dict(zip(header, cells))
        rows.append(row)
    return pd.DataFrame(rows)


def markdown_section(path: Path, heading: str) -> list[str]:
    lines = path.read_text().splitlines()
    start = None
    for index, line in enumerate(lines):
        if line.strip() == heading:
            start = index + 1
            break
    if start is None:
        return []
    section: list[str] = []
    for line in lines[start:]:
        if line.startswith("## "):
            break
        section.append(line)
    return section


def parse_library_match_markdown(path: Path) -> pd.DataFrame:
    hits = parse_markdown_table(markdown_section(path, "## Metabolon Library Hits"))
    if hits.empty:
        return pd.DataFrame()
    rename = {
        "rank": "md_rank",
        "dossier predicted match": "md_predicted_match",
        "candidate InChIKey": "md_candidate_inchikey",
        "DB chemical_id": "db_chemical_id",
        "DB canonical_name": "db_canonical_name",
        "match basis": "md_match_basis",
        "in_global_library": "md_in_global_library",
        "global_library_status": "md_global_library_status",
        "pathway": "md_pathway",
    }
    hits = hits.rename(columns=rename)
    for column in rename.values():
        if column not in hits.columns:
            hits[column] = ""
    hits["md_rank"] = pd.to_numeric(hits["md_rank"].astype(str).str.replace(":", "", regex=False), errors="coerce")
    hits["md_source_path"] = str(path)
    return hits[list(rename.values()) + ["md_source_path"]]


def load_local_metabolon_catalog(catalog_paths: list[Path]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for path in catalog_paths:
        if not path.exists():
            continue
        df = pd.read_csv(path, low_memory=False)
        if "canonical_smiles" not in df.columns:
            continue
        for _, row in df.iterrows():
            canonical = canonicalize_smiles(row.get("canonical_smiles"))
            if not canonical:
                continue
            library_ids = set(split_library_ids(row.get("lib_ids")))
            library_ids.update(split_library_ids(row.get("parent_mz_source_lib")))
            rows.append(
                {
                    "catalog_canonical_smiles": canonical,
                    "catalog_inchikey": computed_inchikey(canonical),
                    "catalog_comp_id": clean(row.get("comp_id")),
                    "catalog_compound_name": clean(row.get("compound_name")),
                    "catalog_library_ids": ";".join(sorted(library_ids, key=str)),
                    "catalog_source_path": str(path),
                    "catalog_smiles_source_csv": clean(row.get("smiles_source_csv")),
                }
            )
    return pd.DataFrame(rows)


def build_metabolon_library_audit(
    evidence: pd.DataFrame,
    *,
    library_match_md: Path,
    catalog_paths: list[Path],
    target_library_id: str = "400",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    catalog = load_local_metabolon_catalog(catalog_paths)
    md_hits = parse_library_match_markdown(library_match_md)
    audit_rows: list[dict[str, Any]] = []

    by_rank = {}
    by_inchikey = {}
    if not md_hits.empty:
        for _, row in md_hits.iterrows():
            rank = as_float(row.get("md_rank"))
            if rank is not None:
                by_rank[int(rank)] = row.to_dict()
            if clean(row.get("md_candidate_inchikey")):
                by_inchikey[clean(row.get("md_candidate_inchikey"))] = row.to_dict()

    for _, row in evidence.iterrows():
        canonical = clean(row.get("canonical_smiles"))
        inchikey = computed_inchikey(canonical)
        rank_float = as_float(row.get("rank"))
        rank_int = int(rank_float) if rank_float is not None else None
        md = by_inchikey.get(inchikey) or {}
        if not md and rank_int is not None and clean(row.get("candidate_family")) in {
            "provisional_non_pure",
            "non_pure_deep_dive",
        }:
            md = by_rank.get(rank_int) or {}
        matches = pd.DataFrame()
        if not catalog.empty:
            matches = catalog.loc[
                (catalog["catalog_canonical_smiles"].eq(canonical))
                | (catalog["catalog_inchikey"].astype(str).eq(inchikey))
            ].copy()
        library_ids: set[str] = set()
        for value in matches.get("catalog_library_ids", pd.Series(dtype=str)).tolist():
            library_ids.update(split_library_ids(value))
        same_library = target_library_id in library_ids
        has_catalog = not matches.empty
        has_md = bool(md)
        if same_library:
            action = f"filtered_same_library_{target_library_id}"
        elif library_ids:
            action = "listed_cross_library_only"
        elif has_md:
            action = "listed_global_library_hit_no_local_library_id"
        else:
            action = "kept_no_metabolon_library_hit"
        if has_catalog and has_md:
            status = "exact_local_catalog_and_global_library_match"
        elif has_catalog:
            status = "exact_local_catalog_match"
        elif has_md:
            status = "global_library_hit_from_match_md"
        else:
            status = "no_metabolon_library_hit"
        audit_rows.append(
            {
                "candidate_id": row.get("candidate_id", ""),
                "canonical_smiles": canonical,
                "display_name": row.get("display_name", ""),
                "candidate_family": row.get("candidate_family", ""),
                "candidate_rank": row.get("rank", ""),
                "candidate_inchikey": inchikey,
                "metabolon_match_status": status,
                "detected_library_ids": ";".join(sorted(library_ids, key=str)),
                "same_library_400": same_library,
                "cross_library_only": bool(library_ids) and not same_library,
                "filter_action": action,
                "catalog_comp_ids": "|".join(sorted(set(matches.get("catalog_comp_id", pd.Series(dtype=str)).astype(str).tolist()))),
                "catalog_compound_names": "|".join(
                    sorted({name for name in matches.get("catalog_compound_name", pd.Series(dtype=str)).astype(str).tolist() if name and name != "nan"})
                ),
                "catalog_source_paths": "|".join(sorted(set(matches.get("catalog_source_path", pd.Series(dtype=str)).astype(str).tolist()))),
                "md_rank": md.get("md_rank", ""),
                "db_chemical_id": md.get("db_chemical_id", ""),
                "db_canonical_name": md.get("db_canonical_name", ""),
                "match_basis": md.get("md_match_basis", "exact local catalog canonical SMILES/InChIKey" if has_catalog else ""),
                "in_global_library": md.get("md_in_global_library", ""),
                "global_library_status": md.get("md_global_library_status", ""),
                "pathway": md.get("md_pathway", ""),
                "notes": (
                    f"Excluded from candidate evidence/network because exact local catalog match includes library {target_library_id}."
                    if same_library
                    else (
                        f"Detected in Metabolon local catalog only on library ids {', '.join(sorted(library_ids, key=str))}; listed separately."
                        if library_ids
                        else (
                            "Global-library hit in match markdown, but no exact local library-id catalog row was found."
                            if has_md
                            else ""
                        )
                    )
                ),
            }
        )

    evidence_ids = set(evidence["candidate_id"].astype(str))
    for _, md in md_hits.iterrows():
        md_rank = as_float(md.get("md_rank"))
        matched = False
        if md_rank is not None:
            matched = any(
                (as_float(row.get("rank")) == md_rank)
                and clean(row.get("candidate_family")) in {"provisional_non_pure", "non_pure_deep_dive"}
                for _, row in evidence.iterrows()
            )
        if matched:
            continue
        audit_rows.append(
            {
                "candidate_id": "",
                "canonical_smiles": "",
                "display_name": md.get("md_predicted_match", ""),
                "candidate_family": "not_in_report_candidate_set",
                "candidate_rank": md.get("md_rank", ""),
                "candidate_inchikey": md.get("md_candidate_inchikey", ""),
                "metabolon_match_status": "global_library_hit_from_match_md",
                "detected_library_ids": "",
                "same_library_400": False,
                "cross_library_only": False,
                "filter_action": "not_in_report_candidate_set",
                "catalog_comp_ids": "",
                "catalog_compound_names": "",
                "catalog_source_paths": "",
                "md_rank": md.get("md_rank", ""),
                "db_chemical_id": md.get("db_chemical_id", ""),
                "db_canonical_name": md.get("db_canonical_name", ""),
                "match_basis": md.get("md_match_basis", ""),
                "in_global_library": md.get("md_in_global_library", ""),
                "global_library_status": md.get("md_global_library_status", ""),
                "pathway": md.get("md_pathway", ""),
                "notes": "Listed from x12798_metabolon_library_matches.md; not part of this report's best-candidate network.",
            }
        )

    audit = pd.DataFrame(audit_rows)
    filtered = audit.loc[audit["filter_action"].astype(str).str.startswith("filtered_same_library_")].copy()
    cross_library = audit.loc[audit["filter_action"].eq("listed_cross_library_only")].copy()
    keep_ids = set(evidence["candidate_id"].astype(str)) - set(filtered["candidate_id"].astype(str))
    filtered_evidence = evidence.loc[evidence["candidate_id"].astype(str).isin(keep_ids)].copy()
    status_map = audit.drop_duplicates(subset=["candidate_id"], keep="first").set_index("candidate_id").to_dict(orient="index")
    for index, row in filtered_evidence.iterrows():
        status = status_map.get(str(row["candidate_id"]), {})
        filtered_evidence.at[index, "metabolon_library_filter_status"] = status.get("filter_action", "")
        filtered_evidence.at[index, "metabolon_detected_library_ids"] = status.get("detected_library_ids", "")
        filtered_evidence.at[index, "metabolon_library_match_basis"] = status.get("match_basis", "")
    return audit, filtered, cross_library, filtered_evidence


def mol_and_fp(smiles: str):
    if not HAS_RDKIT:
        return None, None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None, None
    return mol, MORGAN.GetFingerprint(mol)


def similarity_bucket(value: float) -> str:
    if value >= 0.70:
        return "very_high"
    if value >= 0.50:
        return "high"
    if value >= 0.25:
        return "moderate"
    return "low"


def known_similarity_bucket(value: float) -> str:
    if value >= 0.50:
        return "medium"
    if value >= 0.35:
        return "context"
    return "low"


def known_color(value: float) -> str:
    if value >= 0.50:
        return "#D73A49"
    if value >= 0.35:
        return "#F28E2B"
    return "#59A14F"


def edge_color(value: float) -> str:
    if value >= 0.70:
        return "#8B1A1A"
    if value >= 0.50:
        return "#D73A49"
    if value >= 0.25:
        return "#F28E2B"
    return "#B7C0C7"


def heatmap_color(value: float) -> str:
    value = max(0.0, min(1.0, float(value or 0.0)))
    stops = [
        (0.0, (44, 123, 182)),
        (0.5, (253, 174, 97)),
        (1.0, (215, 25, 28)),
    ]
    for (left_v, left_rgb), (right_v, right_rgb) in zip(stops, stops[1:]):
        if value <= right_v:
            frac = (value - left_v) / max(right_v - left_v, 1e-9)
            rgb = tuple(round(left_rgb[i] + frac * (right_rgb[i] - left_rgb[i])) for i in range(3))
            return f"#{rgb[0]:02X}{rgb[1]:02X}{rgb[2]:02X}"
    return "#D7191C"


def slugify(value: Any) -> str:
    text = re.sub(r"[^a-z0-9]+", "-", clean(value).lower()).strip("-")
    return text[:80] or "unknown"


def split_structure_values(value: Any) -> list[str]:
    if is_missing(value):
        return []
    text = clean(value)
    for separator in ["|", ";"]:
        text = text.replace(separator, ",")
    return [part.strip() for part in text.split(",") if part.strip() and part.strip().lower() not in {"nan", "none", "null"}]


def load_global_library_structures(reference_db: Path) -> pd.DataFrame:
    if not reference_db.exists():
        raise FileNotFoundError(f"Reference DB not found: {reference_db}")
    with sqlite3.connect(reference_db) as conn:
        raw = pd.read_sql_query(
            """
            SELECT chemical_id, canonical_name, smiles, inchikeys, super_pathway, pathway,
                   metabolon_super, metabolon_sub, in_global_library, global_library_status
            FROM metabolites
            WHERE smiles IS NOT NULL AND trim(smiles) != ''
              AND (
                lower(coalesce(in_global_library, '')) LIKE '%true%'
                OR lower(coalesce(global_library_status, '')) LIKE 'yes%'
              )
            """,
            conn,
        )
    rows: list[dict[str, Any]] = []
    for _, row in raw.iterrows():
        for smiles in split_structure_values(row.get("smiles")):
            canonical = canonicalize_smiles(smiles)
            if not canonical:
                continue
            _, fingerprint = mol_and_fp(canonical)
            if fingerprint is None:
                continue
            rows.append(
                {
                    "chemical_id": clean(row.get("chemical_id")),
                    "canonical_name": clean(row.get("canonical_name")),
                    "canonical_smiles": canonical,
                    "inchikeys": clean(row.get("inchikeys")),
                    "super_pathway": clean(row.get("super_pathway")),
                    "pathway": clean(row.get("pathway")),
                    "metabolon_super": clean(row.get("metabolon_super")),
                    "metabolon_sub": clean(row.get("metabolon_sub")),
                    "in_global_library": clean(row.get("in_global_library")),
                    "global_library_status": clean(row.get("global_library_status")),
                    "_fingerprint": fingerprint,
                }
            )
    return pd.DataFrame(rows)


def build_anchor_known_similarity(
    anchor_smiles: str,
    reference_db: Path,
    top_n: int = 30,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    if not HAS_RDKIT:
        empty_neighbors = pd.DataFrame(
            columns=[
                "chemical_id",
                "canonical_name",
                "canonical_smiles",
                "inchikeys",
                "similarity_to_anchor",
                "similarity_bucket",
                "super_pathway",
                "pathway",
                "metabolon_super",
                "metabolon_sub",
                "in_global_library",
                "global_library_status",
                "asset_path",
            ]
        )
        empty_nodes = pd.DataFrame()
        empty_edges = pd.DataFrame()
        network_json = {
            "nodes": [],
            "edges": [],
            "rendered_edge_similarity_floor": 0.35,
            "layout": "deterministic_known_neighbor_layout_unavailable_without_rdkit",
        }
        return empty_neighbors, empty_neighbors.copy(), empty_nodes, empty_edges, network_json
    anchor_canonical = canonicalize_smiles(anchor_smiles)
    _, anchor_fp = mol_and_fp(anchor_canonical)
    if anchor_fp is None:
        raise ValueError(f"Could not fingerprint anchor SMILES: {anchor_smiles}")
    known = load_global_library_structures(reference_db)
    all_rows: list[dict[str, Any]] = []
    for _, row in known.iterrows():
        similarity = float(TanimotoSimilarity(anchor_fp, row["_fingerprint"]))
        bucket = known_similarity_bucket(similarity)
        all_rows.append(
            {
                "chemical_id": row["chemical_id"],
                "canonical_name": row["canonical_name"],
                "canonical_smiles": row["canonical_smiles"],
                "inchikeys": row["inchikeys"],
                "similarity_to_anchor": round(similarity, 6),
                "similarity_bucket": bucket,
                "super_pathway": row["super_pathway"],
                "pathway": row["pathway"],
                "metabolon_super": row["metabolon_super"],
                "metabolon_sub": row["metabolon_sub"],
                "in_global_library": row["in_global_library"],
                "global_library_status": row["global_library_status"],
                "asset_path": "",
            }
        )
    all_known = pd.DataFrame(all_rows).sort_values(
        ["similarity_to_anchor", "chemical_id", "canonical_name"],
        ascending=[False, True, True],
        kind="mergesort",
    )
    neighbors = (
        all_known.drop_duplicates(subset=["canonical_smiles"], keep="first")
        .sort_values(["similarity_to_anchor", "chemical_id", "canonical_name"], ascending=[False, True, True], kind="mergesort")
        .head(top_n)
        .copy()
        .reset_index(drop=True)
    )
    neighbors.insert(0, "rank", range(1, len(neighbors) + 1))

    node_rows: list[dict[str, Any]] = [
        {
            "id": "anchor_best_candidate",
            "label": f"X-12798 best candidate: {ANCHOR_CANDIDATE_NAME}",
            "node_type": "anchor_best_candidate",
            "chemical_id": "",
            "canonical_name": ANCHOR_CANDIDATE_NAME,
            "canonical_smiles": anchor_canonical,
            "similarity_to_anchor": 1.0,
            "similarity_bucket": "anchor",
            "super_pathway": "predicted",
            "pathway": "benzamide-like xenobiotic/benzoate-family",
            "x": 500.0,
            "y": 310.0,
            "radius": 22.0,
            "color": "#1F4E79",
            "asset_path": "assets/x12798_anchor_best_candidate.svg",
        }
    ]
    edge_rows: list[dict[str, Any]] = []
    n = max(len(neighbors), 1)
    for index, row in neighbors.iterrows():
        angle = (2.0 * math.pi * index / n) - math.pi / 2.0
        similarity = float(row["similarity_to_anchor"])
        x_pos = 500.0 + 280.0 * math.cos(angle)
        y_pos = 310.0 + 245.0 * math.sin(angle)
        node_id = f"known_{clean(row['chemical_id'])}_{int(row['rank'])}"
        node_rows.append(
            {
                "id": node_id,
                "label": clean(row["canonical_name"]),
                "node_type": "known_metabolite",
                "chemical_id": clean(row["chemical_id"]),
                "canonical_name": clean(row["canonical_name"]),
                "canonical_smiles": clean(row["canonical_smiles"]),
                "similarity_to_anchor": round(similarity, 6),
                "similarity_bucket": clean(row["similarity_bucket"]),
                "super_pathway": clean(row["super_pathway"]),
                "pathway": clean(row["pathway"]),
                "x": round(x_pos, 2),
                "y": round(y_pos, 2),
                "radius": round(9.0 + similarity * 16.0, 2),
                "color": known_color(similarity),
                "asset_path": clean(row.get("asset_path")),
            }
        )
        edge_rows.append(
            {
                "source": "anchor_best_candidate",
                "target": node_id,
                "similarity": round(similarity, 6),
                "similarity_bucket": clean(row["similarity_bucket"]),
                "label": f"{similarity:.2f}",
                "edge_type": "anchor_to_known",
            }
        )

    nodes = pd.DataFrame(node_rows)
    edges = pd.DataFrame(edge_rows)
    network_json = {
        "nodes": safe_json_records(nodes),
        "edges": safe_json_records(edges),
        "anchor_candidate_name": ANCHOR_CANDIDATE_NAME,
        "anchor_candidate_smiles": anchor_canonical,
        "layout": "best_candidate_centered_known_similarity",
        "rendered_edge_similarity_floor": 0.35,
    }
    return all_known.reset_index(drop=True), neighbors, nodes, edges, network_json


def mcs_atom_mapping(ref: Any, cand: Any) -> dict[int, int]:
    try:
        mcs = rdFMCS.FindMCS([ref, cand], timeout=2, matchValences=False, ringMatchesRingOnly=True)
    except Exception:
        return {}
    if not mcs.smartsString:
        return {}
    patt = Chem.MolFromSmarts(mcs.smartsString)
    if patt is None:
        return {}
    ref_match = ref.GetSubstructMatch(patt)
    cand_match = cand.GetSubstructMatch(patt)
    if not ref_match or not cand_match or len(ref_match) != len(cand_match):
        return {}
    return {int(ref_idx): int(cand_idx) for ref_idx, cand_idx in zip(ref_match, cand_match)}


def atom_token(mol: Any, atom_idx: int, include_bond: bool = False, origin_idx: int | None = None) -> str:
    atom = mol.GetAtomWithIdx(int(atom_idx))
    token = "|".join(
        [
            atom.GetSymbol(),
            "arom" if atom.GetIsAromatic() else "aliph",
            str(atom.GetHybridization()),
            f"chg{atom.GetFormalCharge()}",
            "ring" if atom.IsInRing() else "chain",
            f"deg{atom.GetDegree()}",
            f"h{atom.GetTotalNumHs()}",
        ]
    )
    if include_bond and origin_idx is not None:
        bond = mol.GetBondBetweenAtoms(int(origin_idx), int(atom_idx))
        bond_token = str(bond.GetBondType()) if bond is not None else "none"
        return f"{bond_token}:{token}"
    return token


def adjacent_carbonyl(mol: Any, atom_idx: int) -> bool:
    atom = mol.GetAtomWithIdx(int(atom_idx))
    atoms_to_check = [atom, *atom.GetNeighbors()]
    for candidate in atoms_to_check:
        if candidate.GetSymbol() != "C":
            continue
        for bond in candidate.GetBonds():
            other = bond.GetOtherAtom(candidate)
            if other.GetSymbol() == "O" and bond.GetBondTypeAsDouble() == 2.0:
                return True
    return False


def adjacent_amide(mol: Any, atom_idx: int) -> bool:
    atom = mol.GetAtomWithIdx(int(atom_idx))
    atoms_to_check = [atom, *atom.GetNeighbors()]
    for candidate in atoms_to_check:
        if candidate.GetSymbol() != "C":
            continue
        has_double_o = False
        has_n = False
        for bond in candidate.GetBonds():
            other = bond.GetOtherAtom(candidate)
            if other.GetSymbol() == "O" and bond.GetBondTypeAsDouble() == 2.0:
                has_double_o = True
            if other.GetSymbol() == "N":
                has_n = True
        if has_double_o and has_n:
            return True
    return False


def radius_environment_tokens(mol: Any, atom_idx: int, radius: int) -> set[str]:
    current = {int(atom_idx)}
    visited = {int(atom_idx)}
    tokens: set[str] = set()
    for distance in range(1, radius + 1):
        next_atoms: set[int] = set()
        for idx in current:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = int(neighbor.GetIdx())
                if neighbor_idx in visited:
                    continue
                next_atoms.add(neighbor_idx)
                tokens.add(f"r{distance}:{atom_token(mol, neighbor_idx, include_bond=True, origin_idx=idx)}")
        visited |= next_atoms
        current = next_atoms
    return tokens


def jaccard(left: set[str], right: set[str]) -> float:
    if not left and not right:
        return 1.0
    union = left | right
    return len(left & right) / max(len(union), 1)


def atom_environment_signature(mol: Any, atom_idx: int) -> dict[str, Any]:
    atom = mol.GetAtomWithIdx(int(atom_idx))
    ring_info = mol.GetRingInfo()
    ring_substituent_count = 0
    if atom.IsInRing():
        for neighbor in atom.GetNeighbors():
            if not neighbor.IsInRing():
                ring_substituent_count += 1
    exocyclic_bond_count = 0
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if atom.IsInRing() and not other.IsInRing():
            exocyclic_bond_count += 1
    return {
        "center_symbol": atom.GetSymbol(),
        "center_aromatic": atom.GetIsAromatic(),
        "center_hybridization": str(atom.GetHybridization()),
        "center_formal_charge": atom.GetFormalCharge(),
        "center_ring": atom.IsInRing(),
        "center_degree": atom.GetDegree(),
        "center_total_h": atom.GetTotalNumHs(),
        "neighbor_signature_radius1": radius_environment_tokens(mol, atom_idx, 1),
        "neighbor_signature_radius2": radius_environment_tokens(mol, atom_idx, 2),
        "adjacent_hetero_count": sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() not in {"C", "H"}),
        "adjacent_carbonyl": adjacent_carbonyl(mol, atom_idx),
        "adjacent_amide": adjacent_amide(mol, atom_idx),
        "ring_substituent_count": ring_substituent_count,
        "exocyclic_bond_count": exocyclic_bond_count,
        "ring_count": ring_info.NumRings(),
    }


def atom_environment_similarity(ref: Any, ref_idx: int, cand: Any, cand_idx: int) -> float:
    ref_sig = atom_environment_signature(ref, ref_idx)
    cand_sig = atom_environment_signature(cand, cand_idx)
    center_keys = [
        "center_symbol",
        "center_aromatic",
        "center_hybridization",
        "center_formal_charge",
        "center_ring",
        "center_degree",
        "center_total_h",
        "adjacent_hetero_count",
        "adjacent_carbonyl",
        "adjacent_amide",
        "ring_substituent_count",
        "exocyclic_bond_count",
    ]
    center_score = sum(1.0 for key in center_keys if ref_sig[key] == cand_sig[key]) / len(center_keys)
    radius1_score = jaccard(ref_sig["neighbor_signature_radius1"], cand_sig["neighbor_signature_radius1"])
    radius2_score = jaccard(ref_sig["neighbor_signature_radius2"], cand_sig["neighbor_signature_radius2"])
    return 0.45 * center_score + 0.35 * radius1_score + 0.20 * radius2_score


def consensus_atom_support(reference_smiles: str, candidates: pd.DataFrame) -> tuple[list[dict[str, Any]], str]:
    if not HAS_RDKIT or candidates.empty:
        return [], "<p class='note'>Scaffold position consensus unavailable because RDKit or candidate rows are missing.</p>"
    ref = Chem.MolFromSmiles(reference_smiles)
    if ref is None:
        return [], "<p class='note'>Scaffold position consensus unavailable because the orientation SMILES could not be parsed.</p>"
    rdDepictor.Compute2DCoords(ref)
    n_atoms = ref.GetNumAtoms()
    support = [0.0] * n_atoms
    weighted_support = [0.0] * n_atoms
    valid = 0
    total_weight = 0.0

    for _, candidate in candidates.iterrows():
        smiles = clean(candidate.get("canonical_smiles"))
        if not smiles:
            continue
        cand = Chem.MolFromSmiles(smiles)
        if cand is None:
            continue
        try:
            mcs = rdFMCS.FindMCS([ref, cand], timeout=2, matchValences=False, ringMatchesRingOnly=True)
        except Exception:
            continue
        if not mcs.smartsString:
            continue
        patt = Chem.MolFromSmarts(mcs.smartsString)
        if patt is None:
            continue
        match = ref.GetSubstructMatch(patt)
        if not match:
            continue
        valid += 1
        score = as_float(candidate.get("score"))
        weight = max(0.05, min(float(score) if score is not None else 0.25, 1.0))
        total_weight += weight
        for atom_idx in match:
            support[int(atom_idx)] += 1.0
            weighted_support[int(atom_idx)] += weight

    confidence = [value / valid if valid else 0.0 for value in support]
    weighted_confidence = [value / total_weight if total_weight else 0.0 for value in weighted_support]

    conf = ref.GetConformer()
    coords = []
    for idx in range(n_atoms):
        pos = conf.GetAtomPosition(idx)
        coords.append((pos.x, pos.y))
    min_x = min(x for x, _ in coords)
    max_x = max(x for x, _ in coords)
    min_y = min(y for _, y in coords)
    max_y = max(y for _, y in coords)
    width = max(max_x - min_x, 1e-6)
    height = max(max_y - min_y, 1e-6)
    canvas_w, canvas_h, margin = 760.0, 450.0, 62.0
    scale = min((canvas_w - 2 * margin) / width, (canvas_h - 2 * margin) / height)

    def xy(atom_idx: int) -> tuple[float, float]:
        x, y = coords[atom_idx]
        return margin + (x - min_x) * scale, canvas_h - (margin + (y - min_y) * scale)

    bond_parts = []
    for bond in ref.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        x1, y1 = xy(a)
        x2, y2 = xy(b)
        bond_parts.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            'stroke="#24445F" stroke-opacity=".42" stroke-width="4" stroke-linecap="round" />'
        )

    rows: list[dict[str, Any]] = []
    atom_parts = []
    for atom in ref.GetAtoms():
        idx = atom.GetIdx()
        x_pos, y_pos = xy(idx)
        value = weighted_confidence[idx]
        color = heatmap_color(value)
        symbol = html.escape(atom.GetSymbol())
        atom_parts.append(
            f'<g class="consensus-atom" data-confidence="{value:.3f}" transform="translate({x_pos:.1f},{y_pos:.1f})">'
            f'<circle r="24" fill="{color}" stroke="white" stroke-width="3">'
            f"<title>Atom {idx + 1} {symbol}: weighted consensus {value:.2f}</title></circle>"
            f'<text text-anchor="middle" y="5">{symbol}</text>'
            f'<text class="atom-confidence" text-anchor="middle" y="39">{value:.2f}</text></g>'
        )
        rows.append(
            {
                "atom_index": idx + 1,
                "element": atom.GetSymbol(),
                "consensus_support_fraction": round(confidence[idx], 6),
                "weighted_consensus_support_fraction": round(value, 6),
                "support_count": round(support[idx], 6),
                "weighted_support": round(weighted_support[idx], 6),
                "valid_candidate_count": valid,
                "total_candidate_weight": round(total_weight, 6),
                "heatmap_color": color,
            }
        )

    svg = "\n".join(
        [
            "<svg class='consensus-svg' viewBox='0 0 760 450' role='img' aria-label='Scaffold position consensus heatmap'>",
            *bond_parts,
            *atom_parts,
            "</svg>",
        ]
    )
    return rows, svg


def environment_aware_consensus_atom_support(reference_smiles: str, candidates: pd.DataFrame) -> tuple[list[dict[str, Any]], str]:
    if not HAS_RDKIT or candidates.empty:
        return [], "<p class='note'>Chemical environment consensus unavailable because RDKit or candidate rows are missing.</p>"
    ref = Chem.MolFromSmiles(reference_smiles)
    if ref is None:
        return [], "<p class='note'>Chemical environment consensus unavailable because the orientation SMILES could not be parsed.</p>"
    rdDepictor.Compute2DCoords(ref)
    n_atoms = ref.GetNumAtoms()
    graph_support = [0.0] * n_atoms
    weighted_graph_support = [0.0] * n_atoms
    environment_support = [0.0] * n_atoms
    weighted_environment_support = [0.0] * n_atoms
    valid = 0
    total_weight = 0.0

    for _, candidate in candidates.iterrows():
        smiles = clean(candidate.get("canonical_smiles"))
        if not smiles:
            continue
        cand = Chem.MolFromSmiles(smiles)
        if cand is None:
            continue
        mapping = mcs_atom_mapping(ref, cand)
        if not mapping:
            continue
        valid += 1
        score = as_float(candidate.get("score"))
        weight = max(0.05, min(float(score) if score is not None else 0.25, 1.0))
        total_weight += weight
        for ref_idx, cand_idx in mapping.items():
            graph_support[int(ref_idx)] += 1.0
            weighted_graph_support[int(ref_idx)] += weight
            env_similarity = atom_environment_similarity(ref, ref_idx, cand, cand_idx)
            if env_similarity >= 0.70:
                environment_support[int(ref_idx)] += 1.0
                weighted_environment_support[int(ref_idx)] += weight

    graph_confidence = [value / valid if valid else 0.0 for value in graph_support]
    weighted_graph_confidence = [value / total_weight if total_weight else 0.0 for value in weighted_graph_support]
    environment_confidence = [value / valid if valid else 0.0 for value in environment_support]
    weighted_environment_confidence = [value / total_weight if total_weight else 0.0 for value in weighted_environment_support]

    conf = ref.GetConformer()
    coords = []
    for idx in range(n_atoms):
        pos = conf.GetAtomPosition(idx)
        coords.append((pos.x, pos.y))
    min_x = min(x for x, _ in coords)
    max_x = max(x for x, _ in coords)
    min_y = min(y for _, y in coords)
    max_y = max(y for _, y in coords)
    width = max(max_x - min_x, 1e-6)
    height = max(max_y - min_y, 1e-6)
    canvas_w, canvas_h, margin = 760.0, 450.0, 62.0
    scale = min((canvas_w - 2 * margin) / width, (canvas_h - 2 * margin) / height)

    def xy(atom_idx: int) -> tuple[float, float]:
        x, y = coords[atom_idx]
        return margin + (x - min_x) * scale, canvas_h - (margin + (y - min_y) * scale)

    bond_parts = []
    for bond in ref.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        x1, y1 = xy(a)
        x2, y2 = xy(b)
        bond_parts.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            'stroke="#24445F" stroke-opacity=".42" stroke-width="4" stroke-linecap="round" />'
        )

    rows: list[dict[str, Any]] = []
    atom_parts = []
    for atom in ref.GetAtoms():
        idx = atom.GetIdx()
        x_pos, y_pos = xy(idx)
        graph_value = weighted_graph_confidence[idx]
        environment_value = weighted_environment_confidence[idx]
        gap = max(graph_value - environment_value, 0.0)
        color = heatmap_color(environment_value)
        symbol = html.escape(atom.GetSymbol())
        atom_parts.append(
            f'<g class="consensus-atom environment-consensus-atom" data-confidence="{environment_value:.3f}" '
            f'data-graph-confidence="{graph_value:.3f}" data-environment-gap="{gap:.3f}" transform="translate({x_pos:.1f},{y_pos:.1f})">'
            f'<circle r="24" fill="{color}" stroke="white" stroke-width="3">'
            f"<title>Atom {idx + 1} {symbol}: environment consensus {environment_value:.2f}; graph consensus {graph_value:.2f}</title></circle>"
            f'<text text-anchor="middle" y="5">{symbol}</text>'
            f'<text class="atom-confidence" text-anchor="middle" y="39">{environment_value:.2f}</text></g>'
        )
        rows.append(
            {
                "atom_index": idx + 1,
                "element": atom.GetSymbol(),
                "graph_consensus_support_fraction": round(graph_confidence[idx], 6),
                "weighted_graph_consensus_support_fraction": round(graph_value, 6),
                "environment_consensus_support_fraction": round(environment_confidence[idx], 6),
                "weighted_environment_consensus_support_fraction": round(environment_value, 6),
                "environment_gap_fraction": round(gap, 6),
                "support_count": round(graph_support[idx], 6),
                "environment_support_count": round(environment_support[idx], 6),
                "weighted_support": round(weighted_graph_support[idx], 6),
                "weighted_environment_support": round(weighted_environment_support[idx], 6),
                "valid_candidate_count": valid,
                "total_candidate_weight": round(total_weight, 6),
                "heatmap_color": color,
            }
        )

    svg = "\n".join(
        [
            "<svg class='consensus-svg environment-consensus-svg' viewBox='0 0 760 450' role='img' aria-label='Chemical environment consensus heatmap'>",
            *bond_parts,
            *atom_parts,
            "</svg>",
        ]
    )
    return rows, svg


def mcs_reference_match(ref: Any, cand: Any) -> set[int]:
    try:
        mcs = rdFMCS.FindMCS([ref, cand], timeout=2, matchValences=False, ringMatchesRingOnly=True)
    except Exception:
        return set()
    if not mcs.smartsString:
        return set()
    patt = Chem.MolFromSmarts(mcs.smartsString)
    if patt is None:
        return set()
    return {int(atom_idx) for atom_idx in ref.GetSubstructMatch(patt)}


def reference_ring_atoms_for_side(ref: Any, side: str) -> set[int]:
    patt = Chem.MolFromSmarts("c1ccccc1")
    if patt is None:
        return set()
    matches = [set(match) for match in ref.GetSubstructMatches(patt)]
    for match in matches:
        for atom_idx in match:
            atom = ref.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in match or neighbor.GetSymbol() != "C":
                    continue
                neighbor_neighbors = list(neighbor.GetNeighbors())
                has_double_o = any(
                    bond.GetBondTypeAsDouble() == 2.0 and bond.GetOtherAtom(neighbor).GetSymbol() == "O"
                    for bond in neighbor.GetBonds()
                )
                has_n = any(n.GetSymbol() == "N" for n in neighbor_neighbors)
                has_methyl = any(
                    n.GetSymbol() == "C"
                    and n.GetIdx() != atom_idx
                    and not n.GetIsAromatic()
                    and n.GetTotalNumHs() >= 2
                    for n in neighbor_neighbors
                )
                if side == "left" and has_double_o and has_methyl:
                    return match
                if side == "right" and has_double_o and has_n:
                    return match
    return matches[0] if matches else set()


def reference_atoms_for_motif(ref: Any, motif: dict[str, str]) -> set[int]:
    motif_id = motif["motif_id"]
    if motif_id == "left_aryl_ring":
        return reference_ring_atoms_for_side(ref, "left")
    if motif_id == "right_phenyl":
        return reference_ring_atoms_for_side(ref, "right")
    patt = Chem.MolFromSmarts(motif["motif_smarts"])
    if patt is None:
        return set()
    matches = ref.GetSubstructMatches(patt)
    if not matches:
        return set()
    return {int(atom_idx) for atom_idx in matches[0]}


def motif_present(candidate_mol: Any, motif: dict[str, str]) -> bool:
    patt = Chem.MolFromSmarts(motif["motif_smarts"])
    if patt is not None and candidate_mol.HasSubstructMatch(patt):
        return True
    if motif["motif_id"] == "acetyl_enol_substituent":
        for fallback in ("[CH2]=[C]([OH])", "[CH2]=[C]([O])", "[CX3H1]=[CX3]([OX2H1])"):
            fallback_patt = Chem.MolFromSmarts(fallback)
            if fallback_patt is not None and candidate_mol.HasSubstructMatch(fallback_patt):
                return True
    return False


def motif_interpretation(present_fraction: float, same_position_fraction: float) -> str:
    if present_fraction >= 0.75 and same_position_fraction >= 0.60:
        return "strongly conserved"
    if present_fraction >= 0.60 and same_position_fraction < 0.60:
        return "present but placement varies"
    if present_fraction < 0.40:
        return "weak or absent"
    return "mixed evidence"


def motif_consensus(reference_smiles: str, candidates: pd.DataFrame, atom_support: pd.DataFrame) -> pd.DataFrame:
    del atom_support
    if not HAS_RDKIT or candidates.empty:
        return pd.DataFrame()
    ref = Chem.MolFromSmiles(reference_smiles)
    if ref is None:
        return pd.DataFrame()
    candidate_records: list[dict[str, Any]] = []
    for _, row in candidates.iterrows():
        smiles = clean(row.get("canonical_smiles"))
        if not smiles:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        candidate_records.append({"smiles": smiles, "mol": mol, "mcs_ref_atoms": mcs_reference_match(ref, mol)})
    denominator = max(len(candidate_records), 1)
    rows: list[dict[str, Any]] = []
    for motif in MOTIF_DEFINITIONS:
        ref_atoms = reference_atoms_for_motif(ref, motif)
        present_count = 0
        same_position_count = 0
        example = ""
        for candidate in candidate_records:
            present = motif_present(candidate["mol"], motif)
            if present:
                present_count += 1
                if not example:
                    example = candidate["smiles"]
            same_position = bool(ref_atoms) and len(ref_atoms & candidate["mcs_ref_atoms"]) / max(len(ref_atoms), 1) >= 0.60
            if present and same_position:
                same_position_count += 1
        present_fraction = present_count / denominator
        same_fraction = same_position_count / denominator
        variant_count = max(present_count - same_position_count, 0)
        variant_fraction = variant_count / denominator
        rows.append(
            {
                "motif_id": motif["motif_id"],
                "motif_label": motif["motif_label"],
                "reader_label": motif["reader_label"],
                "motif_present_count": int(present_count),
                "motif_present_fraction": round(float(present_fraction), 6),
                "same_position_count": int(same_position_count),
                "same_position_fraction": round(float(same_fraction), 6),
                "variant_or_repositioned_count": int(variant_count),
                "variant_or_repositioned_fraction": round(float(variant_fraction), 6),
                "interpretation": motif_interpretation(present_fraction, same_fraction),
                "example_candidate_smiles": example,
            }
        )
    return pd.DataFrame(rows)


def is_fused_or_cyclized_pattern(mol: Any) -> bool:
    ring_info = mol.GetRingInfo()
    atom_rings = [set(ring) for ring in ring_info.AtomRings()]
    if len(atom_rings) < 2:
        return False
    for left_index, left in enumerate(atom_rings):
        for right in atom_rings[left_index + 1 :]:
            if left & right:
                return True
    return len(atom_rings) >= 3


def pure_denovo_source_text(row: pd.Series) -> str:
    return "|".join(
        [
            clean(row.get("candidate_family")),
            clean(row.get("source_section")),
            clean(row.get("source_branch")),
            clean(row.get("provenance")),
        ]
    ).lower()


def structural_possibility_from_flags(row: pd.Series, fused_or_cyclized: bool) -> str:
    benzamide = bool(row.get("benzamide_core_present"))
    amide = bool(row.get("amide_linker_present"))
    acetyl_enol = bool(row.get("acetyl_enol_substituent_present"))
    para = bool(row.get("para_placement_present"))
    aryl = bool(row.get("left_aryl_ring_present")) or bool(row.get("right_phenyl_present"))
    pure_denovo = "pure_denovo" in clean(row.get("pure_denovo_source_text"))
    if benzamide and amide and acetyl_enol and para:
        return "benzamide_like_acetyl_para"
    if benzamide and amide and acetyl_enol and not para:
        return "benzamide_like_enol_or_tautomer"
    if benzamide and amide:
        return "benzamide_like_substitution_variant"
    if amide and aryl:
        return "aryl_amide_analog"
    if pure_denovo and fused_or_cyclized and not amide:
        return "fused_or_cyclized_analog"
    return "non_benzamide_outlier"


def candidate_motif_flags(candidates: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, candidate in candidates.iterrows():
        smiles = clean(candidate.get("canonical_smiles"))
        mol = Chem.MolFromSmiles(smiles) if HAS_RDKIT and smiles else None
        motif_flags = {f"{motif['motif_id']}_present": False for motif in MOTIF_DEFINITIONS}
        fused_or_cyclized = False
        if mol is not None:
            motif_flags = {f"{motif['motif_id']}_present": motif_present(mol, motif) for motif in MOTIF_DEFINITIONS}
            fused_or_cyclized = is_fused_or_cyclized_pattern(mol)
        row = {
            "candidate_id": clean(candidate.get("candidate_id")),
            "canonical_smiles": smiles,
            "display_name": clean(candidate.get("display_name")),
            "candidate_family": clean(candidate.get("candidate_family")),
            "rank": clean(candidate.get("rank")),
            "score": clean(candidate.get("score")),
            **motif_flags,
            "pure_denovo_source_text": pure_denovo_source_text(candidate),
        }
        possibility_class = structural_possibility_from_flags(pd.Series(row), fused_or_cyclized)
        row["structural_possibility_class"] = possibility_class
        row["structural_possibility_label"] = POSSIBILITY_BY_ID[possibility_class]["label"]
        rows.append(row)
    columns = [
        "candidate_id",
        "canonical_smiles",
        "display_name",
        "candidate_family",
        "rank",
        "score",
        "benzamide_core_present",
        "left_aryl_ring_present",
        "amide_linker_present",
        "right_phenyl_present",
        "acetyl_enol_substituent_present",
        "para_placement_present",
        "structural_possibility_class",
        "structural_possibility_label",
    ]
    return pd.DataFrame(rows, columns=columns)


def dominant_text(values: pd.Series) -> str:
    cleaned = [clean(value) for value in values if clean(value)]
    if not cleaned:
        return ""
    return pd.Series(cleaned).mode().iloc[0]


def possibility_summary(candidate_flags: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    merged = candidate_flags.merge(
        evidence[
            [
                "candidate_id",
                "ms2_embedding_cosine",
                "gate_status",
                "canonical_smiles",
                "display_name",
            ]
        ],
        on="candidate_id",
        how="left",
        suffixes=("", "_evidence"),
    )
    merged["_score_num"] = pd.to_numeric(merged["score"], errors="coerce").fillna(0.0)
    merged["_rank_num"] = pd.to_numeric(merged["rank"], errors="coerce").fillna(999999.0)
    total = max(len(candidate_flags), 1)
    rows: list[dict[str, Any]] = []
    for possibility in POSSIBILITY_CLASSES:
        group = merged.loc[merged["structural_possibility_class"].eq(possibility["id"])].copy()
        if group.empty:
            top = {}
        else:
            top = group.sort_values(["_score_num", "_rank_num"], ascending=[False, True], kind="mergesort").iloc[0].to_dict()
        rows.append(
            {
                "structural_possibility_class": possibility["id"],
                "structural_possibility_label": possibility["label"],
                "candidate_count": int(len(group)),
                "candidate_fraction": round(float(len(group) / total), 6),
                "top_candidate_id": clean(top.get("candidate_id")),
                "top_candidate_label": clean(top.get("display_name")),
                "top_candidate_smiles": clean(top.get("canonical_smiles")),
                "top_score": clean(top.get("score")),
                "top_ms2_embedding_cosine": clean(top.get("ms2_embedding_cosine")),
                "dominant_gate_status": dominant_text(group.get("gate_status", pd.Series(dtype=str))),
                "reader_interpretation": possibility["reader_interpretation"],
            }
        )
    return pd.DataFrame(rows)


def motif_lookup(motif_df: pd.DataFrame, motif_id: str) -> dict[str, Any]:
    if motif_df.empty or "motif_id" not in motif_df.columns:
        return {}
    rows = motif_df.loc[motif_df["motif_id"].astype(str).eq(motif_id)]
    return rows.iloc[0].to_dict() if not rows.empty else {}


def tree_edge_style(edge_type: str) -> tuple[str, str]:
    if edge_type == "strong_support":
        return "#4F8F61", "solid"
    if edge_type == "mixed_variant":
        return "#F28E2B", "solid"
    if edge_type == "weak_absent":
        return "#9AA7B2", "dashed"
    return "#6B7C8A", "solid"


def add_tree_edge(rows: list[dict[str, Any]], source: str, target: str, label: str, edge_type: str, support: Any) -> None:
    color, stroke_style = tree_edge_style(edge_type)
    support_value = as_float(support)
    rows.append(
        {
            "source": source,
            "target": target,
            "label": label,
            "edge_type": edge_type,
            "support_fraction": round(float(support_value if support_value is not None else 0.0), 6),
            "color": color,
            "stroke_style": stroke_style,
        }
    )


def build_motif_decision_tree(
    motif_df: pd.DataFrame,
    candidate_flags: pd.DataFrame,
    evidence: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary = possibility_summary(candidate_flags, evidence)
    summary_by_class = {row["structural_possibility_class"]: row for _, row in summary.iterrows()}
    node_rows: list[dict[str, Any]] = [
        {
            "id": "root",
            "label": "X-12798 structural possibilities",
            "node_type": "root",
            "x": 60,
            "y": 330,
            "color": "#1F4E79",
            "motif_id": "",
            "motif_present_fraction": "",
            "same_position_fraction": "",
            "candidate_count": int(len(candidate_flags)),
            "top_candidate_id": "",
            "top_candidate_label": "",
            "top_candidate_smiles": "",
            "interpretation": "Retained predicted-candidate set after same-library 400 filtering",
        }
    ]
    decision_specs = [
        ("benzamide_core_decision", "Benzamide-like scaffold?", "benzamide_core", 260, 140),
        ("amide_linker_decision", "Amide bridge retained?", "amide_linker", 460, 100),
        ("aryl_position_decision", "Aryl/phenyl same placement?", "right_phenyl", 460, 300),
        ("acetyl_enol_decision", "Acetyl/enol-like substituent?", "acetyl_enol_substituent", 660, 80),
        ("para_placement_decision", "Same-ring para placement?", "para_placement", 660, 220),
    ]
    for node_id, label, motif_id, x_pos, y_pos in decision_specs:
        motif = motif_lookup(motif_df, motif_id)
        node_rows.append(
            {
                "id": node_id,
                "label": label,
                "node_type": "motif_decision",
                "x": x_pos,
                "y": y_pos,
                "color": "#FFFFFF",
                "motif_id": motif_id,
                "motif_present_fraction": clean(motif.get("motif_present_fraction")),
                "same_position_fraction": clean(motif.get("same_position_fraction")),
                "candidate_count": clean(motif.get("motif_present_count")),
                "top_candidate_id": "",
                "top_candidate_label": "",
                "top_candidate_smiles": "",
                "interpretation": clean(motif.get("interpretation")),
            }
        )
    leaf_y = {
        "benzamide_like_acetyl_para": 35,
        "benzamide_like_enol_or_tautomer": 135,
        "benzamide_like_substitution_variant": 235,
        "aryl_amide_analog": 335,
        "fused_or_cyclized_analog": 435,
        "non_benzamide_outlier": 535,
    }
    for possibility in POSSIBILITY_CLASSES:
        summary_row = summary_by_class.get(possibility["id"], {})
        node_rows.append(
            {
                "id": possibility["id"],
                "label": possibility["label"],
                "node_type": "possibility_leaf",
                "x": 900,
                "y": leaf_y[possibility["id"]],
                "color": possibility["color"],
                "motif_id": "",
                "motif_present_fraction": "",
                "same_position_fraction": "",
                "candidate_count": clean(summary_row.get("candidate_count"), "0"),
                "top_candidate_id": clean(summary_row.get("top_candidate_id")),
                "top_candidate_label": clean(summary_row.get("top_candidate_label")),
                "top_candidate_smiles": clean(summary_row.get("top_candidate_smiles")),
                "interpretation": possibility["reader_interpretation"],
            }
        )

    benzamide = motif_lookup(motif_df, "benzamide_core")
    amide = motif_lookup(motif_df, "amide_linker")
    right_phenyl = motif_lookup(motif_df, "right_phenyl")
    acetyl = motif_lookup(motif_df, "acetyl_enol_substituent")
    para = motif_lookup(motif_df, "para_placement")
    edge_rows: list[dict[str, Any]] = []
    add_tree_edge(edge_rows, "root", "benzamide_core_decision", "start", "neutral", 1.0)
    add_tree_edge(edge_rows, "benzamide_core_decision", "amide_linker_decision", "present", "weak_absent", benzamide.get("motif_present_fraction"))
    add_tree_edge(edge_rows, "benzamide_core_decision", "aryl_position_decision", "weak/absent", "weak_absent", 1.0 - (as_float(benzamide.get("motif_present_fraction")) or 0.0))
    add_tree_edge(edge_rows, "amide_linker_decision", "acetyl_enol_decision", "retained", "mixed_variant", amide.get("motif_present_fraction"))
    add_tree_edge(edge_rows, "amide_linker_decision", "fused_or_cyclized_analog", "lost/changed", "weak_absent", 1.0 - (as_float(amide.get("motif_present_fraction")) or 0.0))
    add_tree_edge(edge_rows, "aryl_position_decision", "aryl_amide_analog", "placement varies", "mixed_variant", right_phenyl.get("variant_or_repositioned_fraction"))
    add_tree_edge(edge_rows, "aryl_position_decision", "non_benzamide_outlier", "weak/absent", "weak_absent", 1.0 - (as_float(right_phenyl.get("motif_present_fraction")) or 0.0))
    add_tree_edge(edge_rows, "acetyl_enol_decision", "para_placement_decision", "present", "weak_absent", acetyl.get("motif_present_fraction"))
    add_tree_edge(edge_rows, "acetyl_enol_decision", "benzamide_like_substitution_variant", "weak/absent", "weak_absent", 1.0 - (as_float(acetyl.get("motif_present_fraction")) or 0.0))
    add_tree_edge(edge_rows, "para_placement_decision", "benzamide_like_acetyl_para", "same placement", "weak_absent", para.get("same_position_fraction"))
    add_tree_edge(edge_rows, "para_placement_decision", "benzamide_like_enol_or_tautomer", "present, placement varies", "mixed_variant", para.get("variant_or_repositioned_fraction"))

    return (
        pd.DataFrame(node_rows, columns=TREE_NODE_COLUMNS),
        pd.DataFrame(edge_rows, columns=TREE_EDGE_COLUMNS),
        summary,
    )


def build_network(evidence: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    valid = evidence.loc[evidence["canonical_smiles"].astype(str).ne("")].copy().reset_index(drop=True)
    fps: dict[str, Any] = {}
    for _, row in valid.iterrows():
        _, fp = mol_and_fp(row["canonical_smiles"])
        fps[row["candidate_id"]] = fp

    node_rows: list[dict[str, Any]] = []
    n = len(valid)
    ring_index = 0
    lib400_winners = set(
        valid.loc[
            valid["provenance"].astype(str).str.contains("encoder_branch_comparison:lib400_only", regex=False),
            "candidate_id",
        ]
    )
    merged_winners = set(
        valid.loc[
            valid["provenance"].astype(str).str.contains("encoder_branch_comparison:merged_libcond", regex=False),
            "candidate_id",
        ]
    )
    branch_winners = lib400_winners | merged_winners
    provisional_ids = set(valid.loc[valid["candidate_family"].eq("provisional_non_pure"), "candidate_id"])
    non_center_count = max(n - len(provisional_ids) - len(branch_winners), 1)

    for _, row in valid.iterrows():
        candidate_id = row["candidate_id"]
        family = row["candidate_family"]
        if candidate_id in provisional_ids:
            x_pos, y_pos = 500.0, 320.0
            radius = 18.0
        elif candidate_id in branch_winners:
            angle = math.radians(210 if candidate_id in lib400_winners else 330)
            x_pos = 500.0 + 135.0 * math.cos(angle)
            y_pos = 320.0 + 135.0 * math.sin(angle)
            radius = 16.0
        else:
            angle = (2.0 * math.pi * ring_index / non_center_count) - math.pi / 2.0
            ring_index += 1
            x_pos = 500.0 + 270.0 * math.cos(angle)
            y_pos = 320.0 + 235.0 * math.sin(angle)
            score = as_float(row.get("score")) or 0.0
            radius = 10.0 + min(max(score, 0.0), 1.0) * 8.0
        label = row["display_name"] or candidate_id
        if len(label) > 34:
            label = label[:31] + "..."
        node_rows.append(
            {
                "id": candidate_id,
                "label": label,
                "canonical_smiles": row["canonical_smiles"],
                "candidate_family": family,
                "source_branch": row.get("source_branch", ""),
                "rank": row.get("rank", ""),
                "score": row.get("score", ""),
                "ms2_embedding_cosine": row.get("ms2_embedding_cosine", ""),
                "predicted_subpathway": row.get("predicted_subpathway", ""),
                "gate_status": row.get("gate_status", ""),
                "x": round(x_pos, 2),
                "y": round(y_pos, 2),
                "radius": round(radius, 2),
                "color": FAMILY_COLORS.get(family, "#4E79A7"),
            }
        )

    matrix_rows: list[dict[str, Any]] = []
    edge_rows: list[dict[str, Any]] = []
    for left_index in range(n):
        left = valid.iloc[left_index]
        for right_index in range(left_index + 1, n):
            right = valid.iloc[right_index]
            left_fp = fps.get(left["candidate_id"])
            right_fp = fps.get(right["candidate_id"])
            tanimoto = float(TanimotoSimilarity(left_fp, right_fp)) if left_fp is not None and right_fp is not None else 0.0
            matrix_rows.append(
                {
                    "source_candidate_id": left["candidate_id"],
                    "source_smiles": left["canonical_smiles"],
                    "target_candidate_id": right["candidate_id"],
                    "target_smiles": right["canonical_smiles"],
                    "tanimoto": round(tanimoto, 6),
                    "source_family": left["candidate_family"],
                    "target_family": right["candidate_family"],
                    "source_label": left["display_name"],
                    "target_label": right["display_name"],
                }
            )
            edge_rows.append(
                {
                    "source": left["candidate_id"],
                    "target": right["candidate_id"],
                    "similarity": round(tanimoto, 6),
                    "similarity_bucket": similarity_bucket(tanimoto),
                    "label": f"{tanimoto:.2f}",
                }
            )

    nodes = pd.DataFrame(node_rows, columns=NODE_COLUMNS)
    edges = pd.DataFrame(edge_rows, columns=EDGE_COLUMNS)
    matrix = pd.DataFrame(matrix_rows)
    network_json = {
        "nodes": safe_json_records(nodes),
        "edges": safe_json_records(edges),
        "rendered_edge_similarity_floor": 0.25,
        "layout": "deterministic_radial_no_reference_anchor",
    }
    return nodes, edges, matrix, network_json


def render_structure_svg(smiles: str, path: Path, legend: str = "") -> bool:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not HAS_RDKIT:
        return False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    try:
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(520, 360)
        options = drawer.drawOptions()
        options.legendFontSize = 18
        options.padding = 0.08
        drawer.DrawMolecule(mol, legend=legend[:80])
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("svg:", "")
        path.write_text(svg)
        return True
    except Exception:
        return False


def render_structure_png(smiles: str, path: Path, legend: str = "") -> bool:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not HAS_RDKIT:
        return False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    try:
        Draw.MolToFile(mol, str(path), size=(520, 360), legend=legend[:80])
        return path.exists()
    except Exception:
        return False


def render_assets(output_dir: Path, evidence: pd.DataFrame, nodes: pd.DataFrame, args: argparse.Namespace) -> dict[str, Any]:
    assets_dir = output_dir / "assets"
    assets_dir.mkdir(parents=True, exist_ok=True)
    rendered: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []

    def render_one(smiles: str, asset_name: str, legend: str) -> None:
        svg_path = assets_dir / f"{asset_name}.svg"
        if render_structure_svg(smiles, svg_path, legend):
            rendered.append({"asset": str(svg_path.relative_to(output_dir)), "smiles": smiles, "format": "svg"})
            return
        png_path = assets_dir / f"{asset_name}.png"
        if render_structure_png(smiles, png_path, legend):
            rendered.append({"asset": str(png_path.relative_to(output_dir)), "smiles": smiles, "format": "png"})
            return
        failures.append({"asset_name": asset_name, "smiles": smiles, "reason": "rdkit_render_failed"})

    for _, row in nodes.iterrows():
        render_one(row["canonical_smiles"], f"candidate-{row['id']}", row["label"])

    provisional = evidence.loc[evidence["candidate_family"].eq("provisional_non_pure")]
    if not provisional.empty:
        row = provisional.iloc[0]
        render_one(row["canonical_smiles"], "x12798_provisional_non_pure", row["display_name"])

    for final_path, asset_name, default_label in [
        (args.lib400_final_json, "x12798_best_pure_denovo_lib400_only", "Best pure de novo lib400-only"),
        (args.merged_final_json, "x12798_best_pure_denovo_merged_libcond", "Best pure de novo merged"),
    ]:
        final = read_json(final_path)
        smiles = clean(final.get("final_smiles"))
        if not smiles:
            continue
        render_one(smiles, asset_name, default_label)

    return {"rendered_assets": rendered, "render_failures": failures}


def render_known_similarity_assets(
    output_dir: Path,
    neighbors: pd.DataFrame,
    known_nodes: pd.DataFrame,
    *,
    max_neighbor_assets: int = 12,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, Any]], list[dict[str, Any]]]:
    assets_dir = output_dir / "assets"
    assets_dir.mkdir(parents=True, exist_ok=True)
    rendered: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []

    anchor_asset = assets_dir / "x12798_anchor_best_candidate.svg"
    if render_structure_svg(ANCHOR_CANDIDATE_SMILES, anchor_asset, ANCHOR_CANDIDATE_NAME):
        rendered.append(
            {
                "asset": str(anchor_asset.relative_to(output_dir)),
                "smiles": canonicalize_smiles(ANCHOR_CANDIDATE_SMILES),
                "format": "svg",
                "kind": "anchor_best_candidate",
            }
        )
    else:
        failures.append(
            {
                "asset": str(anchor_asset.relative_to(output_dir)),
                "smiles": canonicalize_smiles(ANCHOR_CANDIDATE_SMILES),
                "reason": "rdkit_render_failed",
                "kind": "anchor_best_candidate",
            }
        )

    neighbors = neighbors.copy()
    known_nodes = known_nodes.copy()
    for index, row in neighbors.head(max_neighbor_assets).iterrows():
        asset_name = f"known-neighbor-{clean(row.get('chemical_id'))}-{slugify(row.get('canonical_name'))}.svg"
        asset_path = assets_dir / asset_name
        smiles = clean(row.get("canonical_smiles"))
        if render_structure_svg(smiles, asset_path, clean(row.get("canonical_name"))):
            rel_path = str(asset_path.relative_to(output_dir))
            neighbors.at[index, "asset_path"] = rel_path
            if "chemical_id" in known_nodes.columns:
                known_nodes.loc[known_nodes["chemical_id"].astype(str).eq(clean(row.get("chemical_id"))), "asset_path"] = rel_path
            rendered.append({"asset": rel_path, "smiles": smiles, "format": "svg", "kind": "known_neighbor"})
        else:
            failures.append(
                {
                    "asset": str(asset_path.relative_to(output_dir)),
                    "smiles": smiles,
                    "reason": "rdkit_render_failed",
                    "kind": "known_neighbor",
                }
            )

    if "id" in known_nodes.columns:
        known_nodes.loc[known_nodes["id"].eq("anchor_best_candidate"), "asset_path"] = "assets/x12798_anchor_best_candidate.svg"
    return neighbors, known_nodes, rendered, failures


def table_html(df: pd.DataFrame, columns: list[str], limit: int | None = None) -> str:
    view = df.copy()
    if limit is not None:
        view = view.head(limit)
    present = [column for column in columns if column in view.columns]
    if not present or view.empty:
        return "<p class='note'>No rows available.</p>"
    header = "".join(f"<th>{html.escape(column)}</th>" for column in present)
    body_rows: list[str] = []
    for _, row in view.iterrows():
        cells = "".join(f"<td>{html.escape(clean(row.get(column)))}</td>" for column in present)
        body_rows.append(f"<tr>{cells}</tr>")
    return f"<div class='table-wrap'><table><thead><tr>{header}</tr></thead><tbody>{''.join(body_rows)}</tbody></table></div>"


def metric_card(label: str, value: Any, note: str = "") -> str:
    return (
        "<div class='metric-card'>"
        f"<div class='metric-label'>{html.escape(label)}</div>"
        f"<div class='metric-value'>{html.escape(clean(value))}</div>"
        f"<div class='metric-note'>{html.escape(note)}</div>"
        "</div>"
    )


def known_neighbor_cards_html(neighbors: pd.DataFrame, output_dir: Path, limit: int = 4) -> str:
    if neighbors.empty:
        return "<p class='note'>No known-neighbor structures available.</p>"
    cards: list[str] = []
    for _, row in neighbors.head(limit).iterrows():
        asset = clean(row.get("asset_path"))
        image = ""
        if asset and (output_dir / asset).exists():
            image = f"<img src='{html.escape(asset)}' alt='{html.escape(clean(row.get('canonical_name')))} structure' />"
        cards.append(
            "<div class='known-card'>"
            f"{image}"
            f"<h3>{html.escape(clean(row.get('canonical_name')))}</h3>"
            f"<p><strong>Tanimoto {as_float(row.get('similarity_to_anchor')) or 0.0:.3f}</strong></p>"
            f"<p class='note'>{html.escape(clean(row.get('super_pathway')))} / {html.escape(clean(row.get('pathway')))}</p>"
            "</div>"
        )
    return f"<div class='known-card-grid'>{''.join(cards)}</div>"


def motif_badge_class(interpretation: str) -> str:
    if interpretation == "strongly conserved":
        return "strong"
    if interpretation in {"present but placement varies", "mixed evidence"}:
        return "mixed"
    return "weak"


def motif_summary_html(motif_df: pd.DataFrame) -> str:
    if motif_df.empty:
        return "<p class='note'>Motif consensus unavailable.</p>"
    cards: list[str] = []
    for _, row in motif_df.iterrows():
        interpretation = clean(row.get("interpretation"))
        present_pct = 100.0 * (as_float(row.get("motif_present_fraction")) or 0.0)
        same_pct = 100.0 * (as_float(row.get("same_position_fraction")) or 0.0)
        cards.append(
            "<div class='motif-card'>"
            f"<div class='motif-title'>{html.escape(clean(row.get('reader_label')))}</div>"
            f"<div class='motif-badge {motif_badge_class(interpretation)}'>{html.escape(interpretation)}</div>"
            "<div class='motif-measures'>"
            f"<span><strong>{present_pct:.0f}%</strong><small>Present in candidates</small></span>"
            f"<span><strong>{same_pct:.0f}%</strong><small>Same placement on scaffold</small></span>"
            "</div>"
            f"<p>{html.escape(motif_plain_language(row))}</p>"
            "</div>"
        )
    return "<div class='motif-grid'>" + "".join(cards) + "</div>"


def motif_plain_language(row: pd.Series) -> str:
    interpretation = clean(row.get("interpretation"))
    label = clean(row.get("reader_label"))
    if interpretation == "strongly conserved":
        return f"{label} appears in most candidates and usually maps to the same part of the provisional scaffold."
    if interpretation == "present but placement varies":
        return f"{label} appears in many candidates, but the candidates do not agree on the same placement in the displayed scaffold."
    if interpretation == "weak or absent":
        return f"{label} is uncommon in the retained candidates; low atom-position color here can reflect motif absence."
    return f"{label} has mixed support; compare presence and same-placement values before interpreting the atom heatmap."


def network_svg(nodes: pd.DataFrame, edges: pd.DataFrame) -> str:
    node_by_id = {row["id"]: row for _, row in nodes.iterrows()}
    edge_parts: list[str] = []
    for _, edge in edges.iterrows():
        source = node_by_id.get(edge["source"])
        target = node_by_id.get(edge["target"])
        if source is None or target is None:
            continue
        similarity = float(edge["similarity"])
        width = 0.8 + similarity * 5.0
        opacity = 0.16 + similarity * 0.55
        visible = "1" if similarity >= 0.25 else "0"
        edge_parts.append(
            "<line class='net-edge' "
            f"x1='{source['x']}' y1='{source['y']}' x2='{target['x']}' y2='{target['y']}' "
            f"stroke='{edge_color(similarity)}' stroke-width='{width:.2f}' opacity='{opacity:.3f}' "
            f"data-similarity='{similarity:.6f}' data-visible-default='{visible}' "
            f"data-source='{html.escape(edge['source'])}' data-target='{html.escape(edge['target'])}' />"
        )
    node_parts: list[str] = []
    for _, node in nodes.iterrows():
        node_parts.append(
            "<g class='net-node' "
            f"data-family='{html.escape(node['candidate_family'])}' data-id='{html.escape(node['id'])}'>"
            f"<circle cx='{node['x']}' cy='{node['y']}' r='{node['radius']}' fill='{node['color']}' />"
            f"<text x='{node['x']}' y='{float(node['y']) + float(node['radius']) + 16:.2f}'>"
            f"{html.escape(node['label'])}</text>"
            f"<title>{html.escape(node['label'])} | {html.escape(node['canonical_smiles'])}</title>"
            "</g>"
        )
    return (
        "<svg class='network-svg' viewBox='0 0 1000 650' role='img' "
        "aria-label='Pairwise predicted-candidate chemical similarity network'>"
        "<rect width='1000' height='650' rx='18' fill='#F7F9FB' />"
        f"<g class='edges'>{''.join(edge_parts)}</g>"
        f"<g class='nodes'>{''.join(node_parts)}</g>"
        "</svg>"
    )


def known_similarity_network_svg(nodes: pd.DataFrame, edges: pd.DataFrame) -> str:
    node_by_id = {row["id"]: row for _, row in nodes.iterrows()}
    edge_parts: list[str] = []
    for _, edge in edges.iterrows():
        source = node_by_id.get(edge["source"])
        target = node_by_id.get(edge["target"])
        if source is None or target is None:
            continue
        similarity = float(edge["similarity"])
        width = 0.8 + similarity * 5.0
        opacity = 0.18 + similarity * 0.62
        visible = "1" if similarity >= 0.35 else "0"
        edge_parts.append(
            "<line class='known-net-edge' "
            f"x1='{source['x']}' y1='{source['y']}' x2='{target['x']}' y2='{target['y']}' "
            f"stroke='{known_color(similarity)}' stroke-width='{width:.2f}' opacity='{opacity:.3f}' "
            f"data-known-similarity='{similarity:.6f}' data-visible-default='{visible}' "
            f"data-source='{html.escape(edge['source'])}' data-target='{html.escape(edge['target'])}'>"
            f"<title>Tanimoto {similarity:.3f}</title></line>"
        )
    node_parts: list[str] = []
    for _, node in nodes.iterrows():
        label = clean(node.get("label"))
        if len(label) > 32:
            label = label[:29] + "..."
        similarity = as_float(node.get("similarity_to_anchor"))
        subtitle = "anchor" if node["node_type"] == "anchor_best_candidate" else f"{similarity or 0.0:.2f}"
        node_parts.append(
            "<g class='known-net-node' "
            f"data-node-type='{html.escape(node['node_type'])}' data-bucket='{html.escape(node['similarity_bucket'])}' "
            f"data-id='{html.escape(node['id'])}'>"
            f"<circle cx='{node['x']}' cy='{node['y']}' r='{node['radius']}' fill='{node['color']}' />"
            f"<text x='{node['x']}' y='{float(node['y']) + float(node['radius']) + 16:.2f}'>{html.escape(label)}</text>"
            f"<text class='known-node-score' x='{node['x']}' y='{float(node['y']) + float(node['radius']) + 31:.2f}'>{html.escape(subtitle)}</text>"
            f"<title>{html.escape(clean(node['canonical_name']))} | {html.escape(clean(node['canonical_smiles']))}</title>"
            "</g>"
        )
    return (
        "<svg class='known-network-svg' viewBox='0 0 1000 650' role='img' "
        "aria-label='Best-candidate similarity to known Metabolon metabolites'>"
        "<rect width='1000' height='650' rx='18' fill='#F7F9FB' />"
        f"<g class='edges'>{''.join(edge_parts)}</g>"
        f"<g class='nodes'>{''.join(node_parts)}</g>"
        "</svg>"
    )


def svg_text_lines(text: Any, x_pos: float, y_pos: float, max_chars: int, line_height: int = 15) -> str:
    words = html.escape(clean(text)).split()
    if not words:
        return ""
    lines: list[str] = []
    current = ""
    for word in words:
        trial = f"{current} {word}".strip()
        if len(trial) > max_chars and current:
            lines.append(current)
            current = word
        else:
            current = trial
    if current:
        lines.append(current)
    parts = []
    for index, line in enumerate(lines[:3]):
        dy = 0 if index == 0 else line_height
        parts.append(f"<tspan x='{x_pos:.1f}' dy='{dy}'>{line}</tspan>")
    return "".join(parts)


def decision_tree_svg(tree_nodes: pd.DataFrame, tree_edges: pd.DataFrame) -> str:
    node_by_id = {row["id"]: row for _, row in tree_nodes.iterrows()}

    def box_size(node_type: str) -> tuple[float, float]:
        if node_type == "root":
            return 175.0, 64.0
        if node_type == "motif_decision":
            return 188.0, 88.0
        return 225.0, 82.0

    edge_parts: list[str] = []
    for _, edge in tree_edges.iterrows():
        source = node_by_id.get(edge["source"])
        target = node_by_id.get(edge["target"])
        if source is None or target is None:
            continue
        source_w, source_h = box_size(source["node_type"])
        target_w, target_h = box_size(target["node_type"])
        x1 = float(source["x"]) + source_w
        y1 = float(source["y"]) + source_h / 2.0
        x2 = float(target["x"])
        y2 = float(target["y"]) + target_h / 2.0
        mid = (x1 + x2) / 2.0
        dash = " stroke-dasharray='7 6'" if edge["stroke_style"] == "dashed" else ""
        edge_parts.append(
            f"<path d='M{x1:.1f},{y1:.1f} C{mid:.1f},{y1:.1f} {mid:.1f},{y2:.1f} {x2:.1f},{y2:.1f}' "
            f"fill='none' stroke='{html.escape(edge['color'])}' stroke-width='3'{dash} marker-end='url(#tree-arrow)' />"
        )
        label_x = (mid + x2) / 2.0
        label_y = (y1 + y2) / 2.0 - 6.0
        edge_parts.append(
            f"<text class='tree-edge-label' x='{label_x:.1f}' y='{label_y:.1f}'>{html.escape(clean(edge['label']))}</text>"
        )

    node_parts: list[str] = []
    for _, node in tree_nodes.iterrows():
        node_type = clean(node["node_type"])
        width, height = box_size(node_type)
        x_pos = float(node["x"])
        y_pos = float(node["y"])
        fill = html.escape(clean(node["color"]))
        stroke = "#24445F" if node_type == "root" else "#C8D3DC"
        text_color = "#FFFFFF" if node_type in {"root", "possibility_leaf"} else "#17212B"
        rx = 12 if node_type != "motif_decision" else 8
        node_parts.append(
            f"<g class='tree-node tree-node-{html.escape(node_type)}' data-id='{html.escape(clean(node['id']))}'>"
            f"<rect x='{x_pos:.1f}' y='{y_pos:.1f}' width='{width:.1f}' height='{height:.1f}' rx='{rx}' "
            f"fill='{fill}' stroke='{stroke}' stroke-width='2' />"
            f"<text class='tree-node-title' fill='{text_color}' x='{x_pos + width / 2.0:.1f}' y='{y_pos + 20:.1f}'>"
            f"{svg_text_lines(node['label'], x_pos + width / 2.0, y_pos + 20.0, 25)}</text>"
        )
        if node_type == "motif_decision":
            present = as_float(node.get("motif_present_fraction")) or 0.0
            same = as_float(node.get("same_position_fraction")) or 0.0
            node_parts.append(
                f"<text class='tree-node-metric' x='{x_pos + width / 2.0:.1f}' y='{y_pos + height - 24:.1f}'>"
                f"present {present:.0%} / same {same:.0%}</text>"
                f"<text class='tree-node-note' x='{x_pos + width / 2.0:.1f}' y='{y_pos + height - 9:.1f}'>"
                f"{html.escape(clean(node.get('interpretation'))[:30])}</text>"
            )
        elif node_type == "possibility_leaf":
            count = clean(node.get("candidate_count"), "0")
            top_label = clean(node.get("top_candidate_label")) or "No retained candidate"
            if len(top_label) > 32:
                top_label = top_label[:29] + "..."
            node_parts.append(
                f"<text class='tree-leaf-count' x='{x_pos + width / 2.0:.1f}' y='{y_pos + height - 28:.1f}'>n={html.escape(count)}</text>"
                f"<text class='tree-leaf-top' x='{x_pos + width / 2.0:.1f}' y='{y_pos + height - 11:.1f}'>"
                f"{html.escape(top_label)}</text>"
            )
        node_parts.append(
            f"<title>{html.escape(clean(node['label']))}: {html.escape(clean(node.get('interpretation')))}</title></g>"
        )

    return (
        "<svg class='decision-tree-svg' viewBox='0 0 1160 700' role='img' "
        "aria-label='X-12798 structural possibility decision tree'>"
        "<defs><marker id='tree-arrow' viewBox='0 0 10 10' refX='9' refY='5' markerWidth='6' markerHeight='6' orient='auto-start-reverse'>"
        "<path d='M 0 0 L 10 5 L 0 10 z' fill='#6B7C8A' /></marker></defs>"
        "<rect width='1160' height='700' rx='18' fill='#F7F9FB' />"
        "<text class='tree-title' x='60' y='36'>Structural possibility branches from motif consensus</text>"
        f"<g class='tree-edges'>{''.join(edge_parts)}</g>"
        f"<g class='tree-nodes'>{''.join(node_parts)}</g>"
        "</svg>"
    )


def build_candidate_specific_consensus(evidence: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    assets_dir = output_dir / "assets"
    atom_dir = output_dir / "data/candidate_specific_consensus_atom_support"
    environment_atom_dir = output_dir / "data/candidate_environment_consensus_atom_support"
    assets_dir.mkdir(parents=True, exist_ok=True)
    atom_dir.mkdir(parents=True, exist_ok=True)
    environment_atom_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for _, candidate in evidence.iterrows():
        candidate_id = clean(candidate.get("candidate_id"))
        smiles = clean(candidate.get("canonical_smiles"))
        if not candidate_id or not smiles:
            continue
        comparison_pool = evidence.loc[~evidence["candidate_id"].astype(str).eq(candidate_id)].copy()
        atom_rows, svg = consensus_atom_support(smiles, comparison_pool)
        environment_rows, environment_svg = environment_aware_consensus_atom_support(smiles, comparison_pool)
        atom_df = pd.DataFrame(atom_rows)
        environment_df = pd.DataFrame(environment_rows)
        atom_csv_rel = f"data/candidate_specific_consensus_atom_support/{candidate_id}.csv"
        environment_atom_csv_rel = f"data/candidate_environment_consensus_atom_support/{candidate_id}.csv"
        asset_rel = f"assets/candidate_consensus_{candidate_id}.svg"
        environment_asset_rel = f"assets/candidate_environment_consensus_{candidate_id}.svg"
        atom_df.to_csv(output_dir / atom_csv_rel, index=False)
        environment_df.to_csv(output_dir / environment_atom_csv_rel, index=False)
        (output_dir / asset_rel).write_text(svg)
        (output_dir / environment_asset_rel).write_text(environment_svg)
        if atom_df.empty:
            mean_weighted = median_weighted = max_weighted = atoms_ge_half = 0.0
            atom_count = valid_count = 0
        else:
            weighted = pd.to_numeric(atom_df["weighted_consensus_support_fraction"], errors="coerce").fillna(0.0)
            mean_weighted = float(weighted.mean())
            median_weighted = float(weighted.median())
            max_weighted = float(weighted.max())
            atoms_ge_half = float((weighted >= 0.5).mean())
            atom_count = int(len(atom_df))
            valid_counts = pd.to_numeric(atom_df["valid_candidate_count"], errors="coerce").dropna()
            valid_count = int(valid_counts.max()) if not valid_counts.empty else 0
        if environment_df.empty:
            env_mean = env_median = env_max = env_atoms_ge_half = env_gap_mean = 0.0
        else:
            env_weighted = pd.to_numeric(
                environment_df["weighted_environment_consensus_support_fraction"], errors="coerce"
            ).fillna(0.0)
            env_gap = pd.to_numeric(environment_df["environment_gap_fraction"], errors="coerce").fillna(0.0)
            env_mean = float(env_weighted.mean())
            env_median = float(env_weighted.median())
            env_max = float(env_weighted.max())
            env_atoms_ge_half = float((env_weighted >= 0.5).mean())
            env_gap_mean = float(env_gap.mean())
        rows.append(
            {
                "candidate_id": candidate_id,
                "display_name": clean(candidate.get("display_name")),
                "canonical_smiles": smiles,
                "candidate_family": clean(candidate.get("candidate_family")),
                "rank": clean(candidate.get("rank")),
                "score": clean(candidate.get("score")),
                "candidate_consensus_asset": asset_rel,
                "candidate_consensus_atom_support_csv": atom_csv_rel,
                "candidate_consensus_atom_count": atom_count,
                "candidate_consensus_valid_comparison_count": valid_count,
                "candidate_consensus_mean_weighted_fraction": round(mean_weighted, 6),
                "candidate_consensus_median_weighted_fraction": round(median_weighted, 6),
                "candidate_consensus_max_weighted_fraction": round(max_weighted, 6),
                "candidate_consensus_atoms_ge_0p5_fraction": round(atoms_ge_half, 6),
                "candidate_environment_consensus_asset": environment_asset_rel,
                "candidate_environment_consensus_atom_support_csv": environment_atom_csv_rel,
                "candidate_environment_mean_weighted_fraction": round(env_mean, 6),
                "candidate_environment_median_weighted_fraction": round(env_median, 6),
                "candidate_environment_max_weighted_fraction": round(env_max, 6),
                "candidate_environment_atoms_ge_0p5_fraction": round(env_atoms_ge_half, 6),
                "candidate_environment_gap_mean_fraction": round(env_gap_mean, 6),
            }
        )
    return pd.DataFrame(rows)


def candidate_consensus_gallery_html(candidate_consensus: pd.DataFrame, output_dir: Path) -> str:
    if candidate_consensus.empty:
        return "<p class='note'>Candidate-specific consensus images are unavailable.</p>"
    cards: list[str] = []
    for _, row in candidate_consensus.iterrows():
        label = clean(row.get("display_name")) or clean(row.get("candidate_id"))
        if len(label) > 56:
            label = label[:53] + "..."
        mean_pct = 100.0 * (as_float(row.get("candidate_consensus_mean_weighted_fraction")) or 0.0)
        ge_half_pct = 100.0 * (as_float(row.get("candidate_consensus_atoms_ge_0p5_fraction")) or 0.0)
        env_mean_pct = 100.0 * (as_float(row.get("candidate_environment_mean_weighted_fraction")) or 0.0)
        env_ge_half_pct = 100.0 * (as_float(row.get("candidate_environment_atoms_ge_0p5_fraction")) or 0.0)
        asset = clean(row.get("candidate_consensus_asset"))
        asset_path = output_dir / asset
        if asset_path.exists():
            structure_visual = (
                "<div class='candidate-consensus-figure candidate-specific-inline-svg'>"
                "<div class='candidate-consensus-panel-label'>Graph overlap</div>"
                f"{asset_path.read_text()}"
                "</div>"
            )
        else:
            structure_visual = (
                "<p class='note'>Candidate-specific consensus SVG missing: "
                f"{html.escape(asset)}</p>"
            )
        env_asset = clean(row.get("candidate_environment_consensus_asset"))
        env_asset_path = output_dir / env_asset
        if env_asset_path.exists():
            environment_visual = (
                "<div class='candidate-consensus-figure candidate-specific-inline-svg environment-inline-svg'>"
                "<div class='candidate-consensus-panel-label'>Environment-aware</div>"
                f"{env_asset_path.read_text()}"
                "</div>"
            )
        else:
            environment_visual = (
                "<p class='note'>Candidate environment consensus SVG missing: "
                f"{html.escape(env_asset)}</p>"
            )
        cards.append(
            "<div class='candidate-consensus-card'>"
            f"<h3>{html.escape(label)}</h3>"
            "<div class='candidate-consensus-visual-pair'>"
            f"{structure_visual}{environment_visual}"
            "</div>"
            "<div class='candidate-consensus-metrics'>"
            f"<span><strong>{mean_pct:.0f}%</strong><small>mean graph-position support</small></span>"
            f"<span><strong>{env_mean_pct:.0f}%</strong><small>mean environment-aware support</small></span>"
            f"<span><strong>{ge_half_pct:.0f}%</strong><small>graph atoms with support >= 0.5</small></span>"
            f"<span><strong>{env_ge_half_pct:.0f}%</strong><small>environment atoms with support >= 0.5</small></span>"
            "</div>"
            f"<p><code>{html.escape(clean(row.get('canonical_smiles')))}</code></p>"
            "</div>"
        )
    return "<div class='candidate-consensus-grid'>" + "".join(cards) + "</div>"


def build_final_consensus_guided_visuals(consensus_guided_root: Path, output_dir: Path) -> dict[str, Any]:
    top20_path = consensus_guided_root / "consensus_guided_top20.csv"
    final_call_path = consensus_guided_root / "final_consensus_guided_call.json"
    if not top20_path.exists():
        return {
            "available": False,
            "reason": f"missing {top20_path}",
            "summary": pd.DataFrame(),
            "graph_rows": pd.DataFrame(),
            "environment_rows": pd.DataFrame(),
            "graph_svg": "",
            "environment_svg": "",
        }
    top20 = pd.read_csv(top20_path, low_memory=False)
    if top20.empty or "canonical_smiles" not in top20.columns:
        return {
            "available": False,
            "reason": "consensus_guided_top20_empty",
            "summary": pd.DataFrame(),
            "graph_rows": pd.DataFrame(),
            "environment_rows": pd.DataFrame(),
            "graph_svg": "",
            "environment_svg": "",
        }
    top20 = top20.copy()
    top20["rank"] = pd.to_numeric(top20.get("rank"), errors="coerce")
    top20 = top20.sort_values("rank", kind="mergesort").reset_index(drop=True)
    top = top20.iloc[0].to_dict()
    top_smiles = clean(top.get("canonical_smiles"))
    comparison_pool = top20.loc[~top20["canonical_smiles"].astype(str).eq(top_smiles)].copy()
    if "consensus_guided_final_score" in comparison_pool.columns:
        comparison_pool["score"] = pd.to_numeric(comparison_pool["consensus_guided_final_score"], errors="coerce").fillna(0.25)
    else:
        comparison_pool["score"] = 0.25
    graph_rows, graph_svg = consensus_atom_support(top_smiles, comparison_pool)
    environment_rows, environment_svg = environment_aware_consensus_atom_support(top_smiles, comparison_pool)
    graph_df = pd.DataFrame(graph_rows)
    environment_df = pd.DataFrame(environment_rows)
    graph_asset = "assets/x12798_final_consensus_guided_graph_position.svg"
    environment_asset = "assets/x12798_final_consensus_guided_environment_position.svg"
    graph_csv = "data/x12798_final_consensus_guided_graph_position_atom_support.csv"
    environment_csv = "data/x12798_final_consensus_guided_environment_position_atom_support.csv"
    summary_csv = "data/x12798_final_consensus_guided_visual_summary.csv"
    summary_json = "data/x12798_final_consensus_guided_visual_summary.json"
    graph_df.to_csv(output_dir / graph_csv, index=False)
    environment_df.to_csv(output_dir / environment_csv, index=False)
    (output_dir / graph_asset).write_text(graph_svg)
    (output_dir / environment_asset).write_text(environment_svg)
    final_call = read_json(final_call_path) if final_call_path.exists() else {}
    graph_weighted = (
        pd.to_numeric(graph_df.get("weighted_consensus_support_fraction"), errors="coerce").fillna(0.0)
        if not graph_df.empty
        else pd.Series(dtype=float)
    )
    env_weighted = (
        pd.to_numeric(environment_df.get("weighted_environment_consensus_support_fraction"), errors="coerce").fillna(0.0)
        if not environment_df.empty
        else pd.Series(dtype=float)
    )
    gap = (
        pd.to_numeric(environment_df.get("environment_gap_fraction"), errors="coerce").fillna(0.0)
        if not environment_df.empty
        else pd.Series(dtype=float)
    )
    summary = pd.DataFrame(
        [
            {
                "top_rank": clean(top.get("rank")),
                "top_smiles": top_smiles,
                "top_display_name": clean(top.get("display_name")),
                "top_final_score": clean(top.get("consensus_guided_final_score")),
                "top_ms2_embedding_cosine": clean(top.get("ms2_embedding_cosine")),
                "top_final_gate_status": clean(top.get("final_gate_status")),
                "final_status": clean(final_call.get("final_status")),
                "ms2_embedding_gate_pass": str(bool(final_call.get("ms2_embedding_gate_pass"))).lower(),
                "comparison_candidate_count": int(len(comparison_pool)),
                "graph_position_mean_weighted_fraction": round(float(graph_weighted.mean()) if len(graph_weighted) else 0.0, 6),
                "graph_position_atoms_ge_0p5_fraction": round(float((graph_weighted >= 0.5).mean()) if len(graph_weighted) else 0.0, 6),
                "environment_position_mean_weighted_fraction": round(float(env_weighted.mean()) if len(env_weighted) else 0.0, 6),
                "environment_position_atoms_ge_0p5_fraction": round(float((env_weighted >= 0.5).mean()) if len(env_weighted) else 0.0, 6),
                "environment_position_gap_mean_fraction": round(float(gap.mean()) if len(gap) else 0.0, 6),
                "graph_position_asset": graph_asset,
                "environment_position_asset": environment_asset,
                "graph_position_atom_support_csv": graph_csv,
                "environment_position_atom_support_csv": environment_csv,
            }
        ]
    )
    summary.to_csv(output_dir / summary_csv, index=False)
    (output_dir / summary_json).write_text(json.dumps(safe_json_records(summary), indent=2) + "\n")
    return {
        "available": True,
        "reason": "",
        "summary": summary,
        "graph_rows": graph_df,
        "environment_rows": environment_df,
        "graph_svg": graph_svg,
        "environment_svg": environment_svg,
        "graph_asset": graph_asset,
        "environment_asset": environment_asset,
        "graph_csv": graph_csv,
        "environment_csv": environment_csv,
        "summary_csv": summary_csv,
        "summary_json": summary_json,
    }


def final_consensus_guided_html(final_consensus: dict[str, Any]) -> str:
    if not final_consensus.get("available"):
        return (
            "<p class='note'>Consensus-guided final candidate visuals are unavailable: "
            f"{html.escape(clean(final_consensus.get('reason')))}</p>"
        )
    summary = final_consensus["summary"].iloc[0].to_dict()
    graph_mean = 100.0 * (as_float(summary.get("graph_position_mean_weighted_fraction")) or 0.0)
    graph_ge = 100.0 * (as_float(summary.get("graph_position_atoms_ge_0p5_fraction")) or 0.0)
    env_mean = 100.0 * (as_float(summary.get("environment_position_mean_weighted_fraction")) or 0.0)
    env_ge = 100.0 * (as_float(summary.get("environment_position_atoms_ge_0p5_fraction")) or 0.0)
    env_gap = 100.0 * (as_float(summary.get("environment_position_gap_mean_fraction")) or 0.0)
    return f"""
      <div class="call-grid">
        <div class="structure-card">
          <h3>Final Consensus Graph-Position Confidence</h3>
          <div class="consensus-box">{final_consensus['graph_svg']}</div>
          <p class="note">Mean graph-position confidence {graph_mean:.0f}%; atoms >= 0.5 support {graph_ge:.0f}%.</p>
        </div>
        <div class="structure-card">
          <h3>Final Consensus Environment-Position Confidence</h3>
          <div class="consensus-box">{final_consensus['environment_svg']}</div>
          <p class="note">Mean environment-position confidence {env_mean:.0f}%; atoms >= 0.5 support {env_ge:.0f}%; mean graph/environment gap {env_gap:.0f}%.</p>
        </div>
      </div>
      <div class="candidate-consensus-metrics">
        <span><strong>{graph_mean:.0f}%</strong><small>mean graph-position confidence</small></span>
        <span><strong>{env_mean:.0f}%</strong><small>mean environment-position confidence</small></span>
        <span><strong>{graph_ge:.0f}%</strong><small>graph atoms with support >= 0.5</small></span>
        <span><strong>{env_ge:.0f}%</strong><small>environment atoms with support >= 0.5</small></span>
      </div>
      <p><strong>Top consensus-guided SMILES:</strong> <code>{html.escape(clean(summary.get('top_smiles')))}</code></p>
      <p class="note">
        Final consensus status: <code>{html.escape(clean(summary.get('final_status')))}</code>.
        Real MS2 gate pass: <code>{html.escape(clean(summary.get('ms2_embedding_gate_pass')))}</code>.
        These images use the same atom-confidence rendering as the scaffold consensus views, but the displayed molecule is the consensus-guided top candidate.
      </p>
      <p class="note">
        Sidecars: <a href="data/x12798_final_consensus_guided_visual_summary.csv">summary</a>,
        <a href="data/x12798_final_consensus_guided_graph_position_atom_support.csv">graph-position atom support</a>,
        <a href="data/x12798_final_consensus_guided_environment_position_atom_support.csv">environment-position atom support</a>.
      </p>
    """


def build_real_latent_calibration_summary(consensus_guided_root: Path, calibration_root: Path, output_dir: Path) -> dict[str, Any]:
    candidates_path = consensus_guided_root / "consensus_guided_candidates_all.csv"
    calibration_summary_path = calibration_root / "latent_gate_calibration_summary.json"
    reliability_path = calibration_root / "latent_gate_reliability_curve.csv"
    summary_csv = "data/x12798_real_latent_calibration_candidate_summary.csv"
    summary_json = "data/x12798_real_latent_calibration_candidate_summary.json"
    target_smiles = [
        "CC(=O)c1ccc(NC(=O)c2ccccc2)cc1",
        "C=C(O)c1ccc(NC(=O)c2ccccc2)cc1",
    ]
    if not candidates_path.exists():
        return {"available": False, "reason": f"missing {candidates_path}", "summary": pd.DataFrame()}
    candidates = pd.read_csv(candidates_path, low_memory=False)
    if candidates.empty or "canonical_smiles" not in candidates.columns:
        return {"available": False, "reason": "consensus_guided_candidates_empty", "summary": pd.DataFrame()}
    sub = candidates.loc[candidates["canonical_smiles"].astype(str).isin(target_smiles)].copy()
    for col in [
        "rank",
        "canonical_smiles",
        "candidate_family",
        "pure_de_novo_eligible",
        "novelty_status",
        "consensus_guided_final_score",
        "ms2_embedding_match_status",
        "ms2_embedding_cosine",
        "ms2_embedding_rank",
        "latent_calibrated_probability",
        "latent_calibrated_gate_pass",
        "tautomer_cluster_id",
        "tautomer_cluster_representative_smiles",
        "tautomer_cluster_size",
        "tautomer_cluster_best_cosine",
        "tautomer_cluster_best_smiles",
        "tautomer_cluster_evidence_note",
        "final_gate_status",
        "gate_failures",
    ]:
        if col not in sub.columns:
            sub[col] = ""
    sub["candidate_drawn_form"] = sub["canonical_smiles"].map(
        {
            "CC(=O)c1ccc(NC(=O)c2ccccc2)cc1": "keto benzamide best-candidate anchor",
            "C=C(O)c1ccc(NC(=O)c2ccccc2)cc1": "enol/tautomer-like pure de novo candidate",
        }
    )
    sub = sub.sort_values(["rank", "canonical_smiles"], kind="mergesort").reset_index(drop=True)
    out_cols = [
        "candidate_drawn_form",
        "rank",
        "canonical_smiles",
        "candidate_family",
        "pure_de_novo_eligible",
        "novelty_status",
        "consensus_guided_final_score",
        "ms2_embedding_match_status",
        "ms2_embedding_cosine",
        "ms2_embedding_rank",
        "latent_calibrated_probability",
        "latent_calibrated_gate_pass",
        "tautomer_cluster_id",
        "tautomer_cluster_representative_smiles",
        "tautomer_cluster_size",
        "tautomer_cluster_best_cosine",
        "tautomer_cluster_best_smiles",
        "tautomer_cluster_evidence_note",
        "final_gate_status",
        "gate_failures",
    ]
    sub[out_cols].to_csv(output_dir / summary_csv, index=False)
    (output_dir / summary_json).write_text(json.dumps(safe_json_records(sub[out_cols]), indent=2) + "\n")
    calibration_summary = read_json(calibration_summary_path) if calibration_summary_path.exists() else {}
    reliability = pd.read_csv(reliability_path) if reliability_path.exists() else pd.DataFrame()
    return {
        "available": True,
        "reason": "",
        "summary": sub[out_cols],
        "summary_csv": summary_csv,
        "summary_json": summary_json,
        "calibration_summary": calibration_summary,
        "reliability": reliability,
        "calibration_summary_path": str(calibration_summary_path),
        "reliability_path": str(reliability_path),
    }


def real_latent_calibration_html(real_latent: dict[str, Any]) -> str:
    if not real_latent.get("available"):
        return f"<p class='note'>Real-latent calibration summary unavailable: {html.escape(clean(real_latent.get('reason')))}</p>"
    summary_df = real_latent["summary"]
    calibration_summary = real_latent.get("calibration_summary", {})
    n_controls = clean(calibration_summary.get("n_known_controls"))
    n_rows = clean(calibration_summary.get("n_rows"))
    auc = as_float(calibration_summary.get("roc_auc"))
    auc_text = f"{auc:.3f}" if auc is not None else "not available"
    return f"""
      <p>
        Raw Spec2Mol cosine is an encoder-space similarity score, not a direct probability of correctness.
        For X-12798, current raw cosines are low, so structure calls remain provisional unless calibrated known-control evidence supports promotion.
      </p>
      <div class="metrics">
        {metric_card("Known-control rows", n_rows or "0", f"{n_controls or '0'} controls with truth labels")}
        {metric_card("Calibration ROC AUC", auc_text, "Using available V14 latent-like score columns")}
        {metric_card("High-confidence gate", "not passed", "Requires raw cosine >= 0.70 and calibrated support")}
      </div>
      {table_html(summary_df, ['candidate_drawn_form','rank','canonical_smiles','pure_de_novo_eligible','novelty_status','ms2_embedding_match_status','ms2_embedding_cosine','ms2_embedding_rank','latent_calibrated_probability','latent_calibrated_gate_pass','tautomer_cluster_size','tautomer_cluster_best_cosine','tautomer_cluster_best_smiles','tautomer_cluster_evidence_note','final_gate_status'], None)}
      <p class="note">
        Keto/enol forms are compared as a tautomer cluster. Cluster evidence can explain why related drawn forms rank differently,
        but it does not convert a low raw cosine into a high-confidence structure.
        Sidecar: <a href="data/x12798_real_latent_calibration_candidate_summary.csv">real-latent calibration candidate summary</a>.
      </p>
    """


def build_v17_spectrum_latent_summary(v17_root: Path, output_dir: Path) -> dict[str, Any]:
    scores_root = v17_root / "x12798_scores"
    top_path = scores_root / "x12798_v17_candidate_scores_top50.csv"
    cluster_path = scores_root / "x12798_v17_tautomer_cluster_scores.csv"
    final_path = scores_root / "x12798_v17_final_call.json"
    if not top_path.exists():
        return {"available": False, "reason": f"missing {top_path}", "summary": pd.DataFrame(), "clusters": pd.DataFrame()}
    top = pd.read_csv(top_path, low_memory=False)
    clusters = pd.read_csv(cluster_path, low_memory=False) if cluster_path.exists() else pd.DataFrame()
    final_call = read_json(final_path) if final_path.exists() else {}
    target_smiles = {
        "CC(=O)c1ccc(NC(=O)c2ccccc2)cc1",
        "C=C(O)c1ccc(NC(=O)c2ccccc2)cc1",
        "C=C(O)c1ccc(N=C(O)c2ccccc2)cc1",
    }
    target = top.loc[top["canonical_smiles"].astype(str).isin(target_smiles)].copy()
    target["v17_drawn_form"] = target["canonical_smiles"].map(
        {
            "CC(=O)c1ccc(NC(=O)c2ccccc2)cc1": "keto benzamide best-candidate anchor",
            "C=C(O)c1ccc(NC(=O)c2ccccc2)cc1": "enol-like pure de novo candidate",
            "C=C(O)c1ccc(N=C(O)c2ccccc2)cc1": "imidic/enol tautomer-cluster member",
        }
    )
    summary_cols = [
        "v17_drawn_form",
        "rank",
        "canonical_smiles",
        "pure_de_novo_eligible",
        "novelty_status",
        "tautomer_cluster_id",
        "v17_spectrum_quality_bin",
        "v17_raw_original_cosine",
        "v17_best_variant_name",
        "v17_best_variant_raw_cosine",
        "v17_projected_cosine",
        "v17_projected_rank",
        "v17_calibrated_probability",
        "v17_calibration_support_n",
        "v17_calibration_positive_n",
        "v17_test_precision_at_threshold",
        "v17_gate_status",
        "v17_gate_failures",
        "final_status",
    ]
    for col in summary_cols:
        if col not in target.columns:
            target[col] = ""
    summary_csv = "data/x12798_v17_spectrum_latent_candidate_summary.csv"
    summary_json = "data/x12798_v17_spectrum_latent_candidate_summary.json"
    target[summary_cols].to_csv(output_dir / summary_csv, index=False)
    (output_dir / summary_json).write_text(json.dumps(safe_json_records(target[summary_cols]), indent=2) + "\n")
    cluster_csv = "data/x12798_v17_tautomer_cluster_summary.csv"
    cluster_json = "data/x12798_v17_tautomer_cluster_summary.json"
    if not clusters.empty:
        clusters.to_csv(output_dir / cluster_csv, index=False)
        (output_dir / cluster_json).write_text(json.dumps(safe_json_records(clusters), indent=2) + "\n")
    return {
        "available": True,
        "reason": "",
        "summary": target[summary_cols],
        "clusters": clusters,
        "summary_csv": summary_csv,
        "summary_json": summary_json,
        "cluster_csv": cluster_csv if not clusters.empty else "",
        "cluster_json": cluster_json if not clusters.empty else "",
        "final_call": final_call,
        "top50_path": str(top_path),
        "cluster_path": str(cluster_path),
        "final_path": str(final_path),
    }


def v17_spectrum_latent_html(v17: dict[str, Any]) -> str:
    if not v17.get("available"):
        return f"<p class='note'>V17 spectrum-latent alignment unavailable: {html.escape(clean(v17.get('reason')))}</p>"
    final_call = v17.get("final_call", {})
    final_status = clean(final_call.get("final_status"), "provisional_consensus_guided_structure")
    top_smiles = clean(final_call.get("top_smiles"), "")
    gate_status = clean(final_call.get("v17_gate_status"), "")
    return f"""
      <p>
        Raw Spec2Mol cosine compares the observed-spectrum encoder latent with a candidate structure latent.
        V17 separates raw cosine, precursor-suppressed/fragment-weighted cosine, projected cosine, and calibrated probability
        because X-12798 is precursor-dominated and the original raw score is not a direct probability of correctness.
      </p>
      <div class="metrics">
        {metric_card("V17 final status", final_status, "Additive exploratory evidence")}
        {metric_card("V17 top SMILES", top_smiles or "not available", "Not a confirmed structure")}
        {metric_card("V17 gate", gate_status or "not available", "Requires calibrated support")}
      </div>
      {table_html(v17["summary"], ['v17_drawn_form','rank','canonical_smiles','pure_de_novo_eligible','tautomer_cluster_id','v17_raw_original_cosine','v17_best_variant_name','v17_best_variant_raw_cosine','v17_projected_cosine','v17_projected_rank','v17_calibrated_probability','v17_calibration_support_n','v17_calibration_positive_n','v17_test_precision_at_threshold','v17_gate_status'], None)}
      <p class="note">
        V17 keeps raw, variant, projected, and calibrated latent evidence separate. Low raw cosine can be explored by projection and
        spectrum-quality calibration, but no candidate is promoted unless the held-out calibration gate passes.
        Sidecars: <a href="data/x12798_v17_spectrum_latent_candidate_summary.csv">candidate summary</a>
        {', <a href="data/x12798_v17_tautomer_cluster_summary.csv">tautomer-cluster summary</a>' if v17.get('cluster_csv') else ''}.
      </p>
    """


def top_rows(df: pd.DataFrame, rank_column: str, n: int = 10) -> pd.DataFrame:
    out = df.copy()
    out["_rank"] = pd.to_numeric(out.get(rank_column), errors="coerce")
    return out.sort_values("_rank").head(n).drop(columns=["_rank"], errors="ignore")


def render_html(
    output_dir: Path,
    evidence: pd.DataFrame,
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    known_neighbors: pd.DataFrame,
    known_nodes: pd.DataFrame,
    known_edges: pd.DataFrame,
    args: argparse.Namespace,
    library_audit: pd.DataFrame,
    same_library_filtered: pd.DataFrame,
    cross_library_matches: pd.DataFrame,
    consensus_svg: str,
    consensus_rows: pd.DataFrame,
    environment_consensus_svg: str,
    environment_consensus_rows: pd.DataFrame,
    motif_df: pd.DataFrame,
    tree_svg: str,
    tree_summary: pd.DataFrame,
    candidate_consensus: pd.DataFrame,
    final_consensus: dict[str, Any],
    real_latent: dict[str, Any],
    v17_latent: dict[str, Any],
) -> str:
    deep_final = read_json(args.deep_dive_final_json)
    lib400_final = read_json(args.lib400_final_json)
    merged_final = read_json(args.merged_final_json)
    comparison = pd.read_csv(args.encoder_comparison_csv, low_memory=False)
    deep = pd.read_csv(args.deep_dive_csv, low_memory=False)
    lib400_top10 = pd.read_csv(args.lib400_pure_top10, low_memory=False)
    merged_top10 = pd.read_csv(args.merged_pure_top10, low_memory=False)
    network = network_svg(nodes, edges)
    known_network = known_similarity_network_svg(known_nodes, known_edges)

    nonpure_name = deep_final.get("final_smiles_name", "N-(4-acetylphenyl)benzamide")
    nonpure_smiles = deep_final.get("final_smiles", "")
    lib400_smiles = lib400_final.get("final_smiles", "")
    lib400_cosine = lib400_final.get("final_ms2_embedding_cosine", "")
    merged_smiles = merged_final.get("final_smiles", "")
    merged_cosine = merged_final.get("final_ms2_embedding_cosine", "")
    gate_pass = bool(lib400_final.get("ms2_embedding_gate_pass")) or bool(merged_final.get("ms2_embedding_gate_pass"))
    final_status = "provisional_non_pure_structure"

    deep_top10 = top_rows(deep, "deep_dive_rank")
    lib400_top10 = top_rows(lib400_top10, "pure_denovo_rank")
    merged_top10 = top_rows(merged_top10, "pure_denovo_rank")

    metric_cards = "".join(
        [
            metric_card("Final status", final_status, "No high-confidence SMILES promotion"),
            metric_card("Candidate nodes", len(nodes), "Unique canonical predicted structures"),
            metric_card("Pairwise edges", len(edges), "All-vs-all Morgan Tanimoto pairs"),
            metric_card("Filtered lib400 hits", len(same_library_filtered), "Exact same-library catalog matches excluded"),
            metric_card("Other-library hits", len(cross_library_matches), "Listed outside prediction network"),
            metric_card("Best real-latent cosine", f"{as_float(lib400_cosine) or 0.0:.4f}", "Threshold is 0.70"),
            metric_card("MS2 gate pass", str(gate_pass), "Both encoder branches failed"),
        ]
    )

    evidence_notes = [
        f"Target formula: {clean(deep_final.get('target_formula'), 'C15H13NO2')}",
        "Observed feature: feature_X_400_mz240.1017_ri2011_X-12798",
        "Corrected isotope evidence supports the formula-level C15H13NO2 assignment.",
        "Positive-mode adduct support is used as a constraint after structural evidence.",
        "Phase-I relationship evidence includes a link to X-12127, but this is not structural confirmation.",
        "Pure de novo candidates did not pass the real MS2 embedding gate.",
    ]
    evidence_list = "".join(f"<li>{html.escape(note)}</li>" for note in evidence_notes)
    cross_library_note = (
        "No cross-library exact local catalog matches were found."
        if cross_library_matches.empty
        else (
            "Cross-library-only matches were found and listed separately; they remain excluded from same-library filtering "
            "because none include library 400."
        )
    )
    closest_knowns = ", ".join(known_neighbors.head(2)["canonical_name"].astype(str).tolist()) if not known_neighbors.empty else ""
    known_similarity_values = pd.to_numeric(known_neighbors.get("similarity_to_anchor"), errors="coerce")
    knowns_ge_0p50 = int((known_similarity_values >= 0.50).sum())
    knowns_ge_0p35 = int((known_similarity_values >= 0.35).sum())
    known_metric_cards = "".join(
        [
            metric_card("Anchor", ANCHOR_CANDIDATE_NAME, "Current best candidate"),
            metric_card("Closest knowns", closest_knowns, "Morgan radius-2 Tanimoto"),
            metric_card("Knowns >=0.50", knowns_ge_0p50, "Top-30 neighbor set"),
            metric_card("Knowns >=0.35", knowns_ge_0p35, "Context-level analogs"),
        ]
    )

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <meta name="x12798-ms2-gate" content="ms2_embedding_cosine >= 0.70" />
  <title>X-12798 Structure Resolution</title>
  <style>
    :root {{
      color-scheme: light;
      --navy:#001D37; --blue:#28496E; --coral:#FF4F57; --purple:#855DA3; --amber:#FCA750; --mist:#D3DEE6;
      --paper:#FFFFFF; --wash:#F1FAFB; --ink:#15304B; --muted:#5C6B7A; --line:rgba(40,73,110,.18);
      --green:#0E8A68; --shadow:0 24px 58px rgba(0,29,55,.14);
    }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; font-family:"Montserrat","Avenir Next","Segoe UI",Arial,sans-serif; color:var(--ink); background:linear-gradient(180deg,#fff 0,#f7fbfc 22%,#eef6f9 100%); line-height:1.5; }}
    header {{ position:relative; padding:54px 54px 34px; background:linear-gradient(135deg,rgba(0,29,55,.98),rgba(40,73,110,.92) 56%,rgba(255,79,87,.88)); color:#fff; overflow:hidden; }}
    header::after {{ content:""; position:absolute; inset:auto 0 0 0; height:10px; background:linear-gradient(90deg,var(--amber),var(--coral) 36%,var(--purple) 68%,var(--mist)); }}
    main {{ max-width:1380px; margin:0 auto; padding:28px 28px 64px; }}
    h1 {{ margin:0; font-size:46px; line-height:1; letter-spacing:0; }}
    h2 {{ margin:0 0 16px; font-size:26px; line-height:1.08; letter-spacing:0; color:var(--navy); }}
    h3 {{ margin:20px 0 10px; font-size:18px; letter-spacing:0; color:var(--navy); }}
    p {{ margin:8px 0; }}
    .eyebrow {{ margin:0 0 10px; color:rgba(255,255,255,.74); font-size:12px; font-weight:800; letter-spacing:.12em; text-transform:uppercase; }}
    .subtitle {{ margin-top:10px; max-width:1040px; color:rgba(255,255,255,.88); font-size:17px; }}
    .hero-grid {{ display:grid; grid-template-columns:minmax(0,1.2fr) minmax(280px,.8fr); gap:24px; align-items:start; }}
    .section {{ position:relative; background:var(--paper); border:1px solid var(--line); border-radius:8px; padding:24px; margin:20px 0; box-shadow:var(--shadow); overflow:hidden; }}
    .section::before {{ content:""; display:block; width:88px; height:5px; border-radius:999px; background:linear-gradient(90deg,var(--navy),var(--coral) 58%,var(--amber)); margin-bottom:16px; }}
    .metrics {{ display:grid; grid-template-columns:repeat(5,minmax(0,1fr)); gap:12px; margin-top:18px; }}
    .metric-card {{ background:linear-gradient(180deg,rgba(241,250,251,.98),#fff); border:1px solid var(--line); border-radius:8px; padding:13px 14px; box-shadow:0 10px 24px rgba(0,29,55,.06); }}
    .metric-label {{ color:var(--blue); font-size:12px; text-transform:uppercase; letter-spacing:.08em; font-weight:800; }}
    .metric-value {{ margin-top:4px; font-size:23px; font-weight:700; color:var(--navy); }}
    .metric-note {{ color:var(--muted); font-size:12px; margin-top:2px; }}
    .call-grid {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:18px; }}
    .structure-card {{ border:1px solid var(--line); border-radius:8px; padding:16px; background:linear-gradient(180deg,rgba(241,250,251,.9),#fff 38%); }}
    .structure-card img {{ width:100%; max-height:260px; object-fit:contain; background:#fff; border:1px solid #e8edf1; border-radius:6px; }}
    code {{ background:#eef3f6; border-radius:4px; padding:1px 4px; }}
    .badge {{ display:inline-block; border-radius:999px; padding:5px 10px; font-size:12px; font-weight:800; border:1px solid rgba(255,255,255,.16); background:rgba(255,255,255,.09); color:#fff; margin-right:6px; }}
    .badge.warn {{ color:#fff; background:rgba(252,167,80,.2); border-color:rgba(252,167,80,.4); }}
    .badge.fail {{ color:#fff; background:rgba(255,79,87,.2); border-color:rgba(255,79,87,.4); }}
    .note {{ color:var(--muted); }}
    .table-wrap {{ overflow:auto; border:1px solid var(--line); border-radius:8px; }}
    table {{ border-collapse:collapse; width:100%; min-width:860px; background:#fff; }}
    th,td {{ padding:8px 10px; border-bottom:1px solid #e8edf1; text-align:left; vertical-align:top; font-size:13px; }}
    th {{ position:sticky; top:0; background:linear-gradient(180deg,#f7fbfc,#eef6f9); z-index:1; font-size:12px; color:var(--blue); }}
    .network-controls {{ display:flex; gap:10px; flex-wrap:wrap; align-items:center; margin-bottom:12px; }}
	    .network-controls button {{ border:1px solid var(--line); background:#fff; border-radius:6px; padding:7px 10px; cursor:pointer; }}
	    .network-controls button.active {{ background:#17324d; color:#fff; border-color:#17324d; }}
	    .network-svg {{ width:100%; height:auto; border:1px solid var(--line); border-radius:8px; background:#f7f9fb; }}
	    .known-network-svg {{ width:100%; height:auto; border:1px solid var(--line); border-radius:8px; background:#f7f9fb; }}
	    .net-node text {{ font-size:12px; text-anchor:middle; fill:#17212b; paint-order:stroke; stroke:#fff; stroke-width:4px; stroke-linejoin:round; }}
	    .net-node circle {{ stroke:#fff; stroke-width:2.2px; filter:drop-shadow(0 1px 2px rgba(10,20,30,.22)); }}
	    .known-net-node text {{ font-size:12px; text-anchor:middle; fill:#17212b; paint-order:stroke; stroke:#fff; stroke-width:4px; stroke-linejoin:round; }}
	    .known-net-node .known-node-score {{ font-size:11px; font-weight:800; fill:#526372; }}
	    .known-net-node circle {{ stroke:#fff; stroke-width:2.2px; filter:drop-shadow(0 1px 2px rgba(10,20,30,.22)); }}
	    .known-card-grid {{ display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:12px; margin:14px 0; }}
	    .known-card {{ background:#fbfcfd; border:1px solid var(--line); border-radius:8px; padding:12px; }}
	    .known-card img {{ width:100%; max-height:160px; object-fit:contain; background:#fff; border:1px solid #e8edf1; border-radius:6px; }}
	    .known-card h3 {{ margin:8px 0 4px; font-size:14px; }}
    .motif-grid {{ display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:14px; }}
    .motif-card {{ background:#fbfcfd; border:1px solid var(--line); border-radius:8px; padding:14px; }}
    .motif-title {{ font-weight:800; margin-bottom:8px; }}
    .motif-badge {{ display:inline-block; border-radius:999px; padding:3px 9px; font-size:12px; font-weight:800; margin-bottom:10px; }}
    .motif-badge.strong {{ background:#edf8ef; color:#276738; border:1px solid #bfe5c7; }}
    .motif-badge.mixed {{ background:#fff6e6; color:#8a4b00; border:1px solid #f2d3a2; }}
    .motif-badge.weak {{ background:#eef3f6; color:#465764; border:1px solid #cbd6df; }}
    .motif-measures {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:8px; margin:4px 0 8px; }}
    .motif-measures span {{ background:#fff; border:1px solid #e8edf1; border-radius:6px; padding:8px; }}
    .motif-measures strong {{ display:block; font-size:20px; }}
    .motif-measures small {{ color:var(--muted); font-size:11px; }}
    .decision-tree-box {{ background:#fff; border:1px solid var(--line); border-radius:8px; padding:12px; overflow:auto; }}
    .decision-tree-svg {{ width:100%; min-width:1040px; height:auto; display:block; }}
    .tree-title {{ font-size:20px; font-weight:800; fill:#17212b; }}
    .tree-node-title {{ font-size:13px; font-weight:800; text-anchor:middle; }}
    .tree-node-metric {{ font-size:12px; font-weight:800; text-anchor:middle; fill:#394956; }}
    .tree-node-note,.tree-leaf-count,.tree-leaf-top {{ font-size:11px; text-anchor:middle; fill:#394956; }}
    .tree-leaf-count,.tree-leaf-top {{ fill:#fff; font-weight:800; }}
    .tree-edge-label {{ font-size:11px; font-weight:800; fill:#394956; paint-order:stroke; stroke:#f7f9fb; stroke-width:4px; }}
    .candidate-consensus-grid {{ display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:14px; }}
    .candidate-consensus-card {{ background:#fbfcfd; border:1px solid var(--line); border-radius:8px; padding:13px; }}
    .candidate-consensus-card h3 {{ margin:0 0 8px; font-size:14px; }}
    .candidate-consensus-visual-pair {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:8px; }}
    .candidate-consensus-figure {{ width:100%; min-height:210px; background:#fff; border:1px solid #e8edf1; border-radius:6px; overflow:hidden; }}
    .candidate-consensus-panel-label {{ color:#394956; background:#f7fafc; border-bottom:1px solid #e8edf1; padding:5px 8px; font-size:11px; font-weight:800; }}
    .candidate-consensus-figure .consensus-svg {{ width:100%; height:auto; display:block; }}
    .candidate-consensus-card code {{ display:block; overflow-wrap:anywhere; font-size:11px; }}
    .candidate-consensus-metrics {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:8px; margin:8px 0; }}
    .candidate-consensus-metrics span {{ background:#fff; border:1px solid #e8edf1; border-radius:6px; padding:8px; }}
    .candidate-consensus-metrics strong {{ display:block; font-size:18px; }}
    .candidate-consensus-metrics small {{ color:var(--muted); font-size:11px; }}
    .consensus-grid {{ display:grid; grid-template-columns:minmax(360px,.9fr) minmax(0,1.1fr); gap:18px; align-items:start; }}
    .consensus-box {{ background:#fff; border:1px solid var(--line); border-radius:8px; padding:12px; overflow:hidden; }}
    .consensus-svg {{ width:100%; height:auto; display:block; }}
    .consensus-atom text {{ fill:#fff; font-size:16px; font-weight:800; pointer-events:none; }}
    .consensus-atom .atom-confidence {{ fill:#17212b; font-size:11px; font-weight:800; }}
    .consensus-scale {{ margin-top:10px; }}
    .consensus-scale-bar {{ height:16px; border-radius:999px; background:linear-gradient(90deg,#2C7BB6 0%,#FDAE61 50%,#D7191C 100%); border:1px solid var(--line); }}
    .consensus-scale-labels {{ display:grid; grid-template-columns:repeat(3,1fr); gap:8px; margin-top:5px; color:var(--muted); font-size:11px; font-weight:700; }}
    .consensus-scale-labels span:nth-child(2) {{ text-align:center; }}
    .consensus-scale-labels span:nth-child(3) {{ text-align:right; }}
    footer {{ color:var(--muted); font-size:12px; padding:0 28px 34px; max-width:1320px; margin:0 auto; }}
	    @media (max-width:980px) {{ .hero-grid,.call-grid,.metrics,.consensus-grid,.motif-grid,.candidate-consensus-grid,.candidate-consensus-visual-pair,.known-card-grid {{ grid-template-columns:1fr; }} header {{ padding:34px 24px 26px; }} h1 {{ font-size:38px; }} }}
  </style>
</head>
<body>
  <header>
    <div class="hero-grid">
      <div>
	        <p class="eyebrow">X-Compound Elucidation Dashboard</p>
	        <h1>X-12798</h1>
	        <p class="subtitle">
	          lib_400_positive Level 4 structure-resolution report for formula <code>C15H13NO2</code>.
	          {html.escape(PREDICTION_SUMMARY_SENTENCE)}
	        </p>
	        <p>
	          <span class="badge warn">{html.escape(final_status)}</span>
	          <span class="badge fail">pure de novo MS2 gate failed</span>
	          <span class="badge">Best-candidate known-similarity anchor</span>
	        </p>
      </div>
      <div>{metric_cards}</div>
    </div>
  </header>
  <main>
    <section class="section">
	      <h2>Top Candidate Compound</h2>
	      <div class="call-grid">
	        <div class="structure-card">
	          <h3>Most Likely Provisional Structure</h3>
	          <img src="assets/x12798_provisional_non_pure.svg" alt="Provisional non-pure structure" />
	          <p><strong>{html.escape(clean(nonpure_name))}</strong></p>
	          <p><code>{html.escape(clean(nonpure_smiles))}</code></p>
	          <p class="note">Status: {html.escape(final_status)}. The anchor is the current best candidate, not a confirmed structure.</p>
        </div>
        <div class="structure-card">
          <h3>Best Pure De Novo Structure Candidate</h3>
          <img src="assets/x12798_best_pure_denovo_lib400_only.svg" alt="Best pure de novo structure" />
          <p><strong>lib400-only branch rank 1</strong></p>
          <p><code>{html.escape(clean(lib400_smiles))}</code></p>
          <p class="note">Real latent cosine {html.escape(f"{as_float(lib400_cosine) or 0.0:.4f}")}; required gate is ms2_embedding_cosine &gt;= 0.70.</p>
        </div>
      </div>
      <p>
	        {html.escape(PREDICTION_SUMMARY_SENTENCE)}
	        The current best candidate remains provisional because the real-latent pure de novo branches do not pass
	        <code>ms2_embedding_cosine &gt;= 0.70</code>. Pure de novo structures are retained as nearby predicted-candidate hypotheses, not promoted calls.
      </p>
    </section>

    <section class="section">
      <h2>Parameters Contributing To Identification Confidence</h2>
      <ul>{evidence_list}</ul>
	      <p class="note">The sections below separate known-metabolite context from predicted-candidate alternatives.</p>
	    </section>

    <section class="section">
      <h2>Real-Latent Evidence Calibration</h2>
      {real_latent_calibration_html(real_latent)}
    </section>

    <section class="section">
      <h2>V17 Spectrum-Latent Alignment</h2>
      {v17_spectrum_latent_html(v17_latent)}
    </section>

    <section class="section">
      <h2>Structure Confidence Heat Map</h2>
      <p>
        These images show the final consensus-guided top candidate with atom-position confidence from the consensus-guided top-20 set.
        Graph-position confidence asks whether each atom maps into the other top consensus-guided candidates; environment-position confidence asks whether
        the mapped atom also has similar local substituent, heteroatom, and ring context.
      </p>
      {final_consensus_guided_html(final_consensus)}
	    </section>

	    <section class="section">
	      <h2>Similarity To Known Metabolites</h2>
	      <p>
	        This network compares the best X-12798 candidate, <strong>{html.escape(ANCHOR_CANDIDATE_NAME)}</strong>, to structured known metabolites
	        in the local Metabolon global library. It is a chemical-context view, not a confirmation claim.
	      </p>
	      <div class="metrics">{known_metric_cards}</div>
	      <div class="network-controls">
	        <span class="note">Known-similarity floor:</span>
	        <button type="button" data-known-threshold="0.35" class="known-threshold active">0.35</button>
	        <button type="button" data-known-threshold="0.50" class="known-threshold">0.50</button>
	        <button type="button" data-known-threshold="0.00" class="known-threshold">All</button>
	      </div>
	      {known_network}
	      {known_neighbor_cards_html(known_neighbors, output_dir, 4)}
	      <h3>Most Similar Structured Knowns</h3>
	      {table_html(known_neighbors, ['rank','chemical_id','canonical_name','similarity_to_anchor','similarity_bucket','super_pathway','pathway','canonical_smiles'], 30)}
	      <p class="note">
	        Sidecars: <a href="data/x12798_known_similarity_neighbors.csv">known neighbors</a>,
	        <a href="data/x12798_known_similarity_network_nodes.csv">known network nodes</a>,
	        <a href="data/x12798_known_similarity_network_edges.csv">known network edges</a>,
	        <a href="data/x12798_all_known_similarity_to_anchor.csv">all known similarities</a>.
	      </p>
	    </section>

	    <section class="section">
	      <h2>Metabolon Library Filtering</h2>
      <p>
        Candidate evidence and the similarity network are filtered before rendering: any prediction with an exact local Metabolon catalog
        match on library <code>400</code> is removed from the predicted-candidate set. This report found
        <strong>{len(same_library_filtered)}</strong> same-library <code>400</code> matches.
      </p>
      <p class="note">{html.escape(cross_library_note)}</p>
      <h3>Cross-Library Metabolon Matches</h3>
      {table_html(cross_library_matches, ['candidate_rank','display_name','canonical_smiles','detected_library_ids','catalog_comp_ids','catalog_compound_names','metabolon_match_status','match_basis','notes'], None)}
      <h3>Same-Library 400 Filtered Predictions</h3>
      {table_html(same_library_filtered, ['candidate_rank','display_name','canonical_smiles','detected_library_ids','catalog_comp_ids','catalog_compound_names','metabolon_match_status','match_basis','notes'], None)}
      <p class="note">
        Full audit sidecar:
        <a href="data/x12798_metabolon_library_match_audit.csv">x12798_metabolon_library_match_audit.csv</a>.
      </p>
    </section>

	    <section class="section">
	      <h2>Similarity Network</h2>
	      <p>
	        This view compares predicted structural possibilities to each other. Use the known-metabolite network above to see how the best candidate
	        relates to structured known metabolites. Edges use Morgan radius-2, 2048-bit Tanimoto similarity.
	      </p>
      <div class="network-controls">
        <span class="note">Similarity floor:</span>
        <button type="button" data-threshold="0.25" class="threshold active">0.25</button>
        <button type="button" data-threshold="0.50" class="threshold">0.50</button>
        <button type="button" data-threshold="0.70" class="threshold">0.70</button>
        <button type="button" data-threshold="0.00" class="threshold">All</button>
      </div>
      {network}
      <p class="note">
        Sidecars: <a href="data/x12798_candidate_similarity_network_nodes.csv">nodes CSV</a>,
        <a href="data/x12798_candidate_similarity_network_edges.csv">edges CSV</a>,
        <a href="data/x12798_candidate_similarity_matrix.csv">similarity matrix</a>.
      </p>
    </section>

    <section class="section">
      <h2>Motif Consensus Summary</h2>
      <p>
        Present in candidates asks whether the motif appears anywhere; same placement asks whether it maps to the same region of the displayed provisional scaffold.
      </p>
      {motif_summary_html(motif_df)}
      <p class="note">
        Sidecars: <a href="data/x12798_motif_consensus.csv">x12798_motif_consensus.csv</a>,
        <a href="data/x12798_motif_consensus.json">x12798_motif_consensus.json</a>.
      </p>
    </section>

    <section class="section">
      <h2>Structural Possibility Decision Tree</h2>
      <p>
        This tree is a hypothesis map, not a confirmation workflow. Each branch asks whether a motif is present anywhere,
        whether it appears in the same scaffold region, or whether candidates support a moved/variant analog.
      </p>
      <div class="decision-tree-box">
        {tree_svg}
      </div>
      <p class="note">
        Low-confidence branches mean either motif absence or motif relocation depending on the edge label.
        Sidecars: <a href="data/x12798_structural_possibility_tree_nodes.csv">tree nodes</a>,
        <a href="data/x12798_structural_possibility_tree_edges.csv">tree edges</a>,
        <a href="data/x12798_structural_possibility_summary.csv">possibility summary</a>.
      </p>
      {table_html(tree_summary, ['structural_possibility_class','candidate_count','top_candidate_label','top_candidate_smiles','reader_interpretation','dominant_gate_status'], None)}
    </section>

    <section class="section">
      <h2>Scaffold Position Consensus</h2>
      <p>
        This view colors atom positions on the provisional non-pure scaffold by how often candidate structures map onto those positions.
        Blue does not necessarily mean the chemical group is absent; it can also mean the group appears elsewhere, as a tautomer,
        or as a scaffold analog. This scaffold is an orientation structure only.
      </p>
      <p>
        Graph-position consensus asks whether an atom is part of a shared MCS graph. Chemical-environment consensus asks whether
        that mapped atom also has similar local substituent and heteroatom context. Aromatic-ring atoms can be graph-consistent while environment-inconsistent.
      </p>
      <div class="consensus-grid">
        <div class="consensus-box">
          <h3>Graph-Position Consensus</h3>
          {consensus_svg}
          <div class="consensus-scale" aria-label="Consensus heatmap scale">
            <div class="consensus-scale-bar"></div>
            <div class="consensus-scale-labels"><span>0.00 low</span><span>0.50</span><span>1.00 high</span></div>
          </div>
        </div>
        <div class="consensus-box">
          <h3>Chemical Environment Consensus</h3>
          {environment_consensus_svg}
          <div class="consensus-scale" aria-label="Environment consensus heatmap scale">
            <div class="consensus-scale-bar"></div>
            <div class="consensus-scale-labels"><span>0.00 low</span><span>0.50</span><span>1.00 high</span></div>
          </div>
        </div>
      </div>
      <div class="consensus-grid">
        <div>
          {table_html(consensus_rows, ['atom_index','element','consensus_support_fraction','weighted_consensus_support_fraction','support_count','valid_candidate_count'], None)}
          <p class="note">
            Sidecars: <a href="data/x12798_consensus_atom_support.csv">x12798_consensus_atom_support.csv</a>,
            <a href="data/x12798_scaffold_position_consensus_atom_support.csv">x12798_scaffold_position_consensus_atom_support.csv</a>.
          </p>
        </div>
        <div>
          {table_html(environment_consensus_rows, ['atom_index','element','weighted_graph_consensus_support_fraction','weighted_environment_consensus_support_fraction','environment_gap_fraction','valid_candidate_count'], None)}
          <p class="note">
            Sidecar: <a href="data/x12798_scaffold_environment_consensus_atom_support.csv">x12798_scaffold_environment_consensus_atom_support.csv</a>.
          </p>
        </div>
      </div>
    </section>

    <section class="section">
      <h2>Best-Candidate Evidence</h2>
      {table_html(evidence, CANDIDATE_COLUMNS, None)}
    </section>

    <section class="section">
      <h2>Best-Candidate Evidence: Candidate-Specific Consensus Images</h2>
      <p>
        This section includes one visual structure for every row in Best-Candidate Evidence. Each image uses the same
        maximum-common-substructure atom-coloring approach as Scaffold Position Consensus, but uses that specific
        candidate structure as the displayed reference. Red atoms are consistently mapped across the other retained
        top structures; blue atoms are uncommon, repositioned, or absent in the other candidates.
        If a benzene ring is red in graph overlap but blue/amber in environment-aware consensus, the ring scaffold is common
        but the local chemical environment is not.
      </p>
      {candidate_consensus_gallery_html(candidate_consensus, output_dir)}
      <p class="note">
        Summary sidecar: <a href="data/x12798_candidate_specific_consensus_summary.csv">x12798_candidate_specific_consensus_summary.csv</a>.
      </p>
    </section>

    <section class="section">
      <h2>Non-Pure Top 10</h2>
      {table_html(deep_top10, ['deep_dive_rank','canonical_smiles','predicted_name','deep_dive_score','isotope_status','ms2_proxy_score','predicted_pathway','predicted_subpathway','source_provenance'], 10)}
    </section>

    <section class="section">
      <h2>Pure De Novo Top 10 By Real-Latent Branch</h2>
      <h3>Lib400-Only Encoder</h3>
      {table_html(lib400_top10, ['pure_denovo_rank','canonical_smiles','constrained_rank_score','ms2_embedding_cosine','ms2_embedding_rank','ms2_embedding_match_status','isotope_adduct_promoted','promotion_block_reason','predicted_pathway','predicted_subpathway'], 10)}
      <h3>Merged Lib-Conditioned Encoder</h3>
      {table_html(merged_top10, ['pure_denovo_rank','canonical_smiles','constrained_rank_score','ms2_embedding_cosine','ms2_embedding_rank','ms2_embedding_match_status','isotope_adduct_promoted','promotion_block_reason','predicted_pathway','predicted_subpathway'], 10)}
    </section>

    <section class="section">
      <h2>Encoder Branch Comparison</h2>
      {table_html(comparison, ['branch_label','top_smiles','top_score','top_margin','top_ms2_embedding_cosine','top_ms2_embedding_rank','top_ms2_embedding_status','top_isotope_adduct_promoted','top_promotion_block_reason','ms2_embedding_gate_pass','final_call_status'], None)}
      <p class="note">Merged branch top SMILES: <code>{html.escape(clean(merged_smiles))}</code>; cosine {html.escape(f"{as_float(merged_cosine) or 0.0:.4f}")}.</p>
    </section>

    <section class="section">
      <h2>Interpretation</h2>
      <p>
        {html.escape(clean(nonpure_name))} remains the most likely provisional structure for X-12798.
        The strongest pure de novo candidate, <code>{html.escape(clean(lib400_smiles))}</code>, is structurally close enough to explore
        scaffold and tautomer alternatives, but its real observed-latent cosine is far below the promotion threshold.
      </p>
      <p>
        Since no pure de novo candidate passes the real MS2 embedding gate, isotope/adduct support is not allowed to promote a high-confidence pure de novo SMILES.
        The appropriate status is provisional, not confirmed.
      </p>
    </section>
  </main>
  <footer>
    Generated from local X-12798 artifacts. Data files live under <code>data/</code>; rendered structures live under <code>assets/</code>.
  </footer>
  <script>
    function updateEdges(threshold) {{
      document.querySelectorAll('.threshold').forEach(btn => {{
        btn.classList.toggle('active', Number(btn.dataset.threshold) === threshold);
      }});
      document.querySelectorAll('.net-edge').forEach(edge => {{
        const sim = Number(edge.dataset.similarity || 0);
        edge.style.display = sim >= threshold ? '' : 'none';
      }});
    }}
    function updateKnownEdges(threshold) {{
      document.querySelectorAll('.known-threshold').forEach(btn => {{
        btn.classList.toggle('active', Number(btn.dataset.knownThreshold) === threshold);
      }});
      document.querySelectorAll('.known-net-edge').forEach(edge => {{
        const sim = Number(edge.dataset.knownSimilarity || 0);
        edge.style.display = sim >= threshold ? '' : 'none';
      }});
      document.querySelectorAll('.known-net-node').forEach(node => {{
        if (node.dataset.nodeType === 'anchor_best_candidate') {{
          node.style.display = '';
          return;
        }}
        const edge = document.querySelector(`.known-net-edge[data-target="${{node.dataset.id}}"]`);
        const sim = edge ? Number(edge.dataset.knownSimilarity || 0) : 0;
        node.style.display = sim >= threshold ? '' : 'none';
      }});
    }}
    document.querySelectorAll('.threshold').forEach(btn => {{
      btn.addEventListener('click', () => updateEdges(Number(btn.dataset.threshold || 0.25)));
    }});
    document.querySelectorAll('.known-threshold').forEach(btn => {{
      btn.addEventListener('click', () => updateKnownEdges(Number(btn.dataset.knownThreshold || 0.35)));
    }});
    updateEdges(0.25);
    updateKnownEdges(0.35);
  </script>
</body>
</html>
"""
    return html_text


def validate_inputs(args: argparse.Namespace) -> None:
    required = [
        args.resolution_md,
        args.encoder_comparison_csv,
        args.deep_dive_csv,
        args.deep_dive_final_json,
        args.lib400_pure_top10,
        args.merged_pure_top10,
        args.lib400_final_json,
        args.merged_final_json,
        args.v15_names_csv,
        args.library_match_md,
        args.reference_db,
        *args.metabolon_library_catalogs,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        message = "Missing required artifacts:\n" + "\n".join(f"- {path}" for path in missing) + "\n\n" + RERUN_COMMANDS
        raise FileNotFoundError(message)


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    validate_inputs(args)
    output_dir = args.output_dir
    data_dir = output_dir / "data"
    assets_dir = output_dir / "assets"
    data_dir.mkdir(parents=True, exist_ok=True)
    assets_dir.mkdir(parents=True, exist_ok=True)

    raw_evidence = collect_candidates(args)
    library_audit, same_library_filtered, cross_library_matches, evidence = build_metabolon_library_audit(
        raw_evidence,
        library_match_md=args.library_match_md,
        catalog_paths=args.metabolon_library_catalogs,
        target_library_id="400",
    )
    evidence_csv = data_dir / "x12798_best_candidate_evidence.csv"
    evidence_json = data_dir / "x12798_best_candidate_evidence.json"
    library_audit.to_csv(data_dir / "x12798_metabolon_library_match_audit.csv", index=False)
    (data_dir / "x12798_metabolon_library_match_audit.json").write_text(
        json.dumps(safe_json_records(library_audit), indent=2) + "\n"
    )
    same_library_filtered.to_csv(data_dir / "x12798_same_library_400_filtered_candidates.csv", index=False)
    cross_library_matches.to_csv(data_dir / "x12798_cross_library_metabolon_matches.csv", index=False)

    nodes, edges, matrix, network_data = build_network(evidence)
    nodes.to_csv(data_dir / "x12798_candidate_similarity_network_nodes.csv", index=False)
    edges.to_csv(data_dir / "x12798_candidate_similarity_network_edges.csv", index=False)
    matrix.to_csv(data_dir / "x12798_candidate_similarity_matrix.csv", index=False)
    (data_dir / "x12798_candidate_similarity_network.json").write_text(json.dumps(network_data, indent=2) + "\n")

    provisional = evidence.loc[evidence["candidate_family"].eq("provisional_non_pure")]
    reference_smiles = clean(provisional.iloc[0]["canonical_smiles"]) if not provisional.empty else ""
    all_known, known_neighbors, known_nodes, known_edges, known_network_data = build_anchor_known_similarity(
        reference_smiles or ANCHOR_CANDIDATE_SMILES,
        args.reference_db,
        top_n=30,
    )
    known_neighbors, known_nodes, known_rendered_assets, known_render_failures = render_known_similarity_assets(
        output_dir,
        known_neighbors,
        known_nodes,
        max_neighbor_assets=12,
    )
    known_network_data = {
        "nodes": safe_json_records(known_nodes),
        "edges": safe_json_records(known_edges),
        "anchor_candidate_name": ANCHOR_CANDIDATE_NAME,
        "anchor_candidate_smiles": reference_smiles or ANCHOR_CANDIDATE_SMILES,
        "layout": "best_candidate_centered_known_similarity",
        "rendered_edge_similarity_floor": 0.35,
    }
    all_known.to_csv(data_dir / "x12798_all_known_similarity_to_anchor.csv", index=False)
    known_neighbors.to_csv(data_dir / "x12798_known_similarity_neighbors.csv", index=False)
    known_nodes.to_csv(data_dir / "x12798_known_similarity_network_nodes.csv", index=False)
    known_edges.to_csv(data_dir / "x12798_known_similarity_network_edges.csv", index=False)
    (data_dir / "x12798_known_similarity_network.json").write_text(json.dumps(known_network_data, indent=2) + "\n")
    consensus_input = evidence.loc[~evidence["canonical_smiles"].eq(reference_smiles)].head(23).copy()
    consensus_rows, consensus_svg = consensus_atom_support(reference_smiles, consensus_input)
    environment_consensus_rows, environment_consensus_svg = environment_aware_consensus_atom_support(
        reference_smiles, consensus_input
    )
    consensus_df = pd.DataFrame(consensus_rows)
    environment_consensus_df = pd.DataFrame(environment_consensus_rows)
    consensus_df.to_csv(data_dir / "x12798_consensus_atom_support.csv", index=False)
    consensus_df.to_csv(data_dir / "x12798_scaffold_position_consensus_atom_support.csv", index=False)
    environment_consensus_df.to_csv(data_dir / "x12798_scaffold_environment_consensus_atom_support.csv", index=False)
    motif_df = motif_consensus(reference_smiles, consensus_input, consensus_df)
    motif_df.to_csv(data_dir / "x12798_motif_consensus.csv", index=False)
    (data_dir / "x12798_motif_consensus.json").write_text(json.dumps(safe_json_records(motif_df), indent=2) + "\n")

    candidate_flags = candidate_motif_flags(evidence)
    tree_nodes, tree_edges, tree_summary = build_motif_decision_tree(motif_df, candidate_flags, evidence)
    candidate_flags.to_csv(data_dir / "x12798_candidate_motif_flags.csv", index=False)
    tree_nodes.to_csv(data_dir / "x12798_structural_possibility_tree_nodes.csv", index=False)
    tree_edges.to_csv(data_dir / "x12798_structural_possibility_tree_edges.csv", index=False)
    tree_summary.to_csv(data_dir / "x12798_structural_possibility_summary.csv", index=False)
    (data_dir / "x12798_structural_possibility_tree.json").write_text(
        json.dumps(
            {
                "nodes": safe_json_records(tree_nodes),
                "edges": safe_json_records(tree_edges),
                "summary": safe_json_records(tree_summary),
                "layout": "deterministic_left_to_right_motif_decision_tree",
            },
            indent=2,
        )
        + "\n"
    )

    candidate_consensus = build_candidate_specific_consensus(evidence, output_dir)
    final_consensus = build_final_consensus_guided_visuals(args.consensus_guided_root, output_dir)
    real_latent = build_real_latent_calibration_summary(args.consensus_guided_root, args.real_latent_calibration_root, output_dir)
    v17_latent = build_v17_spectrum_latent_summary(args.v17_spectrum_latent_root, output_dir)
    candidate_consensus.to_csv(data_dir / "x12798_candidate_specific_consensus_summary.csv", index=False)
    (data_dir / "x12798_candidate_specific_consensus_summary.json").write_text(
        json.dumps(safe_json_records(candidate_consensus), indent=2) + "\n"
    )
    consensus_columns = [
        "candidate_id",
        "candidate_consensus_asset",
        "candidate_consensus_mean_weighted_fraction",
        "candidate_consensus_atoms_ge_0p5_fraction",
        "candidate_environment_consensus_asset",
        "candidate_environment_mean_weighted_fraction",
        "candidate_environment_atoms_ge_0p5_fraction",
    ]
    for column in consensus_columns:
        if column not in candidate_consensus.columns:
            candidate_consensus[column] = ""
    evidence = evidence.drop(columns=[column for column in consensus_columns if column in evidence.columns and column != "candidate_id"])
    evidence = evidence.merge(candidate_consensus[consensus_columns], on="candidate_id", how="left")
    evidence.to_csv(evidence_csv, index=False)
    evidence_json.write_text(json.dumps(safe_json_records(evidence), indent=2) + "\n")

    asset_manifest = render_assets(output_dir, evidence, nodes, args)
    (output_dir / "assets/x12798_consensus_structure_heatmap.svg").write_text(consensus_svg)
    (output_dir / "assets/x12798_scaffold_position_consensus.svg").write_text(consensus_svg)
    (output_dir / "assets/x12798_scaffold_environment_consensus.svg").write_text(environment_consensus_svg)
    tree_svg = decision_tree_svg(tree_nodes, tree_edges)
    (output_dir / "assets/x12798_structural_possibility_decision_tree.svg").write_text(tree_svg)
    html_text = render_html(
        output_dir,
        evidence,
        nodes,
        edges,
        known_neighbors,
        known_nodes,
        known_edges,
        args,
        library_audit,
        same_library_filtered,
        cross_library_matches,
        consensus_svg,
        consensus_df,
        environment_consensus_svg,
        environment_consensus_df,
        motif_df,
        tree_svg,
        tree_summary,
        candidate_consensus,
        final_consensus,
        real_latent,
        v17_latent,
    )
    html_path = output_dir / "index.html"
    html_path.write_text(html_text)

    manifest = {
        "out_dir": str(output_dir),
        "html": str(html_path),
        "control_id": "X-12798",
        "scope": "lib_400_positive",
        "final_status": "provisional_non_pure_structure",
        "n_candidate_nodes": int(len(nodes)),
        "n_pairwise_edges": int(len(edges)),
        "has_confirmed_structure": False,
        "source_artifacts": [
            str(args.resolution_md),
            str(args.encoder_comparison_csv),
            str(args.deep_dive_csv),
            str(args.deep_dive_final_json),
            str(args.lib400_pure_top10),
            str(args.merged_pure_top10),
            str(args.lib400_final_json),
            str(args.merged_final_json),
            str(args.v15_names_csv),
            str(args.library_match_md),
            str(args.consensus_guided_root),
            str(args.real_latent_calibration_root),
            str(args.v17_spectrum_latent_root),
            str(args.reference_db),
            *[str(path) for path in args.metabolon_library_catalogs],
        ],
        "anchor_candidate_name": ANCHOR_CANDIDATE_NAME,
        "anchor_candidate_smiles": reference_smiles or ANCHOR_CANDIDATE_SMILES,
        "all_known_similarity_to_anchor_csv": "data/x12798_all_known_similarity_to_anchor.csv",
        "known_similarity_neighbors_csv": "data/x12798_known_similarity_neighbors.csv",
        "known_similarity_network_nodes_csv": "data/x12798_known_similarity_network_nodes.csv",
        "known_similarity_network_edges_csv": "data/x12798_known_similarity_network_edges.csv",
        "known_similarity_network_json": "data/x12798_known_similarity_network.json",
        "known_similarity_network_asset_kind": "inline_svg",
        "n_known_similarity_nodes": int(len(known_nodes)),
        "n_known_similarity_neighbors": int(len(known_neighbors)),
        "n_all_known_similarity_rows": int(len(all_known)),
        "closest_known_name": clean(known_neighbors.iloc[0]["canonical_name"]) if not known_neighbors.empty else "",
        "closest_known_similarity": (
            float(known_neighbors.iloc[0]["similarity_to_anchor"]) if not known_neighbors.empty else 0.0
        ),
        "known_similarity_rendered_assets": known_rendered_assets,
        "known_similarity_render_failures": known_render_failures,
        "n_filtered_same_library_400": int(len(same_library_filtered)),
        "n_cross_library_metabolon_matches": int(len(cross_library_matches)),
        "consensus_heatmap_asset": "assets/x12798_consensus_structure_heatmap.svg",
        "motif_consensus_csv": "data/x12798_motif_consensus.csv",
        "scaffold_position_consensus_csv": "data/x12798_scaffold_position_consensus_atom_support.csv",
        "scaffold_position_consensus_asset": "assets/x12798_scaffold_position_consensus.svg",
        "scaffold_environment_consensus_csv": "data/x12798_scaffold_environment_consensus_atom_support.csv",
        "scaffold_environment_consensus_asset": "assets/x12798_scaffold_environment_consensus.svg",
        "structural_possibility_tree_nodes_csv": "data/x12798_structural_possibility_tree_nodes.csv",
        "structural_possibility_tree_edges_csv": "data/x12798_structural_possibility_tree_edges.csv",
        "structural_possibility_summary_csv": "data/x12798_structural_possibility_summary.csv",
        "structural_possibility_tree_asset": "assets/x12798_structural_possibility_decision_tree.svg",
        "candidate_motif_flags_csv": "data/x12798_candidate_motif_flags.csv",
        "candidate_specific_consensus_summary_csv": "data/x12798_candidate_specific_consensus_summary.csv",
        "candidate_specific_consensus_atom_support_dir": "data/candidate_specific_consensus_atom_support",
        "candidate_environment_consensus_atom_support_dir": "data/candidate_environment_consensus_atom_support",
        "n_candidate_specific_consensus_assets": int(len(candidate_consensus)),
        "n_candidate_environment_consensus_assets": int(len(candidate_consensus)),
        "final_consensus_guided_visual_summary_csv": (
            final_consensus.get("summary_csv", "") if final_consensus.get("available") else ""
        ),
        "final_consensus_guided_graph_position_csv": (
            final_consensus.get("graph_csv", "") if final_consensus.get("available") else ""
        ),
        "final_consensus_guided_environment_position_csv": (
            final_consensus.get("environment_csv", "") if final_consensus.get("available") else ""
        ),
        "final_consensus_guided_graph_position_asset": (
            final_consensus.get("graph_asset", "") if final_consensus.get("available") else ""
        ),
        "final_consensus_guided_environment_position_asset": (
            final_consensus.get("environment_asset", "") if final_consensus.get("available") else ""
        ),
        "real_latent_calibration_candidate_summary_csv": (
            real_latent.get("summary_csv", "") if real_latent.get("available") else ""
        ),
        "real_latent_calibration_candidate_summary_json": (
            real_latent.get("summary_json", "") if real_latent.get("available") else ""
        ),
        "real_latent_calibration_available": bool(real_latent.get("available")),
        "v17_spectrum_latent_candidate_summary_csv": (
            v17_latent.get("summary_csv", "") if v17_latent.get("available") else ""
        ),
        "v17_tautomer_cluster_summary_csv": (
            v17_latent.get("cluster_csv", "") if v17_latent.get("available") else ""
        ),
        "v17_spectrum_latent_available": bool(v17_latent.get("available")),
        "n_consensus_atoms": int(len(consensus_df)),
        "n_motif_consensus_rows": int(len(motif_df)),
        **asset_manifest,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        manifest = build_report(args)
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        return 1
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
