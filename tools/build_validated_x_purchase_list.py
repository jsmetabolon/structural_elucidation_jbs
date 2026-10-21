#!/usr/bin/env python3
"""Build a validated purchase shortlist for frequently detected X compounds.

The X annotation, predicted structure, and PubChem identity are deliberately
kept separate. Formula agreement is a gate, never proof of structure identity.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors

JBS_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = JBS_ROOT.parent
DEFAULT_RERANK_ROOT = JBS_ROOT / "output/pipeline_equalization_peptide_corrected"
DEFAULT_MATRIX = (
    WORKSPACE_ROOT
    / "structural_elucidation_pathway_runs/docs/paper/data/all_metabolites_matrix_detection.xlsx"
)
DEFAULT_KNOWN_X = JBS_ROOT / "data/known_X.csv"
DEFAULT_STRUCTURES = JBS_ROOT / "data/known_X_all_available_structures_simple.csv"
DEFAULT_CATALOG = (
    WORKSPACE_ROOT
    / "spec2mol_merged_fullrun_metabolon_only/models/chemical_gap_analysis_v1/catalogs/metabolon_catalog.csv"
)
DEFAULT_OUTPUT_DIR = (
    WORKSPACE_ROOT / "spec2mol_merged_fullrun_metabolon_only/models/chemical_gap_analysis_v1"
)
DEFAULT_CACHE = DEFAULT_OUTPUT_DIR / "purchase_vendor_lookup_cache.json"
VENDOR_MODULE = WORKSPACE_ROOT / "spec2mol_merged_fullrun_metabolon_only/build_purchase_candidate_list_v1.py"

FINAL_NAME = "validated_x_purchase_candidates.csv"
AUDIT_NAME = "validated_x_purchase_audit.csv"
SUMMARY_JSON_NAME = "validated_x_purchase_summary.json"
SUMMARY_MD_NAME = "validated_x_purchase_summary.md"
DEPRECATION_NAME = "predicted_x_high_detection_purchasable_rank1_vendor_confirmed.DEPRECATED.md"

FREQUENCY_COLUMNS = ("Human Plasma", "Human Serum", "Human Urine")
RERANK_COLUMNS = {
    "x_compound_id",
    "control_id",
    "feature_label",
    "rank",
    "predicted_smiles",
    "candidate_smiles",
    "smiles",
    "predicted_formula",
    "candidate_formula",
    "target_formula",
    "rerank_score",
    "rerank_score_calibrated",
    "unified_final_score",
    "evidence_confidence",
    "mass_match_score",
    "candidate_source",
    "source_label",
    "library_name",
    "modification",
    "anchor_source",
    "anchor_smiles",
    "is_valid",
    "strict_mass_pass",
}
ARTIFACT_RE = re.compile(
    r"modified[_ -]?smiles|mass[_ -]?match[_ -]?only|formula[_ -]?only|malformed|decoy|isotopologue",
    re.IGNORECASE,
)
GENERIC_ANNOTATION_RE = re.compile(r"^(unknown|x-\d+|unnamed|unidentified)$", re.IGNORECASE)
MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def text(value: Any) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    value = str(value).strip()
    return "" if value.lower() in {"nan", "none", "<na>"} else value


def truthy(value: Any) -> bool:
    return text(value).lower() in {"1", "true", "yes", "y"}


def number(value: Any, default: float = 0.0) -> float:
    try:
        result = float(value)
        return result if math.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def canonicalize(smiles: Any) -> str:
    mol = Chem.MolFromSmiles(text(smiles))
    return Chem.MolToSmiles(mol, isomericSmiles=True) if mol is not None else ""


def formula_from_smiles(smiles: Any) -> str:
    mol = Chem.MolFromSmiles(text(smiles))
    return rdMolDescriptors.CalcMolFormula(mol) if mol is not None else ""


def inchikey_from_smiles(smiles: Any) -> str:
    mol = Chem.MolFromSmiles(text(smiles))
    return Chem.MolToInchiKey(mol) if mol is not None else ""


def tanimoto(smiles_a: Any, smiles_b: Any) -> float:
    a, b = Chem.MolFromSmiles(text(smiles_a)), Chem.MolFromSmiles(text(smiles_b))
    if a is None or b is None:
        return 0.0
    return float(DataStructs.TanimotoSimilarity(MORGAN.GetFingerprint(a), MORGAN.GetFingerprint(b)))


def connectivity_key(inchikey: Any) -> str:
    return text(inchikey).split("-")[0]


def normalize_formula(value: Any) -> str:
    return re.sub(r"\s+", "", text(value))


def annotation_classes(name: Any) -> set[str]:
    name = text(name).lower()
    classes: set[str] = set()
    rules = {
        "amino_acid": r"alanine|glycine|valine|leucine|isoleucine|serine|threonine|tryptophan|"
        r"tyrosine|phenylalanine|glutam|aspart|amino acid|oxindolyl",
        "fatty_acid": r"fatty acid|palmit|stear|oleic|linole|arachidon|enoic acid|dioic acid",
        "peptide": r"dipeptide|tripeptide|tetrapeptide|peptide",
        "steroid": r"steroid|sterol|cholest|androst|pregnen|cortisol|bile acid",
        "sulfate": r"sulfate|sulfonate",
        "phosphate": r"phosphat|nucleotide",
        "carnitine": r"carnitine",
        "purine": r"purine|guanine|guanosine|adenine|adenosine|urate|uric acid|xanthine",
        "sugar": r"glucuron|glycoside|glucoside|galactoside|riboside",
        "oxindole": r"oxindol",
    }
    for label, pattern in rules.items():
        if re.search(pattern, name):
            classes.add(label)
    return classes


def structure_classes(smiles: Any) -> set[str]:
    mol = Chem.MolFromSmiles(text(smiles))
    if mol is None:
        return set()
    classes: set[str] = set()
    patterns = {
        "amino_acid": "[N;!$(N-C=O);!$([N+](C)(C)C)]-[C;X4]-[C](=O)[O;H1,-1]",
        "peptide": "C(=O)N[C;X4]C(=O)",
        "sulfate": "[C,c]-O-S(=O)(=O)[O;H1,-1]",
        "phosphate": "P(=O)(O)O",
        "carnitine": "[N+](C)(C)C",
        "purine": "n1cnc2ncnc12",
        "sugar": "O[C;R][C;R](O)[C;R](O)",
        "oxindole": "O=C1N[c,C]2[c,C][c,C][c,C][c,C]2C1",
    }
    for label, smarts in patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is not None and mol.HasSubstructMatch(pattern):
            classes.add(label)
    carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    aromatic = sum(atom.GetIsAromatic() for atom in mol.GetAtoms())
    carboxyl = Chem.MolFromSmarts("[C](=O)[O;H1,-1]")
    if carbons >= 8 and aromatic <= 2 and carboxyl is not None and mol.HasSubstructMatch(carboxyl):
        classes.add("fatty_acid")
    if carbons >= 17 and mol.GetRingInfo().NumRings() >= 3:
        classes.add("steroid")
    return classes


RESIDUE_SMARTS = {
    "gly": "[NX3][CH2][CX3](=O)",
    "ala": "[NX3][CH]([CH3])[CX3](=O)",
    "pro": "O=[CX3][CH]1CCCN1",
    "asn": "[NX3][CH](C[CX3]([NX3H2])=O)[CX3](=O)",
    "asp": "[NX3][CH](C[CX3](=O)[OX2H1,OX1-])[CX3](=O)",
    "gln": "[NX3][CH](CC[CX3]([NX3H2])=O)[CX3](=O)",
    "glu": "[NX3][CH](CC[CX3](=O)[OX2H1,OX1-])[CX3](=O)",
    "ser": "[NX3][CH]([CH2][OX2H1])[CX3](=O)",
    "thr": "[NX3][CH]([CH1]([CH3])[OX2H1])[CX3](=O)",
    "val": "[NX3][CH]([CH1]([CH3])[CH3])[CX3](=O)",
    "leu": "[NX3][CH]([CH2][CH1]([CH3])[CH3])[CX3](=O)",
    "ile": "[NX3][CH]([CH1]([CH3])[CH2][CH3])[CX3](=O)",
    "phe": "[NX3][CH](Cc1ccccc1)[CX3](=O)",
    "tyr": "[NX3][CH](Cc1ccc([OX2H1])cc1)[CX3](=O)",
    "trp": "[NX3][CH](Cc1c[nH]c2ccccc12)[CX3](=O)",
    "his": "[NX3][CH](Cc1c[nH]cn1)[CX3](=O)",
    "met": "[NX3][CH](CC[SX2][CH3])[CX3](=O)",
    "cys": "[NX3][CH]([CH2][SX2H1])[CX3](=O)",
    "lys": "[NX3][CH](CCC[CH2][NX3])[CX3](=O)",
    "arg": "[NX3][CH](CCCNC(=N)N)[CX3](=O)",
}
RESIDUE_ALIASES = {
    "glycine": "gly",
    "alanine": "ala",
    "proline": "pro",
    "asparagine": "asn",
    "aspartate": "asp",
    "aspartic acid": "asp",
    "glutamine": "gln",
    "glutamate": "glu",
    "glutamic acid": "glu",
    "serine": "ser",
    "threonine": "thr",
    "valine": "val",
    "leucine": "leu",
    "isoleucine": "ile",
    "phenylalanine": "phe",
    "tyrosine": "tyr",
    "tryptophan": "trp",
    "histidine": "his",
    "methionine": "met",
    "cysteine": "cys",
    "lysine": "lys",
    "arginine": "arg",
}
PEPTIDE_LENGTHS = {"dipeptide": 2, "tripeptide": 3, "tetrapeptide": 4}
PEPTIDE_BOND_SMARTS = "[CX3](=O)[NX3][CX4][CX3](=O)"
URONIC_ACID_SMARTS = "[OX2;R]1[CH1]([CX3](=O)[OX2H1,OX1-])[CH1][CH1][CH1][CH1]1"
GLYCINE_CONJUGATE_SMARTS = "[CX3](=O)[NX3][CH2][CX3](=O)[OX2H1,OX1-]"
QUATERNARY_AMMONIUM_SMARTS = "[N+]([CH3])([CH3])[CH3]"
HETEROGENEOUS_ANNOTATION_RE = re.compile(r"\bor\b", re.IGNORECASE)


def has_substructure(smiles: Any, smarts: str) -> bool:
    mol = Chem.MolFromSmiles(text(smiles))
    pattern = Chem.MolFromSmarts(smarts)
    return bool(mol is not None and pattern is not None and mol.HasSubstructMatch(pattern))


def declared_residues(name: Any) -> set[str]:
    lowered = text(name).lower()
    residues = {code for alias, code in RESIDUE_ALIASES.items() if re.search(rf"\b{alias}\b", lowered)}
    residues.update(code for code in RESIDUE_SMARTS if re.search(rf"\b{code}\b", lowered))
    return residues


def residues_present(smiles: Any) -> set[str]:
    return {code for code, smarts in RESIDUE_SMARTS.items() if has_substructure(smiles, smarts)}


def peptide_residue_count(smiles: Any) -> int:
    mol = Chem.MolFromSmiles(text(smiles))
    pattern = Chem.MolFromSmarts(PEPTIDE_BOND_SMARTS)
    if mol is None or pattern is None:
        return 0
    matches = mol.GetSubstructMatches(pattern)
    return len(matches) + 1 if matches else 0


def peptide_consistency(name: Any, smiles: Any) -> tuple[bool | None, str]:
    """Check declared peptide length and named residues against the structure."""
    lowered = text(name).lower()
    declared_length = next((size for label, size in PEPTIDE_LENGTHS.items() if label in lowered), None)
    expected_residues = declared_residues(name)
    if declared_length is None and not expected_residues:
        return None, ""
    if declared_length is not None:
        observed_length = peptide_residue_count(smiles)
        if observed_length and observed_length != declared_length:
            return False, f"peptide_length_mismatch:declared={declared_length};observed={observed_length}"
    if expected_residues:
        missing = sorted(expected_residues - residues_present(smiles))
        if missing:
            return False, "missing_declared_residues:" + "|".join(missing)
    return True, "peptide_composition_match"


def unexpected_motifs(name: Any, smiles: Any) -> list[str]:
    """Flag strong structural motifs the annotation does not account for."""
    lowered = text(name).lower()
    flagged: list[str] = []
    if has_substructure(smiles, QUATERNARY_AMMONIUM_SMARTS) and "carnitine" not in lowered:
        flagged.append("unexpected_quaternary_ammonium")
    glycosidic = re.search(r"glucoside|galactoside|riboside|glycoside", lowered)
    if glycosidic and "glucuron" not in lowered and has_substructure(smiles, URONIC_ACID_SMARTS):
        flagged.append("uronic_acid_but_annotation_says_glycoside")
    if "glucuron" in lowered and not has_substructure(smiles, URONIC_ACID_SMARTS):
        flagged.append("missing_uronic_acid")
    if re.search(r"glycine conjugate|glycine of|\bglycine\b", lowered) and not has_substructure(
        smiles, GLYCINE_CONJUGATE_SMARTS
    ):
        flagged.append("missing_glycine_conjugate")
    return flagged


def class_consistency(name: Any, smiles: Any) -> tuple[bool | None, str]:
    unexpected = unexpected_motifs(name, smiles)
    if unexpected:
        return False, "motif_conflict:" + "|".join(unexpected)
    peptide_ok, peptide_reason = peptide_consistency(name, smiles)
    if peptide_ok is False:
        return False, peptide_reason
    expected = annotation_classes(name)
    if not expected:
        return (None, "no_reliable_annotation_class") if peptide_ok is None else (True, peptide_reason)
    observed = structure_classes(smiles)
    missing = sorted(expected - observed)
    if missing:
        return False, "missing_annotation_class:" + "|".join(missing)
    matched = "annotation_class_match:" + "|".join(sorted(expected))
    return True, f"{matched};{peptide_reason}" if peptide_ok else matched


def candidate_consistency(
    annotation_name: Any,
    annotation_formula: Any,
    annotation_smiles: Any,
    candidate_smiles: Any,
) -> dict[str, Any]:
    canonical = canonicalize(candidate_smiles)
    candidate_formula = formula_from_smiles(canonical)
    if not canonical:
        return {
            "annotation_structure_consistency": "mismatch",
            "consistency_reason": "invalid_candidate_smiles",
            "formula_agreement": False,
            "annotation_tanimoto": 0.0,
            "class_agreement": False,
        }
    formula = normalize_formula(annotation_formula)
    formula_agreement = bool(formula and formula == candidate_formula)
    truth = canonicalize(annotation_smiles)
    similarity = tanimoto(truth, canonical) if truth else 0.0
    undecided_annotation = bool(HETEROGENEOUS_ANNOTATION_RE.search(text(annotation_name)))
    if truth and inchikey_from_smiles(truth) == inchikey_from_smiles(canonical):
        outcome, reason, class_agreement = "exact", "exact_annotation_inchikey", True
    elif truth:
        class_ok, class_reason = class_consistency(annotation_name, canonical)
        class_agreement = class_ok
        if class_ok is not False and undecided_annotation:
            outcome, reason = "ambiguous", "heterogeneous_annotation_alternatives"
        elif formula_agreement and similarity >= 0.35 and class_ok is not False:
            outcome, reason = "plausible_analog", f"truth_tanimoto={similarity:.3f};{class_reason}"
        elif formula_agreement and similarity >= 0.65:
            outcome, reason = "plausible_analog", f"truth_tanimoto={similarity:.3f}"
        else:
            outcome, reason = "mismatch", f"truth_tanimoto={similarity:.3f};{class_reason}"
    else:
        class_ok, class_reason = class_consistency(annotation_name, canonical)
        class_agreement = class_ok
        generic = not text(annotation_name) or bool(GENERIC_ANNOTATION_RE.match(text(annotation_name)))
        if class_ok is not False and undecided_annotation:
            outcome, reason = "ambiguous", "heterogeneous_annotation_alternatives"
        elif formula_agreement and class_ok is True:
            outcome, reason = "plausible_analog", class_reason
        elif generic or class_ok is None:
            outcome, reason = "ambiguous", class_reason
        else:
            outcome, reason = "mismatch", class_reason
    if formula and not formula_agreement:
        outcome, reason = "mismatch", f"formula_mismatch:{formula}!={candidate_formula}"
    return {
        "annotation_structure_consistency": outcome,
        "consistency_reason": reason,
        "formula_agreement": formula_agreement,
        "annotation_tanimoto": round(similarity, 6),
        "class_agreement": class_agreement,
    }


def is_artifact(row: dict[str, Any]) -> tuple[bool, str]:
    fields = " ".join(
        text(row.get(column))
        for column in ("candidate_source", "source_label", "library_name", "modification", "anchor_source")
    )
    match = ARTIFACT_RE.search(fields)
    if match:
        return True, f"artifact_label:{match.group(0).lower().replace(' ', '_')}"
    if "anchor" in fields.lower() and text(row.get("anchor_source")).lower() in {"nan", "none", "malformed"}:
        return True, "malformed_anchor_label"
    if text(row.get("is_valid")) and not truthy(row.get("is_valid")):
        return True, "upstream_invalid"
    if text(row.get("strict_mass_pass")) and not truthy(row.get("strict_mass_pass")):
        return True, "strict_mass_gate_failed"
    return False, ""


def load_detection(path: Path, sheet_name: str, blood_threshold: float, urine_threshold: float) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=sheet_name)
    df = df.rename(
        columns={
            "X Compound ID": "x_compound_id",
            "Name": "x_annotation_name",
            "Top Match Formula": "x_annotation_formula",
            "MSI Sub-level": "msi_sublevel",
            "MSI Level": "msi_level",
            "Human Plasma": "human_plasma",
            "Human Serum": "human_serum",
            "Human Urine": "human_urine",
        }
    )
    for column in ("human_plasma", "human_serum", "human_urine"):
        df[column] = pd.to_numeric(df[column], errors="coerce").fillna(0.0)
        if df[column].max() <= 1.0:
            df[column] *= 100.0
    df["max_blood_detection"] = df[["human_plasma", "human_serum"]].max(axis=1)
    df["urine_detection"] = df["human_urine"]
    df["frequency_eligible"] = (df["max_blood_detection"] >= blood_threshold) | (
        df["urine_detection"] >= urine_threshold
    )
    df["msi_sublevel_1"] = pd.to_numeric(df["msi_sublevel"], errors="coerce").eq(1)
    keep = [
        "x_compound_id",
        "x_annotation_name",
        "x_annotation_formula",
        "msi_level",
        "msi_sublevel",
        "human_plasma",
        "human_serum",
        "human_urine",
        "max_blood_detection",
        "urine_detection",
        "frequency_eligible",
        "msi_sublevel_1",
    ]
    return df[keep].dropna(subset=["x_compound_id"]).drop_duplicates("x_compound_id")


def load_truth(known_x: Path, structures: Path) -> pd.DataFrame:
    known = pd.read_csv(known_x, low_memory=False)
    known = known.rename(
        columns={
            "Unknown": "x_compound_id",
            "KNOWN": "known_name",
            "Tentative information": "tentative_name",
            "Formula": "known_formula",
        }
    )
    known["x_compound_id"] = known["x_compound_id"].map(text)
    known = known.drop_duplicates("x_compound_id")
    available = pd.read_csv(structures)
    available["x_compound_id"] = available["x_compound_id"].map(text)
    available = available.drop_duplicates("x_compound_id").rename(columns={"smiles": "x_annotation_smiles"})
    merged = known.merge(available, how="outer", on="x_compound_id", suffixes=("", "_structure"))
    merged["resolved_known_name"] = merged.get("tentative_name", pd.Series(index=merged.index, dtype=object)).map(text)
    if "known_name" in merged:
        merged["resolved_known_name"] = merged["resolved_known_name"].where(
            merged["resolved_known_name"].ne(""), merged["known_name"].map(text)
        )
    if "known_name_structure" in merged:
        merged["resolved_known_name"] = merged["resolved_known_name"].where(
            merged["resolved_known_name"].ne(""), merged["known_name_structure"].map(text)
        )
    return merged[
        ["x_compound_id", "resolved_known_name", "known_formula", "x_annotation_smiles"]
    ].drop_duplicates("x_compound_id")


def infer_lib_id(path: Path) -> str:
    match = re.search(r"rerank_lib(\d+)", str(path))
    return match.group(1) if match else ""


def load_candidates(rerank_root: Path, top_k: int) -> pd.DataFrame:
    paths = sorted(rerank_root.glob("rerank_lib*/unified_rerank_all_top10_ppm10.csv"))
    if not paths:
        raise FileNotFoundError(f"No rerank files found under {rerank_root}")
    frames: list[pd.DataFrame] = []
    for path in paths:
        frame = pd.read_csv(path, usecols=lambda column: column in RERANK_COLUMNS, low_memory=False)
        frame["recommended_lib_id"] = infer_lib_id(path)
        frames.append(frame)
    df = pd.concat(frames, ignore_index=True)
    if "control_id" in df:
        df["x_compound_id"] = df["control_id"].map(text)
    elif "x_compound_id" in df:
        df["x_compound_id"] = df["x_compound_id"].map(text)
    elif "feature_label" in df:
        df["x_compound_id"] = df["feature_label"].map(
            lambda value: (re.search(r"(X-\d+)", text(value)) or [""])[0]
        )
    else:
        raise ValueError("Rerank inputs do not contain a usable X compound identifier")
    df["rank"] = pd.to_numeric(df.get("rank"), errors="coerce")
    df = df[df["rank"].between(1, top_k)].copy()
    smile_columns = [column for column in ("predicted_smiles", "candidate_smiles", "smiles") if column in df]
    if not smile_columns:
        raise ValueError("Rerank inputs do not contain a candidate SMILES column")
    df["recommended_predicted_smiles"] = ""
    for column in smile_columns:
        df["recommended_predicted_smiles"] = df["recommended_predicted_smiles"].where(
            df["recommended_predicted_smiles"].ne(""), df[column].map(text)
        )
    score_columns = [
        column
        for column in ("rerank_score_calibrated", "rerank_score", "unified_final_score", "evidence_confidence")
        if column in df
    ]
    df["selection_evidence_score"] = 0.0
    for column in score_columns:
        values = pd.to_numeric(df[column], errors="coerce")
        df["selection_evidence_score"] = df["selection_evidence_score"].where(
            df["selection_evidence_score"].ne(0), values.fillna(0.0)
        )
    return df


def load_metabolon_keys(path: Path) -> set[str]:
    catalog = pd.read_csv(path, usecols=lambda column: column in {"inchikey14", "canonical_smiles"}, low_memory=False)
    keys = {connectivity_key(value) for value in catalog.get("inchikey14", []) if text(value)}
    keys.update(
        connectivity_key(inchikey_from_smiles(smiles))
        for smiles in catalog.get("canonical_smiles", [])
        if canonicalize(smiles)
    )
    return {key for key in keys if key}


def load_cache(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text()) if path.exists() else {}
    except (OSError, json.JSONDecodeError):
        return {}


def save_cache(path: Path, cache: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cache, indent=2, sort_keys=True) + "\n")


def load_vendor_fetcher():
    spec = importlib.util.spec_from_file_location("_purchase_vendor_lookup", VENDOR_MODULE)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load vendor lookup module {VENDOR_MODULE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.fetch_pubchem_vendor_record_cached


def normalize_vendor_record(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "pubchem_name": text(
            record.get("pubchem_name") or record.get("pubchem_compound_name") or record.get("compound_name")
        ),
        "pubchem_cid": text(record.get("pubchem_cid")),
        "pubchem_url": text(record.get("pubchem_url")),
        "pubchem_inchikey": text(record.get("pubchem_inchikey")),
        "pubchem_molecular_formula": text(record.get("pubchem_molecular_formula")),
        "vendor_sources": text(record.get("vendor_sources") or record.get("vendor")),
        "vendor_count": int(number(record.get("vendor_count"), 0)),
        "pubchem_vendor_evidence": text(record.get("pubchem_vendor_evidence")),
        "availability_tier": text(record.get("availability_tier")),
    }


def validate_vendor_identity(candidate: dict[str, Any], record: dict[str, Any]) -> tuple[bool, str]:
    if not record["pubchem_cid"]:
        return False, "no_pubchem_cid"
    if not record["pubchem_name"]:
        return False, "blank_pubchem_name"
    candidate_key = text(candidate["recommended_predicted_inchikey"])
    if not record["pubchem_inchikey"]:
        return False, "cached_pubchem_inchikey_missing"
    if candidate_key != record["pubchem_inchikey"]:
        return False, "pubchem_inchikey_mismatch"
    formula = text(candidate["recommended_predicted_formula"])
    if not record["pubchem_molecular_formula"]:
        return False, "cached_pubchem_formula_missing"
    if normalize_formula(formula) != normalize_formula(record["pubchem_molecular_formula"]):
        return False, "pubchem_formula_mismatch"
    if not record["vendor_sources"] or record["vendor_count"] < 1:
        return False, "no_vendor_listing"
    if re.search(r"(?:^|[- ])d\d+\b|isotope|labelled|labeled", record["pubchem_name"], re.IGNORECASE):
        return False, "isotope_labeled_pubchem_identity"
    return True, "exact_pubchem_structure_vendor_listed"


def rank_candidates(
    candidates: pd.DataFrame,
    detection: pd.DataFrame,
    truth: pd.DataFrame,
    metabolon_keys: set[str],
    min_evidence_score: float,
) -> pd.DataFrame:
    work = candidates.merge(detection, on="x_compound_id", how="inner").merge(truth, on="x_compound_id", how="left")
    work = work[work["frequency_eligible"]].copy()
    rows: list[dict[str, Any]] = []
    for raw in work.to_dict("records"):
        row = dict(raw)
        row["x_annotation_name"] = text(row.get("x_annotation_name")) or text(row.get("resolved_known_name"))
        row["x_annotation_formula"] = text(row.get("x_annotation_formula")) or text(row.get("known_formula"))
        smiles = canonicalize(row["recommended_predicted_smiles"])
        row["recommended_predicted_smiles"] = smiles
        row["recommended_predicted_formula"] = formula_from_smiles(smiles)
        row["recommended_predicted_inchikey"] = inchikey_from_smiles(smiles)
        row["recommended_predicted_inchikey14"] = connectivity_key(row["recommended_predicted_inchikey"])
        consistency = candidate_consistency(
            row["x_annotation_name"],
            row["x_annotation_formula"],
            row.get("x_annotation_smiles"),
            smiles,
        )
        row.update(consistency)
        artifact, artifact_reason = is_artifact(row)
        row["artifact_candidate"] = artifact
        row["already_in_metabolon_library"] = row["recommended_predicted_inchikey14"] in metabolon_keys
        reasons: list[str] = []
        if row["msi_sublevel_1"]:
            reasons.append("msi_sublevel_1")
        if artifact:
            reasons.append(artifact_reason)
        if row["annotation_structure_consistency"] not in {"exact", "plausible_analog"}:
            reasons.append(f"annotation_{row['annotation_structure_consistency']}")
        if row["already_in_metabolon_library"]:
            reasons.append("already_in_metabolon_library")
        if number(row["selection_evidence_score"]) < min_evidence_score:
            reasons.append("weak_prediction_score")
        row["pre_vendor_exclusion_reason"] = "|".join(reasons)
        consistency_weight = {"exact": 100.0, "plausible_analog": 60.0, "ambiguous": 10.0, "mismatch": 0.0}
        row["selection_score"] = round(
            consistency_weight[row["annotation_structure_consistency"]]
            + 20.0 * number(row["annotation_tanimoto"])
            + 5.0 * int(bool(row["formula_agreement"]))
            + 10.0 * number(row["selection_evidence_score"])
            + 2.0 / max(number(row.get("rank"), 10.0), 1.0),
            6,
        )
        rows.append(row)
    audit = pd.DataFrame(rows)
    if audit.empty:
        return audit
    duplicate_map = (
        audit[audit["recommended_predicted_inchikey"].ne("")]
        .groupby("recommended_predicted_inchikey")["x_compound_id"]
        .agg(lambda values: "|".join(sorted(set(map(str, values)))))
        .to_dict()
    )
    audit["prediction_reused_by_x_ids"] = audit["recommended_predicted_inchikey"].map(duplicate_map).fillna("")
    audit["prediction_reused"] = audit["prediction_reused_by_x_ids"].str.contains(r"\|", regex=True)
    audit = audit.sort_values(
        ["x_compound_id", "pre_vendor_exclusion_reason", "selection_score", "rank"],
        ascending=[True, True, False, True],
        kind="stable",
    )
    audit["selected_for_vendor_validation"] = False
    eligible = audit[audit["pre_vendor_exclusion_reason"].eq("")].drop_duplicates("x_compound_id")
    audit.loc[eligible.index, "selected_for_vendor_validation"] = True
    alternatives: dict[str, str] = {}
    for x_id, group in audit.groupby("x_compound_id", sort=False):
        payload = [
            {
                "lib": text(item.get("recommended_lib_id")),
                "rank": int(number(item.get("rank"))),
                "smiles": item["recommended_predicted_smiles"],
                "consistency": item["annotation_structure_consistency"],
                "selection_score": item["selection_score"],
                "rejection": item["pre_vendor_exclusion_reason"],
            }
            for item in group.to_dict("records")
        ]
        alternatives[x_id] = json.dumps(payload, separators=(",", ":"))
    audit["per_lib_top_k_alternatives"] = audit["x_compound_id"].map(alternatives)
    audit["selection_reason"] = audit.apply(
        lambda row: (
            f"best_valid_across_libraries;consistency={row['annotation_structure_consistency']};"
            f"lib={row['recommended_lib_id']};rank={int(number(row['rank']))}"
            if row["selected_for_vendor_validation"]
            else row["pre_vendor_exclusion_reason"]
        ),
        axis=1,
    )
    return audit


def apply_vendor_validation(
    audit: pd.DataFrame,
    cache: dict[str, Any],
    vendor_mode: str,
    timeout: float,
    delay_seconds: float,
) -> pd.DataFrame:
    audit = audit.copy()
    for column, default in {
        "pubchem_name": "",
        "pubchem_cid": "",
        "pubchem_url": "",
        "pubchem_inchikey": "",
        "pubchem_molecular_formula": "",
        "vendor_sources": "",
        "vendor_count": 0,
        "pubchem_vendor_evidence": "",
        "availability_tier": "",
        "pubchem_vendor_listed": False,
        "pubchem_identity_validated": False,
        "vendor_validation_reason": "not_selected",
    }.items():
        audit[column] = default
    fetcher = load_vendor_fetcher() if vendor_mode == "network" else None
    selected = audit.index[audit["selected_for_vendor_validation"]]
    for index in selected:
        candidate = audit.loc[index].to_dict()
        full_key = candidate["recommended_predicted_inchikey"]
        key14 = candidate["recommended_predicted_inchikey14"]
        raw = cache.get(full_key) or cache.get(key14) or {}
        if vendor_mode == "network" and (not raw or not raw.get("pubchem_inchikey")):
            raw = fetcher(
                full_key,
                candidate["recommended_predicted_smiles"],
                timeout=timeout,
                delay_seconds=delay_seconds,
                cache=cache,
            )
        record = normalize_vendor_record(dict(raw))
        valid, reason = validate_vendor_identity(candidate, record)
        for column, value in record.items():
            audit.at[index, column] = value
        audit.at[index, "pubchem_vendor_listed"] = bool(record["vendor_count"] and record["vendor_sources"])
        audit.at[index, "pubchem_identity_validated"] = valid
        audit.at[index, "vendor_validation_reason"] = reason
    audit["purchase_recommended"] = audit["selected_for_vendor_validation"] & audit["pubchem_identity_validated"]
    audit["exclusion_reason"] = audit["pre_vendor_exclusion_reason"]
    vendor_rejected = audit["selected_for_vendor_validation"] & ~audit["pubchem_identity_validated"]
    audit.loc[vendor_rejected, "exclusion_reason"] = audit.loc[vendor_rejected, "vendor_validation_reason"]
    not_best = ~audit["selected_for_vendor_validation"] & audit["exclusion_reason"].eq("")
    audit.loc[not_best, "exclusion_reason"] = "lower_ranked_valid_alternative"
    return audit


def build_summary(audit: pd.DataFrame, blood_threshold: float, urine_threshold: float) -> dict[str, Any]:
    eligible_audit = audit[~audit["msi_sublevel_1"]] if not audit.empty else audit
    x_ids = eligible_audit["x_compound_id"].nunique() if not audit.empty else 0
    selected = audit[audit.get("selected_for_vendor_validation", False)] if not audit.empty else audit
    final = audit[audit.get("purchase_recommended", False)] if not audit.empty else audit
    selected_keys = selected.groupby("recommended_predicted_inchikey")["x_compound_id"].nunique()
    duplicated_selected_keys = set(selected_keys[selected_keys.gt(1)].index)
    duplicate_rows = (
        selected[selected["recommended_predicted_inchikey"].isin(duplicated_selected_keys)]
        if not selected.empty
        else selected
    )
    reason_counts = Counter()
    for reasons in audit.get("exclusion_reason", pd.Series(dtype=str)).map(text):
        reason_counts.update(filter(None, reasons.split("|")))
    return {
        "schema_version": 1,
        "thresholds": {
            "blood_plasma_or_serum": blood_threshold,
            "human_urine": urine_threshold,
        },
        "eligible_detected_x": int(x_ids),
        "candidate_rescued_x": int(
            selected[pd.to_numeric(selected.get("rank"), errors="coerce").gt(1)]["x_compound_id"].nunique()
            if not selected.empty
            else 0
        ),
        "vendor_listed_recommendations": int(final["x_compound_id"].nunique()),
        "existing_library_exclusions": int(
            audit[audit.get("already_in_metabolon_library", False)]["x_compound_id"].nunique()
            if not audit.empty
            else 0
        ),
        "mismatch_exclusions": int(
            audit[audit.get("annotation_structure_consistency", "").eq("mismatch")]["x_compound_id"].nunique()
            if not audit.empty
            else 0
        ),
        "unresolved_x": int(x_ids - final["x_compound_id"].nunique()),
        "duplicate_structures": sorted(
            set(duplicate_rows.get("prediction_reused_by_x_ids", pd.Series(dtype=str)).map(text))
        ),
        "exclusion_reason_counts": dict(sorted(reason_counts.items())),
        "known_failure_cases": {
            x_id: sorted(set(group["exclusion_reason"].map(text)))
            for x_id, group in audit[audit["x_compound_id"].isin(
                ["X-11787", "X-12100", "X-11795", "X-17438", "X-18059", "X-18921", "X-23173", "X-24400", "X-25519"]
            )].groupby("x_compound_id")
        },
    }


def write_outputs(audit: pd.DataFrame, output_dir: Path, summary: dict[str, Any]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    ordered = [
        "x_compound_id",
        "x_annotation_name",
        "x_annotation_formula",
        "x_annotation_smiles",
        "msi_level",
        "msi_sublevel",
        "human_plasma",
        "human_serum",
        "human_urine",
        "recommended_predicted_smiles",
        "recommended_predicted_formula",
        "recommended_predicted_inchikey",
        "rank",
        "recommended_lib_id",
        "selection_evidence_score",
        "selection_score",
        "selection_reason",
        "annotation_structure_consistency",
        "consistency_reason",
        "annotation_tanimoto",
        "formula_agreement",
        "pubchem_name",
        "pubchem_cid",
        "pubchem_url",
        "pubchem_inchikey",
        "pubchem_molecular_formula",
        "vendor_sources",
        "vendor_count",
        "pubchem_vendor_evidence",
        "pubchem_vendor_listed",
        "pubchem_identity_validated",
        "prediction_reused",
        "prediction_reused_by_x_ids",
        "already_in_metabolon_library",
        "purchase_recommended",
        "exclusion_reason",
        "per_lib_top_k_alternatives",
    ]
    columns = [column for column in ordered if column in audit] + [
        column for column in audit.columns if column not in ordered
    ]
    audit[columns].to_csv(output_dir / AUDIT_NAME, index=False)
    audit[audit["purchase_recommended"]][columns].sort_values(
        ["max_blood_detection", "urine_detection"], ascending=False
    ).to_csv(output_dir / FINAL_NAME, index=False)
    (output_dir / SUMMARY_JSON_NAME).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    lines = [
        "# Validated X Purchase Candidate Summary",
        "",
        f"- Eligible detected X compounds: {summary['eligible_detected_x']}",
        f"- Rescued from rank >1: {summary['candidate_rescued_x']}",
        f"- Exact-identity vendor-listed recommendations: {summary['vendor_listed_recommendations']}",
        f"- Existing-library exclusions: {summary['existing_library_exclusions']}",
        f"- Annotation mismatch exclusions: {summary['mismatch_exclusions']}",
        f"- Unresolved X compounds: {summary['unresolved_x']}",
        "",
        "Vendor listing confirms availability only after PubChem InChIKey and formula match the predicted structure.",
    ]
    (output_dir / SUMMARY_MD_NAME).write_text("\n".join(lines) + "\n")
    (output_dir / DEPRECATION_NAME).write_text(
        "# Deprecated output\n\n"
        "`predicted_x_high_detection_purchasable_rank1_vendor_confirmed.csv` is retained only for "
        "historical traceability. It conflates annotation and purchase identities and must not be used "
        f"for purchasing. Use `{FINAL_NAME}` and `{AUDIT_NAME}` instead.\n"
    )


def build(args: argparse.Namespace) -> tuple[pd.DataFrame, dict[str, Any]]:
    detection = load_detection(
        args.matrix,
        args.matrix_sheet,
        args.blood_threshold,
        args.urine_threshold,
    )
    truth = load_truth(args.known_x, args.available_structures)
    candidates = load_candidates(args.rerank_root, args.top_k)
    metabolon_keys = load_metabolon_keys(args.metabolon_catalog)
    audit = rank_candidates(candidates, detection, truth, metabolon_keys, args.min_evidence_score)
    cache = load_cache(args.vendor_cache)
    audit = apply_vendor_validation(audit, cache, args.vendor_mode, args.timeout, args.delay_seconds)
    if args.vendor_mode == "network":
        save_cache(args.vendor_cache, cache)
    summary = build_summary(audit, args.blood_threshold, args.urine_threshold)
    write_outputs(audit, args.output_dir, summary)
    return audit, summary


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rerank-root", type=Path, default=DEFAULT_RERANK_ROOT)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--matrix-sheet", default="Latest_Unique_X")
    parser.add_argument("--known-x", type=Path, default=DEFAULT_KNOWN_X)
    parser.add_argument("--available-structures", type=Path, default=DEFAULT_STRUCTURES)
    parser.add_argument("--metabolon-catalog", type=Path, default=DEFAULT_CATALOG)
    parser.add_argument("--vendor-cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--blood-threshold", type=float, default=10.0)
    parser.add_argument("--urine-threshold", type=float, default=10.0)
    parser.add_argument("--min-evidence-score", type=float, default=0.05)
    parser.add_argument("--vendor-mode", choices=("cached", "network"), default="cached")
    parser.add_argument("--timeout", type=float, default=20.0)
    parser.add_argument("--delay-seconds", type=float, default=0.2)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    _, summary = build(args)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
