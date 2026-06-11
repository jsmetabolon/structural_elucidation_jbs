#!/usr/bin/env python3
"""Assign Metabolon pathway labels to published dual-lane structure predictions."""

from __future__ import annotations

import argparse
import json
import math
import re
import sqlite3
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
from pandas.errors import EmptyDataError

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem, Descriptors

    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

from phase0_utils import canonicalize_smiles, extract_xid, read_csv_clean

try:
    from lactoyl_amino_acid_v1_utils import THEORETICAL_N_LACTOYL_SMILES, infer_lactoyl_amino_acid
except ImportError:
    THEORETICAL_N_LACTOYL_SMILES = {}

    def infer_lactoyl_amino_acid(value: Any) -> str:
        return ""


DEFAULT_DELIVERY_MANIFEST = Path("models/x209_full_list_dual_lane_manifest.json")
DEFAULT_REFERENCE_DB = Path("Metabolite Synonym Database/metabolite_synonyms_merged_with_mtx_giga.sqlite")
DEFAULT_OUTPUT_DIR = Path("models/x209_predicted_pathway_annotations")
DEFAULT_KNOWN_X_CSV = Path("known_X.csv")
DEFAULT_KNOWN_STRUCTURES_CSV = Path("known_X_all_available_structures_simple.csv")
DEFAULT_POWER_SWEEP_CSV = Path("models/master_209_pathway_predictions/master_209_pathway_predictions.csv")
DEFAULT_TOP_STRUCTURES_PER_FEATURE = 10
DEFAULT_TOP_K_NEIGHBORS = 10
DEFAULT_NEIGHBOR_SIMILARITY_FLOOR = 0.40
DEFAULT_VOTE_SIMILARITY_POWER = 3.0
POWER_SWEEP_VALUES = (2.0, 3.0, 4.0, 5.0)
POWER_SWEEP_MIN_SAMPLES = 20
LOW_SIMILARITY_PENALTY_THRESHOLD = 0.50
LOW_SIMILARITY_PENALTY_FLOOR = 0.10
MW_PRIOR_STD_MULTIPLIER = 2.0
MW_PRIOR_FLOOR = 0.30
SUBPATHWAY_TIE_MARGIN = 0.08
SUBPATHWAY_NAME_BOOST = 0.15
XENOBIOTIC_CLASSIFIER_THRESHOLD = 0.80
XENOBIOTIC_BOOST_FACTOR = 1.35
ENDOGENOUS_PENALTY_FACTOR = 0.70
MIN_POWER_TUNING_IMPROVEMENT = 0.02
DEFAULT_REGRESSION_BASELINE_JSON = Path("models/x209_predicted_pathway_annotations/baseline_accuracy_before_changes.json")
DEFAULT_REGRESSION_FAIL_DROP_PP = 0.5
OFFICIAL_LANES = ("coverage", "hitmax")
SIMILARITY_BUCKET_ORDER = (">=0.90", "0.80-0.899", "0.70-0.799", "<0.70", "missing")
ENDOGENOUS_SUPER_PATHWAYS = {
    "Amino Acid",
    "Carbohydrate",
    "Cofactors and Vitamins",
    "Energy",
    "Lipid",
    "Nucleotide",
    "Partially Characterized Molecules",
    "Peptide",
}
TOKEN_STOPWORDS = {
    "acid",
    "and",
    "beta",
    "compound",
    "containing",
    "derivative",
    "drug",
    "ester",
    "form",
    "isomer",
    "metabolism",
    "mono",
    "none",
    "other",
    "pathway",
    "with",
}
PEPTIDE_LIKE_SUBPATHWAYS = {
    "Dipeptide",
    "Gamma-glutamyl Amino Acid",
    "Polypeptide",
}
METHIONINE_CYSTEINE_TAURINE_SUBPATHWAY = "Methionine, Cysteine, SAM and Taurine Metabolism"
SULFUR_AMINO_NAME_TOKENS = (
    "taurine",
    "cysteine",
    "cystine",
    "cysteinyl",
    "cysteamine",
    "homocysteine",
    "homocystine",
    "methionine",
    "thiazolidine",
    "adenosylcysteine",
    "s-adenosyl",
    "sam ",
    "sam-",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Assign pathway and sub-pathway labels to the published dual-lane structure predictions "
            "using only existing Metabolon library metadata."
        )
    )
    parser.add_argument(
        "--delivery-manifest",
        type=Path,
        default=DEFAULT_DELIVERY_MANIFEST,
        help="Published dual-lane delivery manifest with ranked policy row paths.",
    )
    parser.add_argument(
        "--reference-db",
        type=Path,
        default=DEFAULT_REFERENCE_DB,
        help="SQLite synonym database containing Metabolon metadata columns.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory to write the pathway assignment bundle.",
    )
    parser.add_argument(
        "--top-structures-per-feature",
        type=int,
        default=DEFAULT_TOP_STRUCTURES_PER_FEATURE,
        help="Number of ranked predicted structures to keep per feature per lane.",
    )
    parser.add_argument(
        "--top-k-neighbors",
        type=int,
        default=DEFAULT_TOP_K_NEIGHBORS,
        help="Number of nearest reference metabolites to keep per predicted structure.",
    )
    parser.add_argument(
        "--neighbor-similarity-floor",
        type=float,
        default=DEFAULT_NEIGHBOR_SIMILARITY_FLOOR,
        help="Similarity floor used only to label strong vs weak evidence rows.",
    )
    parser.add_argument(
        "--vote-similarity-power",
        type=float,
        default=DEFAULT_VOTE_SIMILARITY_POWER,
        help="Exponent applied to similarity in vote weighting (vote = structure_weight * similarity^power).",
    )
    parser.add_argument(
        "--evaluate",
        action="store_true",
        help="Evaluate pathway/subpathway accuracy against known_X.csv labels.",
    )
    parser.add_argument(
        "--known-x-csv",
        type=Path,
        default=DEFAULT_KNOWN_X_CSV,
        help="Path to known_X.csv file used for evaluation.",
    )
    parser.add_argument(
        "--known-structures-csv",
        type=Path,
        default=DEFAULT_KNOWN_STRUCTURES_CSV,
        help="Path to known_X_all_available_structures_simple.csv used for evaluation diagnostics.",
    )
    parser.add_argument(
        "--power-sweep-csv",
        type=Path,
        default=DEFAULT_POWER_SWEEP_CSV,
        help=(
            "Optional labeled CSV used to tune vote_similarity_power by super-pathway "
            "(expects master_canonical_smiles and existing_pathway columns)."
        ),
    )
    parser.add_argument(
        "--disable-pathway-power-sweep",
        action="store_true",
        help="Disable per-pathway vote_similarity_power tuning and use a single global exponent.",
    )
    parser.add_argument(
        "--regression-check",
        action="store_true",
        help="When used with --evaluate, compare current accuracy vs a baseline JSON and fail on regressions.",
    )
    parser.add_argument(
        "--regression-baseline-json",
        type=Path,
        default=DEFAULT_REGRESSION_BASELINE_JSON,
        help="Baseline JSON used by --regression-check.",
    )
    parser.add_argument(
        "--regression-fail-drop-pp",
        type=float,
        default=DEFAULT_REGRESSION_FAIL_DROP_PP,
        help="Fail threshold in percentage points for --regression-check.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def read_csv_or_empty(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return read_csv_clean(path)
    except EmptyDataError:
        return pd.DataFrame()


def parse_truthy(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False
    text = str(value).strip().lower()
    return text in {"1", "true", "t", "yes", "y"}


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def coalesce_smiles(row: pd.Series) -> str:
    for column in ("smiles", "predicted_smiles"):
        if column not in row.index:
            continue
        value = clean_text(row[column])
        if value:
            return value
    return ""


def safe_numeric(value: object, default: float = 0.0) -> float:
    numeric = pd.to_numeric(value, errors="coerce")
    if pd.isna(numeric):
        return float(default)
    return float(numeric)


def resolve_lane_inputs(manifest_payload: dict[str, Any]) -> dict[str, dict[str, Path | str]]:
    official_lanes = manifest_payload.get("official_lanes")
    if not isinstance(official_lanes, dict):
        raise KeyError("Delivery manifest missing 'official_lanes'.")

    resolved: dict[str, dict[str, Path | str]] = {}
    for lane in OFFICIAL_LANES:
        lane_payload = official_lanes.get(lane)
        if not isinstance(lane_payload, dict):
            raise KeyError(f"Delivery manifest missing official lane '{lane}'.")
        ranked_csv = lane_payload.get("ranked_policy_rows_csv")
        if not ranked_csv:
            raise KeyError(f"Official lane '{lane}' missing ranked_policy_rows_csv.")
        resolved[lane] = {
            "policy": str(lane_payload.get("policy", "")),
            "ranked_policy_rows_csv": Path(str(ranked_csv)).expanduser().resolve(),
        }
    return resolved


@lru_cache(maxsize=65536)
def morgan_fingerprint(smiles: str):
    if not HAS_RDKIT:
        return None
    normalized = str(smiles).strip()
    if not normalized:
        return None
    mol = Chem.MolFromSmiles(normalized)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


@lru_cache(maxsize=65536)
def molecular_weight(smiles: str) -> float | None:
    if not HAS_RDKIT:
        return None
    normalized = str(smiles).strip()
    if not normalized:
        return None
    mol = Chem.MolFromSmiles(normalized)
    if mol is None:
        return None
    return float(Descriptors.MolWt(mol))


@lru_cache(maxsize=65536)
def molecular_descriptor_tuple(smiles: str) -> tuple[float, float, float, float, float] | None:
    if not HAS_RDKIT:
        return None
    normalized = str(smiles).strip()
    if not normalized:
        return None
    mol = Chem.MolFromSmiles(normalized)
    if mol is None:
        return None
    atoms = list(mol.GetAtoms())
    atom_count = float(len(atoms))
    if atom_count <= 0:
        return None
    hetero_count = float(sum(1 for atom in atoms if atom.GetAtomicNum() not in {1, 6}))
    halogen_count = float(sum(1 for atom in atoms if atom.GetAtomicNum() in {9, 17, 35, 53}))
    ring_count = float(mol.GetRingInfo().NumRings())
    sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)N")
    has_sulfonamide = 1.0 if mol.HasSubstructMatch(sulfonamide_pattern) else 0.0
    return (
        float(Descriptors.MolWt(mol)),
        ring_count,
        halogen_count,
        has_sulfonamide,
        hetero_count / atom_count,
    )


def build_pathway_mw_stats(reference_df: pd.DataFrame) -> dict[str, dict[str, float]]:
    if "metabolon_super" not in reference_df.columns or "molecular_weight" not in reference_df.columns:
        return {}
    stats: dict[str, dict[str, float]] = {}
    for pathway, group in reference_df.groupby("metabolon_super", sort=True):
        pathway = clean_text(pathway)
        if not pathway:
            continue
        values = pd.to_numeric(group["molecular_weight"], errors="coerce").dropna()
        if values.empty:
            continue
        mean_value = float(values.mean())
        std_value = float(values.std(ddof=0) if len(values) > 1 else 0.0)
        if not math.isfinite(std_value):
            std_value = 0.0
        spread = max(1.0, std_value * MW_PRIOR_STD_MULTIPLIER)
        lower = mean_value - (MW_PRIOR_STD_MULTIPLIER * std_value)
        upper = mean_value + (MW_PRIOR_STD_MULTIPLIER * std_value)
        stats[pathway] = {
            "mean": mean_value,
            "std": std_value,
            "lower": lower,
            "upper": upper,
            "spread": spread,
            "count": float(len(values)),
        }
    return stats


def compute_pathway_mw_prior(
    predicted_mw: float | None,
    metabolon_super: str,
    pathway_mw_stats: dict[str, dict[str, float]],
) -> float:
    if predicted_mw is None or not math.isfinite(float(predicted_mw)):
        return 1.0
    pathway = clean_text(metabolon_super)
    stats = pathway_mw_stats.get(pathway)
    if not stats:
        return 1.0
    lower = float(stats.get("lower", float("-inf")))
    upper = float(stats.get("upper", float("inf")))
    if lower <= float(predicted_mw) <= upper:
        return 1.0
    if float(predicted_mw) < lower:
        distance = lower - float(predicted_mw)
    else:
        distance = float(predicted_mw) - upper
    spread = max(1.0, float(stats.get("spread", 1.0)))
    prior = math.exp(-distance / spread)
    return max(MW_PRIOR_FLOOR, min(1.0, float(prior)))


def tokenize_name(text: str) -> set[str]:
    normalized = clean_text(text).lower()
    if not normalized:
        return set()
    raw_tokens = re.findall(r"[a-z0-9][a-z0-9'\-/]{2,}", normalized)
    tokens: set[str] = set()
    for token in raw_tokens:
        token = token.strip("-_/'")
        if not token or token in TOKEN_STOPWORDS:
            continue
        if token.isdigit():
            continue
        tokens.add(token)
    return tokens


def build_subpathway_name_prior_map(
    reference_df: pd.DataFrame,
    min_count: int = 4,
    min_precision: float = 0.90,
) -> dict[tuple[str, str], dict[str, float | str]]:
    if "canonical_name" not in reference_df.columns:
        return {}
    token_counts: dict[tuple[str, str], Counter[str]] = {}
    for row in reference_df[["canonical_name", "metabolon_super", "metabolon_sub"]].itertuples(index=False):
        pathway = clean_text(getattr(row, "metabolon_super", ""))
        subpathway = clean_text(getattr(row, "metabolon_sub", ""))
        if not pathway or not subpathway:
            continue
        tokens = tokenize_name(clean_text(getattr(row, "canonical_name", "")))
        for token in tokens:
            key = (pathway, token)
            if key not in token_counts:
                token_counts[key] = Counter()
            token_counts[key][subpathway] += 1

    prior_map: dict[tuple[str, str], dict[str, float | str]] = {}
    for key, counter in token_counts.items():
        total = int(sum(counter.values()))
        if total < int(min_count):
            continue
        subpathway, count = sorted(counter.items(), key=lambda item: (-item[1], item[0]))[0]
        precision = float(count) / float(total)
        if precision < float(min_precision):
            continue
        prior_map[key] = {
            "subpathway": str(subpathway),
            "precision": precision,
            "count": float(total),
        }
    return prior_map


def subpathway_name_boost_scores(
    query_name: str,
    winning_pathway: str,
    subpathway_name_priors: dict[tuple[str, str], dict[str, float | str]],
) -> dict[str, float]:
    boosts: dict[str, float] = {}
    pathway = clean_text(winning_pathway)
    if not pathway or not subpathway_name_priors:
        return boosts
    for token in tokenize_name(query_name):
        payload = subpathway_name_priors.get((pathway, token))
        if not payload:
            continue
        subpathway = clean_text(payload.get("subpathway", ""))
        precision = float(payload.get("precision", 0.0))
        if not subpathway:
            continue
        boosts[subpathway] = boosts.get(subpathway, 0.0) + precision
    return boosts


def apply_subpathway_name_boost(
    subpathway_votes: pd.DataFrame,
    query_name: str,
    winning_pathway: str,
    winning_pathway_vote_weight: float,
    subpathway_name_priors: dict[tuple[str, str], dict[str, float | str]],
) -> tuple[pd.DataFrame, bool]:
    if subpathway_votes.empty:
        return subpathway_votes, False
    if not query_name:
        return subpathway_votes, False
    scores = subpathway_name_boost_scores(
        query_name=query_name,
        winning_pathway=winning_pathway,
        subpathway_name_priors=subpathway_name_priors,
    )
    if not scores:
        return subpathway_votes, False
    adjusted = subpathway_votes.copy()
    applied = False
    for idx, row in adjusted.iterrows():
        subpathway = clean_text(row["matched_metabolon_sub"])
        if subpathway not in scores:
            continue
        raw_boost = float(scores[subpathway]) * SUBPATHWAY_NAME_BOOST * float(winning_pathway_vote_weight)
        capped_boost = min(0.30 * float(winning_pathway_vote_weight), raw_boost)
        adjusted.at[idx, "vote_weight"] = float(adjusted.at[idx, "vote_weight"]) + capped_boost
        applied = True
    if not applied:
        return subpathway_votes, False
    adjusted = adjusted.sort_values(["vote_weight", "matched_metabolon_sub"], ascending=[False, True]).reset_index(
        drop=True
    )
    return adjusted, True


def train_xenobiotic_classifier(reference_df: pd.DataFrame) -> dict[str, Any] | None:
    feature_rows: list[dict[str, float]] = []
    for row in reference_df[["canonical_smiles", "metabolon_super"]].itertuples(index=False):
        smiles = clean_text(getattr(row, "canonical_smiles", ""))
        descriptors = molecular_descriptor_tuple(smiles)
        if descriptors is None:
            continue
        feature_rows.append(
            {
                "molecular_weight": float(descriptors[0]),
                "ring_count": float(descriptors[1]),
                "halogen_count": float(descriptors[2]),
                "sulfonamide": float(descriptors[3]),
                "hetero_ratio": float(descriptors[4]),
                "is_xenobiotic": 1.0 if clean_text(getattr(row, "metabolon_super", "")) == "Xenobiotics" else 0.0,
            }
        )
    if not feature_rows:
        return None
    frame = pd.DataFrame(feature_rows)
    positives = frame.loc[frame["is_xenobiotic"].eq(1.0)].copy()
    negatives = frame.loc[frame["is_xenobiotic"].eq(0.0)].copy()
    if positives.empty or negatives.empty:
        return None

    model: dict[str, Any] = {
        "features": ["molecular_weight", "ring_count", "halogen_count", "sulfonamide", "hetero_ratio"],
        "priors": {
            "xenobiotic": float(len(positives)) / float(len(frame)),
            "non_xenobiotic": float(len(negatives)) / float(len(frame)),
        },
        "stats": {},
    }
    for feature in model["features"]:
        pos_series = positives[feature].astype(float)
        neg_series = negatives[feature].astype(float)
        model["stats"][feature] = {
            "xenobiotic_mean": float(pos_series.mean()),
            "xenobiotic_std": float(pos_series.std(ddof=0) if len(pos_series) > 1 else 1.0),
            "non_xenobiotic_mean": float(neg_series.mean()),
            "non_xenobiotic_std": float(neg_series.std(ddof=0) if len(neg_series) > 1 else 1.0),
        }
    return model


def gaussian_logpdf(value: float, mean: float, std: float) -> float:
    bounded_std = max(float(std), 1e-3)
    variance = bounded_std * bounded_std
    return -0.5 * math.log(2.0 * math.pi * variance) - ((float(value) - float(mean)) ** 2) / (2.0 * variance)


def predict_xenobiotic_probability(smiles: str, classifier_model: dict[str, Any] | None) -> float:
    if classifier_model is None:
        return 0.5
    descriptors = molecular_descriptor_tuple(smiles)
    if descriptors is None:
        return 0.5
    descriptor_map = {
        "molecular_weight": float(descriptors[0]),
        "ring_count": float(descriptors[1]),
        "halogen_count": float(descriptors[2]),
        "sulfonamide": float(descriptors[3]),
        "hetero_ratio": float(descriptors[4]),
    }
    priors = classifier_model.get("priors", {})
    stats = classifier_model.get("stats", {})
    xen_prior = max(float(priors.get("xenobiotic", 0.5)), 1e-6)
    non_prior = max(float(priors.get("non_xenobiotic", 0.5)), 1e-6)
    xen_log = math.log(xen_prior)
    non_log = math.log(non_prior)
    for feature in classifier_model.get("features", []):
        value = float(descriptor_map.get(feature, 0.0))
        feature_stats = stats.get(feature, {})
        xen_log += gaussian_logpdf(
            value,
            float(feature_stats.get("xenobiotic_mean", 0.0)),
            float(feature_stats.get("xenobiotic_std", 1.0)),
        )
        non_log += gaussian_logpdf(
            value,
            float(feature_stats.get("non_xenobiotic_mean", 0.0)),
            float(feature_stats.get("non_xenobiotic_std", 1.0)),
        )
    max_log = max(xen_log, non_log)
    xen_prob = math.exp(xen_log - max_log)
    non_prob = math.exp(non_log - max_log)
    denom = xen_prob + non_prob
    if denom <= 0:
        return 0.5
    return float(xen_prob / denom)


def aggregate_feature_xenobiotic_probability(
    candidate_group: pd.DataFrame,
    classifier_model: dict[str, Any] | None,
) -> float:
    if classifier_model is None or candidate_group.empty:
        return 0.5
    weighted_sum = 0.0
    total_weight = 0.0
    for row in candidate_group[["canonical_smiles", "structure_weight"]].itertuples(index=False):
        smiles = clean_text(getattr(row, "canonical_smiles", ""))
        if not smiles:
            continue
        structure_weight = safe_numeric(getattr(row, "structure_weight", 0.0), default=0.0)
        if structure_weight <= 0:
            continue
        probability = predict_xenobiotic_probability(smiles, classifier_model)
        weighted_sum += structure_weight * probability
        total_weight += structure_weight
    if total_weight <= 0:
        return 0.5
    return weighted_sum / total_weight


def apply_xenobiotic_vote_adjustment(
    pathway_votes: pd.DataFrame,
    xenobiotic_probability: float,
) -> tuple[pd.DataFrame, bool]:
    if pathway_votes.empty:
        return pathway_votes, False
    if xenobiotic_probability < XENOBIOTIC_CLASSIFIER_THRESHOLD:
        return pathway_votes, False
    adjusted = pathway_votes.copy()
    xen_mask = adjusted["matched_metabolon_super"].astype(str).eq("Xenobiotics")
    if not bool(xen_mask.any()):
        return pathway_votes, False
    endo_mask = adjusted["matched_metabolon_super"].map(
        lambda pathway: clean_text(pathway) in ENDOGENOUS_SUPER_PATHWAYS
    )
    if not bool(endo_mask.any()):
        return pathway_votes, False
    adjusted.loc[endo_mask, "vote_weight"] *= ENDOGENOUS_PENALTY_FACTOR
    adjusted.loc[xen_mask, "vote_weight"] *= XENOBIOTIC_BOOST_FACTOR
    adjusted = adjusted.sort_values(["vote_weight", "matched_metabolon_super"], ascending=[False, True]).reset_index(
        drop=True
    )
    return adjusted, True


def choose_query_name(candidate_group: pd.DataFrame, fallback_name: str = "") -> str:
    if "query_name" in candidate_group.columns:
        for value in candidate_group["query_name"].tolist():
            text = clean_text(value)
            if text:
                return text
    return clean_text(fallback_name)


@lru_cache(maxsize=128)
def compiled_smarts(pattern: str):
    if not HAS_RDKIT:
        return None
    return Chem.MolFromSmarts(pattern)


def smiles_has_smarts(smiles: str, pattern: str) -> bool:
    if not HAS_RDKIT:
        return False
    normalized = clean_text(smiles)
    if not normalized:
        return False
    mol = Chem.MolFromSmiles(normalized)
    query = compiled_smarts(pattern)
    if mol is None or query is None:
        return False
    return bool(mol.HasSubstructMatch(query))


@lru_cache(maxsize=1)
def lactoyl_smiles_to_aa_map() -> dict[str, str]:
    out: dict[str, str] = {}
    for aa, smiles in THEORETICAL_N_LACTOYL_SMILES.items():
        canonical = canonicalize_smiles(clean_text(smiles))
        if canonical:
            out[canonical] = str(aa)
    return out


def infer_lactoyl_amino_acid_from_smiles(smiles: str) -> str:
    canonical = canonicalize_smiles(clean_text(smiles))
    if not canonical:
        return ""
    return clean_text(lactoyl_smiles_to_aa_map().get(canonical, ""))


def infer_sulfur_amino_subpathway(
    *,
    smiles: str = "",
    query_name: str = "",
    candidate_name: str = "",
) -> dict[str, str]:
    combined_text = " ".join(
        part.lower() for part in [clean_text(query_name), clean_text(candidate_name)] if clean_text(part)
    )
    matched_name = any(token in combined_text for token in SULFUR_AMINO_NAME_TOKENS)
    canonical = canonicalize_smiles(clean_text(smiles)) or clean_text(smiles)
    mw = molecular_weight(canonical)

    taurine_like = smiles_has_smarts(canonical, "[NX3,NX4+][CH2][CH2]S(=O)(=O)[O-,O]")
    cysteine_like = smiles_has_smarts(canonical, "[NX3,NX4+][CH]([CH2]S)[C](=O)[O-,O]")
    methionine_like = smiles_has_smarts(canonical, "[NX3,NX4+][CH]([CH2][CH2]S[#6,#8,#16])[C](=O)[O-,O]")
    thiazolidine_like = smiles_has_smarts(canonical, "O=C([O-,O])C1CSCN1")
    adenosyl_homocysteine_like = "n1cnc2" in canonical.lower() and "sccc(" in canonical.lower() and "c(=o)o" in canonical.lower()

    short_chain_taurine_like = taurine_like and (mw is not None and float(mw) <= 260.0)
    sulfur_like = any(
        [
            matched_name,
            short_chain_taurine_like,
            cysteine_like,
            methionine_like,
            thiazolidine_like,
            adenosyl_homocysteine_like,
        ]
    )
    if not sulfur_like:
        return {}
    return {
        "predicted_pathway": "Amino Acid",
        "predicted_subpathway": METHIONINE_CYSTEINE_TAURINE_SUBPATHWAY,
    }


def choose_representative_smiles(candidate_group: pd.DataFrame, fallback_smiles: str = "") -> str:
    if not candidate_group.empty:
        ordered = candidate_group.copy()
        for column in ("structure_weight", "policy_rank", "canonical_smiles"):
            if column not in ordered.columns:
                ordered[column] = ""
        ordered["structure_weight"] = pd.to_numeric(ordered["structure_weight"], errors="coerce").fillna(0.0)
        ordered["policy_rank"] = pd.to_numeric(ordered["policy_rank"], errors="coerce").fillna(9999.0)
        ordered = ordered.sort_values(
            ["structure_weight", "policy_rank", "canonical_smiles"],
            ascending=[False, True, True],
            kind="mergesort",
        )
        for value in ordered["canonical_smiles"].tolist():
            text = clean_text(value)
            if text:
                return text
    return clean_text(fallback_smiles)


def infer_pathway_override(
    *,
    smiles: str = "",
    query_name: str = "",
    candidate_name: str = "",
    current_pathway: str = "",
    current_subpathway: str = "",
) -> dict[str, str]:
    sulfur_override = infer_sulfur_amino_subpathway(
        smiles=smiles,
        query_name=query_name,
        candidate_name=candidate_name,
    )
    if sulfur_override:
        return sulfur_override
    aa = (
        infer_lactoyl_amino_acid(query_name)
        or infer_lactoyl_amino_acid(candidate_name)
        or infer_lactoyl_amino_acid_from_smiles(smiles)
    )
    if not aa:
        return {}
    pathway = clean_text(current_pathway)
    subpathway = clean_text(current_subpathway)
    override: dict[str, str] = {}
    if pathway != "Amino Acid":
        override["predicted_pathway"] = "Amino Acid"
    if (
        subpathway == ""
        or subpathway in PEPTIDE_LIKE_SUBPATHWAYS
        or pathway in {"", "Peptide"}
    ):
        override["predicted_subpathway"] = "Lactoyl Amino Acid"
    return override


def predict_pathway_from_neighbors(neighbor_rows: pd.DataFrame, similarity_power: float) -> str:
    if neighbor_rows.empty:
        return ""
    power = float(similarity_power)
    votes = (
        neighbor_rows.assign(vote_weight=neighbor_rows["matched_similarity"].astype(float).pow(power))
        .groupby("matched_metabolon_super", as_index=False)["vote_weight"]
        .sum()
        .sort_values(["vote_weight", "matched_metabolon_super"], ascending=[False, True])
        .reset_index(drop=True)
    )
    if votes.empty:
        return ""
    return clean_text(votes.iloc[0]["matched_metabolon_super"])


def load_power_sweep_labels(power_sweep_csv: Path) -> pd.DataFrame:
    sweep_path = power_sweep_csv.expanduser().resolve()
    if not sweep_path.exists():
        return pd.DataFrame(columns=["canonical_smiles", "true_pathway"])
    try:
        labeled_df = read_csv_clean(sweep_path, dtype=str).fillna("")
    except EmptyDataError:
        return pd.DataFrame(columns=["canonical_smiles", "true_pathway"])

    required = {"master_canonical_smiles", "existing_pathway"}
    if not required.issubset(set(labeled_df.columns)):
        return pd.DataFrame(columns=["canonical_smiles", "true_pathway"])
    out = labeled_df.rename(
        columns={
            "master_canonical_smiles": "canonical_smiles",
            "existing_pathway": "true_pathway",
        }
    )[["canonical_smiles", "true_pathway"]].copy()
    out["canonical_smiles"] = out["canonical_smiles"].map(clean_text)
    out["true_pathway"] = out["true_pathway"].map(clean_text)
    out = out.loc[out["canonical_smiles"].ne("") & out["true_pathway"].ne("")].drop_duplicates().reset_index(drop=True)
    return out


def tune_similarity_power_by_pathway(
    reference_df: pd.DataFrame,
    *,
    top_k_neighbors: int,
    power_sweep_csv: Path | None,
) -> tuple[dict[str, float], dict[str, Any]]:
    if power_sweep_csv is None:
        return {}, {"enabled": False, "reason": "no_calibration_csv"}
    labeled_df = load_power_sweep_labels(power_sweep_csv)
    if labeled_df.empty:
        return {}, {"enabled": False, "reason": "missing_or_empty_calibration_labels"}
    available_pathways = {clean_text(value) for value in reference_df["metabolon_super"].tolist() if clean_text(value)}
    label_pathways = {clean_text(value) for value in labeled_df["true_pathway"].tolist() if clean_text(value)}
    if not available_pathways.intersection(label_pathways):
        return {}, {"enabled": False, "reason": "no_pathway_overlap_between_labels_and_reference"}

    neighbor_df = build_neighbor_table(
        predicted_smiles=labeled_df["canonical_smiles"].tolist(),
        reference_df=reference_df,
        top_k_neighbors=max(int(top_k_neighbors), 1),
        similarity_floor=0.0,
    )
    if neighbor_df.empty:
        return {}, {"enabled": False, "reason": "no_neighbors_for_calibration_labels"}

    merged = labeled_df.merge(
        neighbor_df,
        left_on="canonical_smiles",
        right_on="predicted_canonical_smiles",
        how="left",
    )
    merged = merged.loc[merged["matched_metabolon_super"].map(clean_text).ne("")].copy()
    if merged.empty:
        return {}, {"enabled": False, "reason": "no_labeled_neighbor_overlap"}

    grouped = {smiles: group.copy() for smiles, group in merged.groupby("canonical_smiles", sort=False)}
    pathway_rows: list[dict[str, Any]] = []
    for smiles, group in grouped.items():
        true_pathway = clean_text(group["true_pathway"].iloc[0])
        if not true_pathway:
            continue
        row_payload: dict[str, Any] = {
            "canonical_smiles": smiles,
            "true_pathway": true_pathway,
        }
        for power in POWER_SWEEP_VALUES:
            row_payload[f"pred_{power}"] = predict_pathway_from_neighbors(group, power)
        pathway_rows.append(row_payload)
    if not pathway_rows:
        return {}, {"enabled": False, "reason": "no_pathway_rows_after_sweep"}

    sweep_df = pd.DataFrame(pathway_rows)
    tuned: dict[str, float] = {}
    diagnostics: dict[str, Any] = {
        "enabled": True,
        "source_csv": str(power_sweep_csv.expanduser().resolve()),
        "n_labeled_rows": int(sweep_df.shape[0]),
        "n_unique_pathways": int(sweep_df["true_pathway"].nunique()),
        "pathway_power_accuracy": {},
        "selected_power_by_pathway": {},
    }

    for pathway, subset in sweep_df.groupby("true_pathway", sort=True):
        n_rows = int(subset.shape[0])
        if n_rows < POWER_SWEEP_MIN_SAMPLES:
            continue
        power_scores: list[tuple[float, float]] = []
        for power in POWER_SWEEP_VALUES:
            predicted_col = f"pred_{power}"
            accuracy = float(subset[predicted_col].eq(pathway).mean())
            power_scores.append((power, accuracy))
        best_power, best_accuracy = sorted(
            power_scores,
            key=lambda item: (
                -item[1],
                abs(item[0] - DEFAULT_VOTE_SIMILARITY_POWER),
                item[0],
            ),
        )[0]
        default_accuracy = float(
            subset[f"pred_{DEFAULT_VOTE_SIMILARITY_POWER}"].eq(pathway).mean()
        )
        if best_accuracy - default_accuracy < MIN_POWER_TUNING_IMPROVEMENT:
            best_power = DEFAULT_VOTE_SIMILARITY_POWER
            best_accuracy = default_accuracy
        tuned[pathway] = float(best_power)
        diagnostics["selected_power_by_pathway"][pathway] = float(best_power)
        diagnostics["pathway_power_accuracy"][pathway] = {
            str(power): round(accuracy, 6)
            for power, accuracy in sorted(power_scores, key=lambda item: item[0])
        }
        diagnostics["pathway_power_accuracy"][pathway]["n_rows"] = n_rows
        diagnostics["pathway_power_accuracy"][pathway]["default_accuracy"] = round(default_accuracy, 6)
        diagnostics["pathway_power_accuracy"][pathway]["best_accuracy"] = round(best_accuracy, 6)
    return tuned, diagnostics


def load_reference_library(db_path: Path) -> pd.DataFrame:
    query = """
        SELECT
            chemical_id,
            canonical_name,
            smiles,
            super_pathway,
            pathway,
            metabolon_super,
            metabolon_sub
        FROM metabolites
        WHERE LOWER(TRIM(COALESCE(on_hplasma_209, ''))) LIKE '%true%'
          AND TRIM(COALESCE(smiles, '')) <> ''
          AND TRIM(COALESCE(metabolon_super, '')) <> ''
          AND TRIM(COALESCE(metabolon_sub, '')) <> ''
        ORDER BY chemical_id
    """
    connection = sqlite3.connect(str(db_path))
    try:
        reference_df = pd.read_sql_query(query, connection)
    finally:
        connection.close()

    if reference_df.empty:
        raise ValueError(f"No eligible reference metabolites found in {db_path}.")

    reference_df["chemical_id"] = reference_df["chemical_id"].astype(str)
    reference_df["canonical_name"] = reference_df["canonical_name"].map(clean_text)
    reference_df["super_pathway"] = reference_df["super_pathway"].map(clean_text)
    reference_df["pathway"] = reference_df["pathway"].map(clean_text)
    reference_df["metabolon_super"] = reference_df["metabolon_super"].map(clean_text)
    reference_df["metabolon_sub"] = reference_df["metabolon_sub"].map(clean_text)
    reference_df["canonical_smiles"] = reference_df["smiles"].map(canonicalize_smiles)
    reference_df = reference_df.dropna(subset=["canonical_smiles"]).copy()
    reference_df["molecular_weight"] = reference_df["canonical_smiles"].map(molecular_weight)
    reference_df = reference_df.dropna(subset=["molecular_weight"]).copy()
    reference_df["fingerprint"] = reference_df["canonical_smiles"].map(morgan_fingerprint)
    reference_df = reference_df.dropna(subset=["fingerprint"]).reset_index(drop=True)

    if reference_df.empty:
        raise ValueError(f"Eligible reference metabolites in {db_path} did not yield any fingerprints.")

    return reference_df[
        [
            "chemical_id",
            "canonical_name",
            "canonical_smiles",
            "super_pathway",
            "pathway",
            "metabolon_super",
            "metabolon_sub",
            "molecular_weight",
            "fingerprint",
        ]
    ].copy()


def compute_raw_structure_weight(row: pd.Series) -> float:
    normalized = pd.to_numeric(row.get("normalized_merge_score"), errors="coerce")
    if pd.notna(normalized) and float(normalized) > 0:
        return float(normalized)
    policy_rank = pd.to_numeric(row.get("policy_rank"), errors="coerce")
    if pd.notna(policy_rank) and float(policy_rank) > 0:
        return 1.0 / float(policy_rank)
    return 1.0


def prepare_lane_candidates(
    ranked_df: pd.DataFrame,
    lane: str,
    top_structures_per_feature: int = DEFAULT_TOP_STRUCTURES_PER_FEATURE,
) -> tuple[pd.DataFrame, list[str]]:
    if "feature_key" not in ranked_df.columns:
        raise KeyError("Ranked policy rows are missing feature_key.")

    out = ranked_df.copy()
    out["feature_key"] = out["feature_key"].astype(str)
    feature_keys = sorted(out["feature_key"].dropna().astype(str).unique().tolist())
    out["include_topk_bool"] = out["include_topk"].apply(parse_truthy) if "include_topk" in out.columns else True
    out["query_name"] = ""
    for name_col in ("query_name", "name", "compound_name", "short_name", "chemical_name", "feature_name"):
        if name_col in out.columns:
            candidate_names = out[name_col].map(clean_text)
            out["query_name"] = out["query_name"].mask(out["query_name"].eq(""), candidate_names)

    rank_source = "policy_rank" if "policy_rank" in out.columns else "rank"
    if rank_source not in out.columns:
        raise KeyError("Ranked policy rows are missing policy_rank/rank.")
    out["policy_rank"] = pd.to_numeric(out[rank_source], errors="coerce")
    out = out.loc[out["include_topk_bool"]].copy()
    out = out.loc[out["policy_rank"].le(float(top_structures_per_feature))].copy()
    out["source_smiles"] = out.apply(coalesce_smiles, axis=1)
    out["canonical_smiles"] = out["source_smiles"].map(canonicalize_smiles)
    out = out.dropna(subset=["canonical_smiles"]).copy()

    if out.empty:
        empty_columns = [
            "feature_key",
            "lane",
            "policy_rank",
            "canonical_smiles",
            "source_smiles",
            "query_name",
            "track",
            "normalized_merge_score",
            "raw_structure_weight",
            "structure_weight",
        ]
        return pd.DataFrame(columns=empty_columns), feature_keys

    if "track" not in out.columns:
        out["track"] = ""
    if "normalized_merge_score" not in out.columns:
        out["normalized_merge_score"] = pd.NA
    out["normalized_merge_score"] = pd.to_numeric(out["normalized_merge_score"], errors="coerce")
    out["raw_structure_weight"] = out.apply(compute_raw_structure_weight, axis=1)
    out["lane"] = lane

    out = out.sort_values(
        ["feature_key", "policy_rank", "raw_structure_weight", "source_smiles"],
        ascending=[True, True, False, True],
    )
    out = out.drop_duplicates(subset=["feature_key", "canonical_smiles"], keep="first").copy()
    out["candidate_rank"] = out.groupby("feature_key").cumcount() + 1
    group_totals = out.groupby("feature_key")["raw_structure_weight"].transform("sum")
    out["structure_weight"] = out["raw_structure_weight"] / group_totals

    return out[
        [
            "feature_key",
            "lane",
            "policy_rank",
            "candidate_rank",
            "canonical_smiles",
            "source_smiles",
            "query_name",
            "track",
            "normalized_merge_score",
            "raw_structure_weight",
            "structure_weight",
        ]
    ].reset_index(drop=True), feature_keys


def build_neighbor_table(
    predicted_smiles: list[str],
    reference_df: pd.DataFrame,
    top_k_neighbors: int = DEFAULT_TOP_K_NEIGHBORS,
    similarity_floor: float = DEFAULT_NEIGHBOR_SIMILARITY_FLOOR,
) -> pd.DataFrame:
    if int(top_k_neighbors) < 1:
        raise ValueError("top_k_neighbors must be >= 1.")

    if not predicted_smiles or reference_df.empty:
        return pd.DataFrame(
            columns=[
                "predicted_canonical_smiles",
                "matched_chemical_id",
                "matched_name",
                "matched_similarity",
                "matched_super_pathway",
                "matched_pathway",
                "matched_metabolon_super",
                "matched_metabolon_sub",
                "match_strength",
                "is_exact_match",
                "neighbor_rank",
            ]
        )

    reference_records = reference_df.to_dict("records")
    reference_fps = reference_df["fingerprint"].tolist()
    records: list[dict[str, Any]] = []

    for predicted in sorted(set(predicted_smiles)):
        predicted_fp = morgan_fingerprint(predicted)
        if predicted_fp is None:
            continue
        similarities = list(DataStructs.BulkTanimotoSimilarity(predicted_fp, reference_fps))
        if not similarities:
            continue

        selected_indices = sorted(
            range(len(similarities)),
            key=lambda idx: (-float(similarities[idx]), str(reference_records[idx]["chemical_id"])),
        )[: min(int(top_k_neighbors), len(similarities))]

        for neighbor_rank, idx in enumerate(selected_indices, start=1):
            reference = reference_records[idx]
            similarity = float(similarities[idx])
            records.append(
                {
                    "predicted_canonical_smiles": predicted,
                    "matched_chemical_id": str(reference["chemical_id"]),
                    "matched_name": str(reference["canonical_name"]),
                    "matched_similarity": similarity,
                    "matched_super_pathway": str(reference["super_pathway"]),
                    "matched_pathway": str(reference["pathway"]),
                    "matched_metabolon_super": str(reference["metabolon_super"]),
                    "matched_metabolon_sub": str(reference["metabolon_sub"]),
                    "match_strength": (
                        "strong_similarity_vote"
                        if similarity >= float(similarity_floor)
                        else "weak_similarity_backfill"
                    ),
                    "is_exact_match": bool(similarity >= 0.999999),
                    "neighbor_rank": neighbor_rank,
                }
            )

    return pd.DataFrame(records)


def build_feature_evidence(
    candidates_df: pd.DataFrame,
    neighbor_df: pd.DataFrame,
    vote_similarity_power: float = DEFAULT_VOTE_SIMILARITY_POWER,
    pathway_similarity_power_map: dict[str, float] | None = None,
    pathway_mw_stats: dict[str, dict[str, float]] | None = None,
    low_similarity_penalty_threshold: float = LOW_SIMILARITY_PENALTY_THRESHOLD,
    low_similarity_penalty_floor: float = LOW_SIMILARITY_PENALTY_FLOOR,
) -> pd.DataFrame:
    if candidates_df.empty or neighbor_df.empty:
        return pd.DataFrame()

    evidence_df = candidates_df.merge(
        neighbor_df,
        left_on="canonical_smiles",
        right_on="predicted_canonical_smiles",
        how="inner",
    )
    if evidence_df.empty:
        return evidence_df

    evidence_df["matched_similarity"] = pd.to_numeric(evidence_df["matched_similarity"], errors="coerce").fillna(0.0)
    evidence_df["structure_weight"] = pd.to_numeric(evidence_df["structure_weight"], errors="coerce").fillna(0.0)
    if pathway_similarity_power_map:
        evidence_df["similarity_power"] = evidence_df["matched_metabolon_super"].map(
            lambda value: float(pathway_similarity_power_map.get(clean_text(value), vote_similarity_power))
        )
    else:
        evidence_df["similarity_power"] = float(vote_similarity_power)
    evidence_df["predicted_molecular_weight"] = evidence_df["canonical_smiles"].map(molecular_weight)
    if pathway_mw_stats:
        evidence_df["matched_pathway_mw_prior"] = evidence_df.apply(
            lambda row: compute_pathway_mw_prior(
                predicted_mw=row["predicted_molecular_weight"],
                metabolon_super=clean_text(row["matched_metabolon_super"]),
                pathway_mw_stats=pathway_mw_stats,
            ),
            axis=1,
        )
    else:
        evidence_df["matched_pathway_mw_prior"] = 1.0
    evidence_df["base_vote_weight"] = evidence_df["structure_weight"] * evidence_df["matched_similarity"].pow(
        evidence_df["similarity_power"]
    )
    evidence_df["vote_weight"] = evidence_df["base_vote_weight"] * evidence_df["matched_pathway_mw_prior"]
    feature_best_similarity = (
        evidence_df.groupby("feature_key", as_index=False)["matched_similarity"].max().rename(
            columns={"matched_similarity": "feature_best_similarity"}
        )
    )
    penalty_threshold = max(1e-6, float(low_similarity_penalty_threshold))
    penalty_floor = max(0.0, min(1.0, float(low_similarity_penalty_floor)))
    feature_best_similarity["feature_similarity_penalty"] = feature_best_similarity["feature_best_similarity"].map(
        lambda value: 1.0
        if float(value) >= penalty_threshold
        else max(penalty_floor, float(value) / penalty_threshold)
    )
    evidence_df = evidence_df.merge(feature_best_similarity, on="feature_key", how="left")
    evidence_df["feature_similarity_penalty"] = pd.to_numeric(
        evidence_df["feature_similarity_penalty"], errors="coerce"
    ).fillna(1.0)
    evidence_df["vote_weight"] = evidence_df["vote_weight"] * evidence_df["feature_similarity_penalty"]
    evidence_df = evidence_df.rename(columns={"lane": "assignment_lane"})
    return evidence_df.sort_values(
        ["feature_key", "assignment_lane", "policy_rank", "matched_similarity", "matched_chemical_id"],
        ascending=[True, True, True, False, True],
    ).reset_index(drop=True)


def collapse_pooled_candidates(pooled_candidates_df: pd.DataFrame) -> pd.DataFrame:
    if pooled_candidates_df.empty:
        return pooled_candidates_df

    collapsed = (
        pooled_candidates_df.sort_values(
            ["feature_key", "policy_rank", "structure_weight", "canonical_smiles"],
            ascending=[True, True, False, True],
        )
        .groupby(["feature_key", "canonical_smiles"], as_index=False)
        .agg(
            {
                "policy_rank": "min",
                "candidate_rank": "min",
                "source_smiles": "first",
                "query_name": "first",
                "track": "first",
                "normalized_merge_score": "max",
                "raw_structure_weight": "sum",
                "structure_weight": "sum",
            }
        )
    )
    collapsed["lane"] = "consensus"
    collapsed = collapsed.sort_values(
        ["feature_key", "policy_rank", "structure_weight", "canonical_smiles"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)
    collapsed["candidate_rank"] = collapsed.groupby("feature_key").cumcount() + 1
    return collapsed[
        [
            "feature_key",
            "lane",
            "policy_rank",
            "candidate_rank",
            "canonical_smiles",
            "source_smiles",
            "query_name",
            "track",
            "normalized_merge_score",
            "raw_structure_weight",
            "structure_weight",
        ]
    ].copy()


def top_predicted_smiles_json(candidate_group: pd.DataFrame) -> str:
    if candidate_group.empty:
        return "[]"
    ordered = (
        candidate_group.sort_values(
            ["structure_weight", "policy_rank", "canonical_smiles"],
            ascending=[False, True, True],
        )["canonical_smiles"]
        .astype(str)
        .tolist()
    )
    return json.dumps(ordered)


def classify_confidence_tier(
    assignment_status: str,
    vote_share: float,
    margin_to_runner_up: float,
    best_match_similarity: float,
    exact_match_found: bool,
) -> str:
    if assignment_status == "no_valid_predicted_structure":
        return "low"
    if exact_match_found:
        return "high"
    if assignment_status == "weak_similarity_backfill":
        return "low"
    if vote_share >= 0.70 and margin_to_runner_up >= 0.25 and best_match_similarity >= 0.75:
        return "high"
    if vote_share >= 0.45 and margin_to_runner_up >= 0.10 and best_match_similarity >= 0.60:
        return "medium"
    return "low"


def summarize_feature_assignment(
    feature_key: str,
    lane: str,
    candidate_group: pd.DataFrame,
    evidence_group: pd.DataFrame,
    query_name: str = "",
    subpathway_name_priors: dict[tuple[str, str], dict[str, float | str]] | None = None,
    xenobiotic_classifier_model: dict[str, Any] | None = None,
) -> dict[str, Any]:
    resolved_query_name = choose_query_name(candidate_group, fallback_name=query_name)
    representative_smiles = choose_representative_smiles(candidate_group)
    xenobiotic_probability = aggregate_feature_xenobiotic_probability(candidate_group, xenobiotic_classifier_model)
    summary: dict[str, Any] = {
        "feature_key": feature_key,
        "lane": lane,
        "query_name": resolved_query_name,
        "predicted_pathway": "",
        "predicted_subpathway": "",
        "confidence_tier": "low",
        "assignment_status": "no_valid_predicted_structure",
        "pathway_confidence_score": 0.0,
        "subpathway_confidence_score": 0.0,
        "vote_share": 0.0,
        "margin_to_runner_up": 0.0,
        "exact_match_found": False,
        "best_match_chemical_id": "",
        "best_match_name": "",
        "best_match_similarity": 0.0,
        "best_match_super_pathway": "",
        "best_match_pathway": "",
        "top_predicted_smiles_used": "[]",
        "n_structures_used": 0,
        "xenobiotic_probability": round(xenobiotic_probability, 6),
        "xenobiotic_adjustment_applied": False,
        "name_pattern_boost_applied": False,
    }
    if candidate_group.empty:
        return summary

    summary["top_predicted_smiles_used"] = top_predicted_smiles_json(candidate_group)
    summary["n_structures_used"] = int(candidate_group.shape[0])
    if evidence_group.empty:
        return summary

    pathway_votes = (
        evidence_group.groupby(["matched_metabolon_super"], dropna=False, as_index=False)["vote_weight"]
        .sum()
        .sort_values(["vote_weight", "matched_metabolon_super"], ascending=[False, True])
        .reset_index(drop=True)
    )
    pathway_votes, xenobiotic_adjustment_applied = apply_xenobiotic_vote_adjustment(
        pathway_votes=pathway_votes,
        xenobiotic_probability=xenobiotic_probability,
    )
    summary["xenobiotic_adjustment_applied"] = bool(xenobiotic_adjustment_applied)

    if not pathway_votes.empty and float(pathway_votes["vote_weight"].sum()) > 0:
        total_vote_weight = float(pathway_votes["vote_weight"].sum())
        winning_pathway_vote = pathway_votes.iloc[0]
        winning_pathway = str(winning_pathway_vote["matched_metabolon_super"])
        winning_pathway_vote_weight = float(winning_pathway_vote["vote_weight"])
        pathway_runner_up_weight = float(pathway_votes.iloc[1]["vote_weight"]) if len(pathway_votes) > 1 else 0.0
        vote_share = winning_pathway_vote_weight / total_vote_weight
        margin_to_runner_up = (winning_pathway_vote_weight - pathway_runner_up_weight) / total_vote_weight
        pathway_confidence_score = vote_share

        winning_pathway_rows = evidence_group.loc[
            evidence_group["matched_metabolon_super"].astype(str).eq(winning_pathway)
        ].copy()
        subpathway_votes = (
            winning_pathway_rows.groupby(["matched_metabolon_sub"], dropna=False, as_index=False)["vote_weight"]
            .sum()
            .sort_values(["vote_weight", "matched_metabolon_sub"], ascending=[False, True])
            .reset_index(drop=True)
        )
        name_pattern_boost_applied = False
        if not subpathway_votes.empty:
            if len(subpathway_votes) > 1 and winning_pathway_vote_weight > 0:
                sub_runner_up_weight = float(subpathway_votes.iloc[1]["vote_weight"])
                sub_margin = (float(subpathway_votes.iloc[0]["vote_weight"]) - sub_runner_up_weight) / float(
                    winning_pathway_vote_weight
                )
            else:
                sub_margin = 1.0
            if sub_margin < SUBPATHWAY_TIE_MARGIN and subpathway_name_priors:
                subpathway_votes, name_pattern_boost_applied = apply_subpathway_name_boost(
                    subpathway_votes=subpathway_votes,
                    query_name=resolved_query_name,
                    winning_pathway=winning_pathway,
                    winning_pathway_vote_weight=winning_pathway_vote_weight,
                    subpathway_name_priors=subpathway_name_priors,
                )
            winning_sub_vote = subpathway_votes.iloc[0]
            winning_subpathway = str(winning_sub_vote["matched_metabolon_sub"])
            winning_subpathway_vote_weight = float(winning_sub_vote["vote_weight"])
            subpathway_confidence_score = (
                winning_subpathway_vote_weight / winning_pathway_vote_weight
                if winning_pathway_vote_weight > 0
                else 0.0
            )
            summary["name_pattern_boost_applied"] = bool(name_pattern_boost_applied)
        else:
            fallback_row = winning_pathway_rows.sort_values(
                ["matched_similarity", "structure_weight", "matched_chemical_id"],
                ascending=[False, False, True],
            ).iloc[0]
            winning_subpathway = str(fallback_row["matched_metabolon_sub"])
            subpathway_confidence_score = 1.0
    else:
        winning_row = evidence_group.sort_values(
            ["matched_similarity", "structure_weight", "matched_chemical_id"],
            ascending=[False, False, True],
        ).iloc[0]
        winning_pathway = str(winning_row["matched_metabolon_super"])
        winning_subpathway = str(winning_row["matched_metabolon_sub"])
        vote_share = 1.0
        margin_to_runner_up = 1.0
        pathway_confidence_score = 1.0
        subpathway_confidence_score = 1.0

    winning_rows = evidence_group.loc[
        evidence_group["matched_metabolon_super"].astype(str).eq(winning_pathway)
        & evidence_group["matched_metabolon_sub"].astype(str).eq(winning_subpathway)
    ].copy()
    winning_rows = winning_rows.sort_values(
        ["matched_similarity", "vote_weight", "matched_chemical_id"],
        ascending=[False, False, True],
    )
    best_match = winning_rows.iloc[0]
    assignment_status = (
        "strong_similarity_vote"
        if bool(winning_rows["match_strength"].eq("strong_similarity_vote").any())
        else "weak_similarity_backfill"
    )
    exact_match_found = bool(winning_rows["is_exact_match"].any())
    best_match_similarity = float(best_match["matched_similarity"])

    summary.update(
        {
            "predicted_pathway": winning_pathway,
            "predicted_subpathway": winning_subpathway,
            "assignment_status": assignment_status,
            "pathway_confidence_score": round(pathway_confidence_score, 6),
            "subpathway_confidence_score": round(subpathway_confidence_score, 6),
            "vote_share": round(vote_share, 6),
            "margin_to_runner_up": round(margin_to_runner_up, 6),
            "exact_match_found": exact_match_found,
            "best_match_chemical_id": str(best_match["matched_chemical_id"]),
            "best_match_name": str(best_match["matched_name"]),
            "best_match_similarity": round(best_match_similarity, 6),
            "best_match_super_pathway": str(best_match["matched_super_pathway"]),
            "best_match_pathway": str(best_match["matched_pathway"]),
        }
    )
    override = infer_pathway_override(
        smiles=representative_smiles,
        query_name=resolved_query_name,
        current_pathway=summary["predicted_pathway"],
        current_subpathway=summary["predicted_subpathway"],
    )
    summary.update(override)
    summary["confidence_tier"] = classify_confidence_tier(
        assignment_status=assignment_status,
        vote_share=float(summary["vote_share"]),
        margin_to_runner_up=float(summary["margin_to_runner_up"]),
        best_match_similarity=best_match_similarity,
        exact_match_found=exact_match_found,
    )
    return summary


def build_lane_summary(
    feature_keys: list[str],
    lane: str,
    lane_candidates_df: pd.DataFrame,
    neighbor_df: pd.DataFrame,
    vote_similarity_power: float = DEFAULT_VOTE_SIMILARITY_POWER,
    pathway_similarity_power_map: dict[str, float] | None = None,
    pathway_mw_stats: dict[str, dict[str, float]] | None = None,
    feature_query_name_map: dict[str, str] | None = None,
    subpathway_name_priors: dict[tuple[str, str], dict[str, float | str]] | None = None,
    xenobiotic_classifier_model: dict[str, Any] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    evidence_df = build_feature_evidence(
        lane_candidates_df,
        neighbor_df,
        vote_similarity_power=vote_similarity_power,
        pathway_similarity_power_map=pathway_similarity_power_map,
        pathway_mw_stats=pathway_mw_stats,
    )
    candidate_groups = {key: group.copy() for key, group in lane_candidates_df.groupby("feature_key", sort=False)}
    evidence_groups = {key: group.copy() for key, group in evidence_df.groupby("feature_key", sort=False)}

    records = [
        summarize_feature_assignment(
            feature_key=feature_key,
            lane=lane,
            candidate_group=candidate_groups.get(feature_key, pd.DataFrame()),
            evidence_group=evidence_groups.get(feature_key, pd.DataFrame()),
            query_name=clean_text((feature_query_name_map or {}).get(feature_key, "")),
            subpathway_name_priors=subpathway_name_priors,
            xenobiotic_classifier_model=xenobiotic_classifier_model,
        )
        for feature_key in feature_keys
    ]
    summary_df = pd.DataFrame(records).sort_values("feature_key").reset_index(drop=True)
    return summary_df, evidence_df


def prefixed_summary_fields(row: pd.Series, prefix: str) -> dict[str, Any]:
    wanted = [
        "predicted_pathway",
        "predicted_subpathway",
        "confidence_tier",
        "assignment_status",
        "pathway_confidence_score",
        "subpathway_confidence_score",
        "vote_share",
        "margin_to_runner_up",
        "exact_match_found",
        "best_match_chemical_id",
        "best_match_name",
        "best_match_similarity",
    ]
    return {f"{prefix}_{field}": row.get(field) for field in wanted}


def build_consensus_row_exact(feature_key: str, coverage_row: pd.Series, hitmax_row: pd.Series) -> dict[str, Any]:
    coverage_vote_share = safe_numeric(coverage_row.get("vote_share"), default=0.0)
    hitmax_vote_share = safe_numeric(hitmax_row.get("vote_share"), default=0.0)
    coverage_similarity = safe_numeric(coverage_row.get("best_match_similarity"), default=0.0)
    hitmax_similarity = safe_numeric(hitmax_row.get("best_match_similarity"), default=0.0)
    representative = coverage_row
    if (hitmax_vote_share, hitmax_similarity) > (coverage_vote_share, coverage_similarity):
        representative = hitmax_row

    record: dict[str, Any] = {
        "feature_key": feature_key,
        "lane_agreement": "exact",
        "consensus_source": str(representative.get("lane", "")),
        "consensus_predicted_pathway": representative.get("predicted_pathway", ""),
        "consensus_predicted_subpathway": representative.get("predicted_subpathway", ""),
        "consensus_confidence_tier": representative.get("confidence_tier", "low"),
        "consensus_assignment_status": representative.get("assignment_status", "no_valid_predicted_structure"),
        "consensus_pathway_confidence_score": representative.get("pathway_confidence_score", 0.0),
        "consensus_subpathway_confidence_score": representative.get("subpathway_confidence_score", 0.0),
        "consensus_vote_share": representative.get("vote_share", 0.0),
        "consensus_margin_to_runner_up": representative.get("margin_to_runner_up", 0.0),
        "consensus_exact_match_found": representative.get("exact_match_found", False),
        "consensus_best_match_chemical_id": representative.get("best_match_chemical_id", ""),
        "consensus_best_match_name": representative.get("best_match_name", ""),
        "consensus_best_match_similarity": representative.get("best_match_similarity", 0.0),
        "consensus_top_predicted_smiles_used": representative.get("top_predicted_smiles_used", "[]"),
    }
    record.update(prefixed_summary_fields(coverage_row, "coverage"))
    record.update(prefixed_summary_fields(hitmax_row, "hitmax"))
    return record


def build_consensus_row_pooled(
    feature_key: str,
    coverage_row: pd.Series,
    hitmax_row: pd.Series,
    pooled_summary: dict[str, Any],
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "feature_key": feature_key,
        "lane_agreement": "disagree",
        "consensus_source": "pooled",
        "consensus_predicted_pathway": pooled_summary.get("predicted_pathway", ""),
        "consensus_predicted_subpathway": pooled_summary.get("predicted_subpathway", ""),
        "consensus_confidence_tier": pooled_summary.get("confidence_tier", "low"),
        "consensus_assignment_status": pooled_summary.get("assignment_status", "no_valid_predicted_structure"),
        "consensus_pathway_confidence_score": pooled_summary.get("pathway_confidence_score", 0.0),
        "consensus_subpathway_confidence_score": pooled_summary.get("subpathway_confidence_score", 0.0),
        "consensus_vote_share": pooled_summary.get("vote_share", 0.0),
        "consensus_margin_to_runner_up": pooled_summary.get("margin_to_runner_up", 0.0),
        "consensus_exact_match_found": pooled_summary.get("exact_match_found", False),
        "consensus_best_match_chemical_id": pooled_summary.get("best_match_chemical_id", ""),
        "consensus_best_match_name": pooled_summary.get("best_match_name", ""),
        "consensus_best_match_similarity": pooled_summary.get("best_match_similarity", 0.0),
        "consensus_top_predicted_smiles_used": pooled_summary.get("top_predicted_smiles_used", "[]"),
    }
    record.update(prefixed_summary_fields(coverage_row, "coverage"))
    record.update(prefixed_summary_fields(hitmax_row, "hitmax"))
    return record


def build_consensus_outputs(
    feature_keys: list[str],
    lane_candidate_map: dict[str, pd.DataFrame],
    lane_summary_map: dict[str, pd.DataFrame],
    neighbor_df: pd.DataFrame,
    vote_similarity_power: float = DEFAULT_VOTE_SIMILARITY_POWER,
    pathway_similarity_power_map: dict[str, float] | None = None,
    pathway_mw_stats: dict[str, dict[str, float]] | None = None,
    feature_query_name_map: dict[str, str] | None = None,
    subpathway_name_priors: dict[tuple[str, str], dict[str, float | str]] | None = None,
    xenobiotic_classifier_model: dict[str, Any] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    lane_row_maps = {
        lane: lane_summary_map[lane].set_index("feature_key", drop=False) for lane in OFFICIAL_LANES
    }
    lane_candidate_groups = {
        lane: {key: group.copy() for key, group in lane_candidate_map[lane].groupby("feature_key", sort=False)}
        for lane in OFFICIAL_LANES
    }

    consensus_records: list[dict[str, Any]] = []
    consensus_evidence_frames: list[pd.DataFrame] = []
    for feature_key in feature_keys:
        coverage_row = lane_row_maps["coverage"].loc[feature_key]
        hitmax_row = lane_row_maps["hitmax"].loc[feature_key]
        coverage_pair = (
            str(coverage_row.get("predicted_pathway", "")),
            str(coverage_row.get("predicted_subpathway", "")),
        )
        hitmax_pair = (
            str(hitmax_row.get("predicted_pathway", "")),
            str(hitmax_row.get("predicted_subpathway", "")),
        )

        if coverage_pair == hitmax_pair:
            consensus_records.append(build_consensus_row_exact(feature_key, coverage_row, hitmax_row))
            continue

        pooled_parts: list[pd.DataFrame] = []
        for lane in OFFICIAL_LANES:
            lane_group = lane_candidate_groups[lane].get(feature_key)
            if lane_group is None or lane_group.empty:
                continue
            pooled_group = lane_group.copy()
            pooled_group["structure_weight"] = pooled_group["structure_weight"] * 0.5
            pooled_parts.append(pooled_group)
        pooled_candidates_df = (
            collapse_pooled_candidates(pd.concat(pooled_parts, ignore_index=True))
            if pooled_parts
            else pd.DataFrame()
        )
        pooled_evidence_df = build_feature_evidence(
            pooled_candidates_df,
            neighbor_df,
            vote_similarity_power=vote_similarity_power,
            pathway_similarity_power_map=pathway_similarity_power_map,
            pathway_mw_stats=pathway_mw_stats,
        )
        pooled_summary = summarize_feature_assignment(
            feature_key=feature_key,
            lane="consensus",
            candidate_group=pooled_candidates_df,
            evidence_group=pooled_evidence_df,
            query_name=clean_text((feature_query_name_map or {}).get(feature_key, "")),
            subpathway_name_priors=subpathway_name_priors,
            xenobiotic_classifier_model=xenobiotic_classifier_model,
        )
        if not pooled_evidence_df.empty:
            pooled_evidence_df = pooled_evidence_df.copy()
            pooled_evidence_df["assignment_lane"] = "consensus"
            pooled_evidence_df["consensus_source"] = "pooled"
            consensus_evidence_frames.append(pooled_evidence_df)
        consensus_records.append(build_consensus_row_pooled(feature_key, coverage_row, hitmax_row, pooled_summary))

    consensus_df = pd.DataFrame(consensus_records).sort_values("feature_key").reset_index(drop=True)
    consensus_evidence_df = (
        pd.concat(consensus_evidence_frames, ignore_index=True) if consensus_evidence_frames else pd.DataFrame()
    )
    return consensus_df, consensus_evidence_df


def validate_assigned_labels(
    summary_df: pd.DataFrame,
    pathway_col: str,
    subpathway_col: str,
    allowed_pairs: set[tuple[str, str]],
) -> None:
    invalid_pairs = []
    for row in summary_df[[pathway_col, subpathway_col]].itertuples(index=False):
        pathway = clean_text(row[0])
        subpathway = clean_text(row[1])
        if not pathway and not subpathway:
            continue
        if (pathway, subpathway) not in allowed_pairs:
            invalid_pairs.append((pathway, subpathway))
    if invalid_pairs:
        raise ValueError(f"Found labels outside the allowed Metabolon library set: {sorted(set(invalid_pairs))[:5]}")


def count_label_pairs(summary_df: pd.DataFrame, pathway_col: str, subpathway_col: str) -> dict[str, int]:
    labels = summary_df[[pathway_col, subpathway_col]].copy()
    labels[pathway_col] = labels[pathway_col].map(clean_text)
    labels[subpathway_col] = labels[subpathway_col].map(clean_text)
    labels = labels.loc[labels[pathway_col].ne("") & labels[subpathway_col].ne("")].copy()
    counter = Counter(
        f"{row[0]} | {row[1]}"
        for row in labels[[pathway_col, subpathway_col]].itertuples(index=False, name=None)
    )
    return {key: counter[key] for key in sorted(counter, key=lambda item: (-counter[item], item))}


def summarize_output_table(
    summary_df: pd.DataFrame,
    pathway_col: str,
    subpathway_col: str,
    extra_counts: dict[str, Any] | None = None,
) -> dict[str, Any]:
    pathway_series = summary_df[pathway_col].map(clean_text)
    subpathway_series = summary_df[subpathway_col].map(clean_text)
    assigned_mask = pathway_series.ne("") & subpathway_series.ne("")
    payload: dict[str, Any] = {
        "n_features": int(summary_df.shape[0]),
        "n_features_with_assignment": int(assigned_mask.sum()),
        "n_features_without_assignment": int((~assigned_mask).sum()),
        "label_pair_counts": count_label_pairs(summary_df, pathway_col, subpathway_col),
    }
    for column in ("assignment_status", "confidence_tier", "lane_agreement"):
        if column in summary_df.columns:
            counter = Counter(summary_df[column].map(clean_text))
            payload[f"{column}_counts"] = {
                key: counter[key]
                for key in sorted(counter, key=lambda item: (-counter[item], item))
                if key
            }
    if extra_counts:
        payload.update(extra_counts)
    return payload


def serialize_json_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for column in columns:
        if column in out.columns:
            out[column] = out[column].apply(
                lambda value: json.dumps(value) if isinstance(value, (list, dict)) else value
            )
    return out


def load_known_truth_labels(known_x_csv: Path) -> pd.DataFrame:
    known_x_path = known_x_csv.expanduser().resolve()
    if not known_x_path.exists():
        return pd.DataFrame(columns=["x_compound_id", "true_pathway", "true_subpathway"])
    try:
        known_df = read_csv_clean(known_x_path, dtype=str, encoding="utf-8-sig")
    except EmptyDataError:
        return pd.DataFrame(columns=["x_compound_id", "true_pathway", "true_subpathway"])
    required = {"LIB_ID", "Unknown", "SuperPathway", "Draft Subpathway"}
    missing = sorted(required.difference(known_df.columns))
    if missing:
        raise ValueError(
            f"Missing required columns in {known_x_csv}: {', '.join(missing)}"
        )

    labels = known_df.rename(
        columns={
            "LIB_ID": "library_id",
            "Unknown": "x_compound_id",
            "SuperPathway": "true_pathway",
            "Draft Subpathway": "true_subpathway",
        }
    ).copy()
    labels["library_id"] = labels["library_id"].map(clean_text)
    labels["x_compound_id"] = labels["x_compound_id"].map(clean_text)
    labels["true_pathway"] = labels["true_pathway"].map(clean_text)
    labels["true_subpathway"] = labels["true_subpathway"].map(clean_text)
    labels = labels.loc[
        labels["library_id"].eq("209")
        & labels["x_compound_id"].ne("")
        & labels["true_pathway"].ne("")
        & labels["true_subpathway"].ne("")
    ].copy()

    pair_counts = (
        labels.groupby("x_compound_id")[["true_pathway", "true_subpathway"]]
        .apply(lambda frame: len(set(zip(frame["true_pathway"], frame["true_subpathway"]))))
    )
    ambiguous = pair_counts[pair_counts.gt(1)]
    if not ambiguous.empty:
        examples = ambiguous.index.tolist()[:5]
        raise ValueError(
            "Known labels are ambiguous for some X-compound IDs "
            f"(multiple pathway/subpathway pairs): {examples}"
        )

    labels = labels.drop_duplicates(subset=["x_compound_id"], keep="first").reset_index(drop=True)
    return labels[["x_compound_id", "true_pathway", "true_subpathway"]].copy()


def load_known_structure_ids(known_structures_csv: Path) -> set[str]:
    known_structures_path = known_structures_csv.expanduser().resolve()
    if not known_structures_path.exists():
        return set()
    try:
        structure_df = read_csv_clean(known_structures_path, dtype=str)
    except EmptyDataError:
        return set()
    if structure_df.empty:
        return set()
    if "x_compound_id" not in structure_df.columns:
        raise ValueError(
            f"Missing required column x_compound_id in {known_structures_csv}."
        )
    return {
        clean_text(value)
        for value in structure_df["x_compound_id"].tolist()
        if clean_text(value)
    }


def similarity_bucket(value: object) -> str:
    similarity = pd.to_numeric(value, errors="coerce")
    if pd.isna(similarity):
        return "missing"
    score = float(similarity)
    if score >= 0.90:
        return ">=0.90"
    if score >= 0.80:
        return "0.80-0.899"
    if score >= 0.70:
        return "0.70-0.799"
    return "<0.70"


def compute_accuracy_counts(
    scored_df: pd.DataFrame,
    *,
    predicted_pathway_col: str,
    predicted_subpathway_col: str,
) -> dict[str, Any]:
    total = int(scored_df.shape[0])
    if total == 0:
        return {
            "n_features": 0,
            "pathway_correct": 0,
            "subpathway_correct": 0,
            "pathway_accuracy": 0.0,
            "subpathway_accuracy": 0.0,
        }

    predicted_pathway = scored_df[predicted_pathway_col].map(clean_text)
    predicted_subpathway = scored_df[predicted_subpathway_col].map(clean_text)
    truth_pathway = scored_df["true_pathway"].map(clean_text)
    truth_subpathway = scored_df["true_subpathway"].map(clean_text)

    pathway_correct_mask = predicted_pathway.eq(truth_pathway)
    subpathway_correct_mask = pathway_correct_mask & predicted_subpathway.eq(truth_subpathway)
    pathway_correct = int(pathway_correct_mask.sum())
    subpathway_correct = int(subpathway_correct_mask.sum())

    return {
        "n_features": total,
        "pathway_correct": pathway_correct,
        "subpathway_correct": subpathway_correct,
        "pathway_accuracy": round(pathway_correct / total, 6),
        "subpathway_accuracy": round(subpathway_correct / total, 6),
    }


def evaluate_summary_against_known(
    summary_df: pd.DataFrame,
    *,
    predicted_pathway_col: str,
    predicted_subpathway_col: str,
    confidence_col: str,
    similarity_col: str,
    known_labels_df: pd.DataFrame,
    known_structure_ids: set[str],
) -> dict[str, Any]:
    if "feature_key" not in summary_df.columns:
        raise KeyError("Summary table is missing feature_key.")

    scored = summary_df.copy()
    scored["feature_key"] = scored["feature_key"].astype(str)
    scored["x_compound_id"] = scored["feature_key"].map(extract_xid).map(clean_text)
    if predicted_pathway_col not in scored.columns or predicted_subpathway_col not in scored.columns:
        raise KeyError(
            f"Summary table missing required prediction columns: {predicted_pathway_col}, {predicted_subpathway_col}"
        )
    scored[predicted_pathway_col] = scored[predicted_pathway_col].map(clean_text)
    scored[predicted_subpathway_col] = scored[predicted_subpathway_col].map(clean_text)
    if confidence_col in scored.columns:
        scored[confidence_col] = scored[confidence_col].map(clean_text)

    scored = scored.merge(known_labels_df, on="x_compound_id", how="left")
    scored["true_pathway"] = scored["true_pathway"].map(clean_text)
    scored["true_subpathway"] = scored["true_subpathway"].map(clean_text)

    known_label_mask = scored["true_pathway"].ne("") & scored["true_subpathway"].ne("")
    scored_with_labels = scored.loc[known_label_mask].copy()
    scored_with_labels["similarity_bucket"] = scored_with_labels[similarity_col].map(similarity_bucket)

    by_confidence: dict[str, Any] = {}
    if confidence_col in scored_with_labels.columns:
        confidence_labels = sorted(
            {
                clean_text(value)
                for value in scored_with_labels[confidence_col].tolist()
                if clean_text(value)
            }
        )
        for tier in confidence_labels:
            tier_rows = scored_with_labels.loc[scored_with_labels[confidence_col].eq(tier)].copy()
            by_confidence[tier] = compute_accuracy_counts(
                tier_rows,
                predicted_pathway_col=predicted_pathway_col,
                predicted_subpathway_col=predicted_subpathway_col,
            )

    by_similarity_bucket: dict[str, Any] = {}
    for bucket in SIMILARITY_BUCKET_ORDER:
        bucket_rows = scored_with_labels.loc[scored_with_labels["similarity_bucket"].eq(bucket)].copy()
        if bucket_rows.empty:
            continue
        by_similarity_bucket[bucket] = compute_accuracy_counts(
            bucket_rows,
            predicted_pathway_col=predicted_pathway_col,
            predicted_subpathway_col=predicted_subpathway_col,
        )

    overall = compute_accuracy_counts(
        scored_with_labels,
        predicted_pathway_col=predicted_pathway_col,
        predicted_subpathway_col=predicted_subpathway_col,
    )
    overall.update(
        {
            "n_features_total": int(summary_df.shape[0]),
            "n_features_with_x_compound_id": int(scored["x_compound_id"].ne("").sum()),
            "n_features_with_known_labels": int(scored_with_labels.shape[0]),
            "n_unique_known_x_compound_ids": int(scored_with_labels["x_compound_id"].nunique()),
            "n_features_in_known_structure_catalog": int(
                scored_with_labels["x_compound_id"].isin(known_structure_ids).sum()
            ),
        }
    )
    return {
        "overall": overall,
        "by_confidence_tier": by_confidence,
        "by_similarity_bucket": by_similarity_bucket,
    }


def evaluate_assignment_outputs(
    *,
    lane_summary_map: dict[str, pd.DataFrame],
    consensus_df: pd.DataFrame,
    known_x_csv: Path,
    known_structures_csv: Path,
) -> dict[str, Any]:
    known_labels_df = load_known_truth_labels(known_x_csv)
    known_structure_ids = load_known_structure_ids(known_structures_csv)

    lanes_payload = {
        lane: evaluate_summary_against_known(
            lane_summary_map[lane],
            predicted_pathway_col="predicted_pathway",
            predicted_subpathway_col="predicted_subpathway",
            confidence_col="confidence_tier",
            similarity_col="best_match_similarity",
            known_labels_df=known_labels_df,
            known_structure_ids=known_structure_ids,
        )
        for lane in OFFICIAL_LANES
    }
    consensus_payload = evaluate_summary_against_known(
        consensus_df,
        predicted_pathway_col="consensus_predicted_pathway",
        predicted_subpathway_col="consensus_predicted_subpathway",
        confidence_col="consensus_confidence_tier",
        similarity_col="consensus_best_match_similarity",
        known_labels_df=known_labels_df,
        known_structure_ids=known_structure_ids,
    )

    return {
        "known_label_source": str(known_x_csv.expanduser().resolve()),
        "known_structure_source": str(known_structures_csv.expanduser().resolve()),
        "known_label_pool": {
            "n_x_compound_ids_with_labels": int(known_labels_df["x_compound_id"].nunique()),
            "n_x_compound_ids_with_structures": int(len(known_structure_ids)),
        },
        "lanes": lanes_payload,
        "consensus": consensus_payload,
    }


def format_percent(value: float) -> str:
    return f"{100.0 * float(value):.1f}%"


def print_evaluation_report(evaluation_payload: dict[str, Any]) -> None:
    print("\nEvaluation vs known_X labels (LIB_ID=209)")
    print("lane       n_eval  pathway_acc  subpathway_acc")

    for lane in (*OFFICIAL_LANES, "consensus"):
        lane_payload = (
            evaluation_payload["lanes"][lane]
            if lane in OFFICIAL_LANES
            else evaluation_payload["consensus"]
        )
        overall = lane_payload["overall"]
        print(
            f"{lane:<10} {overall['n_features']:>6d}  "
            f"{format_percent(overall['pathway_accuracy']):>11}  "
            f"{format_percent(overall['subpathway_accuracy']):>14}"
        )

    for lane in OFFICIAL_LANES:
        print(f"\n{lane} by confidence_tier")
        print("tier       n_eval  pathway_acc  subpathway_acc")
        tier_payload = evaluation_payload["lanes"][lane]["by_confidence_tier"]
        for tier in ("high", "medium", "low"):
            if tier not in tier_payload:
                continue
            metrics = tier_payload[tier]
            print(
                f"{tier:<10} {metrics['n_features']:>6d}  "
                f"{format_percent(metrics['pathway_accuracy']):>11}  "
                f"{format_percent(metrics['subpathway_accuracy']):>14}"
            )

        print(f"{lane} by best_match_similarity")
        print("bucket       n_eval  pathway_acc  subpathway_acc")
        bucket_payload = evaluation_payload["lanes"][lane]["by_similarity_bucket"]
        for bucket in SIMILARITY_BUCKET_ORDER:
            if bucket not in bucket_payload:
                continue
            metrics = bucket_payload[bucket]
            print(
                f"{bucket:<12} {metrics['n_features']:>6d}  "
                f"{format_percent(metrics['pathway_accuracy']):>11}  "
                f"{format_percent(metrics['subpathway_accuracy']):>14}"
            )


def collect_regression_metrics(payload: dict[str, Any]) -> dict[str, float]:
    metrics: dict[str, float] = {}
    if not isinstance(payload, dict):
        return metrics

    if "evaluation" in payload and isinstance(payload["evaluation"], dict):
        return collect_regression_metrics(payload["evaluation"])

    if "pathway_accuracy" in payload and "subpathway_accuracy" in payload:
        metrics["consensus.pathway_accuracy"] = float(payload["pathway_accuracy"])
        metrics["consensus.subpathway_accuracy"] = float(payload["subpathway_accuracy"])
        return metrics

    lanes_payload = payload.get("lanes")
    if isinstance(lanes_payload, dict):
        for lane in OFFICIAL_LANES:
            lane_payload = lanes_payload.get(lane, {})
            if not isinstance(lane_payload, dict):
                continue
            overall = lane_payload.get("overall", {})
            if not isinstance(overall, dict):
                continue
            if "pathway_accuracy" in overall and "subpathway_accuracy" in overall:
                metrics[f"{lane}.pathway_accuracy"] = float(overall["pathway_accuracy"])
                metrics[f"{lane}.subpathway_accuracy"] = float(overall["subpathway_accuracy"])

    consensus_payload = payload.get("consensus")
    if isinstance(consensus_payload, dict):
        overall = consensus_payload.get("overall", {})
        if isinstance(overall, dict) and "pathway_accuracy" in overall and "subpathway_accuracy" in overall:
            metrics["consensus.pathway_accuracy"] = float(overall["pathway_accuracy"])
            metrics["consensus.subpathway_accuracy"] = float(overall["subpathway_accuracy"])

    return metrics


def run_regression_check(
    evaluation_payload: dict[str, Any],
    baseline_json_path: Path,
    fail_drop_pp: float,
) -> int:
    baseline_path = baseline_json_path.expanduser().resolve()
    if not baseline_path.exists():
        print(f"\nWARNING: regression baseline not found at {baseline_path}; skipping regression check.")
        return 0
    try:
        baseline_payload = load_json(baseline_path)
    except Exception as exc:  # pragma: no cover - defensive parse path
        print(f"\nWARNING: failed to parse baseline JSON ({baseline_path}): {exc}; skipping regression check.")
        return 0

    current_metrics = collect_regression_metrics(evaluation_payload)
    baseline_metrics = collect_regression_metrics(baseline_payload)
    shared_metrics = sorted(set(current_metrics).intersection(set(baseline_metrics)))
    if not shared_metrics:
        print("\nWARNING: no overlapping accuracy metrics found between baseline and current run.")
        return 0

    fail_threshold = max(0.0, float(fail_drop_pp)) / 100.0
    failing_metrics: list[str] = []
    print("\nRegression check against baseline:")
    for metric in shared_metrics:
        baseline_value = float(baseline_metrics[metric])
        current_value = float(current_metrics[metric])
        drop = baseline_value - current_value
        delta_pp = 100.0 * (current_value - baseline_value)
        status = "OK"
        if drop > 0:
            status = "WARNING"
        if drop > fail_threshold:
            status = "FAIL"
            failing_metrics.append(metric)
        print(
            f"- {metric}: baseline={format_percent(baseline_value)}, "
            f"current={format_percent(current_value)}, delta={delta_pp:+.2f}pp [{status}]"
        )

    if failing_metrics:
        print(
            f"\nERROR: Regression check failed (drop > {float(fail_drop_pp):.2f}pp) "
            f"for {len(failing_metrics)} metric(s): {', '.join(failing_metrics)}"
        )
        return 1
    print(f"\nRegression check passed (no metric dropped by more than {float(fail_drop_pp):.2f}pp).")
    return 0

def run_assignment(
    *,
    delivery_manifest: Path,
    reference_db: Path,
    output_dir: Path,
    top_structures_per_feature: int = DEFAULT_TOP_STRUCTURES_PER_FEATURE,
    top_k_neighbors: int = DEFAULT_TOP_K_NEIGHBORS,
    neighbor_similarity_floor: float = DEFAULT_NEIGHBOR_SIMILARITY_FLOOR,
    vote_similarity_power: float = DEFAULT_VOTE_SIMILARITY_POWER,
    evaluate: bool = False,
    known_x_csv: Path = DEFAULT_KNOWN_X_CSV,
    known_structures_csv: Path = DEFAULT_KNOWN_STRUCTURES_CSV,
    power_sweep_csv: Path | None = DEFAULT_POWER_SWEEP_CSV,
    disable_pathway_power_sweep: bool = False,
) -> dict[str, Any]:
    if not HAS_RDKIT:
        raise RuntimeError("RDKit is required to assign predicted pathways.")
    if int(top_structures_per_feature) < 1:
        raise ValueError("top_structures_per_feature must be >= 1.")
    if int(top_k_neighbors) < 1:
        raise ValueError("top_k_neighbors must be >= 1.")
    if float(vote_similarity_power) <= 0:
        raise ValueError("vote_similarity_power must be > 0.")

    manifest_path = delivery_manifest.expanduser().resolve()
    reference_db_path = reference_db.expanduser().resolve()
    output_dir_path = output_dir.expanduser().resolve()
    resolved_power_sweep_csv = (
        power_sweep_csv.expanduser().resolve() if isinstance(power_sweep_csv, Path) else None
    )
    output_dir_path.mkdir(parents=True, exist_ok=True)

    manifest_payload = load_json(manifest_path)
    lane_inputs = resolve_lane_inputs(manifest_payload)
    reference_df = load_reference_library(reference_db_path)
    pathway_mw_stats = build_pathway_mw_stats(reference_df)
    subpathway_name_priors = build_subpathway_name_prior_map(reference_df)
    xenobiotic_classifier_model = train_xenobiotic_classifier(reference_df)
    if disable_pathway_power_sweep:
        pathway_similarity_power_map: dict[str, float] = {}
        power_sweep_diagnostics: dict[str, Any] = {
            "enabled": False,
            "reason": "disabled_by_flag",
        }
    else:
        pathway_similarity_power_map, power_sweep_diagnostics = tune_similarity_power_by_pathway(
            reference_df,
            top_k_neighbors=int(top_k_neighbors),
            power_sweep_csv=resolved_power_sweep_csv,
        )
    allowed_pairs = {
        (str(row.metabolon_super), str(row.metabolon_sub))
        for row in reference_df[["metabolon_super", "metabolon_sub"]].itertuples(index=False)
    }

    lane_candidates: dict[str, pd.DataFrame] = {}
    feature_key_set: set[str] = set()
    for lane in OFFICIAL_LANES:
        ranked_path = Path(str(lane_inputs[lane]["ranked_policy_rows_csv"]))
        ranked_df = read_csv_or_empty(ranked_path)
        prepared_df, feature_keys = prepare_lane_candidates(
            ranked_df=ranked_df,
            lane=lane,
            top_structures_per_feature=top_structures_per_feature,
        )
        lane_candidates[lane] = prepared_df
        feature_key_set.update(feature_keys)
    feature_keys = sorted(feature_key_set)
    feature_query_name_map: dict[str, str] = {}
    for lane in OFFICIAL_LANES:
        lane_df = lane_candidates[lane]
        if "query_name" not in lane_df.columns or lane_df.empty:
            continue
        for row in lane_df[["feature_key", "query_name"]].itertuples(index=False):
            feature_key = clean_text(getattr(row, "feature_key", ""))
            query_name = clean_text(getattr(row, "query_name", ""))
            if not feature_key or not query_name or feature_key in feature_query_name_map:
                continue
            feature_query_name_map[feature_key] = query_name

    predicted_smiles = sorted(
        {
            str(smiles)
            for lane_df in lane_candidates.values()
            for smiles in lane_df.get("canonical_smiles", pd.Series(dtype=str)).dropna().astype(str).tolist()
        }
    )
    neighbor_df = build_neighbor_table(
        predicted_smiles=predicted_smiles,
        reference_df=reference_df,
        top_k_neighbors=top_k_neighbors,
        similarity_floor=neighbor_similarity_floor,
    )

    lane_summary_map: dict[str, pd.DataFrame] = {}
    evidence_frames: list[pd.DataFrame] = []
    output_files: dict[str, str] = {}
    for lane in OFFICIAL_LANES:
        summary_df, evidence_df = build_lane_summary(
            feature_keys=feature_keys,
            lane=lane,
            lane_candidates_df=lane_candidates[lane],
            neighbor_df=neighbor_df,
            vote_similarity_power=vote_similarity_power,
            pathway_similarity_power_map=pathway_similarity_power_map,
            pathway_mw_stats=pathway_mw_stats,
            feature_query_name_map=feature_query_name_map,
            subpathway_name_priors=subpathway_name_priors,
            xenobiotic_classifier_model=xenobiotic_classifier_model,
        )
        validate_assigned_labels(summary_df, "predicted_pathway", "predicted_subpathway", allowed_pairs)

        lane_summary_map[lane] = summary_df
        if not evidence_df.empty:
            evidence_frames.append(evidence_df)

        summary_path = output_dir_path / f"{lane}_pathway_assignments.csv"
        summary_df.to_csv(summary_path, index=False)
        output_files[f"{lane}_summary_csv"] = str(summary_path)

    consensus_df, consensus_evidence_df = build_consensus_outputs(
        feature_keys=feature_keys,
        lane_candidate_map=lane_candidates,
        lane_summary_map=lane_summary_map,
        neighbor_df=neighbor_df,
        vote_similarity_power=vote_similarity_power,
        pathway_similarity_power_map=pathway_similarity_power_map,
        pathway_mw_stats=pathway_mw_stats,
        feature_query_name_map=feature_query_name_map,
        subpathway_name_priors=subpathway_name_priors,
        xenobiotic_classifier_model=xenobiotic_classifier_model,
    )
    validate_assigned_labels(
        consensus_df,
        "consensus_predicted_pathway",
        "consensus_predicted_subpathway",
        allowed_pairs,
    )
    if not consensus_evidence_df.empty:
        evidence_frames.append(consensus_evidence_df)

    consensus_path = output_dir_path / "consensus_pathway_assignments.csv"
    consensus_df.to_csv(consensus_path, index=False)
    output_files["consensus_summary_csv"] = str(consensus_path)

    evidence_path = output_dir_path / "pathway_assignment_evidence.csv"
    if evidence_frames:
        evidence_df = pd.concat(evidence_frames, ignore_index=True)
    else:
        evidence_df = pd.DataFrame()
    serialize_json_columns(evidence_df, []).to_csv(evidence_path, index=False)
    output_files["evidence_csv"] = str(evidence_path)

    run_summary = {
        "delivery_manifest": str(manifest_path),
        "reference_db": str(reference_db_path),
        "output_dir": str(output_dir_path),
        "config": {
            "top_structures_per_feature": int(top_structures_per_feature),
            "top_k_neighbors": int(top_k_neighbors),
            "neighbor_similarity_floor": float(neighbor_similarity_floor),
            "vote_similarity_power": float(vote_similarity_power),
            "pathway_similarity_power_map": {
                key: float(pathway_similarity_power_map[key]) for key in sorted(pathway_similarity_power_map)
            },
            "low_similarity_penalty_threshold": float(LOW_SIMILARITY_PENALTY_THRESHOLD),
            "low_similarity_penalty_floor": float(LOW_SIMILARITY_PENALTY_FLOOR),
            "mw_prior_floor": float(MW_PRIOR_FLOOR),
            "subpathway_tie_margin": float(SUBPATHWAY_TIE_MARGIN),
            "xenobiotic_classifier_threshold": float(XENOBIOTIC_CLASSIFIER_THRESHOLD),
        },
        "reference_library": {
            "eligible_metabolites": int(reference_df.shape[0]),
            "n_unique_pathways": int(reference_df["metabolon_super"].nunique()),
            "n_unique_subpathways": int(reference_df["metabolon_sub"].nunique()),
        },
        "inputs": {
            lane: {
                "policy": str(lane_inputs[lane]["policy"]),
                "ranked_policy_rows_csv": str(lane_inputs[lane]["ranked_policy_rows_csv"]),
                "n_selected_structures": int(lane_candidates[lane].shape[0]),
            }
            for lane in OFFICIAL_LANES
        },
        "lanes": {
            lane: summarize_output_table(
                lane_summary_map[lane],
                pathway_col="predicted_pathway",
                subpathway_col="predicted_subpathway",
            )
            for lane in OFFICIAL_LANES
        },
        "consensus": summarize_output_table(
            consensus_df,
            pathway_col="consensus_predicted_pathway",
            subpathway_col="consensus_predicted_subpathway",
        ),
        "pathway_power_sweep": power_sweep_diagnostics,
        "xenobiotic_classifier": {
            "enabled": xenobiotic_classifier_model is not None,
            "feature_count": int(len((xenobiotic_classifier_model or {}).get("features", []))),
        },
        "output_files": output_files,
    }
    if evaluate:
        run_summary["evaluation"] = evaluate_assignment_outputs(
            lane_summary_map=lane_summary_map,
            consensus_df=consensus_df,
            known_x_csv=known_x_csv,
            known_structures_csv=known_structures_csv,
        )
    run_summary_path = output_dir_path / "pathway_assignment_run_summary.json"
    output_files["run_summary_json"] = str(run_summary_path)
    run_summary["output_files"] = output_files
    write_json(run_summary_path, run_summary)
    return run_summary


def main() -> None:
    args = parse_args()
    if args.regression_check and not args.evaluate:
        raise ValueError("--regression-check requires --evaluate.")
    run_summary = run_assignment(
        delivery_manifest=args.delivery_manifest,
        reference_db=args.reference_db,
        output_dir=args.output_dir,
        top_structures_per_feature=args.top_structures_per_feature,
        top_k_neighbors=args.top_k_neighbors,
        neighbor_similarity_floor=args.neighbor_similarity_floor,
        vote_similarity_power=args.vote_similarity_power,
        evaluate=bool(args.evaluate),
        known_x_csv=args.known_x_csv,
        known_structures_csv=args.known_structures_csv,
        power_sweep_csv=args.power_sweep_csv,
        disable_pathway_power_sweep=bool(args.disable_pathway_power_sweep),
    )
    if args.evaluate and "evaluation" in run_summary:
        print_evaluation_report(run_summary["evaluation"])
        if args.regression_check:
            exit_code = run_regression_check(
                evaluation_payload=run_summary["evaluation"],
                baseline_json_path=args.regression_baseline_json,
                fail_drop_pp=float(args.regression_fail_drop_pp),
            )
            if exit_code != 0:
                raise SystemExit(exit_code)
    print(json.dumps(run_summary["output_files"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
