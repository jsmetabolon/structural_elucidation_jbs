#!/usr/bin/env python3
"""Add coarse lipid class annotations to equalization rerank top-10 CSVs.

Parses Metabolon-style lipid names and SMILES substructures (when RDKit is
available) to emit lipid class, chain totals, and shorthand notation.

Pipeline context:
    Second enrichment stage; reads pathway-enriched top-10 CSVs from
    ``pipeline_equalization/rerank_lib*/``.

Key I/O:
    Inputs: ``*_with_pathways.csv`` per library.
    Outputs: ``*_lipid_class.csv``, ``lipid_class_enrichment_summary.json``.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from structural_elucidation.paths import REPO_ROOT, resolve

try:
    from rdkit import Chem
    from rdkit.Chem.rdchem import BondType

    HAS_RDKIT = True
except ImportError:  # pragma: no cover - optional in lightweight environments
    Chem = None  # type: ignore[assignment]
    BondType = None  # type: ignore[assignment]
    HAS_RDKIT = False


# ── Defaults and output schema ───────────────────────────────────────
DEFAULT_EQUALIZATION_ROOT: Path = Path("models/pipeline_equalization")
DEFAULT_LIBS: str = "209,305,400,402"
DEFAULT_INPUT_NAME: str = "unified_rerank_all_top10_with_pathways.csv"
DEFAULT_OUTPUT_NAME: str = "unified_rerank_all_top10_with_pathways_lipid_class.csv"
LIPID_COLUMNS: list[str] = [
    "lipid_class",
    "lipid_chains_total",
    "lipid_shorthand",
    "lipid_classification_source",
]

# ── Name-parse dictionaries and regex patterns ───────────────────────
NAME_COLUMNS: list[str] = [
    "library_name",
    "candidate_name",
    "name",
    "compound_name",
    "known_compound_name",
    "annotation_text",
    "feature_label",
]

ACYL_TO_CHAIN: dict[str, str] = {  # common acyl prefix → shorthand chain token
    "palmitoyl": "16:0",
    "palmitoleoyl": "16:1",
    "palmitelaidoyl": "16:1",
    "stearoyl": "18:0",
    "oleoyl": "18:1",
    "linoleoyl": "18:2",
    "linolenoyl": "18:3",
    "arachidonoyl": "20:4",
    "archidonoyl": "20:4",
    "arachidonyl": "20:4",
    "arachidoyl": "20:0",
    "arachidyl": "20:0",
    "eicosatrienoyl": "20:3",
    "eicosenoyl": "20:1",
    "docosahexaenoyl": "22:6",
    "docosapentaenoyl": "22:5",
    "dihomo-linolenoyl": "20:3",
    "dihomo-linoleoyl": "20:2",
    "myristoyl": "14:0",
    "lauroyl": "12:0",
    "pentadecanoyl": "15:0",
    "heptadecanoyl": "17:0",
    "margaroyl": "17:0",
    "nonadecanoyl": "19:0",
    "nervonoyl": "24:1",
    "lignoceroyl": "24:0",
    "behenoyl": "22:0",
    "erucoyl": "22:1",
}

HEADGROUP_NAME_PATTERNS: list[tuple[re.Pattern[str], str, str]] = [  # (regex, full_class, lyso_class)
    (re.compile(r"(?:\bGPC\b|phospho\s*choline|phosphocholine|glycerophospho\s*choline)", re.I), "PC", "LPC"),
    (
        re.compile(r"(?:\bGPE\b|phospho\s*ethanolamine|phosphoethanolamine|glycerophospho\s*ethanolamine)", re.I),
        "PE",
        "LPE",
    ),
    (re.compile(r"(?:\bGPS\b|phospho\s*serine|phosphoserine|glycerophospho\s*serine)", re.I), "PS", "LPS"),
    (re.compile(r"(?:\bGPI\b|phospho\s*inositol|phosphoinositol|glycerophospho\s*inositol)", re.I), "PI", "LPI"),
    (re.compile(r"(?:\bGPG\b|phospho\s*glycerol|phosphoglycerol|glycerophospho\s*glycerol)", re.I), "PG", "LPG"),
]
DIRECT_CLASSES: tuple[str, ...] = (
    "HexCer",
    "Cer",
    "SM",
    "LPC",
    "LPE",
    "LPS",
    "LPI",
    "LPG",
    "PC",
    "PE",
    "PS",
    "PI",
    "PG",
    "DAG",
    "TAG",
    "MAG",
)
CLASS_ALIASES: dict[str, str] = {
    "DG": "DAG",
    "TG": "TAG",
    "MG": "MAG",
    "CER": "Cer",
    "HEXCER": "HexCer",
    "STEROL": "sterol",
}
CHAIN_RE: re.Pattern[str] = re.compile(r"(?P<prefix>O-|P-|[dmte])?(?P<carbons>\d{1,2}):(?P<double_bonds>\d{1,2})", re.I)
DIRECT_PAREN_RE: re.Pattern[str] = re.compile(
    r"\b(?P<class>HexCer|Cer|SM|LPC|LPE|LPS|LPI|LPG|PC|PE|PS|PI|PG|DAG|TAG|MAG|DG|TG|MG)\s*"
    r"\((?P<body>[^)]*\d{1,2}:\d{1,2}[^)]*)\)",
    re.I,
)
DIRECT_SPACE_RE: re.Pattern[str] = re.compile(
    r"\b(?P<class>HexCer|Cer|SM|LPC|LPE|LPS|LPI|LPG|PC|PE|PS|PI|PG|DAG|TAG|MAG|DG|TG|MG)\s+"
    r"(?P<body>(?:O-|P-|[dmte])?\d{1,2}:\d{1,2}(?:\([^)]*\))?"
    r"(?:[/_ -](?:O-|P-|[dmte])?\d{1,2}:\d{1,2}(?:\([^)]*\))?)*)",
    re.I,
)


# ── Data classes ─────────────────────────────────────────────────────
@dataclass(frozen=True)
class Chain:
    """One acyl or sphingoid chain token (e.g. 18:1)."""

    carbons: int
    double_bonds: int
    prefix: str = ""

    @property
    def token(self) -> str:
        return f"{self.prefix}{self.carbons}:{self.double_bonds}"


@dataclass(frozen=True)
class LipidAnnotation:
    """Coarse lipid classification result for one candidate row."""

    lipid_class: str = ""
    lipid_chains_total: str = ""
    lipid_shorthand: str = ""
    lipid_classification_source: str = ""

    @property
    def is_empty(self) -> bool:
        return not self.lipid_class

    def as_dict(self) -> dict[str, str]:
        """Serialize annotation fields for DataFrame assignment."""
        return {
            "lipid_class": self.lipid_class,
            "lipid_chains_total": self.lipid_chains_total,
            "lipid_shorthand": self.lipid_shorthand,
            "lipid_classification_source": self.lipid_classification_source,
        }


def parse_args() -> argparse.Namespace:
    """Build the lipid-class enrichment CLI argument parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--equalization-root", type=Path, default=DEFAULT_EQUALIZATION_ROOT)
    parser.add_argument("--libs", nargs="+", default=[DEFAULT_LIBS])
    parser.add_argument("--input-name", default=DEFAULT_INPUT_NAME)
    parser.add_argument("--output-name", default=DEFAULT_OUTPUT_NAME)
    parser.add_argument("--summary-json", type=Path, default=None)
    parser.add_argument("--allow-missing", action="store_true")
    return parser.parse_args()


def parse_libs(values: list[str]) -> list[int]:
    libs: list[int] = []
    for value in values:
        for token in str(value).split(","):
            token = token.strip()
            if token:
                libs.append(int(token))
    return libs


def clean_text(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return str(value).strip()


def clean_display_name(value: object) -> str:
    text = clean_text(value)
    text = re.sub(r"\s+", " ", text.replace("_", " ")).strip()
    return text


def normalize_class(value: str) -> str:
    text = clean_text(value)
    if not text:
        return ""
    upper = text.upper()
    if upper in CLASS_ALIASES:
        return CLASS_ALIASES[upper]
    for klass in DIRECT_CLASSES:
        if upper == klass.upper():
            return klass
    return text


def parse_chain_token(value: str) -> Chain | None:
    token = clean_text(value)
    token = re.sub(r"\([^)]*\)", "", token)
    match = CHAIN_RE.search(token)
    if not match:
        return None
    prefix = clean_text(match.group("prefix"))
    if prefix:
        prefix = prefix[0].lower() if len(prefix) == 1 else prefix.upper()
    return Chain(int(match.group("carbons")), int(match.group("double_bonds")), prefix)


def extract_chains(value: str) -> list[Chain]:
    chains: list[Chain] = []
    for match in CHAIN_RE.finditer(clean_text(value)):
        chain = parse_chain_token(match.group(0))
        if chain:
            chains.append(chain)
    return chains


def total_chain_notation(chains: list[Chain]) -> str:
    if not chains:
        return ""
    return f"{sum(chain.carbons for chain in chains)}:{sum(chain.double_bonds for chain in chains)}"


def shorthand_for(lipid_class: str, chains: list[Chain], *, source: str) -> LipidAnnotation:
    klass = normalize_class(lipid_class)
    if not klass:
        return LipidAnnotation()
    if klass == "sterol":
        return LipidAnnotation("sterol", "", "sterol", source)
    if not chains:
        return LipidAnnotation(klass, "", klass, source)
    total = total_chain_notation(chains)
    if klass in {"SM", "Cer", "HexCer"} and len(chains) >= 2:
        detail = "/".join(chain.token for chain in chains)
        shorthand = f"{klass} {detail}"
    else:
        shorthand = f"{klass} {total}"
    return LipidAnnotation(klass, total, shorthand, source)


def parse_direct_lipid_token(value: object) -> LipidAnnotation:
    text = clean_display_name(value)
    if not text:
        return LipidAnnotation()
    search_text = re.sub(r"(?<=\d:\d)\([^)]*\)", "", text)
    for pattern in (DIRECT_PAREN_RE, DIRECT_SPACE_RE):
        for match in pattern.finditer(search_text):
            chains = extract_chains(match.group("body"))
            if chains:
                return shorthand_for(match.group("class"), chains, source="name_parse")
    sm_match = re.search(r"\b(?P<c>\d{1,2}):(?P<db>\d{1,2})\s+SM\b", search_text, flags=re.I)
    if sm_match:
        chains = [Chain(int(sm_match.group("c")), int(sm_match.group("db")))]
        return shorthand_for("SM", chains, source="name_parse")
    return LipidAnnotation()


def acyl_prefix_to_chain(prefix: str) -> Chain | None:
    text = clean_display_name(prefix).lower().rstrip("-")
    if not text:
        return None
    candidates = [text]
    if not text.endswith(("oyl", "yl")):
        candidates.extend([f"{text}oyl", f"{text}yl"])
    if text.startswith("di"):
        candidates.append(text[2:])
        if not text[2:].endswith(("oyl", "yl")):
            candidates.append(f"{text[2:]}oyl")
    for candidate in candidates:
        chain_text = ACYL_TO_CHAIN.get(candidate)
        if chain_text:
            return parse_chain_token(chain_text)
    return None


def find_named_acyl_chains(text: str) -> list[Chain]:
    hits: list[tuple[int, Chain]] = []
    lowered = text.lower()
    for acyl, chain_text in ACYL_TO_CHAIN.items():
        pattern = re.compile(rf"\b{re.escape(acyl)}(?=\b|-|glycerophospho|sphingo|ceramide)", re.I)
        for match in pattern.finditer(lowered):
            chain = parse_chain_token(chain_text)
            if chain:
                hits.append((match.start(), chain))
    hits.sort(key=lambda item: item[0])
    out: list[Chain] = []
    seen_positions: set[int] = set()
    for pos, chain in hits:
        if pos in seen_positions:
            continue
        seen_positions.add(pos)
        out.append(chain)
    return out


def detect_named_headgroup(text: str, chain_count: int) -> str:
    for pattern, full_class, lyso_class in HEADGROUP_NAME_PATTERNS:
        if pattern.search(text):
            lyso_hint = bool(re.search(r"\blyso|glycerophospho", text, flags=re.I))
            return lyso_class if chain_count <= 1 or lyso_hint else full_class
    if re.search(r"\b(?:diacylglycerol|DAG)\b", text, flags=re.I):
        return "DAG"
    if re.search(r"\b(?:triacylglycerol|triglyceride|TAG)\b", text, flags=re.I):
        return "TAG"
    return ""


def parse_metabolon_lipid_name(value: object) -> LipidAnnotation:
    text = clean_display_name(value)
    if not text:
        return LipidAnnotation()
    direct = parse_direct_lipid_token(text)
    if not direct.is_empty:
        return direct
    if re.search(r"\b(?:cholesterol|cholesteryl|sterol|steroid)\b", text, flags=re.I):
        return LipidAnnotation("sterol", "", "sterol", "name_parse")
    sphingo_match = re.search(r"\b(?:sphingomyelin|SM)\b", text, flags=re.I)
    if sphingo_match:
        chains = find_named_acyl_chains(text)
        if chains:
            return shorthand_for("SM", [Chain(18, 1, "d"), chains[0]], source="name_parse")
    cer_match = re.search(r"\b(?:ceramide|erythrosphingosine|erythrosphinganine)\b", text, flags=re.I)
    if cer_match:
        chains = find_named_acyl_chains(text)
        if chains:
            base = Chain(18, 0 if re.search(r"sphinganine", text, flags=re.I) else 1, "d")
            return shorthand_for("Cer", [base, chains[0]], source="name_parse")
    chains = find_named_acyl_chains(text)
    if chains:
        headgroup = detect_named_headgroup(text, len(chains))
        if headgroup:
            return shorthand_for(headgroup, chains[:3], source="name_parse")
    return LipidAnnotation()


def annotate_from_name(row: pd.Series) -> LipidAnnotation:
    for column in NAME_COLUMNS:
        if column not in row.index:
            continue
        annotation = parse_metabolon_lipid_name(row.get(column))
        if not annotation.is_empty:
            return annotation
    return LipidAnnotation()


@lru_cache(maxsize=128)
def compiled_smarts(smarts: str):
    if not HAS_RDKIT:
        return None
    return Chem.MolFromSmarts(smarts)


def mol_from_smiles(smiles: object):
    if not HAS_RDKIT:
        return None
    text = clean_text(smiles)
    if not text:
        return None
    try:
        return Chem.MolFromSmiles(text)
    except Exception:
        return None


def has_substructure(mol: Any, smarts: str) -> bool:
    pattern = compiled_smarts(smarts)
    return bool(mol is not None and pattern is not None and mol.HasSubstructMatch(pattern))


def count_substructures(mol: Any, smarts: str) -> int:
    pattern = compiled_smarts(smarts)
    if mol is None or pattern is None:
        return 0
    return len(mol.GetSubstructMatches(pattern, uniquify=True))


def is_carbonyl_carbon(atom: Any) -> bool:
    if atom.GetAtomicNum() != 6:
        return False
    has_double_o = False
    has_single_hetero = False
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetAtomicNum() == 8 and bond.GetBondType() == BondType.DOUBLE:
            has_double_o = True
        if other.GetAtomicNum() in {7, 8} and bond.GetBondType() == BondType.SINGLE:
            has_single_hetero = True
    return has_double_o and has_single_hetero


def carbonyl_hetero_types(atom: Any) -> set[int]:
    out: set[int] = set()
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetAtomicNum() in {7, 8} and bond.GetBondType() == BondType.SINGLE:
            out.add(other.GetAtomicNum())
    return out


def carbon_chain_path_from(mol: Any, start_idx: int, blocked: set[int]) -> list[int]:
    def valid_carbon(idx: int) -> bool:
        atom = mol.GetAtomWithIdx(idx)
        return atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() and idx not in blocked

    def dfs(idx: int, visited: tuple[int, ...]) -> list[int]:
        best = list(visited)
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in visited or not valid_carbon(n_idx):
                continue
            bond = mol.GetBondBetweenAtoms(idx, n_idx)
            if bond is None or bond.IsInRing():
                continue
            candidate = dfs(n_idx, visited + (n_idx,))
            if len(candidate) > len(best):
                best = candidate
        return best

    return dfs(start_idx, (start_idx,)) if valid_carbon(start_idx) else []


def count_double_bonds_on_path(mol: Any, path: list[int]) -> int:
    count = 0
    for left, right in zip(path, path[1:]):
        bond = mol.GetBondBetweenAtoms(int(left), int(right))
        if bond is not None and bond.GetBondType() == BondType.DOUBLE:
            count += 1
    return count


def acyl_chains_from_smiles(mol: Any, hetero_atomic_num: int | None = None) -> list[Chain]:
    chains: list[Chain] = []
    seen_carbonyls: set[int] = set()
    for atom in mol.GetAtoms():
        if not is_carbonyl_carbon(atom):
            continue
        heteros = carbonyl_hetero_types(atom)
        if hetero_atomic_num is not None and hetero_atomic_num not in heteros:
            continue
        carbonyl_idx = atom.GetIdx()
        if carbonyl_idx in seen_carbonyls:
            continue
        best_path: list[int] = []
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() != 6:
                continue
            path = carbon_chain_path_from(mol, n_idx, {carbonyl_idx})
            if len(path) > len(best_path):
                best_path = path
        if not best_path:
            continue
        seen_carbonyls.add(carbonyl_idx)
        chains.append(Chain(1 + len(best_path), count_double_bonds_on_path(mol, best_path)))
    return sorted(chains, key=lambda chain: (chain.carbons, chain.double_bonds), reverse=True)


def has_sugar_like_ring(mol: Any) -> bool:
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) not in {5, 6}:
            continue
        atomic_numbers = [mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring]
        if any(number not in {6, 8} for number in atomic_numbers):
            continue
        ring_set = set(ring)
        exocyclic_oxygen_neighbors: set[int] = set()
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring_set and neighbor.GetAtomicNum() == 8:
                    exocyclic_oxygen_neighbors.add(neighbor.GetIdx())
        if len(exocyclic_oxygen_neighbors) >= 3:
            return True
    return False


def is_sterol_like(mol: Any) -> bool:
    if has_substructure(mol, "C1CCC2(C(C1)CCC1C2CCC2(C)C(CC1)CC2)C"):
        return True
    rings = mol.GetRingInfo().AtomRings()
    carbon_rings = [ring for ring in rings if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring)]
    return len(carbon_rings) >= 4 and sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) >= 17


def has_ceramide_backbone(mol: Any) -> bool:
    return any(
        has_substructure(mol, smarts)
        for smarts in (
            "[CX3](=O)[NX3][CH]([CH2][OX2])[CH]([OX2])",
            "[CX3](=O)[NX3][C]([CH2][OX2])[C]([OX2])",
        )
    )


def detect_smiles_headgroup(mol: Any) -> str:
    ester_count = count_substructures(mol, "[CX3](=O)[OX2H0][#6]")
    amide_count = count_substructures(mol, "[CX3](=O)[NX3]")
    if has_substructure(mol, "[PX4](=[OX1])([OX2H,OX1-])([OX2][#6])[OX2][CH2][CH2][NX4+]([CH3])([CH3])[CH3]") or has_substructure(
        mol, "[PX4](=[OX1])([OX2][#6])([OX2][#6])[OX2][CH2][CH2][NX4+]([CH3])([CH3])[CH3]"
    ):
        return "SM" if amide_count else ("LPC" if ester_count <= 1 else "PC")
    if has_substructure(mol, "[PX4](=[OX1])([OX2H,OX1-])([OX2][#6])[OX2][CH2][CH2][NX3;!$(N=*)]") or has_substructure(
        mol, "[PX4](=[OX1])([OX2][#6])([OX2][#6])[OX2][CH2][CH2][NX3;!$(N=*)]"
    ):
        return "LPE" if ester_count <= 1 else "PE"
    if has_substructure(mol, "[PX4](=[OX1])([OX2H,OX1-])([OX2][#6])[OX2][CH2][CH]([NX3])[CX3](=O)[OX2H,OX1-]"):
        return "LPS" if ester_count <= 1 else "PS"
    if has_substructure(mol, "[PX4](=[OX1])([OX2H,OX1-])([OX2][#6])[OX2][CH2][CH]([OX2H])[CH2][OX2H]"):
        return "LPG" if ester_count <= 1 else "PG"
    if has_substructure(mol, "[PX4](=[OX1])([OX2H,OX1-])([OX2][#6])[OX2][C@H]1[C@H]([OX2H])[C@H]([OX2H])[C@H]([OX2H])[C@H]([OX2H])[CH2]1"):
        return "LPI" if ester_count <= 1 else "PI"
    if amide_count and has_ceramide_backbone(mol):
        return "HexCer" if has_sugar_like_ring(mol) else "Cer"
    if ester_count >= 2 and has_substructure(mol, "[OX2][CH2][CH]([OX2])[CH2][OX2]"):
        return "TAG" if ester_count >= 3 else "DAG"
    if is_sterol_like(mol):
        return "sterol"
    return ""


def annotate_from_smiles(smiles: object) -> LipidAnnotation:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return LipidAnnotation()
    lipid_class = detect_smiles_headgroup(mol)
    if not lipid_class:
        return LipidAnnotation()
    if lipid_class == "sterol":
        return LipidAnnotation("sterol", "", "sterol", "smiles_smarts")
    ester_chains = acyl_chains_from_smiles(mol, hetero_atomic_num=8)
    amide_chains = acyl_chains_from_smiles(mol, hetero_atomic_num=7)
    if lipid_class in {"SM", "Cer", "HexCer"}:
        acyl = amide_chains[0] if amide_chains else None
        chains = [Chain(18, 1, "d")] + ([acyl] if acyl else [])
    elif lipid_class in {"PC", "LPC", "PE", "LPE", "PS", "LPS", "PI", "LPI", "PG", "LPG", "DAG", "TAG"}:
        chains = ester_chains
    else:
        chains = ester_chains or amide_chains
    return shorthand_for(lipid_class, chains, source="smiles_smarts")


def annotate_row(row: pd.Series) -> LipidAnnotation:
    """Classify one row: prefer name parse, fall back to SMILES SMARTS."""
    name_annotation = annotate_from_name(row)
    if not name_annotation.is_empty:
        return name_annotation
    for column in ("predicted_smiles", "canonical_smiles", "smiles", "candidate_smiles"):
        if column in row.index:
            smiles_annotation = annotate_from_smiles(row.get(column))
            if not smiles_annotation.is_empty:
                return smiles_annotation
    return LipidAnnotation()


# ── I/O and enrichment driver ────────────────────────────────────────
def rerank_path(root: Path, lib_id: int, name: str) -> Path:
    """Path to a named CSV under ``rerank_lib{lib_id}/``."""
    return root / f"rerank_lib{int(lib_id)}" / name


def load_inputs(root: Path, libs: list[int], input_name: str, allow_missing: bool) -> dict[int, pd.DataFrame]:
    frames: dict[int, pd.DataFrame] = {}
    missing: list[Path] = []
    for lib_id in libs:
        path = rerank_path(root, lib_id, input_name)
        if not path.exists():
            missing.append(path)
            continue
        frames[int(lib_id)] = pd.read_csv(path, low_memory=False)
    if missing and not allow_missing:
        joined = "\n".join(f"  {path}" for path in missing)
        raise FileNotFoundError(f"Missing rerank input CSVs:\n{joined}")
    return frames


def ordered_columns(original_columns: list[str], enriched: pd.DataFrame) -> list[str]:
    columns = list(original_columns)
    for column in LIPID_COLUMNS:
        if column not in columns and column in enriched.columns:
            columns.append(column)
    return columns + [column for column in enriched.columns if column not in columns]


def enrich_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Add lipid class columns to all rows in a rerank top-10 frame."""
    original_columns = list(frame.columns)
    annotations = [annotate_row(row).as_dict() for _, row in frame.iterrows()]
    out = frame.copy()
    annotation_frame = pd.DataFrame(annotations, index=out.index)
    for column in LIPID_COLUMNS:
        out[column] = annotation_frame[column].fillna("")
    return out.loc[:, ordered_columns(original_columns, out)]


def summarize_lib(lib_id: int, input_path: Path, output_path: Path, frame: pd.DataFrame) -> dict[str, Any]:
    lipid_class = frame["lipid_class"].fillna("").astype(str)
    shorthand = frame["lipid_shorthand"].fillna("").astype(str)
    source = frame["lipid_classification_source"].fillna("").astype(str)
    return {
        "lib_id": int(lib_id),
        "input_csv": str(input_path),
        "output_csv": str(output_path),
        "n_rows": int(len(frame)),
        "n_rows_with_lipid_class": int(lipid_class.ne("").sum()),
        "n_rows_with_lipid_shorthand": int(shorthand.ne("").sum()),
        "lipid_class_counts": {str(k): int(v) for k, v in lipid_class.loc[lipid_class.ne("")].value_counts().sort_index().items()},
        "lipid_classification_source_counts": {str(k): int(v) for k, v in source.loc[source.ne("")].value_counts().sort_index().items()},
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    """Run lipid-class enrichment for all requested libraries."""
    repo_root = REPO_ROOT
    root = resolve(args.equalization_root, repo_root)
    libs = parse_libs(args.libs)
    if not libs:
        raise ValueError("No libraries requested")
    frames = load_inputs(root, libs, args.input_name, bool(args.allow_missing))
    if not frames:
        raise FileNotFoundError(f"No rerank inputs found under {root}")
    lib_summaries: list[dict[str, Any]] = []
    for lib_id in libs:
        if lib_id not in frames:
            continue
        input_path = rerank_path(root, lib_id, args.input_name)
        output_path = rerank_path(root, lib_id, args.output_name)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        enriched = enrich_frame(frames[lib_id])
        enriched.to_csv(output_path, index=False)
        lib_summaries.append(summarize_lib(lib_id, input_path, output_path, enriched))
    summary_path = resolve(args.summary_json, repo_root) if args.summary_json else root / "lipid_class_enrichment_summary.json"
    summary = {
        "equalization_root": str(root),
        "input_name": args.input_name,
        "output_name": args.output_name,
        "has_rdkit": bool(HAS_RDKIT),
        "libs": lib_summaries,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def main() -> int:
    """CLI entry point."""
    run(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
