"""First-class neighbor-consensus pathway vote over the reference DB.

Replaces the legacy EC2-side builder (``build_v14_all_library_prediction_review_v1``)
that used to cast the neighbor votes behind pathway enrichment. Without it the
enrichment stage silently degraded to exact-SMILES matches only, dropping the
mechanism that produced most delivered assignments. This module implements the
documented mechanism in-repo:

- Morgan radius-2, 2048-bit fingerprints (RDKit),
- k=10 nearest reference metabolites at Tanimoto floor 0.30,
- similarity**2 vote weights; subpathway voted only within the winning super.

The tally itself is shared with the report layer
(``structural_elucidation.reports.pathway_explain.pathway_vote``), so the vote
that assigns a pathway and the vote a report explains are one implementation.

Status and tier mapping (vendored policy, documented here because the legacy
builder's exact thresholds were never recorded; the parity tool
``tools/pathway_vote_parity.py`` quantifies any divergence from delivered
assignments):

- ``strong_similarity_vote``: best neighbor similarity >= 0.5,
- ``weak_similarity_backfill``: floor <= best similarity < 0.5,
- ``no_reference_match``: no neighbor at or above the floor.
- Tier ``high``: best similarity >= 0.7 and winning vote share >= 0.5;
  ``medium``: best similarity >= 0.5; ``low`` otherwise.

The production builder additionally applied molecular-weight priors and SMARTS
overrides that were never documented; they are deliberately not reproduced.
"""

from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from structural_elucidation.reports.pathway_explain import (
    VOTE_K,
    VOTE_SIMILARITY_FLOOR,
    VOTE_WEIGHT_POWER,
    PathwayVote,
    pathway_vote,
)
from structural_elucidation.utils.chemistry import canonicalize_smiles

STRONG_SIMILARITY = 0.5
HIGH_TIER_SIMILARITY = 0.7
HIGH_TIER_SHARE = 0.5
MAX_PROVENANCE_VOTERS = 10


@dataclass
class ReferenceIndex:
    """Fingerprinted reference metabolites with Metabolon labels."""

    fingerprints: list[Any]
    rows: list[dict[str, str]]  # smiles, metabolon_super, metabolon_sub

    def __len__(self) -> int:
        return len(self.rows)


def load_reference_rows(reference_db: Path) -> pd.DataFrame:
    """Reference metabolites with a super label (sub may be empty).

    Broader than ``pathway.load_reference_pathways`` (which requires a sub
    label too) because a neighbor with only a super label still carries super
    information; this matches the population the audit's leave-one-out votes
    over.
    """

    query = (
        "SELECT smiles, metabolon_super, metabolon_sub FROM metabolites "
        "WHERE TRIM(COALESCE(smiles,'')) <> '' AND TRIM(COALESCE(metabolon_super,'')) <> ''"
    )
    connection = sqlite3.connect(str(reference_db))
    try:
        reference = pd.read_sql_query(query, connection)
    finally:
        connection.close()
    reference["metabolon_super"] = reference["metabolon_super"].fillna("").astype(str).str.strip()
    reference["metabolon_sub"] = reference["metabolon_sub"].fillna("").astype(str).str.strip()
    return reference


def _fingerprint(smiles: str):
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator

    mol = Chem.MolFromSmiles(str(smiles or "").split("|")[0].strip())
    if mol is None:
        return None
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return generator.GetFingerprint(mol)


def build_reference_index(reference: pd.DataFrame) -> ReferenceIndex:
    """Fingerprint every parseable reference SMILES once, keeping labels aligned."""

    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprints: list[Any] = []
    rows: list[dict[str, str]] = []
    for _, row in reference.iterrows():
        smiles = str(row.get("smiles") or "").split("|")[0].strip()
        if not smiles:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        fingerprints.append(generator.GetFingerprint(mol))
        rows.append(
            {
                "smiles": smiles,
                "metabolon_super": str(row.get("metabolon_super") or "").strip(),
                "metabolon_sub": str(row.get("metabolon_sub") or "").strip(),
            }
        )
    return ReferenceIndex(fingerprints=fingerprints, rows=rows)


def vote_for_smiles(
    smiles: str,
    index: ReferenceIndex,
    *,
    k: int = VOTE_K,
    floor: float = VOTE_SIMILARITY_FLOOR,
    power: float = VOTE_WEIGHT_POWER,
) -> PathwayVote:
    """Score one target against the reference index and run the shared vote."""

    from rdkit.DataStructs import BulkTanimotoSimilarity

    target = _fingerprint(smiles)
    if target is None or not index.fingerprints:
        return PathwayVote(floor=floor, k=k, power=power)
    similarities = BulkTanimotoSimilarity(target, index.fingerprints)
    scored_rows = (
        {
            "similarity_to_target": similarity,
            "smiles": row["smiles"],
            "metabolon_super": row["metabolon_super"],
            "metabolon_sub": row["metabolon_sub"],
        }
        for similarity, row in zip(similarities, index.rows)
        if similarity >= floor
    )
    return pathway_vote(scored_rows, floor=floor, k=k, power=power)


def _no_match_row(smiles: str) -> dict[str, Any]:
    return {
        "canonical_smiles": smiles,
        "predicted_pathway": "",
        "predicted_subpathway": "",
        "pathway_confidence_tier": "low",
        "pathway_assignment_status": "no_reference_match",
        "pathway_vote_share": 0.0,
        "pathway_best_reference_similarity": 0.0,
        "pathway_vote_voters_json": "",
    }


def assignment_row(smiles: str, vote: PathwayVote) -> dict[str, Any]:
    """Map one vote onto the enrichment assignment-map schema."""

    if not vote.available:
        return _no_match_row(smiles)
    best_similarity = float(vote.best_similarity or 0.0)
    if best_similarity >= STRONG_SIMILARITY:
        status = "strong_similarity_vote"
    else:
        status = "weak_similarity_backfill"
    if best_similarity >= HIGH_TIER_SIMILARITY and vote.super_share >= HIGH_TIER_SHARE:
        tier = "high"
    elif best_similarity >= STRONG_SIMILARITY:
        tier = "medium"
    else:
        tier = "low"
    voters = [voter.as_dict() for voter in vote.voters[:MAX_PROVENANCE_VOTERS]]
    return {
        "canonical_smiles": smiles,
        "predicted_pathway": vote.winning_super,
        "predicted_subpathway": vote.winning_sub,
        "pathway_confidence_tier": tier,
        "pathway_assignment_status": status,
        "pathway_vote_share": round(vote.super_share, 4),
        "pathway_best_reference_similarity": round(best_similarity, 4),
        "pathway_vote_voters_json": json.dumps(voters, sort_keys=True),
    }


def vote_assignment_map(
    canonical_smiles: list[str],
    reference: pd.DataFrame | Path,
    *,
    k: int = VOTE_K,
    floor: float = VOTE_SIMILARITY_FLOOR,
    power: float = VOTE_WEIGHT_POWER,
) -> pd.DataFrame:
    """Neighbor-vote assignments for a list of canonical SMILES.

    ``reference`` is either the loaded reference DataFrame (smiles,
    metabolon_super, metabolon_sub) or a path to the synonym SQLite.
    """

    if isinstance(reference, (str, Path)):
        reference = load_reference_rows(Path(reference))
    index = build_reference_index(reference)
    rows: list[dict[str, Any]] = []
    for raw in canonical_smiles:
        smiles = str(raw or "").strip()
        if not smiles:
            continue
        canonical = canonicalize_smiles(smiles) or smiles
        rows.append(assignment_row(canonical, vote_for_smiles(canonical, index, k=k, floor=floor, power=power)))
    return pd.DataFrame(rows)
