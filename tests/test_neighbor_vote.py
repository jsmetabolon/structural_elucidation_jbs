"""Unit tests for the vendored neighbor-consensus pathway vote (T1)."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("rdkit")

from structural_elucidation.pipeline.enrichment import pathway  # noqa: E402
from structural_elucidation.pipeline.enrichment.neighbor_vote import (  # noqa: E402
    assignment_row,
    build_reference_index,
    vote_assignment_map,
    vote_for_smiles,
)
from structural_elucidation.reports.pathway_explain import PathwayVote, Voter  # noqa: E402

GLYCINE = "NCC(=O)O"
ALANINE = "CC(N)C(=O)O"
SERINE = "OCC(N)C(=O)O"
AMINOBUTYRATE = "CCC(N)C(=O)O"  # not in reference; alanine-like
HEXANOATE = "CCCCCC(=O)O"
OCTANOATE = "CCCCCCCC(=O)O"
TETRACHLOROMETHANE = "ClC(Cl)(Cl)Cl"  # shares no bits with the reference set

REFERENCE = pd.DataFrame(
    [
        {"smiles": GLYCINE, "metabolon_super": "Amino Acid", "metabolon_sub": "Glycine Metabolism"},
        {"smiles": ALANINE, "metabolon_super": "Amino Acid", "metabolon_sub": "Alanine Metabolism"},
        {"smiles": SERINE, "metabolon_super": "Amino Acid", "metabolon_sub": "Serine Metabolism"},
        {"smiles": HEXANOATE, "metabolon_super": "Lipid", "metabolon_sub": "Medium Chain Fatty Acid"},
        {"smiles": OCTANOATE, "metabolon_super": "Lipid", "metabolon_sub": "Medium Chain Fatty Acid"},
    ]
)


@pytest.fixture()
def reference_db(tmp_path: Path) -> Path:
    db_path = tmp_path / "metabolite_synonyms.sqlite"
    con = sqlite3.connect(db_path)
    con.execute("CREATE TABLE metabolites (smiles TEXT, metabolon_super TEXT, metabolon_sub TEXT)")
    con.executemany(
        "INSERT INTO metabolites (smiles, metabolon_super, metabolon_sub) VALUES (?, ?, ?)",
        REFERENCE.itertuples(index=False),
    )
    con.commit()
    con.close()
    return db_path


def test_identical_structure_votes_strongly_for_its_class():
    result = vote_assignment_map([ALANINE], REFERENCE)
    row = result.iloc[0]
    assert row["predicted_pathway"] == "Amino Acid"
    assert row["pathway_assignment_status"] == "strong_similarity_vote"
    assert row["pathway_confidence_tier"] == "high"
    assert row["pathway_best_reference_similarity"] == 1.0
    assert "Alanine Metabolism" in row["pathway_vote_voters_json"]


def test_near_neighbor_query_inherits_cluster_super():
    index = build_reference_index(REFERENCE)
    vote = vote_for_smiles(AMINOBUTYRATE, index)
    assert vote.available
    assert vote.winning_super == "Amino Acid"
    # Sub is voted only within the winning super.
    assert vote.winning_sub in {"Glycine Metabolism", "Alanine Metabolism", "Serine Metabolism"}


def test_dissimilar_query_gets_no_reference_match():
    result = vote_assignment_map([TETRACHLOROMETHANE], REFERENCE)
    row = result.iloc[0]
    assert row["pathway_assignment_status"] == "no_reference_match"
    assert row["predicted_pathway"] == ""
    assert row["pathway_confidence_tier"] == "low"


def test_assignment_row_status_and_tier_thresholds():
    def vote_with(similarity: float, share: float) -> PathwayVote:
        voter = Voter(
            chemical_id="",
            name="",
            smiles="C",
            similarity=similarity,
            super_pathway="Lipid",
            sub_pathway="",
            weight=similarity**2,
        )
        return PathwayVote(
            voters=[voter],
            super_tally=[("Lipid", similarity**2)],
            winning_super="Lipid",
            super_share=share,
            total_weight=similarity**2,
        )

    assert assignment_row("C", vote_with(0.75, 0.9))["pathway_confidence_tier"] == "high"
    assert assignment_row("C", vote_with(0.75, 0.4))["pathway_confidence_tier"] == "medium"
    assert assignment_row("C", vote_with(0.55, 0.9))["pathway_confidence_tier"] == "medium"
    weak = assignment_row("C", vote_with(0.35, 0.9))
    assert weak["pathway_confidence_tier"] == "low"
    assert weak["pathway_assignment_status"] == "weak_similarity_backfill"


def test_generate_map_uses_exact_then_neighbor_vote(reference_db: Path):
    result = pathway.generate_pathway_assignment_map([ALANINE, AMINOBUTYRATE, TETRACHLOROMETHANE], reference_db)
    by_smiles = {row["canonical_smiles"]: row for _, row in result.iterrows()}
    canonical = {name: pathway.canonicalize_for_pathway(name) for name in (ALANINE, AMINOBUTYRATE, TETRACHLOROMETHANE)}

    exact = by_smiles[canonical[ALANINE]]
    assert exact["pathway_assignment_status"] == "exact_reference_match"
    assert exact["predicted_pathway"] == "Amino Acid"

    voted = by_smiles[canonical[AMINOBUTYRATE]]
    assert voted["pathway_assignment_status"] in {"strong_similarity_vote", "weak_similarity_backfill"}
    assert voted["predicted_pathway"] == "Amino Acid"

    unmatched = by_smiles[canonical[TETRACHLOROMETHANE]]
    assert unmatched["pathway_assignment_status"] == "no_reference_match"


def test_enriched_rerank_rows_carry_pathways_but_not_voter_provenance(reference_db: Path):
    assignment_map = pathway.generate_pathway_assignment_map([AMINOBUTYRATE], reference_db)
    assert "pathway_vote_voters_json" in assignment_map.columns
    frame = pd.DataFrame([{"control_id": "X-1", "rank": 1, "predicted_smiles": AMINOBUTYRATE}])
    enriched = pathway.enrich_pathway_annotations(frame, assignment_map)
    assert enriched.loc[0, "predicted_pathway"] == "Amino Acid"
    assert "pathway_vote_voters_json" not in enriched.columns
