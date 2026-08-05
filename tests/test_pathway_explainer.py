"""Tests for the reconstructed pathway vote and the section that explains it."""
from __future__ import annotations

import json
import random
import sys
from collections import defaultdict
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "structural_elucidation" / "reports"))

import confidence as C  # noqa: E402

from structural_elucidation.reports import pathway_explain as PE  # noqa: E402


def _reference(name, similarity, super_pathway, sub_pathway="", chemical_id=""):
    return {
        "chemical_id": chemical_id or name.upper(),
        "canonical_name": name,
        "smiles": "C" * (len(name) % 5 + 1),
        "similarity_to_target": similarity,
        "metabolon_super": super_pathway,
        "metabolon_sub": sub_pathway,
    }


# --- the vote ----------------------------------------------------------------


def test_votes_are_weighted_by_similarity_squared():
    vote = PE.pathway_vote(
        [
            _reference("close", 0.9, "Lipid", "Sterol"),
            _reference("mid", 0.6, "Lipid", "Fatty Acid"),
            _reference("other", 0.8, "Amino Acid", "Glycine"),
        ]
    )
    assert [round(voter.weight, 4) for voter in vote.voters] == [0.81, 0.64, 0.36]
    # Lipid 0.81 + 0.36 against Amino Acid 0.64.
    assert vote.winning_super == "Lipid"
    assert vote.super_share == pytest.approx(1.17 / 1.81)


def test_a_single_close_neighbor_outweighs_two_distant_ones():
    """The squared weight is the whole point: distance is discounted hard."""

    vote = PE.pathway_vote(
        [
            _reference("close", 0.90, "Lipid"),
            _reference("far_a", 0.55, "Xenobiotics"),
            _reference("far_b", 0.55, "Xenobiotics"),
        ]
    )
    assert vote.winning_super == "Lipid"
    # Unweighted, Xenobiotics would win two votes to one.
    assert sum(1 for v in vote.voters if v.super_pathway == "Xenobiotics") == 2


def test_neighbors_below_the_floor_do_not_vote():
    vote = PE.pathway_vote(
        [_reference("kept", 0.31, "Lipid"), _reference("dropped", 0.29, "Amino Acid")]
    )
    assert [voter.name for voter in vote.voters] == ["kept"]
    assert vote.n_above_floor == 1


def test_only_the_k_closest_neighbors_vote():
    rows = [_reference(f"n{i}", 0.9 - i * 0.01, "Lipid") for i in range(20)]
    rows.append(_reference("distant", 0.31, "Amino Acid"))
    vote = PE.pathway_vote(rows)
    assert len(vote.voters) == PE.VOTE_K
    assert vote.n_above_floor == 21
    assert "distant" not in {voter.name for voter in vote.voters}


def test_unlabeled_references_are_not_voters():
    """The audit votes over labeled references only, so an unlabeled neighbor is
    not a vote for "no pathway" - it is simply not in the electorate."""

    vote = PE.pathway_vote(
        [_reference("unlabeled", 0.99, ""), _reference("labeled", 0.4, "Lipid")]
    )
    assert [voter.name for voter in vote.voters] == ["labeled"]
    assert vote.winning_super == "Lipid"


def test_subpathway_is_decided_only_among_the_winning_pathways_voters():
    vote = PE.pathway_vote(
        [
            _reference("a", 0.9, "Lipid", "Sterol"),
            _reference("b", 0.7, "Lipid", "Sterol"),
            _reference("c", 0.85, "Amino Acid", "Glycine"),
        ]
    )
    assert vote.winning_super == "Lipid"
    assert vote.winning_sub == "Sterol"
    # Glycine belongs to the losing pathway and never enters the subpathway tally.
    assert "Glycine" not in dict(vote.sub_tally)


def test_a_pathway_whose_voters_carry_no_subpathway_yields_no_subpathway():
    vote = PE.pathway_vote([_reference("a", 0.9, "Lipid", "")])
    assert vote.winning_super == "Lipid"
    assert vote.winning_sub == ""
    assert vote.sub_share == 0.0


def test_ties_resolve_toward_the_closest_neighbor():
    vote = PE.pathway_vote(
        [_reference("closest", 0.8, "Lipid"), _reference("equal", 0.8, "Amino Acid")]
    )
    assert vote.winning_super == "Lipid"
    assert vote.super_share == pytest.approx(0.5)


def test_an_empty_electorate_produces_an_unavailable_vote():
    empty = PE.pathway_vote([])
    assert empty.available is False
    assert empty.winning_super == ""
    assert empty.best_similarity is None
    assert PE.pathway_vote([_reference("far", 0.1, "Lipid")]).available is False


def test_ledger_rows_mark_who_backed_the_winner():
    vote = PE.pathway_vote(
        [_reference("a", 0.9, "Lipid", "Sterol"), _reference("b", 0.85, "Amino Acid")]
    )
    rows = PE.vote_rows(vote)
    assert [row["rank"] for row in rows] == [1, 2]
    assert [row["backed_winning_super"] for row in rows] == [1, 0]


# --- equivalence with the published mechanism --------------------------------


def _audit_vote(rows, k, floor, power):
    """The audit's inner loop, transcribed.

    Copied from ``leave_one_out_validation`` in tools/audit_pathway_predictions.py
    so a divergence in floor, k, weighting or tie-breaking fails here rather than
    silently making the report's reconstruction something other than the
    mechanism the manuscript describes.
    """
    sims = [row["similarity_to_target"] for row in rows]
    supers = [row["metabolon_super"] for row in rows]
    subs = [row["metabolon_sub"] for row in rows]
    cand = sorted(
        ((sims[j], j) for j in range(len(rows)) if sims[j] >= floor and supers[j]), reverse=True
    )[:k]
    if not cand:
        return "", ""
    votes: dict[str, float] = defaultdict(float)
    for sim, j in cand:
        votes[supers[j]] += sim**power
    pred_super = max(votes, key=votes.get)
    sub_votes: dict[str, float] = defaultdict(float)
    for sim, j in cand:
        if supers[j] == pred_super and subs[j]:
            sub_votes[subs[j]] += sim**power
    pred_sub = max(sub_votes, key=sub_votes.get) if sub_votes else ""
    return pred_super, pred_sub


def test_reconstruction_matches_the_published_mechanism_on_random_electorates():
    """Randomized equivalence against a transcription of the audit's own loop.

    The committed leave-one-out accuracy is the published measurement of how well
    this vote performs, and the report borrows it to describe itself. That is only
    legitimate if the report runs the same vote.
    """

    rng = random.Random(20260805)
    pathways = ["Lipid", "Amino Acid", "Xenobiotics", "Nucleotide", "Cofactor"]
    subs = ["Sterol", "Fatty Acid", "Glycine", "Purine", ""]
    for _ in range(500):
        rows = [
            _reference(
                f"c{index}",
                round(rng.uniform(0.05, 0.99), 4),
                rng.choice(pathways) if rng.random() > 0.1 else "",
                rng.choice(subs),
            )
            for index in range(rng.randint(0, 25))
        ]
        # Distinct similarities keep the comparison about the algorithm rather
        # than about two equally valid tie-breaks; ties are covered separately.
        seen: set[float] = set()
        rows = [r for r in rows if not (r["similarity_to_target"] in seen or seen.add(r["similarity_to_target"]))]
        mine = PE.pathway_vote(rows)
        theirs = _audit_vote(rows, PE.VOTE_K, PE.VOTE_SIMILARITY_FLOOR, PE.VOTE_WEIGHT_POWER)
        assert (mine.winning_super, mine.winning_sub) == theirs


def test_the_vote_parameters_are_the_published_ones():
    assert (PE.VOTE_K, PE.VOTE_SIMILARITY_FLOOR, PE.VOTE_WEIGHT_POWER) == (10, 0.30, 2.0)


@pytest.mark.slow
def test_reconstruction_reproduces_the_committed_leave_one_out_accuracy():
    """End-to-end fidelity against the committed cohort result.

    Needs the reference database, which is fetched rather than committed, so this
    runs where the bundle exists and skips elsewhere.
    """

    database = ROOT / "vendor" / "metabolite_synonyms.sqlite"
    if not database.exists():
        pytest.skip("reference metabolite database is not present in this checkout")
    pytest.importorskip("rdkit")

    import sqlite3

    import pandas as pd
    from rdkit import Chem, RDLogger, rdBase  # noqa: F401
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.DataStructs import BulkTanimotoSimilarity

    RDLogger.DisableLog("rdApp.*")
    with sqlite3.connect(str(database)) as connection:
        reference = pd.read_sql_query(
            "SELECT smiles, metabolon_super, metabolon_sub FROM metabolites "
            "WHERE TRIM(COALESCE(smiles,'')) <> '' AND TRIM(COALESCE(metabolon_super,'')) <> ''",
            connection,
        )
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprints, supers, subs = [], [], []
    for _, row in reference.iterrows():
        mol = Chem.MolFromSmiles(str(row["smiles"]).split("|")[0].strip())
        if mol is None:
            continue
        fingerprints.append(generator.GetFingerprint(mol))
        supers.append(str(row["metabolon_super"]).strip())
        subs.append(str(row["metabolon_sub"] or "").strip())

    correct = evaluated = 0
    for index in range(len(fingerprints)):
        similarities = BulkTanimotoSimilarity(fingerprints[index], fingerprints)
        rows = [
            _reference(f"r{other}", similarities[other], supers[other], subs[other])
            for other in range(len(fingerprints))
            if other != index
        ]
        vote = PE.pathway_vote(rows)
        if not vote.available:
            continue
        evaluated += 1
        correct += int(vote.winning_super.lower() == supers[index].lower())

    committed = json.loads((ROOT / "docs/paper/data/pathway/audit_pathway_loo.json").read_text())
    assert abs(correct / evaluated - committed["super_accuracy"]) < 0.02


# --- reliability lookup ------------------------------------------------------


def test_similarity_bins_match_the_audit_edges():
    assert PE.similarity_bin(0.1) == "[0.0,0.3)"
    assert PE.similarity_bin(0.3) == "[0.3,0.5)"
    assert PE.similarity_bin(0.62) == "[0.5,0.7)"
    assert PE.similarity_bin(1.0) == "[0.7,1.0)"
    assert PE.similarity_bin(None) == ""


def test_reliability_is_quoted_from_the_committed_audit():
    measured = PE.measured_reliability(best_similarity=0.62)
    committed = json.loads(
        (ROOT / "docs/paper/data/pathway/audit_pathway_summary.json").read_text()
    )["accuracy"]["by_similarity_bin"]["[0.5,0.7)"]
    assert measured["super_accuracy"] == committed["super"]["accuracy"]
    assert measured["sub_accuracy"] == committed["sub"]["accuracy"]
    assert measured["n"] == committed["super"]["n"]


def test_reliability_falls_back_to_provenance_without_a_similarity():
    measured = PE.measured_reliability(best_similarity=None, source="exact")
    assert "direct catalogue lookup" in measured["basis"]
    assert PE.measured_reliability(best_similarity=None, source="predicted")["basis"].endswith(
        "neighbor consensus"
    )


def test_reliability_is_omitted_rather_than_invented_when_the_audit_is_missing():
    assert PE.measured_reliability(best_similarity=0.62, audit={}) is None


def test_annotation_ceiling_reads_the_harmonization_artifact():
    ceiling = PE.annotation_ceiling()
    committed = json.loads((ROOT / "docs/paper/data/taxonomy_harmonization.json").read_text())
    assert ceiling["super_agreement"] == committed["super_pathway"]["raw"]["reference_ceiling"]["accuracy"]
    assert ceiling["sub_agreement"] == committed["subpathway"]["raw"]["reference_ceiling"]["accuracy"]
    assert PE.annotation_ceiling(Path("/nonexistent/harmonization.json")) is None


# --- the section -------------------------------------------------------------


def _pathway_info(**over):
    base = dict(
        super="Lipid", sub="Sterol", source="predicted", tier="medium",
        similarity=0.62, status="strong_similarity_vote", basis=PE.BASIS_NEIGHBOR_VOTE,
    )
    base.update(over)
    return base


def _vote():
    return PE.pathway_vote(
        [
            _reference("cholate", 0.82, "Lipid", "Sterol"),
            _reference("glycocholate", 0.71, "Lipid", "Sterol"),
            _reference("taurine", 0.44, "Amino Acid", "Methionine"),
        ]
    )


def test_section_shows_the_voters_and_the_tally():
    body = PE.pathway_section_html(_pathway_info(), _vote(), public=False)
    assert "cholate" in body
    assert "glycocholate" in body
    assert "Vote weight" in body
    assert "Subpathway" in body
    # The losing pathway's subpathway must not appear in the subpathway tally.
    assert "Methionine" not in body.split("Subpathway is decided only")[-1]


def test_section_trims_the_voter_list_for_public_reports():
    rows = [_reference(f"voter{index}", 0.9 - index * 0.02, "Lipid", "Sterol") for index in range(10)]
    vote = PE.pathway_vote(rows)
    public = PE.pathway_section_html(_pathway_info(), vote, public=True)
    internal = PE.pathway_section_html(_pathway_info(), vote, public=False)
    assert public.count("<tr>") < internal.count("<tr>")
    assert "further voter" in public
    assert "voter9" in internal


def test_section_names_the_lookup_route_rather_than_calling_it_all_exact():
    structural = PE.pathway_section_html(
        _pathway_info(source="exact", basis=PE.BASIS_INCHIKEY, similarity=1.0),
        PE.PathwayVote(),
        public=True,
    )
    by_name = PE.pathway_section_html(
        _pathway_info(source="exact", basis=PE.BASIS_CANONICAL_NAME, similarity=1.0),
        PE.PathwayVote(),
        public=True,
    )
    assert "InChIKey skeleton matches" in structural
    assert "name match is weaker than a structure match" in by_name
    assert PE.BASIS_IS_STRUCTURAL[PE.BASIS_INCHIKEY] is True
    assert PE.BASIS_IS_STRUCTURAL[PE.BASIS_CANONICAL_NAME] is False


def test_section_warns_when_the_pathway_rests_on_a_rejected_structure():
    withheld = PE.pathway_section_html(_pathway_info(), _vote(), public=True, msi_level=4)
    shown = PE.pathway_section_html(_pathway_info(), _vote(), public=True, msi_level=3)
    assert "does not stand behind" in withheld
    assert "does not stand behind" not in shown


def test_section_discloses_a_reconstruction_that_disagrees():
    agreeing = PE.pathway_section_html(_pathway_info(super="Lipid"), _vote(), public=True)
    disagreeing = PE.pathway_section_html(_pathway_info(super="Nucleotide"), _vote(), public=True)
    assert "reproduces the assigned pathway" in agreeing
    assert "does not reproduce the assigned pathway" in disagreeing
    assert "Nucleotide" in disagreeing


def test_section_survives_having_no_voters_at_all():
    body = PE.pathway_section_html(
        _pathway_info(similarity=0.0), PE.PathwayVote(), public=True
    )
    assert "No catalogued compound was similar enough to vote" in body


def test_section_quotes_measured_reliability_and_its_ceiling():
    body = PE.pathway_section_html(_pathway_info(), _vote(), public=True)
    assert "How often this kind of assignment has been right" in body
    assert "still agree only" in body


def test_section_carries_the_candidate_list_agreement_when_supplied():
    body = PE.pathway_section_html(
        _pathway_info(), _vote(), public=True, candidate_agreement="High Confidence: 8/10 map"
    )
    assert "8/10" in body
    assert "Across the candidate list" in body


# --- rendered page -----------------------------------------------------------

pytest.importorskip("pandas", reason="the report renderer needs pandas")
pytest.importorskip("rdkit", reason="the report renderer needs rdkit")

from structural_elucidation.reports import compound as R  # noqa: E402


def _render(variant: str, *, msi_level: int = 3, vote=None, pathway=None) -> str:
    conf = C.compute_confidence(
        dict(
            source="predicted", structure_similarity=0.60, mass_match_score=0.99,
            ms2_spectrum_present=False, ms2_score=None, pathway_vote_share=0.65,
            observed_mz_available=True,
        )
        if msi_level == 3
        else dict(
            source="predicted", structure_similarity=0.10, mass_match_score=0.0,
            abs_mass_error_da=176.0, observed_mz_available=True, pathway_vote_share=0.65,
        )
    )
    assert conf.msi_level == msi_level, f"fixture drifted: got L{conf.msi_level}"
    return R.render_html(
        out_dir=Path("/tmp/does-not-need-to-exist"),
        target={
            "control_id": "X-99999", "name": "Bile-acid-like conjugate",
            "smiles": "OC(=O)CCC(N)C(=O)O", "formula": "C5H9NO4",
            "exact_mass": 147.053158, "inchikey": "WHUUTDBJXJRKMK-UHFFFAOYSA-N",
            "tentative_descriptor": {"label": "Bile-acid-like conjugate", "confidence": "medium"},
        },
        evidence={"parent_mz": 146.0459, "mass_error_da": 0.0004},
        db_detail={},
        confidence=conf,
        conf_fields={"ms2_detail": {"present": False, "n_pairs": 0}},
        pathway_info=pathway if pathway is not None else _pathway_info(),
        synonyms=[],
        neighbors=[],
        nodes=[],
        edges=[],
        candidate_scores=[
            {"rank": 1, "predicted_formula": "C5H9NO4", "predicted_smiles": "OC(=O)CCC(N)C(=O)O",
             "candidate_report_score": 0.04, "predicted_pathway": "Lipid"},
        ],
        fragment_candidate_rows=[],
        structure_panel_html="<svg id='structure'></svg>",
        consensus_svg="<svg id='consensus'></svg>",
        candidate_consensus_maps={1: "<svg/>"},
        candidate_consensus_summaries={},
        default_consensus_rank=1,
        structural_possibility_tree_svg="",
        structural_possibility_rows=[],
        feature_total_count=0,
        extended_network_svg_html="",
        extended_unknown_rows=[],
        extended_known_rows=[],
        extended_stats={},
        correlation_rows=[],
        correlation_blocker=None,
        variant=variant,
        pathway_vote=vote if vote is not None else _vote(),
    )


@pytest.mark.parametrize("variant", ["public", "internal"])
def test_rendered_page_carries_the_pathway_section(variant):
    page = _render(variant)
    assert R.PATHWAY_SECTION_TITLE in page
    assert "cholate" in page
    assert "Assigned pathway" in page
    assert "Assigned subpathway" in page


def test_rendered_page_warns_at_a_withheld_level():
    assert "does not stand behind" in _render("public", msi_level=4)
    assert "does not stand behind" not in _render("public", msi_level=3)


def test_rendered_page_survives_a_missing_vote():
    page = _render("public", vote=PE.PathwayVote())
    assert R.PATHWAY_SECTION_TITLE in page
    assert "No catalogued compound was similar enough to vote" in page


def test_the_section_heading_has_a_guide_entry():
    """The public info button links to guides.html by slug, and the guide list is
    built separately from the headings, so the two can silently drift apart."""

    page = _render("public")
    slug = R.slugify_label(R.PATHWAY_SECTION_TITLE)
    assert f"guides.html#guide-{slug}" in page
    guides = R.render_guides_html("X-99999", [(slug, R.PATHWAY_SECTION_TITLE)])
    assert f"id='guide-{slug}'" in guides
