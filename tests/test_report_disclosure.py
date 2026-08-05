"""What a report is allowed to disclose at each MSI level.

These render the real page rather than checking the policy helpers in isolation,
because the gap this closes was exactly that: the policy was stated in the docs
and the renderer ignored it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "structural_elucidation" / "reports"))

import confidence as C  # noqa: E402

pytest.importorskip("pandas", reason="the report renderer needs pandas")
pytest.importorskip("rdkit", reason="the report renderer needs rdkit")

from structural_elucidation.reports import compound as R  # noqa: E402

SMILES = "OC(=O)CCC(N)C(=O)O"
INCHIKEY = "WHUUTDBJXJRKMK-UHFFFAOYSA-N"


def _confidence(msi_level: int, **over) -> C.Confidence:
    """Build a Confidence whose MSI level is the one under test."""

    profiles = {
        1: dict(metabolon_level="1", truth_match=True, structure_similarity=0.95,
                mass_match_score=0.99, ms2_spectrum_present=True, ms2_score=0.8),
        2: dict(structure_similarity=0.95, mass_match_score=0.99,
                ms2_spectrum_present=True, ms2_score=0.8),
        3: dict(structure_similarity=0.60, mass_match_score=0.99,
                ms2_spectrum_present=False, ms2_score=None),
        4: dict(structure_similarity=0.10, mass_match_score=0.0,
                abs_mass_error_da=176.0, observed_mz_available=True),
        5: dict(structure_similarity=0.10, mass_match_score=None,
                observed_mz_available=False),
    }
    fields = dict(source="predicted", pathway_vote_share=0.8, observed_mz_available=True)
    fields.update(profiles[msi_level])
    fields.update(over)
    conf = C.compute_confidence(fields)
    assert conf.msi_level == msi_level, f"fixture drifted: got L{conf.msi_level}"
    return conf


def _render(conf: C.Confidence, *, variant: str, candidates: int = 10, **over) -> str:
    kwargs = dict(
        out_dir=Path("/tmp/does-not-need-to-exist"),
        target={
            "control_id": "X-99999",
            "name": "Glutamate-like conjugate",
            "smiles": SMILES,
            "formula": "C5H9NO4",
            "exact_mass": 147.053158,
            "inchikey": INCHIKEY,
            "tentative_descriptor": {"label": "Glutamate-like conjugate", "confidence": "medium"},
        },
        evidence={"parent_mz": 146.0459, "mass_error_da": 0.0004},
        db_detail={},
        confidence=conf,
        conf_fields={"ms2_detail": {"present": True, "n_pairs": 12}},
        pathway_info={"source": "predicted", "similarity": 0.6, "super": "Amino Acid"},
        synonyms=[{"synonym": "L-glutamic acid"}],
        neighbors=[],
        nodes=[],
        edges=[],
        candidate_scores=[
            {
                "rank": rank,
                "predicted_formula": "C5H9NO4",
                "predicted_smiles": SMILES,
                "candidate_report_score": 0.04 / rank,
            }
            for rank in range(1, candidates + 1)
        ],
        fragment_candidate_rows=[
            {
                "rank": rank,
                "predicted_formula": "C5H9NO4",
                "predicted_smiles": SMILES,
                "candidate_report_score": 0.04 / rank,
                "forward_spectrum_score": 0.5,
                "predicted_spectrum_score_v2": 0.6,
            }
            for rank in range(1, candidates + 1)
        ],
        structure_panel_html="<svg id='the-structure-drawing'></svg>",
        consensus_svg="<svg id='the-consensus-map'></svg>",
        candidate_consensus_maps={rank: "<svg/>" for rank in range(1, candidates + 1)},
        candidate_consensus_summaries={},
        default_consensus_rank=1,
        structural_possibility_tree_svg="",
        structural_possibility_rows=[],
        feature_total_count=12,
        extended_network_svg_html="",
        extended_unknown_rows=[],
        extended_known_rows=[],
        extended_stats={},
        correlation_rows=[],
        correlation_blocker=None,
        variant=variant,
    )
    kwargs.update(over)
    return R.render_html(**kwargs)


@pytest.mark.parametrize("level", [4, 5])
@pytest.mark.parametrize("variant", ["public", "internal"])
def test_levels_four_and_five_withhold_the_structure(level, variant):
    page = _render(_confidence(level), variant=variant)
    assert "the-structure-drawing" not in page
    assert "the-consensus-map" not in page
    assert "Structure withheld" in page
    assert '<option value="' not in page
    # The evidence that does survive is the mass, and the formula at Level 4.
    assert "146.0459" in page
    if level == 4:
        assert "C5H9NO4" in page


@pytest.mark.parametrize("level", [4, 5])
def test_a_withheld_public_report_names_no_structure_anywhere(level):
    """The per-candidate fit table used to smuggle the structure back in.

    It lists formulas and SMILES for the top candidates, so withholding the
    structure panel alone left the structure on the page.
    """

    page = _render(_confidence(level), variant="public")
    assert SMILES not in page
    assert INCHIKEY not in page
    assert "L-glutamic acid" not in page
    assert "no per-candidate spectral fits" in page


@pytest.mark.parametrize("level", [1, 2, 3])
@pytest.mark.parametrize("variant", ["public", "internal"])
def test_levels_one_to_three_still_show_the_structure(level, variant):
    page = _render(_confidence(level), variant=variant)
    assert "the-structure-drawing" in page
    assert SMILES in page
    assert "Structure withheld" not in page


def test_public_shows_three_candidates_at_level_three_and_internal_shows_ten():
    conf = _confidence(3)
    public = _render(conf, variant="public")
    internal = _render(conf, variant="internal")
    assert "Top 3 potential structures" in public
    assert public.count('<option value="') == 3
    assert internal.count('<option value="') == 10


def test_public_asserts_a_single_structure_at_level_two():
    conf = _confidence(2)
    public = _render(conf, variant="public")
    assert "Proposed structure" in public
    assert public.count('<option value="') == 1
    # Internal keeps the alternates for triage.
    assert _render(conf, variant="internal").count('<option value="') == 10


def test_public_candidate_count_is_configurable():
    conf = _confidence(3)
    page = _render(conf, variant="public", public_candidate_n=5)
    assert page.count('<option value="') == 5
    assert "Top 5 potential structures" in page


def test_a_refuted_structure_is_withheld_and_badged():
    conf = C.compute_confidence(
        dict(
            source="predicted", structure_similarity=0.95, mass_match_score=0.99,
            ms2_spectrum_present=True, ms2_score=0.8, observed_mz_available=True,
            pathway_vote_share=0.8, assurance=C.ASSURANCE_STANDARD_REFUTED,
        )
    )
    page = _render(conf, variant="public")
    assert "Ruled out by an authentic-standard run" in page
    assert "the-structure-drawing" not in page
    assert SMILES not in page
