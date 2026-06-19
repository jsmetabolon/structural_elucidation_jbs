"""Tests for peptide composition normalization and row-level enrichment."""

from __future__ import annotations

import pandas as pd

from structural_elucidation.pipeline.enrichment.peptide_composition import (
    annotate_row,
    enrich_frame,
    normalize_peptide_composition,
    peptide_composition_label,
)


def test_normalize_peptide_composition_collapses_xle() -> None:
    """Leu/Ile ambiguity collapses to ``Xle`` in canonical composition strings."""
    assert normalize_peptide_composition("Val-Leu-Gly") == "Gly-Val-Xle"
    assert normalize_peptide_composition("isoleucine, alanine") == "Ala-Xle"


def test_peptide_composition_label_generation() -> None:
    """Human-readable labels reflect peptide length and residue order."""
    assert peptide_composition_label("Ala-Val") == "dipeptide: Ala, Val"
    assert peptide_composition_label("Gly-Val-Xle") == "tripeptide: Gly, Val, Xle"
    assert peptide_composition_label("Ala-Ala") == "dipeptide: Ala, Ala"


def test_high_confidence_ordered_tripeptide() -> None:
    """Strong spectrum and ordering evidence yields high composition and order confidence."""
    row = pd.Series(
        {
            "normalized_peptide_composition": "Gly-Val-Xle",
            "candidate_ordered_sequence": "Gly-Val-Leu",
            "is_peptide_candidate": True,
            "composition_confidence": 0.9,
            "composition_support_score": 0.8,
            "immonium_evidence": True,
            "peptide_neutral_mass_fit": True,
            "ordered_sequence_within_composition_rank": 1,
            "spectrum_b_ion_score": 0.8,
            "spectrum_y_ion_score": 0.7,
            "peptide_connectivity_confidence": 0.9,
            "ordered_sequence_support_score": 0.8,
        }
    )
    out = annotate_row(row)
    assert out["peptide_composition_label"] == "tripeptide: Gly, Val, Xle"
    assert out["peptide_composition_confidence"] > 0.85
    assert out["peptide_ordering_confidence"] > 0.85
    assert out["peptide_ordering_label"] == "Gly-Val-Xle"
    assert out["peptide_length"] == 3


def test_low_ordering_confidence_keeps_composition_deliverable() -> None:
    """Weak ordering evidence must not block composition-only deliverables."""
    row = pd.Series(
        {
            "normalized_peptide_composition": "Ala-Val",
            "is_peptide_candidate": True,
            "composition_support_score": 0.7,
            "mass_match_score": 0.9,
            "ordered_sequence_within_composition_rank": 4,
            "spectrum_b_ion_score": 0.0,
            "spectrum_y_ion_score": 0.0,
            "ordered_sequence_support_score": 0.0,
        }
    )
    out = annotate_row(row)
    assert out["peptide_composition_label"] == "dipeptide: Ala, Val"
    assert out["peptide_ordering_label"] == "order uncertain"
    assert out["peptide_ordering_confidence"] < 0.65


def test_identical_residue_dipeptide_has_unambiguous_order() -> None:
    """Homogeneous dipeptides have trivially unique sequence order."""
    row = pd.Series(
        {
            "normalized_peptide_composition": "Ala-Ala",
            "is_peptide_candidate": True,
            "composition_support_score": 0.5,
        }
    )
    out = annotate_row(row)
    assert out["peptide_ordering_confidence"] == 1.0
    assert out["peptide_ordering_label"] == "Ala-Ala"


def test_non_peptide_rows_get_empty_columns() -> None:
    """Non-peptide candidates receive blank peptide enrichment columns."""
    enriched = enrich_frame(pd.DataFrame([{"predicted_smiles": "CCO", "rank": 1}]))
    assert enriched.loc[0, "peptide_composition_label"] == ""
    assert enriched.loc[0, "peptide_ordering_label"] == ""
    assert enriched.loc[0, "peptide_length"] == ""
