#!/usr/bin/env python3
# flake8: noqa
"""Static MS2/MS1 diagnostic tables for truth-blind chemistry classes.

evidence_confidence DOWNSTREAM CONTRACT:
- chemistry_class_evidence_per_control.csv writes both:
  - confidence: str tier in {low, medium, high, very_high}
  - evidence_score: float in [0.0, 1.0]
- Truth-blind enumerators join evidence_score into candidate rows as
  evidence_confidence (float).
- F2 promotion gate build_proposed_unknowns_v1._confidence_high accepts both.
- Numeric-to-tier mapping: <0.2 low, 0.2-0.4 medium, 0.4-0.6 high,
  >=0.6 very_high.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any


DIAGNOSTIC_FRAGMENTS_NEG: dict[str, dict[str, Any]] = {
    "sulfate": {
        "diagnostic_mz": [(79.9574, 0.005), (96.9601, 0.005)],
        "diagnostic_neutral_losses": [(79.9568, 0.010), (97.9674, 0.010)],
        "min_diagnostic_intensity_rel": 0.01,
        "min_evidence_count": 1,
        "deconj_mass_delta": 79.9568,
        "confounders": ["sulfonate(no_NL_80)", "phosphate(check_P_in_MS1)"],
        "enumerator_name": "sulfate_position_v1",
    },
    "glucuronide": {
        "diagnostic_mz": [(175.0241, 0.005), (113.0244, 0.005), (85.0295, 0.005)],
        "diagnostic_neutral_losses": [(176.0321, 0.010)],
        "min_diagnostic_intensity_rel": 0.01,
        "min_evidence_count": 1,
        "deconj_mass_delta": 176.0321,
        "enumerator_name": "glucuronide_v1",
    },
    "taurine": {
        "diagnostic_mz": [(124.0074, 0.005), (79.9574, 0.005), (80.0, 0.010)],
        "diagnostic_neutral_losses": [(125.0152, 0.010)],
        "min_diagnostic_intensity_rel": 0.05,
        "min_evidence_count": 1,
        "deconj_mass_delta": 107.0041,
        "enumerator_name": "taurine_variant_v1",
    },
    "glycine": {
        "diagnostic_mz": [(74.0247, 0.005)],
        "diagnostic_neutral_losses": [(75.0320, 0.010)],
        "min_diagnostic_intensity_rel": 0.02,
        "min_evidence_count": 1,
        "deconj_mass_delta": 57.0215,
        "enumerator_name": "glycine_v1",
    },
    "glutathione": {
        "diagnostic_mz": [(272.0888, 0.010), (254.0782, 0.010), (143.0464, 0.010), (128.0353, 0.010)],
        "diagnostic_neutral_losses": [(305.0681, 0.010), (129.0426, 0.010)],
        "min_diagnostic_intensity_rel": 0.01,
        "min_evidence_count": 2,
        "deconj_mass_delta": 305.0681,
        "enumerator_name": "glutathione_v1",
    },
    "cysteinyl": {
        "diagnostic_mz": [(121.0197, 0.005), (88.0, 0.010)],
        "diagnostic_neutral_losses": [(122.0276, 0.010)],
        "min_diagnostic_intensity_rel": 0.01,
        "min_evidence_count": 1,
        "deconj_mass_delta": 121.0197,
        "enumerator_name": "cysteinyl_v1",
    },
    "n_acetyl": {
        "diagnostic_mz": [],
        "diagnostic_neutral_losses": [(42.0106, 0.010), (59.0371, 0.010)],
        "min_diagnostic_intensity_rel": 0.05,
        "min_evidence_count": 1,
        "deconj_mass_delta": 42.0106,
        "enumerator_name": "n_acetyl_v1",
    },
    "n_lactoyl": {
        "diagnostic_mz": [(89.0244, 0.005), (71.0139, 0.005)],
        "diagnostic_neutral_losses": [(90.0317, 0.010), (72.0211, 0.010)],
        "min_diagnostic_intensity_rel": 0.02,
        "min_evidence_count": 1,
        "deconj_mass_delta": 72.0211,
        "enumerator_name": "lactoyl_amino_acid_v1",
    },
    "methylation": {
        "diagnostic_mz": [],
        "diagnostic_neutral_losses": [],
        "ms1_parent_delta_da": 14.0157,
        "deconj_mass_delta": 14.0157,
        "min_diagnostic_intensity_rel": 0.0,
        "min_evidence_count": 1,
        "enumerator_name": "methylation_v1",
    },
}


DIAGNOSTIC_FRAGMENTS_POS: dict[str, dict[str, Any]] = deepcopy(DIAGNOSTIC_FRAGMENTS_NEG)
DIAGNOSTIC_FRAGMENTS_POS.update(
    {
        "sulfate": {
            **DIAGNOSTIC_FRAGMENTS_NEG["sulfate"],
            "diagnostic_mz": [(80.9652, 0.010), (98.9747, 0.010)],
        },
        "taurine": {
            **DIAGNOSTIC_FRAGMENTS_NEG["taurine"],
            "diagnostic_mz": [(126.0220, 0.010), (108.0114, 0.010), (80.9652, 0.010)],
        },
        "glycine": {
            **DIAGNOSTIC_FRAGMENTS_NEG["glycine"],
            "diagnostic_mz": [(76.0393, 0.010), (30.0344, 0.010)],
        },
        "glutathione": {
            **DIAGNOSTIC_FRAGMENTS_NEG["glutathione"],
            "diagnostic_mz": [(308.0911, 0.010), (179.0485, 0.010), (162.0219, 0.010), (130.0499, 0.010)],
        },
        "cysteinyl": {
            **DIAGNOSTIC_FRAGMENTS_NEG["cysteinyl"],
            "diagnostic_mz": [(122.0270, 0.010), (76.0221, 0.010)],
        },
        "n_lactoyl": {
            **DIAGNOSTIC_FRAGMENTS_NEG["n_lactoyl"],
            "diagnostic_mz": [(90.0550, 0.010), (72.0444, 0.010)],
        },
    }
)


CONFIDENCE_ORDER = {"low": 0, "medium": 1, "high": 2, "very_high": 3}


def diagnostic_table(ion_mode: str) -> dict[str, dict[str, Any]]:
    return DIAGNOSTIC_FRAGMENTS_POS if str(ion_mode).lower().startswith("pos") else DIAGNOSTIC_FRAGMENTS_NEG


def confidence_tier(evidence_score: float, n_evidence: int, min_evidence_count: int) -> str:
    if float(evidence_score) >= 0.6 and int(n_evidence) >= int(min_evidence_count) + 1:
        return "very_high"
    if float(evidence_score) >= 0.4 and int(n_evidence) >= int(min_evidence_count):
        return "high"
    if float(evidence_score) >= 0.2 and int(n_evidence) >= int(min_evidence_count):
        return "medium"
    return "low"


def confidence_at_least(value: Any, minimum: str = "medium") -> bool:
    return CONFIDENCE_ORDER.get(str(value or "").lower(), 0) >= CONFIDENCE_ORDER.get(str(minimum).lower(), 1)


TRUTH_BLIND_COMBO_SOURCES = {
    "sulfate_position_v1",
    "taurine_variant_v1",
    "lactoyl_amino_acid_v1",
    "glucuronide_v1",
    "glycine_v1",
    "glutathione_v1",
    "cysteinyl_v1",
    "n_acetyl_v1",
    "methylation_v1",
    "hydroxyl_v1",
    "dehydro_v1",
    "demethyl_v1",
    "wgcna_seeded_v1",
}
