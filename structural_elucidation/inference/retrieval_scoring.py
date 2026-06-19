#!/usr/bin/env python3
"""Weighted retrieval reranking over latent, mass, fragment, and HTML signals.

Combines cosine similarity with precursor mass fit, fragment explanation,
HTML formula bonus, prevalence, and structural-neighbor scores into a
single ``retrieval_score`` used by Track A mass-gated KNN.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

# ── scoring component registry ────────────────────────────────────────────────

VALID_COMPONENTS: tuple[str, ...] = (
    "latent_cosine",
    "precursor_mass_score",
    "fragment_explained_intensity_fraction",
    "fragment_matched_peak_fraction",
    "fragment_match_count_score",
    "html_formula_bonus",
    "isotope_formula_score",
    "substructure_consensus_score",
    "tanimoto_neighbor_score",
    "prevalence_score",
)

DEFAULT_RETRIEVAL_WEIGHTS: dict[str, float] = {
    "latent_cosine": 0.30,
    "precursor_mass_score": 0.15,
    "fragment_explained_intensity_fraction": 0.18,
    "fragment_matched_peak_fraction": 0.08,
    "fragment_match_count_score": 0.08,
    "html_formula_bonus": 0.04,
    "isotope_formula_score": 0.02,
    "substructure_consensus_score": 0.07,
    "tanimoto_neighbor_score": 0.05,
    "prevalence_score": 0.03,
}


@dataclass(frozen=True)
class RetrievalScoringConfig:
    """Immutable reranking weights and latent-cosine transform mode."""

    weights: dict[str, float]
    latent_cos_transform: str = "raw"

    def to_dict(self) -> dict[str, object]:
        """Serialize config for run summaries and reproducibility logs."""
        return {
            "weights": {k: float(v) for k, v in self.weights.items()},
            "latent_cos_transform": str(self.latent_cos_transform),
        }


def normalize_weights(weights: dict[str, float]) -> dict[str, float]:
    """Restrict to ``VALID_COMPONENTS`` and L1-normalize to sum 1.

    Raises:
        ValueError: When all component weights are zero or negative total.
    """
    cleaned = {name: float(weights.get(name, 0.0)) for name in VALID_COMPONENTS}
    total = float(sum(cleaned.values()))
    if total <= 0:
        raise ValueError("Retrieval weights must sum to a positive value.")
    return {name: float(value / total) for name, value in cleaned.items()}


def resolve_scoring_config(weights_json_path: str | Path | None) -> RetrievalScoringConfig:
    """Load reranking weights from JSON or fall back to ``DEFAULT_RETRIEVAL_WEIGHTS``.

    Args:
        weights_json_path: Optional JSON with ``weights`` and
            ``latent_cos_transform`` (``raw`` or ``scaled01``).

    Returns:
        Normalized ``RetrievalScoringConfig``.

    Raises:
        FileNotFoundError: When a non-empty path does not exist.
        ValueError: On malformed JSON or invalid transform mode.
    """
    if weights_json_path is None or not str(weights_json_path).strip():
        return RetrievalScoringConfig(weights=normalize_weights(DEFAULT_RETRIEVAL_WEIGHTS), latent_cos_transform="raw")

    path = Path(weights_json_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Retrieval weights JSON not found: {path}")

    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"Retrieval weights JSON must be an object: {path}")

    raw_weights = payload.get("weights", payload)
    if not isinstance(raw_weights, dict):
        raise ValueError(f"Retrieval weights JSON 'weights' must be an object: {path}")

    latent_cos_transform = str(payload.get("latent_cos_transform", "raw")).strip().lower()
    if latent_cos_transform not in {"raw", "scaled01"}:
        raise ValueError(f"Invalid latent_cos_transform={latent_cos_transform!r}. Expected 'raw' or 'scaled01'.")

    normalized = normalize_weights({k: float(raw_weights.get(k, 0.0)) for k in VALID_COMPONENTS})
    return RetrievalScoringConfig(weights=normalized, latent_cos_transform=latent_cos_transform)


def latent_cos_feature_value(latent_cosine: float, mode: str) -> float:
    """Map raw cosine [-1, 1] to feature space (optionally scaled to [0, 1])."""
    if mode == "scaled01":
        return float((float(latent_cosine) + 1.0) / 2.0)
    return float(latent_cosine)


def compute_retrieval_score(
    latent_cosine: float,
    precursor_mass_score: float,
    fragment_explained_intensity_fraction: float,
    fragment_matched_peak_fraction: float,
    html_formula_bonus: float = 0.0,
    isotope_formula_score: float = 0.0,
    fragment_match_count_score: float = 0.0,
    substructure_consensus_score: float = 0.0,
    tanimoto_neighbor_score: float = 0.0,
    prevalence_score: float = 0.0,
    config: RetrievalScoringConfig | None = None,
) -> float:
    """Weighted linear combination of retrieval feature components.

    Args:
        latent_cosine: Cosine similarity between query and bank latent.
        precursor_mass_score: Mass-fit score within the active gate window.
        fragment_explained_intensity_fraction: Fraction of query intensity
            explained by candidate fragment spectrum.
        fragment_matched_peak_fraction: Fraction of query peaks matched.
        html_formula_bonus: 1.0 when candidate formula is in HTML usable set.
        isotope_formula_score: Reserved isotope-formula component (often 0).
        fragment_match_count_score: Normalized fragment match count.
        substructure_consensus_score: MCS-based neighbor consensus.
        tanimoto_neighbor_score: Fingerprint similarity among top analogs.
        prevalence_score: FASSTrecords prevalence for candidate InChIKey14.
        config: Scoring weights; defaults to built-in profile when ``None``.

    Returns:
        Scalar retrieval score in approximately [0, 1] when features are
        themselves bounded.
    """
    scoring_cfg = config or resolve_scoring_config(None)
    features = {
        "latent_cosine": latent_cos_feature_value(latent_cosine, scoring_cfg.latent_cos_transform),
        "precursor_mass_score": float(precursor_mass_score),
        "fragment_explained_intensity_fraction": float(fragment_explained_intensity_fraction),
        "fragment_matched_peak_fraction": float(fragment_matched_peak_fraction),
        "fragment_match_count_score": float(fragment_match_count_score),
        "html_formula_bonus": float(html_formula_bonus),
        "isotope_formula_score": float(isotope_formula_score),
        "substructure_consensus_score": float(substructure_consensus_score),
        "tanimoto_neighbor_score": float(tanimoto_neighbor_score),
        "prevalence_score": float(prevalence_score),
    }
    return float(sum(float(scoring_cfg.weights[name]) * float(features[name]) for name in VALID_COMPONENTS))
