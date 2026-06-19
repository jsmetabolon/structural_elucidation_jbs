"""Structural elucidation inference and reranking pipeline.

Mass-spectrometry structural elucidation for Metabolon X controls: candidate
generation from library anchors and truth-blind enumerators, NCE verifier
scoring, cross-source pool merge, gradient-boosted unified rerank, and
post-rerank pathway/lipid/peptide enrichment.

Pipeline context:
    Subpackages under ``structural_elucidation`` feed the equalization
    deliverable tree (``models/pipeline_equalization/``). Typical flow:
    candidates → verify → pool → rerank → enrichment.

Key I/O:
    Inputs: feature rows, neutral-mass hypotheses, diagnosis CSVs, vendor
    libraries (HMDB, NCE snapshot), observed MS2 lookups.
    Outputs: ``candidate_metrics_all.csv``, ``combined_scored_candidates.csv``,
    ``unified_candidates_per_control.csv``, ``unified_rerank_top10.csv``,
    enrichment-augmented top-10 CSVs per library.
"""

__name__ = "structural_elucidation"
__description__ = "Mass-spectrometry structural elucidation inference and reranking"
