"""Post-rerank enrichment utilities for equalization deliverables.

Adds interpretable annotation columns to ``unified_rerank_all_top10.csv``
artifacts produced by the equalization pipeline, one library at a time.

Pipeline context:
    Run after ``rerank`` inside ``models/pipeline_equalization/``. Stages are
    independent CLIs: pathway → lipid_class → peptide_composition.

Key I/O:
    Inputs: per-lib rerank top-10 CSVs under ``rerank_lib{lib_id}/``.
    Outputs: suffixed CSVs (``_with_pathways``, ``_lipid_class``,
    ``_peptide_composition``) plus per-stage summary JSON.
"""
