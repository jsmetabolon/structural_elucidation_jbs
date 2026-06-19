"""Core structural-elucidation pipeline stages.

Stages (in typical execution order):
    1. ``candidates`` — library-anchored candidate generation for missed S0.
    2. ``verify`` — NCE predicted-spectrum verifier scoring across sources.
    3. ``pool_v1`` / ``pool`` — deduplicated cross-source candidate pool.
    4. ``rerank`` — gradient-boosted unified reranker with LOOCV gating.
    5. ``enrichment`` — pathway, lipid-class, and peptide-composition columns.

Key I/O:
    Each stage is a CLI module writing under ``models/<stage>_v*/`` with
    summary JSON and markdown reports for acceptance checks.
"""
