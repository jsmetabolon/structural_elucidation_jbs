"""Vendored third-party packages for spec2mol v21 inference.

Pipeline context:
    Makes the MOSES SMILES decoder and NCE elucidation stack importable from
    ``vendor/`` without separate pip installs. First-party code prepends
    ``vendor/`` to ``sys.path`` and registers ``moses.*`` aliases before
    unpickling decoder checkpoints.

Subpackages:
    ``decoder`` — constrained de novo SMILES VAE (Track B generation).
    ``nce`` — spectrum verification and chemistry utilities (rerank path).
"""
