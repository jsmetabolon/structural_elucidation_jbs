"""Adapters that bridge external pipelines into structural elucidation inputs.

The :mod:`structural_elucidation.adapters.untargeted_inputs` module converts a
neutral, dependency-light handoff produced by the untargeted metabolomics
pipeline (``untargeted_workflow``) into the per-library input tables that the
full elucidation stack (candidates -> verify -> pool -> rerank -> reports)
consumes for brand-new, un-curated features.
"""
