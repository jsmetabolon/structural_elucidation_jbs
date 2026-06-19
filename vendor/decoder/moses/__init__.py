"""MOSES molecular-sequence modeling utilities (vendored).

Pipeline context:
    Upstream MOSES decoder code, aliased as ``moses`` during checkpoint
    unpickling. Inference touches ``models_storage`` and ``trans.model``;
    training helpers remain under ``trans.trainer``.

See:
    ``structural_elucidation.inference.predict.register_moses_aliases``
"""

__version__ = '0.2'
