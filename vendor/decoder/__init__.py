"""MOSES-based SMILES decoder vendored for de novo generation.

Pipeline context:
    ``structural_elucidation.inference.predict.load_decoder`` loads pickled
    vocab/config and ``.pt`` weights from ``models/decoder*``, then instantiates
    ``TranslationModel`` via ``decoder.moses.models_storage.ModelsStorage``.
    Track B sampling conditions encoder latents on precursor mass and RI.

Key I/O:
    Inputs: ``vocab.nb``, ``config.nb``, and ``model*.pt`` under ``models/decoder/``.
    Outputs: sampled SMILES strings and optional per-token mean log-probabilities.
"""
