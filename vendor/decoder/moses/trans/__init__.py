"""Translation VAE subpackage for the vendored MOSES decoder.

Pipeline context:
    Re-exports ``TranslationModel`` (inference), ``TranslationTrainer``
    (offline training), and ``translation_parser`` (CLI hyperparameters).
    ``predict.load_decoder`` constructs ``TranslationModel`` from this module.

Public API:
    ``TranslationModel``, ``TranslationTrainer``, ``translation_parser``.
"""

from .config import get_parser as translation_parser
from .model import TranslationModel
from .trainer import TranslationTrainer

__all__ = ['trans_parser', 'TranslationModel', 'TranslationTrainer']
