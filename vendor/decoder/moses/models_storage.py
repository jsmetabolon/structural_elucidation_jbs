"""Model registry for vendored MOSES decoder architectures.

Pipeline context:
    ``structural_elucidation.inference.predict.load_decoder`` imports
    ``ModelsStorage``, resolves the ``translation`` entry, and constructs
    ``TranslationModel`` from pickled vocab/config plus ``.pt`` weights.

Key I/O:
    ``get_model_class("translation")`` → ``TranslationModel`` constructor.
"""

from moses.trans import (
    TranslationModel,
    TranslationTrainer,
    translation_parser,
)


class ModelsStorage:
    """Registry mapping model names to class, trainer, and CLI parser."""

    def __init__(self):
        self._models = {}
        self.add_model(
            "translation",
            TranslationModel,
            TranslationTrainer,
            translation_parser,
        )

    def add_model(self, name, class_, trainer_, parser_):
        """Register a decoder architecture under ``name``."""
        self._models[name] = {
            "class": class_,
            "trainer": trainer_,
            "parser": parser_,
        }

    def get_model_names(self):
        """Return registered model name strings (e.g. ``"translation"``)."""
        return list(self._models.keys())

    def get_model_trainer(self, name):
        """Return the trainer class registered for ``name``."""
        return self._models[name]["trainer"]

    def get_model_class(self, name):
        """Return the ``nn.Module`` class registered for ``name``."""
        return self._models[name]["class"]

    def get_model_train_parser(self, name):
        """Return the argparse builder registered for ``name``."""
        return self._models[name]["parser"]
