"""Shared trainer base class for vendored MOSES decoder models.

Pipeline context:
    ``TranslationTrainer`` subclasses ``MosesTrainer`` for offline decoder
    training. Inference does not use this module directly; it is kept so
    ``models_storage`` can register trainer and parser entry points.
"""

from abc import ABC, abstractmethod
from torch.utils.data import DataLoader
from moses.utils import set_torch_seed_to_all_gens


class MosesTrainer(ABC):
    """Base trainer providing DataLoader construction for MOSES decoders."""

    @property
    def n_workers(self):
        """Return ``config.n_workers``, mapping 1 → 0 for in-process loading."""
        n_workers = self.config.n_workers
        return n_workers if n_workers != 1 else 0

    def get_collate_device(self, model):
        """Device for batch collation; CPU when worker processes are enabled."""
        n_workers = self.n_workers
        return 'cpu' if n_workers > 0 else model.device

    def get_dataloader(self, model, data, collate_fn=None, shuffle=True):
        """Build a PyTorch ``DataLoader`` with MOSES seeding and collation."""
        if collate_fn is None:
            collate_fn = self.get_collate_fn(model)
        return DataLoader(data, batch_size=self.config.n_batch,
                          shuffle=shuffle,
                          num_workers=self.n_workers, collate_fn=collate_fn,
                          worker_init_fn=set_torch_seed_to_all_gens
                          if self.n_workers > 0 else None)

    def get_collate_fn(self, model):
        """Return batch collation callable; subclasses override for SMILES tuples."""
        return None

    @abstractmethod
    def get_vocabulary(self, data):
        """Build a token vocabulary from training SMILES strings."""
        pass

    @abstractmethod
    def fit(self, model, train_data, val_data=None):
        """Train ``model`` on ``train_data``, optionally validating on ``val_data``."""
        pass
