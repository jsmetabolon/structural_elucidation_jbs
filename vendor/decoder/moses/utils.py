"""SMILES tokenization, vocabulary, and training helpers for the MOSES decoder.

Pipeline context:
  Pickled ``SmilesOneHotVocab`` objects in ``models/decoder/vocab.nb`` supply
  token IDs for ``TranslationModel``. ``structural_elucidation.inference.predict``
  also calls ``vocab.string2ids`` when caching encoder means for alignment.

Key I/O:
  ``SmilesVocab.string2ids`` / ``ids2string`` — SMILES ↔ token sequences.
  ``SmilesOneHotVocab`` — one-hot embedding matrix copied into the model.
"""

import random
import re
import torch
import numpy as np
import pandas as pd
from collections import UserList, defaultdict
from rdkit import Chem


REGEX_SML = r'Cl|Br|[#%\)\(\+\-1032547698:=@CBFIHONPS\[\]cionps]'


def normalize_smiles(smi, canonical, isomeric):
    normalized = Chem.MolToSmiles(
        Chem.MolFromSmiles(smi), canonical=canonical, isomericSmiles=isomeric
    )
    return normalized

def compute_reconstruction_trans(model, train_loader):
    samples = []
    real_smiles = []
    rand_smiles_converted = []
    for i, input_batch in enumerate(train_loader):
        randomized_smiles, canonical_smiles, regression_targets = zip(
                *input_batch
        )
        mu, logvar, z, kl_loss, recon_loss = model(
                randomized_smiles, canonical_smiles
        )
        current_samples = model.sample(len(mu), z=mu)
        samples.extend(current_samples)
        real_smiles.extend([model.tensor2string(i_x) for i_x in canonical_smiles])
        rand_smiles_converted.extend([model.tensor2string(i_x) for i_x in randomized_smiles])
    samples = pd.DataFrame({'REAL_CAN':real_smiles, 'GENERATED': samples, 'RANDOMIZED_CONVERTED':rand_smiles_converted})
    samples['MATCH'] = samples['REAL_CAN'] == samples['GENERATED']
    print(samples.head())
    total = len(samples)
    match = samples['MATCH'].sum()
    pct_match = samples['MATCH'].mean()
    print(f'Total: {total} Matched: {match} Percent Matched {pct_match}')

    return pct_match


# https://pytorch.org/docs/stable/data.html#torch.utils.data.DataLoader
def set_torch_seed_to_all_gens(_):
    seed = torch.initial_seed() % (2**32 - 1)
    random.seed(seed)
    np.random.seed(seed)


class SS:
    bos = '<bos>'
    eos = '<eos>'
    pad = '<pad>'
    unk = '<unk>'


def smiles_tokenize(smiles):
    """Split a SMILES string into regex-defined atom and bond tokens."""
    return re.findall(REGEX_SML, smiles)


class SmilesVocab:
    """Character-level SMILES vocabulary with BOS/EOS/PAD/UNK specials."""

    @classmethod
    def from_data(cls, data, tokenizer=smiles_tokenize, max_smiles=1000000, *args, **kwargs):
        """Build vocabulary from an iterable of SMILES strings."""
        chars = set()
        for string in data[:max_smiles]:
            chars.update(tokenizer(string))

        return cls(chars, *args, **kwargs)

    def __init__(self, chars, tokenizer=smiles_tokenize, ss=SS):
        """Initialize from a set of token characters plus special symbols."""
        if (ss.bos in chars) or (ss.eos in chars) or \
                (ss.pad in chars) or (ss.unk in chars):
            raise ValueError('SS in chars')

        all_syms = sorted(list(chars)) + [ss.bos, ss.eos, ss.pad, ss.unk]

        self.ss = ss
        self.c2i = {c: i for i, c in enumerate(all_syms)}
        self.i2c = {i: c for i, c in enumerate(all_syms)}
        self.tokenizer = tokenizer

    def __len__(self):
        return len(self.c2i)

    @property
    def bos(self):
        return self.c2i[self.ss.bos]

    @property
    def eos(self):
        return self.c2i[self.ss.eos]

    @property
    def pad(self):
        return self.c2i[self.ss.pad]

    @property
    def unk(self):
        return self.c2i[self.ss.unk]

    def char2id(self, char):
        if char not in self.c2i:
            return self.unk

        return self.c2i[char]

    def id2char(self, id):
        if id not in self.i2c:
            return self.ss.unk

        return self.i2c[id]

    def string2ids(self, string, add_bos=False, add_eos=False):
        """Tokenize ``string`` and return integer token IDs."""
        ids = [self.char2id(c) for c in self.tokenizer(string)]

        if add_bos:
            ids = [self.bos] + ids
        if add_eos:
            ids = ids + [self.eos]

        return ids

    def ids2string(self, ids, rem_bos=True, rem_eos=True):
        """Decode token IDs back to a SMILES string."""
        if len(ids) == 0:
            return ''
        if rem_bos and ids[0] == self.bos:
            ids = ids[1:]
        if rem_eos and ids[-1] == self.eos:
            ids = ids[:-1]

        string = ''.join([self.id2char(id) for id in ids])

        return string




class SmilesOneHotVocab(SmilesVocab):
    """``SmilesVocab`` with a fixed one-hot embedding matrix for the decoder."""

    def __init__(self, *args, **kwargs):
        super(SmilesOneHotVocab, self).__init__(*args, **kwargs)
        self.vectors = torch.eye(len(self.c2i))

class Logger(UserList):
    def __init__(self, data=None):
        super().__init__()
        self.sdata = defaultdict(list)
        for step in (data or []):
            self.append(step)

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.data[key]
        elif isinstance(key, slice):
            return Logger(self.data[key])
        else:
            ldata = self.sdata[key]
            if isinstance(ldata[0], dict):
                return Logger(ldata)
            else:
                return ldata

    def append(self, step_dict):
        super().append(step_dict)
        for k, v in step_dict.items():
            self.sdata[k].append(v)

    def save(self, path):
        df = pd.DataFrame(list(self))
        df.to_csv(path, index=None)


