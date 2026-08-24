#!/usr/bin/env python3
"""Persistent on-disk cache for verifier predicted-spectrum scores.

The NCE matcher call is by far the most expensive operation on the verifier
path, and the same (SMILES, observed spectrum, lib) key recurs both across
the many enumerator input CSVs inside one ``verify`` run and across repeat
runs. This cache keys scores by a SHA-1 digest of the full scoring key so
the file stays small regardless of spectrum JSON size, and persists as plain
CSV so it adds no dependencies.

Caching is output-identical: a hit returns the exact float a fresh matcher
call would produce for the same key (the scorer is deterministic).
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path

import pandas as pd

# (smiles, observed_ion_bins_json, observed_loss_bins_json, lib_id)
ScoreKey = tuple[str, str, str, int]


def key_digest(key: ScoreKey) -> str:
    """Stable SHA-1 digest for one scoring key."""
    smiles, ion_json, loss_json, lib_id = key
    raw = "\x1f".join((str(smiles), str(ion_json), str(loss_json), str(int(lib_id))))
    return hashlib.sha1(raw.encode("utf-8")).hexdigest()


class VerifierScoreCache:
    """Score cache with optional CSV persistence (``key_sha1,score`` columns)."""

    def __init__(self, path: Path | str | None = None) -> None:
        """Load existing cache from *path* when present; corrupt files are ignored."""
        self.path = Path(path) if path is not None else None
        self._scores: dict[str, float] = {}
        self._n_new = 0
        if self.path is not None and self.path.exists():
            try:
                frame = pd.read_csv(self.path, dtype={"key_sha1": str})
                self._scores = dict(
                    zip(frame["key_sha1"].astype(str), pd.to_numeric(frame["score"], errors="coerce"))
                )
                self._scores = {k: float(v) for k, v in self._scores.items() if pd.notna(v)}
            except Exception as exc:  # a corrupt cache must never block scoring
                logging.warning("Ignoring unreadable score cache %s (%s)", self.path, exc)
                self._scores = {}

    def __len__(self) -> int:
        """Number of cached scores."""
        return len(self._scores)

    @property
    def n_new(self) -> int:
        """Number of scores added since load/last save."""
        return self._n_new

    def get(self, key: ScoreKey) -> float | None:
        """Return the cached score for *key* or ``None`` on miss."""
        return self._scores.get(key_digest(key))

    def put(self, key: ScoreKey, score: float) -> None:
        """Store *score* for *key*."""
        digest = key_digest(key)
        if digest not in self._scores:
            self._n_new += 1
        self._scores[digest] = float(score)

    def save(self) -> None:
        """Atomically persist the cache when a path is configured and dirty."""
        if self.path is None or self._n_new == 0:
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        frame = pd.DataFrame(
            {"key_sha1": list(self._scores.keys()), "score": list(self._scores.values())}
        )
        tmp = self.path.with_name(self.path.name + ".tmp")
        frame.to_csv(tmp, index=False)
        tmp.replace(self.path)
        self._n_new = 0
