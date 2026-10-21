"""Bit-exact NCE evidence binning for synthesized observed MS2 lookups.

This module reproduces the integer bin ids used by the vendored NCE
``elucidation.modelingAssistant.ModelingAssistant`` (``ion_bin`` / ``loss_bin``)
so that observed fragment/loss evidence built from an external pipeline lines up
exactly with the predicted-spectrum model keys used by
:mod:`structural_elucidation.nce.verifier_v2`.

Only ``math``/``numpy`` are imported here (no pandas/rdkit/torch), so the
functions are cheap to unit-test in any environment.

Conventions (mirrors ``vendor/nce/elucidation/binning.py`` and
``vendor/nce/elucidation/modelingAssistant.py``):

* ``bin_for_mass(m) = 299999 + int32(log(1000/m) / log(1 - 10ppm))``.
* Ion m/z are *neutralized* by the electron-mass convention before binning
  (``+e`` in positive mode, ``-e`` in negative mode).
* An ion at neutral mass ``m`` maps to ``bin_for_mass(m)`` when ``m >= 50`` and
  ``-bin_for_mass(|m - 100|)`` otherwise.
* A neutral loss ``L`` (``parent_mz - fragment_mz``) maps to
  ``-bin_for_mass(L + 100)`` (always negative).
* ``nce_bin_name`` turns a signed int into the model key: ``bin_{n}`` for
  non-negative ids, ``bin_m{n}`` for negative ids.  The verifier re-applies this
  naming to the integers we store, so storing signed ints keeps ion vs loss keys
  distinct and matched against the trained models.
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence

import numpy as np

# ── grid + neutralization constants (must match vendored NCE) ─────────────────
_TEN_PPM_FACTOR: float = 1.0 - 10.0 / 1e6
_LOG_TEN_PPM_FACTOR: float = math.log(_TEN_PPM_FACTOR)
# Must match ``vendor/nce/elucidation/chemFacts.electron_mass`` exactly so the
# neutralized ion m/z (and therefore the bin ids) line up with the NCE models.
ELECTRON_MASS: float = 0.00054858
NEGATIVE_LIBS: frozenset[int] = frozenset({209, 305})
POSITIVE_LIBS: frozenset[int] = frozenset({400, 402})


def ion_mode_for_lib(lib_id: int) -> str:
    """Return the canonical ion mode (``negative``/``positive``) for a library."""
    return "negative" if int(lib_id) in NEGATIVE_LIBS else "positive"


def bin_for_mass(mass: float) -> int:
    """Assign a mass to the legacy NCE logarithmic 10 ppm bin grid.

    Mirrors ``vendor/nce/elucidation/binning.bin_for_mass`` including the
    ``numpy.int32`` truncation-toward-zero behaviour.
    """
    value = float(mass)
    if value <= 0 or not np.isfinite(value):
        return 0
    return int(299999 + np.int32(math.log(1000.0 / value) / _LOG_TEN_PPM_FACTOR))


def neutralize_ion_mz(mz: float, lib_id: int) -> float:
    """Neutralize an observed ion m/z using the NCE electron-mass convention."""
    value = float(mz)
    if int(lib_id) in POSITIVE_LIBS:
        return value + ELECTRON_MASS
    if int(lib_id) in NEGATIVE_LIBS:
        return value - ELECTRON_MASS
    return value


def ion_bin_int(mz: float, lib_id: int) -> int:
    """Signed integer bin id for an observed ion m/z (mirrors ``ion_bin``)."""
    m = neutralize_ion_mz(mz, lib_id)
    if m >= 50:
        return bin_for_mass(m)
    return -bin_for_mass(abs(m - 100.0))


def loss_bin_int(loss: float) -> int:
    """Signed integer bin id for a neutral loss (mirrors ``loss_bin``).

    ``loss`` is ``parent_mz - fragment_mz`` (a positive neutral loss).  The
    returned id is always negative so ``nce_bin_name`` renders it as ``bin_m*``.
    """
    return -bin_for_mass(abs(-float(loss) - 100.0))


def nce_bin_name(bin_id: int) -> str:
    """Format a signed integer bin id the way NCE spectrum dicts name keys."""
    value = int(bin_id)
    return f"bin_{value}" if value >= 0 else f"bin_m{-value}"


def observed_bins_for_spectrum(
    pairs: Sequence[tuple[float, float]] | Iterable[tuple[float, float]],
    parent_mz: float,
    lib_id: int,
    *,
    min_loss_da: float = 0.5,
) -> tuple[list[int], list[int]]:
    """Compute sorted, de-duplicated signed ion and loss bin ids for a spectrum.

    Replicates ``ModelingAssistant.bin_spectrum`` narrow binning: every peak
    contributes one ion bin, and peaks below ``parent_mz - min_loss_da``
    additionally contribute a neutral-loss bin.

    Args:
        pairs: Iterable of ``(mz, intensity)`` peaks (intensity is unused here).
        parent_mz: Precursor m/z used for neutral-loss computation.
        lib_id: Metabolon library id (controls neutralization + ion mode).
        min_loss_da: Minimum fragment-to-precursor gap to record a loss bin.

    Returns:
        ``(ion_bins, loss_bins)`` as sorted lists of unique signed integers.
    """
    parent = float(parent_mz)
    ion_bins: set[int] = set()
    loss_bins: set[int] = set()
    for peak in pairs:
        try:
            mz = float(peak[0])
        except (TypeError, ValueError, IndexError):
            continue
        if not np.isfinite(mz) or mz < 0.5:
            continue
        ion_bins.add(ion_bin_int(mz, lib_id))
        if np.isfinite(parent) and mz < parent - float(min_loss_da):
            loss = parent - mz
            if loss > 0.5:
                loss_bins.add(loss_bin_int(loss))
    return sorted(ion_bins), sorted(loss_bins)
