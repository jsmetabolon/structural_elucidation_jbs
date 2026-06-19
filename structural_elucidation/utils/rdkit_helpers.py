"""Legacy RDKit and spectrum helpers retained for encoder compatibility.

Pipeline context:
    ``spec2vec`` and related binning utilities predate the
    ``structural_elucidation`` package refactor. ``spectrum.py`` delegates
    vectorization here so training-time binning behavior stays unchanged.

Note:
    Prefer ``structural_elucidation.utils.chemistry`` for new SMILES
    canonicalization; this module keeps British spellings
    (``canonicalise_smiles``) for backward compatibility.
"""

from __future__ import annotations

import pickle
from collections.abc import Iterable, Set
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem
from rdkit.Chem import Mol


# ── SMILES validation ────────────────────────────────────────────────
def check_smiles(smi: str) -> bool:
    """Return whether ``smi`` parses as a valid RDKit molecule.

    Args:
        smi: SMILES string.

    Returns:
        ``True`` when RDKit accepts the structure.
    """
    mol = Chem.MolFromSmiles(smi)
    return mol is not None


def canonicalise_smiles(smi: str) -> str:
    """Return isomeric canonical SMILES for ``smi``.

    Args:
        smi: Input SMILES (must already be valid).

    Returns:
        Canonical SMILES string.

    Raises:
        TypeError: When ``smi`` does not parse (``MolFromSmiles`` returns ``None``).
    """
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def normalize_smiles(smi: str, canonical: bool, isomeric: bool) -> str:
    """Normalize SMILES with explicit canonical / isomeric flags.

    Args:
        smi: Input SMILES.
        canonical: Passed to ``MolToSmiles(canonical=...)``.
        isomeric: Passed to ``MolToSmiles(isomericSmiles=...)``.

    Returns:
        Normalized SMILES string.
    """
    normalized = Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=canonical, isomericSmiles=isomeric)
    return normalized


# ── Molecular descriptors ────────────────────────────────────────────
def get_species(mol: Mol) -> set[str]:
    """Return the set of element symbols present in ``mol``.

    Args:
        mol: RDKit molecule.

    Returns:
        Unique element symbols (e.g. ``{"C", "H", "O"}``).
    """
    species: set[str] = set()
    for atom in mol.GetAtoms():
        species.add(atom.GetSymbol())
    return species


# ── Spectrum vectorization ───────────────────────────────────────────
def spec2vec(
    spectrum: Set[tuple[float, float]] | Iterable[tuple[float, float]],
    minmass: float,
    maxmass: float,
    resolution: int,
) -> NDArray[np.float32]:
    """Convert peak tuples into a binned, max-normalized intensity vector.

    Each bin accumulates intensities for peaks rounded to the nearest bin
    index. The vector is divided by its maximum so the top peak is 1.0.

    Args:
        spectrum: ``(m/z, intensity)`` peak set or iterable.
        minmass: Minimum m/z (Da) for bin zero.
        maxmass: Maximum m/z (Da) for the last bin.
        resolution: Bins per Da (e.g. 2 => 0.5 Da bins).

    Returns:
        Float32 vector of length ``int((maxmass - minmass) * resolution) + 1``.
    """
    minmass = float(minmass)
    maxmass = float(maxmass)
    resolution = int(resolution)
    length = int((maxmass - minmass) * resolution) + 1
    vec = np.zeros(length, dtype=np.float32)
    if len(spectrum) == 0:
        return vec

    for mass, abund in spectrum:
        try:
            mz = float(mass)
            intensity = float(abund)
            idx = int(round((mz - minmass) * resolution))
            if 0 <= idx < length:
                vec[idx] = vec[idx] + intensity
        except Exception:
            continue
    max_val = float(np.max(vec)) if len(vec) else 0.0
    if max_val <= 0:
        return vec
    return vec / max_val


def spec2tuples(spectrum: str) -> set[tuple[float, float]]:
    """Parse NIST-style newline peak text into ``(m/z, intensity)`` tuples.

    Args:
        spectrum: Multi-line string with ``mass intensity`` per line.

    Returns:
        Set of peak tuples.
    """
    tuples: set[tuple[float, float]] = set()

    for peak in spectrum.split("\n"):
        mass, abund, *c = peak.split(" ")
        mass = float(mass)
        abund = float(abund)
        tuples.add((mass, abund))
    return tuples


def spec2tuples_fromCSV(spectrum_file: str | Path) -> set[tuple[float, float]]:
    """Load comma-separated m/z,intensity rows from a CSV file.

    Args:
        spectrum_file: Path to two-column CSV (m/z, intensity).

    Returns:
        Set of peak tuples.
    """
    tuples: set[tuple[float, float]] = set()
    entries = Path(spectrum_file).read_text().split("\n")
    for i in range(0, len(entries) - 1):
        entry = entries[i].split(",")
        mass = float(entry[0])
        abund = float(entry[1])
        tuples.add((mass, abund))
    return tuples


def min_max_mass(spectrum: str) -> tuple[float, float]:
    """Return the first and last m/z values in NIST-style peak text.

    Args:
        spectrum: Multi-line NIST peak block.

    Returns:
        ``(start_mz, end_mz)`` from first and last parsed peaks.
    """
    for peak in spectrum.split("\n"):
        mass, abund, *c = peak.split(" ")
        mass = round(float(mass), 2)
        start = mass
        break
    for peak in spectrum.split("\n"):
        mass, abund, *c = peak.split(" ")
        mass = round(float(mass), 2)
    end = mass
    return (start, end)


def remove_large(spectrum: Iterable[tuple[float, float]]) -> bool:
    """Return ``True`` when any peak exceeds 500 Da (legacy training filter)."""
    remove = False
    for mass, abund in spectrum:
        mass = float(mass)
        if mass > 500:
            remove = True
    return remove


def writeout_spec(spectrum: str, outfile: str | Path) -> None:
    """Write NIST peaks as comma-separated m/z, relative-intensity rows.

    Intensities are scaled by the maximum abundance in the block.

    Args:
        spectrum: NIST-style peak text.
        outfile: Output CSV path.
    """
    writeout = open(outfile, "w")
    max_abund = 0
    for peak in spectrum.split("\n"):
        mass, abund, *c = peak.split(" ")
        if float(abund) > max_abund:
            max_abund = float(abund)
    max_abund = round(max_abund, 2)
    vec = np.zeros(197500)
    for peak in spectrum.split("\n"):
        mass, abund, *c = peak.split(" ")
        abund = float(abund)
        rel_abund = abund / max_abund
        writeout.write(mass + ", " + str(rel_abund) + "\n")
    writeout.close()
    return


# ── Value scaling ────────────────────────────────────────────────────
def normalise_value(x: float, min_value: float, max_value: float) -> float:
    """Linearly scale ``x`` from ``[min_value, max_value]`` to ``[0, 1]``."""
    return (x - min_value) / (max_value - min_value)


def retrieve_value(x: float, min_value: float, max_value: float) -> float:
    """Inverse of ``normalise_value``: map ``[0, 1]`` back to the original range."""
    return x * (max_value - min_value) + min_value


# ── Pickle and text I/O ────────────────────────────────────────────────
def save_obj(obj: object, filename: str) -> None:
    """Persist ``obj`` to ``filename.pkl`` with highest pickle protocol.

    Args:
        obj: Pickle-serializable object.
        filename: Base path without ``.pkl`` suffix.
    """
    with open(filename + ".pkl", "wb") as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(filename: str) -> object:
    """Load an object from ``filename.pkl``.

    Args:
        filename: Base path without ``.pkl`` suffix.

    Returns:
        Unpickled object.
    """
    with open(filename + ".pkl", "rb") as f:
        return pickle.load(f)


def writeout_smiles(filename: str | Path, smiles_set: Iterable[str]) -> None:
    """Write one SMILES per line to a text file.

    Args:
        filename: Output path.
        smiles_set: SMILES strings to write.
    """
    outfile = open(filename, "w")
    for smiles in smiles_set:
        outfile.write(smiles + "\n")
    outfile.close()
    return


def read_smiles(filename: str | Path) -> set[str]:
    """Read line-separated SMILES, keeping only RDKit-valid structures.

    Invalid lines are printed and skipped.

    Args:
        filename: Input text file path.

    Returns:
        Set of valid SMILES strings.
    """
    lines = Path(filename).read_text().split("\n")
    smiles: set[str] = set()
    for smi in lines:
        if check_smiles(smi):
            smiles.add(smi)
        else:
            print("Problem with SMILES: ", smi)
    return smiles


def save_csv(directory: str | Path, data: NDArray[Any]) -> None:
    """Save a numpy array as comma-separated text.

    Args:
        directory: Output file path.
        data: Array to write.
    """
    with open(directory, "wb") as f:
        np.savetxt(f, data, fmt="%s", delimiter=",")
    return


def create_valid(datafile: str | Path, valid_size: int) -> set[str]:
    """Collect up to ``valid_size`` SMILES lines from a text file (legacy helper).

    Args:
        datafile: Line-separated SMILES file.
        valid_size: Maximum number of lines to retain.

    Returns:
        Set of SMILES strings (first N non-empty lines).
    """
    validset: set[str] = set()
    lines = Path(datafile).read_text().split("\n")
    for smi in lines:
        if len(validset) < valid_size:
            validset.add(smi)
    return validset


# ── Dictionary attribute access ──────────────────────────────────────
class dotdict(dict):
    """Dict subclass exposing keys as attribute access (``d.key``)."""

    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
