from bisect import bisect_left, bisect_right
from collections import Counter

import elucidation.chemFacts as cf
import numpy as np
import pandas as pd
import regex as re
from elucidation.binning import bin_for_mass, bin_base
from elucidation.chemFacts import (
    fmla_mass,
    canonicalize_fmla,
)
from elucidation.conveniences import lzip, qw, lmap
from elucidation.tableLoaders import fmla_finder_cache_loader


# These are the formula finder's precomputed tables.  They contain every possible formula that could be returned,
# sorted by mass.  Very small masses are in a separate table and queried differently.
# If run as a main program, the tables will be constructed and pickled.
# If not running as main and tables are not found, we will throw an exception, because building takes ~10 minutes.
#
# I'm just gonna grab a section from the Word document:
#
# All possible formulas are precomputed and stored in a table that is loaded from disk when the program starts.  After analyzing the library, I found that formulas satisfying these conditions would include nearly every formula already in the library:
#
# •	0 <= C <= 55
# •	0 <= H <= 94
# •	0 <= N <= 11
# •	0 <= O <= 22
# •	0 <= S <= 2
# •	PLUS at most one option from: Br, Br’ ([81Br]), Cl, Cl2, F, F2, I, P, and P2.  (Both isotopes of bromine are included because they are equally abundant in nature.)
#
# This gives about 44 million formulas, but fewer than 39 million have mass less than 1000 Da.
# The formulas can easily be encoded in 4 bytes (28 bits, in fact) giving a table of about 155 Mb.
# The formulas are sorted by mass and stored in a numpy array.  They are indexed by another array of 300000, one for each bin.
# To query for formulas for mass m, we find m’s bin (in the usual scheme) and go to the formula array at that index.  Then we can scan forward through the array to until we get to the first formula with mass too big to match the query mass.
# Small masses less than 50 Da are handled in a separate table.  We do retain these formulas.
# They can be useful in various calculations even though the instruments do not measure reliably in that range.
#
# Formulas are packed into 32-bit integers, which is not described in the word document.
# The choice of C can be packed into ceiling(log2(55)) = 6 bits, H in 7-bits, etc., for a total of 28 bits.
# We only have to pack when building the tables, and unpack when querying the tables.
# I call the integers "packed 6-tuples".
# To handle the possibilities that we might want additional elements or more of some element, we can execute multiple
# queries after subtracting the "new stuff" from the query mass, but it's not being used right now.


class FormulaFinder(object):
    __instances__ = {}
    __large_bins__ = None
    __large_ints__ = None
    __small_ints__ = None
    __small_masses__ = None
    __options__ = None
    __git_branch__ = None
    __timestamp__ = None

    @staticmethod
    def get_instance(lib_id):
        if lib_id is None:  # No instances have been created, so load tables into class variables.
            lib_id = 0
        if FormulaFinder.__large_bins__ is None:
            d = fmla_finder_cache_loader()
            FormulaFinder.__large_ints__ = d["large_ints"]
            FormulaFinder.__large_bins__ = d["large_bins"]
            FormulaFinder.__small_ints__ = d["small_ints"]
            FormulaFinder.__small_masses__ = d["small_masses"]
            FormulaFinder.__options__ = d["options"]
            FormulaFinder.__timestamp__ = d["timestamp"]
            FormulaFinder.__git_branch__ = d["git_branch"]
        if lib_id not in FormulaFinder.__instances__:
            FormulaFinder.__instances__[lib_id] = FormulaFinder(lib_id)

        return FormulaFinder.__instances__[lib_id]

    def __init__(self, lib_id):
        if lib_id in FormulaFinder.__instances__:
            raise Exception("try again with FormulaFinder.get_instance(lib_id)")
        self.lib_id = lib_id  # TODO: fix this to use cachier again.

    def get_timestamp(self):
        """
        This will return the timestamp of the model file, for inclusion in reports.
        @return: something like 2021-03-04-1609
        @rtype: str
        """
        return FormulaFinder.__timestamp__

    def get_gitbranch(self):
        """
        This will return the git branch of the program that built the model file, for inclusion in reports.
        @return: something like master@wcejc23c
        @rtype: str
        """
        return FormulaFinder.__git_branch__

    @staticmethod
    def formula_from_int(packed):
        """
        Given a PACKED 6-tuple of indices, we compute the molecular formula string.
        :param packed: a packed 6-tuple (32-bit integer)
        :return: formula string
        """
        return FormulaFinder.formula_from_tuple(FormulaFinder.unpack_fmla_tuple(packed))

    @staticmethod
    def mass_from_int(packed):
        """
        Given a PACKED 6-tuple of indices, we compute the mass of the formula.
        :param packed: packed 6-tuple
        :return: a mass
        """
        return FormulaFinder.mass_from_tuple(FormulaFinder.unpack_fmla_tuple(packed))

    @staticmethod
    def fmlas_for_large_mass(mass, tol):
        """
        Suggests formulas for masses >=50.  Normally called via fmlas_for_mass.
        :param mass: neutral mass
        :param tol: absolute tolerance for query (not ppm)
        :return: a list of (formula,mass) pairs.
        """
        lo = mass - tol
        if lo > 1000:
            return []
        if lo < 50:
            lo = 50
        hi = mass + tol
        if hi > 1000:
            hi = 1000
        results = []

        start_at = FormulaFinder.__large_bins__[bin_for_mass(lo)]
        large_ints = FormulaFinder.__large_ints__
        # save_i = 0
        for i in range(start_at, 50000000):  # "forever"
            # save_i = i
            w = FormulaFinder.mass_from_int(large_ints[i])
            if w < lo:
                continue
            if w > hi:
                break
            results.append((FormulaFinder.formula_from_int(large_ints[i]), w))
        # print(start_at, save_i, save_i-start_at)
        return results

    @staticmethod
    def fmlas_for_small_mass(mass, tol):
        """
        Suggests formulas for masses <50u.  Normally called via fmlas_for_mass.
        :param mass: neutral mass
        :param tol: absolute tolerance for query (not ppm)
        :return: a list of (formula,mass) pairs.
        """
        small_masses, small_ints = (
            FormulaFinder.__small_masses__,
            FormulaFinder.__small_ints__,
        )
        lo = bisect_left(small_masses, mass - tol)
        hi = bisect_right(small_masses, mass + tol)
        ints = small_ints[lo:hi]
        fmlas = lmap(FormulaFinder.formula_from_int, ints)
        masses = small_masses[lo:hi]
        return lzip(fmlas, masses)

    def fmlas_for_precursor_mz(self, mz, tol=None, ppm=None):
        return self.fmlas_for_neutral_mass(
            cf.neutralize_precursor_mz(mz, self.lib_id), tol=tol, ppm=ppm
        )

    def fmlas_for_ion_mz(self, mz, tol=None, ppm=None):
        return self.fmlas_for_neutral_mass(cf.neutralize_ion_mz(mz, self.lib_id), tol=tol, ppm=ppm)

    def fmlas_for_neutral_mass(self, neutral_mass, tol=None, ppm=None):
        """
        Suggest formulas for a mass +/- an absolute or relative tolerance.
        :param neutral_mass: query mass
        :param tol: absolute tolerance.  Can't specify both tol and ppm.
        :param ppm: relative tolerance in ppm.
        :return: a dataframe of masses and ppm deviations.
        """
        if tol is not None and ppm is not None:
            raise Exception("Can't specify both tol and ppm.")
        if neutral_mass < 50:
            if tol is None:
                tol = 0.001
            results = FormulaFinder.fmlas_for_small_mass(neutral_mass, tol)
        elif neutral_mass < 100:
            if tol is None:
                tol = 0.001
            results = FormulaFinder.fmlas_for_large_mass(neutral_mass, tol)
        else:
            if ppm is None:
                ppm = 5
            tol = neutral_mass * ppm * 1e-6
            results = FormulaFinder.fmlas_for_large_mass(neutral_mass, tol)

        if not results:
            df = pd.DataFrame(dict(cand_fmla=[], ppm=[]))
            return df.astype({"cand_fmla": str, "ppm": float})
        fmlas, masses = list(zip(*(results)))

        fmla_df = pd.DataFrame(dict(cand_fmla=fmlas, mass=masses))
        fmla_df.cand_fmla = fmla_df.cand_fmla.apply(canonicalize_fmla)

        fmla_df = fmla_df.drop_duplicates()
        fmla_df["ppm"] = ((fmla_df.mass / neutral_mass) - 1) * 1e6
        fmla_df["absppm"] = fmla_df.ppm.apply(np.fabs)
        fmla_df = fmla_df.groupby("cand_fmla").apply(lambda t: t.max())
        fmla_df = fmla_df.sort_values("absppm")[qw("cand_fmla mass ppm")].reset_index(drop=True)
        return fmla_df

    @staticmethod
    def unpack_fmla_tuple(packed):
        """
        unpack a 32-bit int into 6 small integers indexing the options list.
        :param packed: the 32-bit integer to unpack.
        :return: A 6-tuple of small integers.
        """
        in_int = packed
        unpacked = []
        for mask, shift in [
            (0xF, 4),
            (0x3, 2),
            (0x1F, 5),
            (0xF, 4),
            (0x7F, 7),
            (0x3F, 6),
        ]:
            unpacked.append(packed & mask)
            packed >>= shift
        if packed != 0:
            print(f"something's wrong: {in_int:x}, {packed:x}")
        return tuple(unpacked[::-1])

    @staticmethod
    def mass_from_tuple(t):
        """
        compute mass of a molecular formula in the form of an UNPACKED 6-tuple of indices.
        :param t: 6-tuple of indices.
        :return: mass
        """
        # It's not hard because the individual masses are stored in the options list.
        pairs = [FormulaFinder.__options__[i][t[i]] for i in range(len(FormulaFinder.__options__))]
        return sum([w for elt, w in pairs])

    @staticmethod
    def formula_from_tuple(t):
        """
        find string form af a molecular formula expressed as an UNPACKED 6-tuple of indices.
        :param t: 6-tuple of indices.
        :return: mass
        """
        # Easy because the pieces are in the options list.
        pairs = [FormulaFinder.__options__[i][t[i]] for i in range(len(FormulaFinder.__options__))]
        return "".join([elt for elt, w in pairs])


####################################################################################################
# Everything below this is for BUILDING the tables, not querying.
####################################################################################################
def _mass_from_tuple(t, options):
    """
    compute mass of a molecular formula in the form of an UNPACKED 6-tuple of indices.
    :param t: 6-tuple of indices.
    :return: mass
    """
    # It's not hard because the individual masses are stored in the options list.
    pairs = [options[i][t[i]] for i in range(len(options))]
    return sum([w for elt, w in pairs])


def pack_fmla_tuple(t):
    """
    pack a tuple of 6 small integers into a 32-bit integer.  only for BUILDING tables.
    :param t: tuple of 6 indices into the options list.
    :return: a 32-bit integer,
    """
    packed = 0
    for opt, shift in zip(t, [0, 7, 4, 5, 2, 4]):
        packed <<= shift
        packed |= opt
    return np.int32(packed)


def _build_large_mass_tables(pairs):
    """
    Builds table for large masses (>50u).
    :param pairs:  (mass, packed formula) pairs.
    :return: ints = list of packed formulas sorted by mass, and bins = indices into ints for fast querying.
    """
    # Focus on large formulas first.
    masses, ints = lzip(*pairs)
    masses, ints = list(masses), list(ints)

    # add sentinels
    masses.append(bin_base(300000))
    ints.append(np.int32(2 ** 31 - 1))
    masses = np.array(masses)
    ints = np.array(ints)

    # The formulas are sorted by mass, so so are the bins.
    binned_masses = np.array(lmap(bin_for_mass, masses))

    # We are going to save:
    # a) A list of packed formulas, sorted by mass
    # b) A list of indices into (a).
    bins = np.array([np.int32(-1)] * 300001)  # intialize bins to -1
    bins[300000] = len(ints) - 1
    for i in range(
        len(binned_masses) - 1, 0, -1
    ):  # in each OCCUPIED bin, record the index of the first formula in that bin.
        bins[binned_masses[i]] = i

    # When enumerating bin b, the stopping point is going to be the first index in bin b+1.
    # However, unoccupied bins still point to -1.  We have to fix that.
    # so we propagate downward.
    for i in range(len(bins) - 1, 0, -1):
        if bins[i] == -1:
            bins[i] = bins[i + 1]

    return ints, bins


def _count_multiplicity(elt, fmlas):
    """
    Analyzes the counts of element elt in a list of formulas.  This is used for explicating the construction of the tables.
    :param elt: The element to be analyze.
    :param fmlas: List of formulas to be analyzed.
    :return: list of triples from a Counter.  See output.
    """
    # given list of formulas, get multiplicity of element elt in each formula.
    pat = r"(?<=" + elt + r")\d*"  # a pattern to extract element + number (or no number)
    lis = [re.findall(pat, f) for f in fmlas]
    lis = [0 if x == [] else 1 if x[0] == "" else int(x[0]) for x in lis]
    # count multiplicities.
    lis = sorted(Counter(lis).items(), reverse=True)
    # unzip counts from counts of counts, and add cumulative sum.
    # Because of reverse sort, cumulative sum will be number of formulas LOST if count is cutoff.
    ks, ns = zip(*lis)
    cum = np.cumsum(ns)
    lis = lzip(ks, ns, cum)[::-1]
    return lis
