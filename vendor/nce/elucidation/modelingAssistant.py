"""Modeling-bin utilities for spectrum-to-structure ML in NCE elucidation.

Maps neutralized m/z and losses to logarithmic bin pairs for training and
runtime matching. Shared by MS2 matchers, SMILES predictors, and Chemont
models in spec2mol v21 inference.
"""

import logging
from collections import defaultdict, Counter
from functools import reduce
from itertools import zip_longest
from operator import itemgetter, add
from pprint import pprint

import numpy as np
import pandas as pd
from elucidation.binning import bin_for_mass, mass_interval, bin_base
from elucidation.chemFacts import neutralize_ion_mz
from elucidation.libraryService import LibraryService
from elucidation.rdkitUtils import (
    canonicalize_smiles,
    smiles_mass,
    fmla_for_smiles,
    fragment_smiles,
)
from elucidation.conveniences import qw, lmap

"""
Modeling bins.
A modeling bin is a string, e.g. "bin_123456" or "bin_m123456".  Strings are better than numbers because
they are easily recognized as representing bins.
The modeling bin of an m/z is always based on a "neutralized" m/z.
In the negative mode, this means adding back a proton, for example.
The reason for this is that we also want to bin SMILES fragments, and our function
returns neutral masses for those. Possibly it would have been better to ionize the SMILES instead of neutralize
the ions.
We must be able to bin m/z in our normal range 50-1000u as well as losses.
Losses are represented by negative numbers.  Since losses aren't measured by the instrument, we can have any negative loss,
even >-50.  Losses have bin names like "bin_m123456".

For illustration, consider now only m/z in the range 100-1000u and library 209:
- First of all, we add back a proton to neutralize m/z, giving m.
- m are measured +/- 5ppm, and bins are 10ppm wide.
- Since errors are relative, we must bin on a logarithmic scale.  The formula (discovery.binning) is
     bin(m) = 299999 + int( log(1000 / m) / log(1-10e-6) )
  This formula puts 50u near 0 and 1000u at 300000.
- Measurement of the mass could fall in either of two adjacent bins, bin(m-5ppm) or bin(m+5ppm), the "bin-pair" for m.
Now suppose we have two LISTS of masses, and we want to find which are common.  We can construct a bin pairs of
bins for one list, and a list of single bins bin(m) for the other list, and take the intersection. We conclude that the
corresponding masses in the second list are in the first list.
We might be wrong, because parts of the bins in bin_pair(m) are >5ppm from m.  But we might be wrong anyway; even two
exactly equal masses could represent different ions. We are making a compromise in accuracy to lower the cost of computation. In the future,
not buying a 1ppm instrument due to cost might also be a sacrifice Metabolon would make.

The general principle is that TRAINING spectra and ions stored in search tables will be represented by bin pairs (or bin ranges, see below),
while the spectra extracted at runtime to be matched against tables will be represented by single bins.

Next consider an ion in the range 50-100u. These can be measured, we assume the accuracy is only an absolute 0.0005,
which is 5ppm at 100u, dropping to 10ppm at 50u. Therefore, we need a span of up to four bins in a training spectrum
to cover the range of possible measurements.

Next consider losses.  A loss is the difference between an observed ion and the mass of the entire molecule. Here's
why losses are interesting:  During fragmentation, when a bond is broken, one piece becomes positively charged and
one negative. In the negative mode, only the negative piece is measured.  However, it may be that a certain
fragment typically becomes positively charged.  Its presence might be indirectly observed by subtracting the mass
of the negative piece from the mass of the entire molecule.

Losses can obviously be <50u.  Losses will be represented by negative numbers.  The name of a loss bin is something
like bin_m123456.  If all losses were >50u, we could simply assign loss L to bin -bin_for_mass(-L).  However,
to accommodate smaller losses, we have to shift to the left 100u.  (50 would suffice, but this lets us bin 0 to
+50u if needed later.)

Note that a loss is always neutral, whether ion minus ion or neutral minus neutral.

What is the accuracy of measurement of a loss?  The best/safest answer I could come up with was
+/-5ppm OF THE MASS OF THE ENTIRE MOLECULE.


This one is for masses and losses.
The 'ppm' is wrt big_m, which is relevant for a loss. Otherwise, big_m=m (I expect.)
"""


class ModelingAssistant:

    __instances__ = dict()

    @staticmethod
    def get_instance(lib_id):
        """
        Enforces a singleton-per-library discipline.
        @param lib_id: library
        @type lib_id: int
        @return: singleton instance for this library.
        @rtype: ModelingAssistant
        """
        if lib_id not in ModelingAssistant.__instances__:
            ma = ModelingAssistant(lib_id)
        return ModelingAssistant.__instances__[lib_id]

    def __init__(self, lib_id):

        if lib_id in ModelingAssistant.__instances__:
            raise Exception("try again with ModelingAssistant.get_instance(lib_id)")
        self.rediscovery_df = None
        self.ms1_rediscovery_df = None
        self.lib_id = lib_id
        self.ms2_df = None
        self.smiles_dict = None
        self.smiles_name_dict = None
        self.smiles_fragment_dict = None
        self.spectrum_dict = None
        self.parent_mass_dict = None
        self.fragment_compids_dict = None
        self.modeling_df = None
        self.matching_df = None
        ModelingAssistant.__instances__[lib_id] = self

    def ion_bin(self, mz):
        """
        returns a bin name for the given m/z, e.g. "bin_12345"
        @param mz: m/z of ion
        @type mz: float
        @return: e.g. "bin_12345"
        @rtype: str
        """
        if mz < 0.5:  # need a little slack for "positive" here to avoid 0.
            raise ValueError("mz must be positive")
        m = neutralize_ion_mz(mz, self.lib_id)
        if m >= 50:  # The range for bin_to_mass is [50,1000]
            return self.name_bin(bin_for_mass(m))
        else:
            return self.name_bin(-bin_for_mass(abs(m - 100)))

    def loss_bin(self, m):
        """
        returns a bin name for a loss.
        @param m: a neutral mass (not m/z) for a loss.
        @type m: float
        @return: bin name, e.g., "bin_m56789"
        @rtype: str
        """
        if m < 0.5:
            raise ValueError("loss must be positive")
        return self.name_bin(-bin_for_mass(abs(-m - 100)))  # The range for bin_to_mass is [50,1000]

    def fragment_bin(self, frag):
        """
        returns a bin name for a smiles fragment based on its neutral mass.
        @param frag: a smiles fragment
        @type frag: str
        @return: bin name, e.g., "bin_56789"
        @rtype: str
        """
        m = smiles_mass(frag)
        if m < 0.5:
            raise ValueError("mz must be positive")
        if m >= 50:  # The range for bin_to_mass is [50,1000]
            return self.name_bin(bin_for_mass(m))
        else:
            return self.name_bin(-bin_for_mass(abs(m - 100)))

    def unname_bin(self, name):
        """
        "unnames" a bin: turns the name into a representative mass.
        @param name: bin name, e.g., "bin_56789"
        @type name: str
        @return: mass
        @rtype: float
        """
        if not name.startswith("bin_m"):  # the easy case; ion mass >50
            return bin_base(int(name[4:]))
        # The hard cases: ion mass < 50 or loss
        else:
            return -bin_base(int(name[5:])) + 100

    def ion_bin_range(self, mz, m_for_ppm):
        """
        returns a *range* of bin names for an m/z +/- tolerance
        @param mz: the m/z
        @type mz: float
        @param m_for_ppm: a mass for determining tolerance ppm; could be parent mass of compound.
        @type m_for_ppm: float
        @return: a list of bin names.
        @rtype: list of str
        """
        # if a loss, it's measured with the same accuracy as the corresponding ion.
        if mz < 0.5:
            raise ValueError("mz must be positive")
        m = neutralize_ion_mz(mz, self.lib_id)
        half_width = mass_interval(m_for_ppm)
        if m >= 50:
            lo_bin = bin_for_mass(m - half_width)
            hi_bin = bin_for_mass(m + half_width)
            return [self.name_bin(b) for b in range(lo_bin, hi_bin + 1)]
        else:
            m = abs(m - 100)  # The range for bin_to_mass is [50,1000]
            lo_bin = bin_for_mass(m - half_width)
            hi_bin = bin_for_mass(m + half_width)
            return [self.name_bin(-b) for b in range(lo_bin, hi_bin + 1)]

    def loss_bin_range(self, loss, m_for_ppm):
        """
        returns a *range* of bin names for a loss +/- tolerance.
        @param loss: a neutral mass
        @type loss: float
        @param m_for_ppm: a mass for determining ppm for tolerance.
        A loss can't be measured more accurately than the precursor ion.
        @type m_for_ppm: float
        @return: a list of bin names.
        @rtype: list of str
        """

        # if a loss, it's measured with the same accuracy as the corresponding ion.
        if loss <= 0.5:
            # Occasionally something bigger than the precursor does make it into the MS2,
            # but this should be filtered out elsewhere.
            raise ValueError(f"Losses should be positive.  {loss}")
        half_width = mass_interval(m_for_ppm)
        loss = abs(-loss - 100)  # The range for bin_to_mass is [50,1000]
        lo_bin = bin_for_mass(loss - half_width)
        hi_bin = bin_for_mass(loss + half_width)
        return [self.name_bin(-b) for b in range(lo_bin, hi_bin + 1)]

    def lost_fragment_bin_range(
        self, fragment, mol_m
    ):  # if a loss, it's measured with the same accuracy as the corresponding ion.
        try:  # In case the smiles is funky:
            if type(mol_m) is str:
                mol_m = smiles_mass(mol_m)
            loss = abs(smiles_mass(fragment) - mol_m - 100)
            half_width = mass_interval(mol_m)
            lo_bin = bin_for_mass(loss - half_width)
            hi_bin = bin_for_mass(loss + half_width)
            return [self.name_bin(-b) for b in range(lo_bin, hi_bin + 1)]
        except Exception as e:
            return []

    def fragment_bin_range(self, fragment):
        """
        Gives a range of bins for a SMILES fragment based on its neutral mass.
        @param fragment: SMILES fragment.
        @type fragment: str
        @return: list of bin names
        @rtype: list of str
        """
        try:  # In case the smiles is funky:
            m = smiles_mass(fragment, lenient=False)
            half_width = mass_interval(m)
            lo_bin = bin_for_mass(m - half_width)
            hi_bin = bin_for_mass(m + half_width)
            return [self.name_bin(b) for b in range(lo_bin, hi_bin + 1)]
        except Exception as e:
            return []

    def range_spectrum(self, pairs, mol_mz):
        """
        Gives a "wide net" list of bins for an MS2 spectrum.  Each m/z contributes several bins:
        - bins for m/z's from mz-5ppm to mz+5ppm (or wider if m/z<100)
        - bins for the neutral loss: the whole molecule minus m/z.  These also cover a range, but the
            loss is "measured" with accuracy of 5ppm of the WHOLE molecule.
        Used to be called target_spectrum
        @param pairs: a raw spectrum, [(m/z, intensity), (m/z, intensity), ...]
        @type pairs: list of tuples (pairs)
        @param mol_mz: m/z of entire molecule, for computing losses.
        @type mol_mz: float
        @return: list of bin names; no intensities
        @rtype: list of str
        """
        bins = reduce(add, [self.ion_bin_range(mz, mz) for mz, i in pairs], [])
        # Don't allow ~0 losses
        bins += reduce(
            add,
            [self.loss_bin_range(mol_mz - mz, mol_mz) for mz, i in pairs if mz < mol_mz - 0.5],
            [],
        )
        return bins

    def range_spectrum_with_intensities(self, pairs, mol_mz):
        """
        Like above, but the resulting spectrum includes intensities.
        Used to be target_spectrum_with_intensities.
        @param pairs: a raw spectrum, [(m/z, intensity), (m/z, intensity), ...]
        @type pairs: list of tuples (pairs)
        @param mol_mz: m/z of entire molecule, for computing losses.
        @type mol_mz: float
        @return: dict(binname=intensity, binname=intensity, ...)
        @rtype: dict
        """
        # Don't allow ~0 losses
        ranges = [(self.ion_bin_range(mz, mz), i) for mz, i in pairs] + [
            (self.loss_bin_range(mol_mz - mz, mol_mz), i) for mz, i in pairs if mz < mol_mz - 0.5
        ]
        d = dict()
        for bins, i in ranges:
            d.update(zip_longest(bins, [], fillvalue=i))
        return d

    def bin_spectrum(self, pairs, mol_mz):
        """
        Gives a "narrow" list of bins for an MS2 spectrum.  Each m/z contributes only two bins:
        - one bin for m/z
        - one bin for the neutral loss: the whole molecule minus m/z.
        Used to be called target_spectrum
        @param pairs: a raw spectrum, [(m/z, intensity), (m/z, intensity), ...]
        @type pairs: list of tuples (pairs)
        @param mol_mz: m/z of entire molecule, for computing losses.
        @type mol_mz: float
        @return: dict(binname=intensity, binname=intensity, ...)
        @rtype: dict
        """
        bins = [(self.ion_bin(mz), i) for mz, i in pairs]
        bins += [(self.loss_bin(mol_mz - mz), i) for mz, i in pairs if mz < mol_mz - 0.5]
        # Build map from bins back to masses/losses
        d = {self.ion_bin(mz): mz for mz, i in pairs}
        d.update({self.loss_bin(mol_mz - mz): mz - mol_mz for mz, i in pairs if mz < mol_mz - 0.5})
        return dict(bins), d

    def spectrum_magnitude(self, pairs, mol_mz):
        """
        Return the vector magnitude of a raw spectrum.  This is used for matching two OBSERVED spectra.
        The theoretical starting point is simply the sqrt of sum of squares of intensities.
        However, we need to include the losses, which are assigned the same intensity as the corresponding m/z.
        Furthermore, an m/z at molecular mass will not generate a loss, so that is accounted for.
        To do all that easily, we create a bin spectrum.
        @param pairs: a raw spectrum of (m/z, intensity) pairs.
        @type pairs: list of tuples (pairs)
        @param mol_mz: mass of parent molecule, for computing losses.
        @type mol_mz: float
        @return: The magnitude
        @rtype: float
        """
        qs, d = self.bin_spectrum(pairs, mol_mz)
        return np.sqrt(sum([i * i for i in qs.values()]))

    def match_specs(self, sp1, sp2):
        """
        Calculate a match score of two spectra.  One is in range form and the other in bin form.
        (Could be in either order, actually, but it needs one of each.)
        Note that the result is NOT scaled by magnitudes of spectra.
        @param sp1: First spectrum, in range form as return by ModelingAssistant.range_spectrum_with_intensities.
        @type sp1: dict: bin name -> intensity
        @param sp2: Second spectrum, in bin form, as returned by ModelingAssistant.bin_spectrum
        @type sp2: dict: bin name -> intensity
        @return: score in [0,1]
        @rtype: float
        """
        if len(sp1) == 0 or len(sp2) == 0:
            logging.warning(f"matching an empty spectrum {sp1} {sp2}")
            return 0
        return sum([sp1[b] * sp2[b] for b in set(sp1.keys()) & set(sp2.keys())])

    def smiles_range_spectrum(self, smiles):
        """
        Generates a synthetic range spectrum for a SMILES string based on masses of its fragments.
        @param smiles: a SMILES string
        @type smiles: str
        @return: a range spectrum [bin, bin, bin]
        @rtype: list of str
        """
        mol_m = smiles_mass(smiles)
        frags = fragment_smiles(smiles)
        bins = reduce(add, [self.fragment_bin_range(frag) for frag in frags], [])
        bins += reduce(add, [self.lost_fragment_bin_range(frag, mol_m) for frag in frags], [])
        return bins

    def reshape_ion_df(self, df, cols=[]):
        """
        Given a dataframe with one row per MS2 ion, aggregates MS2 spectra for each compound and returns
        a dataframe with one row for each compound.
        @param df: dataframe with one row per MS2 ion
        @type df: dataframe
        @param cols: a list of auxilliary columns that must be preserved (e.g., SMILES)
        @type cols: list of str
        @return: dataframe with one row for each compound
        @rtype: dataframe
        """
        df = df[df.mass_ratio >= 1].copy()
        d = dict(df.groupby("comp_id").apply(lambda tab: list(zip(tab.mass, tab.mass_ratio))))
        df = df[qw("comp_id lib_id ri parent_mass name ") + cols].drop_duplicates()

        df["raw_spectrum"] = df.comp_id.apply(d.get)

        df["range_spectrum"] = [
            self.range_spectrum(spec, parent)
            for spec, parent in zip(df.raw_spectrum, df.parent_mass)
        ]
        df["bin_spectrum"] = [
            self.bin_spectrum(spec, parent) for spec, parent in zip(df.raw_spectrum, df.parent_mass)
        ]
        df["magnitude"] = [
            self.spectrum_magnitude(spec, parent)
            for spec, parent in zip(df.raw_spectrum, df.parent_mass)
        ]

        return df.reset_index(drop=True).copy()

    def get_matching_df(self):
        """
        Gets a dataframe of MS2s used for library matching against named compounds.
        @return:  dataframe with one row per compound.
        @rtype: dataframe
        """
        if self.matching_df is None:
            df = LibraryService.get_instance().get_named_ms2(self.lib_id)
            df = self.reshape_ion_df(df, cols=qw("formula smiles"))
            df.smiles = df.smiles.apply(canonicalize_smiles)
            df["fragments"] = df.smiles.apply(fragment_smiles)
            df["frag_formulas"] = df.fragments.apply(lambda f: sorted(set(map(fmla_for_smiles, f))))
            df["smiles_range_spectrum"] = df.smiles.apply(self.smiles_range_spectrum)
            self.matching_df = df
        return self.matching_df

    def get_modeling_df(self):
        """
        Gets a dataframe of MS2s usable for training ML models. We get confused by charges...
        @return:  dataframe with one row per compound.
        @rtype: dataframe
        """
        df = self.get_matching_df()
        df = df[
            ~df.smiles.str.contains("+]", regex=False) & ~df.smiles.str.contains("-]", regex=False)
        ]
        return df

    def get_rediscovery_df(self):
        """
        Gets a dataframe of MS2s of X- and NMD- compounds for rediscovery check.
        @return:  dataframe with one row per compound.
        @rtype: dataframe
        """
        if self.rediscovery_df is None:
            df = LibraryService.get_instance().get_unnamed_ms2(self.lib_id)
            df = self.reshape_ion_df(df, cols=qw("comments retired_for_comp_id"))
            self.rediscovery_df = df
        return self.rediscovery_df

    def get_MS1_rediscovery_df(self):
        """
        Gets a dataframe of MS2s of X- and NMD- compounds for rediscovery check.
        @return:  dataframe with one row per compound.
        @rtype: dataframe
        """
        if self.ms1_rediscovery_df is None:
            df = LibraryService.get_instance().get_ms1_rediscovery(self.lib_id)
            self.ms1_rediscovery_df = df
        return self.ms1_rediscovery_df

    def _get_matrix(self, col, min_count=5, comp_ids=None):
        """
        builds a 0/1 training matrix based on the lists in column 'col' of the modeling df.
        @param col: The column of feature lists in the modeling df to base the matrix on.
        @type col: str
        @param min_count: A feature must exist in at least this many compounds to be included.
        @type min_count: int
        @param comp_ids: The compounds to be in matrix, if not all.
        @type comp_ids: list of int
        @return: training matrix, (feature->column index) map, (compound ID -> row index) map
        @rtype: (array, dict, list)
        """
        df = self.get_modeling_df()
        if comp_ids is not None:
            df = df[df.comp_id.isin(comp_ids)]
        lists = df[col]
        c = Counter(reduce(add, list(lists)))
        cols = sorted([k for k, v, in c.items() if v >= min_count])
        d = {val: i for i, val in enumerate(cols)}
        a = np.zeros([len(lists), len(d)])
        for i, lis in enumerate(lists):
            for col in lis:
                if col in d:
                    a[i][d[col]] = 1
        return a, d, list(df.comp_id)

    def get_matrix_for_fragments(self, min_count=5, comp_ids=None):
        """
        builds a 0/1 training matrix based on SMILES fragment lists of the modeling df.
        @param min_count: A fragment must exist in at least this many compounds to be included.
        @type min_count: int
        @param comp_ids: The compounds to be in matrix, if not all.
        @type comp_ids: list of int
        @return: training matrix, (feature->column index) map, (compound ID -> row index) map
        @rtype: (array, dict, list)
        """
        return self._get_matrix("fragments", min_count=min_count, comp_ids=comp_ids)

    def get_matrix_for_frag_formulas(self, min_count=5, comp_ids=None):
        """
        builds a 0/1 training matrix based on SMILES fragment formula lists of the modeling df.
        @param min_count: A fragment formulo must exist in at least this many compounds to be included.
        @type min_count: int
        @param comp_ids: The compounds to be in matrix, if not all.
        @type comp_ids: list of int
        @return: training matrix, (feature->column index) map, (compound ID -> row index) map
        @rtype: (array, dict, list)
        """
        return self._get_matrix("frag_formulas", min_count=min_count, comp_ids=comp_ids)

    def get_matrix_for_bin_spectrum(self, min_count=5, comp_ids=None):
        """
        builds a 0/1 training matrix based on spectrum bin lists of the modeling df.
        @param min_count: A bin must exist in at least this many compounds to be included.
        @type min_count: int
        @param comp_ids: The compounds to be in matrix, if not all.
        @type comp_ids: list of int
        @return: training matrix, (feature->column index) map, (compound ID -> row index) map
        @rtype: (array, dict, list)
        """
        df = self.get_modeling_df()
        if comp_ids is not None:
            df = df[df.comp_id.isin(set(comp_ids))]
        dicts = df.bin_spectrum
        lists = [list(d[0].keys()) for d in dicts]
        return self._get_matrix_for_spectrum(lists, min_count=min_count, comp_ids=df.comp_id)

    def get_matrix_for_range_spectrum(self, min_count=5, comp_ids=None):
        """
        builds a 0/1 training matrix based on spectrum bin range lists of the modeling df.
        @param min_count: A bin must exist in at least this many compounds to be included.
        @type min_count: int
        @param comp_ids: The compounds to be in matrix, if not all.
        @type comp_ids: list of int
        @return: training matrix, (feature->column index) map, (compound ID -> row index) map
        @rtype: (array, dict, list)
        """
        df = self.get_modeling_df()
        if comp_ids is not None:
            df = df[df.comp_id.isin(comp_ids)]
        lists = df.range_spectrum
        return self._get_matrix_for_spectrum(lists, min_count=min_count, comp_ids=df.comp_id)

    def _get_matrix_for_spectrum(self, lists, min_count=5, comp_ids=None):
        c = Counter(reduce(add, list(lists)))
        cols = sorted([k for k, v in c.items() if v >= min_count])
        d = {val: i for i, val in enumerate(cols)}
        a = np.zeros([len(lists), len(d)])
        for i, lis in enumerate(lists):
            for col in lis:
                if col in d:
                    a[i][d[col]] = 1
        return a, d, list(comp_ids)

    def get_fragment_matrix(self):
        df = self.get_modeling_df()
        c = Counter()
        for t in df.itertuples():
            c.update(t.bin_spectrum)
        lis = [k for k, v in c.items() if v >= 5]
        d = dict([(ion, i) for i, ion in enumerate(sorted(lis))])

    def name_bin(self, b):
        """
        Given an integer bin number, it turns it into a string.  The advantage is that the number is immediately
        recognizable as a bin number.
        @param b: integer bin number
        @type b: int
        @return: bin name string
        @rtype: str
        """
        return f"bin_{b}" if b >= 0 else f"bin_m{-b}"

    def find_matching_models(self, models, features):
        """
        Given feature values (typically a candidate compound), predicts T/F for every model in list, and returns
        names of models predicted T.
        @param models: A list of models, of which one component is an sklearn model with predict method, e.g., logistic regression.
        @type models:
        @param features:
        @type features:
        @return:
        @rtype:
        """
        # We don't need to derive bins from formulas because we derive formulas from fragments.
        d = models[0]["x_cols"]  # Currently all models are the same.
        a = np.zeros(len(d))
        # for f in fragments | formulas:
        for f in features:
            if f in d:
                a[d[f]] = 1
        a = a.reshape(1, -1)
        # trace_match(models, features, a, "CHEMONTID:0004322")
        # trace_match(models, features, a, "*P(=O)(O)O")
        return [
            (m, m["model"].predict_proba(a)[0][1]) for m in models if m["model"].predict(a)[0] == 1
        ]


def trace_match(models, features, a, smiles):
    """
    Traces matching of model for debugging
    @param models:
    @type models:
    @param features:
    @type features:
    @param a:
    @type a:
    @param smiles:
    @type smiles:
    @return:
    @rtype:
    """
    # Pick out the model from the list.
    model = [m["model"] for m in models if m["y"] == smiles]
    if not model:
        print("No model for", smiles)
        return
    model, *_ = model

    print(smiles)
    pprint(features)
    d = models[0]["x_cols"]
    coefs = model.coef_[0]
    d1 = {k: coefs[v] for k, v in d.items()}
    pairs = sorted(d1.items(), key=itemgetter(1))
    df = pd.DataFrame(pairs, columns=qw("bin coef"))
    df["matched"] = df.bin.isin(set(features))
    for t in df.itertuples():
        print(t.bin, t.coef, t.matched)
    print()
    features = sorted(features)
    pprint(pd.DataFrame(dict(bin=features, coef=lmap(d1.get, features))))
    print()
    print(model.predict(a), model.predict_proba(a))


def invert_dict_of_lists(d0):
    """
    Best explained by example.
    IN: {'a': ['x', 'y', 'z'], 'b': ['w', 'x', 'y'], 'c': ['w', 'z']}
    OUT: {'x': ['a', 'b'], 'y': ['a', 'b'], 'z': ['a', 'c'], 'w': ['b', 'c']}
    """
    d1 = defaultdict(list)
    for k, lis in d0.items():
        for e in lis:
            d1[e].append(k)
    return d1
