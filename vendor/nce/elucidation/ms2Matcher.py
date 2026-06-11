# coding: utf-8
from elucidation.binning import minus_5ppm, plus_5ppm
from elucidation.modelingAssistant import ModelingAssistant
from elucidation.conveniences import qw
import pandas as pd


class Ms2Matcher(object):
    """
    The Ms2Matcher class is for matching an observed spectrum to the library to get a dataframe of the most similar
    spectra in the library.  The parent masses need not be equal.
    The matching considers both ions and losses. The maximum score is theoretically 1, when something is matched to itself.
    But unless precursor masses are equal, we can't expect to match both an ion and its corresponding loss, so in practice,
    0.5 is about as high as it gets.
    """

    __instances__ = {}  # singleton-ish: one per library.

    @staticmethod
    def get_instance(lib_id):
        """
        Gets the singleton Ms2Matcher object for the given library.
        @param lib_id: library
        @type lib_id: int
        @return: the singleton object
        @rtype: Ms2Matcher
        """
        if lib_id not in Ms2Matcher.__instances__:
            _ = Ms2Matcher(lib_id)
        return Ms2Matcher.__instances__[lib_id]

    def __init__(self, lib_id):
        if lib_id in Ms2Matcher.__instances__:
            raise Exception("try again with Ms2Matcher.get_instance(lib_id)")
        self.lib_id = lib_id
        self.ma = ModelingAssistant.get_instance(lib_id)
        Ms2Matcher.__instances__[lib_id] = self

    def get_best_ms2_library_matches(self, masses, intensities, parent_mz, top=10):
        """
        Compares an observed MS2 spectrum with give masses and intensities to entire MLIMS library.
        Returns a dataframe with the most similar entries as defined by cosine similarity.
        @param masses: m/z's actually of spectrum
        @type masses: list of float
        @param intensities: corresponding intensities
        @type intensities: list of float
        @param parent_mz: m/z of precursor
        @type parent_mz: float
        @param top: number of matches to return
        @type top:int
        @return: a dataframe of matches
        @rtype: DataFrame
        """
        pairs = list(zip(masses, intensities))
        ma = self.ma
        df = ma.get_matching_df()
        test_spec = ma.range_spectrum_with_intensities(pairs, parent_mz)
        test_magnitude = ma.spectrum_magnitude(pairs, parent_mz)
        df["match_score"] = df.bin_spectrum.apply(lambda qs: ma.match_specs(qs[0], test_spec))
        df.match_score = df.match_score / (df.magnitude * test_magnitude)
        df = df[qw("match_score lib_id comp_id parent_mass ri formula smiles name")]
        df = df[df.match_score > 0.1].sort_values("match_score", ascending=False).head(top).copy()
        return df

    def get_rediscovery_matches(self, masses, intensities, parent_mz, ri):
        """
        Compares an observed MS2 spectrum with give masses and intensities to all X- and NMD- compounds with same mass.
        Returns a dataframe with the most similar entries as defined by cosine similarity.
        @param masses: m/z's actually of spectrum
        @type masses: list of float
        @param intensities: corresponding intensities
        @type intensities: list of float
        @param parent_mz: m/z of precursor
        @type parent_mz: float
        @param ri: RI of feature
        @type ri: int
        @return: a dataframe of matches
        @rtype: DataFrame
        """
        pairs = list(zip(masses, intensities))
        ma = self.ma

        df = ma.get_rediscovery_df()
        df = df[
            (minus_5ppm(parent_mz) <= df.parent_mass) & (df.parent_mass <= plus_5ppm(parent_mz))
        ]
        df = df[(df.ri - ri).abs() <= 50].copy()
        test_spec = ma.range_spectrum_with_intensities(pairs, parent_mz)
        test_magnitude = ma.spectrum_magnitude(pairs, parent_mz)
        df["match_score"] = df.bin_spectrum.apply(lambda qs: ma.match_specs(qs[0], test_spec))
        df.match_score = df.match_score / (df.magnitude * test_magnitude)
        df = df[qw("match_score lib_id comp_id parent_mass ri name comments retired_for_comp_id")]
        df.columns = qw("match_score lib_id comp_id parent_mass ri name comments retired_for")
        df = (
            df[(df.match_score > 0.9) | ((df.ri - ri).abs() <= 15)]
            .sort_values("match_score", ascending=False)
            .copy()
        )

        ms1_df = ma.get_MS1_rediscovery_df()
        ms1_df = ms1_df[
            (minus_5ppm(parent_mz) <= ms1_df.parent_mass)
            & (ms1_df.parent_mass <= plus_5ppm(parent_mz))
        ]
        ms1_df = ms1_df[(ms1_df.ri - ri).abs() <= 30].copy()
        ms1_df = ms1_df[
            qw("match_score lib_id comp_id parent_mass ri name comments retired_for_comp_id")
        ]
        ms1_df.columns = qw("match_score lib_id comp_id parent_mass ri name comments retired_for")
        ms1_df = ms1_df[~ms1_df["comp_id"].isin(df["comp_id"])]
        df = pd.concat([df, ms1_df])
        return df
