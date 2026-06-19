"""Single-scan MS1 or MS2 spectrum wrapper with binned intensity lookup.

Rehydrates pickled scan dicts from discovery and supports fast m/z bin and
intensity queries used by MS1Analyzer and MS2Consensus.
"""

import logging
from operator import itemgetter

import numpy as np
import pandas as pd
from elucidation.binning import bin_for_mass, within_5ppm
from elucidation.conveniences import lmap
from numpy import searchsorted


# This class represents a single observed spectrum, MS1 or MS2.
class Spectrum(object):
    """One Thermo MS1 or MS2 scan with sentinel-padded mass arrays for lookup."""

    def __init__(self, d):
        """
        Create a Spectrum object based on a single scan from a Thermo file.
        :param d: a dictionary for the spectrum as pickled in a feature by discovery.
        """

        # All of these will be set from the dictionary, and this is just here to shut up pycharm.
        self.masses = None
        self.intensities = None
        self.ri = None
        self.precursor = None
        self.sample_name = None

        for k, v in d.items():
            self.__setattr__(k, v)

        # Add sentinels to masses and intensities for binary searches.
        self.bins = np.array([0] + list(bin_for_mass(self.masses)) + [1000000000])
        self.masses = np.array([0] + list(self.masses) + [9999])
        self.intensities = np.array([0] + list(self.intensities) + [0])
        self.rel_intensities = []  # Must know normalizing mass to compute.
        self.buckets = []
        self.norm_int = None

    def get_bins(self):
        return self.bins

    def get_intensity(self, mz):
        i = searchsorted(self.masses, mz)
        if within_5ppm(self.masses[i], mz):
            return self.intensities[i]
        if within_5ppm(self.masses[i - 1], mz):
            return self.intensities[i - 1]
        return 0

    def __getitem__(
        self, mz
    ):  # I am so cool. For a Spectrum sp, sp[mass] is the corresponding intensity!
        return self.get_intensity(mz)

    def get_mass_and_intensity(self, mz):
        i = searchsorted(self.masses, mz)
        if within_5ppm(self.masses[i], mz):
            return self.masses[i], self.intensities[i]
        if within_5ppm(self.masses[i - 1], mz):
            return self.masses[i - 1], self.intensities[i - 1]
        return 0, 0

    def get_rel_intensity(self, mz):
        i = searchsorted(self.masses, mz)
        if within_5ppm(self.masses[i], mz):
            return self.rel_intensities[i]
        if within_5ppm(self.masses[i - 1], mz):
            return self.rel_intensities[i - 1]
        return 0

    def get_bucket_mz(self, b):
        i = searchsorted(self.buckets, b)
        if self.buckets[i] != b:
            return -1
        return self.masses[i]

    def get_bucket_int(self, b):
        i = searchsorted(self.buckets, b)
        if self.buckets[i] != b:
            return 0
        return self.intensities[i]

    def get_bucket_rel_int(self, b):
        i = searchsorted(self.buckets, b)
        if self.buckets[i] != b:
            return 0
        return self.rel_intensities[i]

    def get_dr_tuples(self, mass):
        feature_m, feature_i = self.get_mass_and_intensity(mass)
        if feature_i == 0:  # Not sure how this happens, but it did once.
            return []
        mzs = self.masses
        ints = self.intensities
        return [
            (
                mzs[j] - feature_m,
                ints[j] / feature_i,
                mzs[j],
                ints[j],
            )  # delta, intensity ratio, buddy_mass, buddy intensity
            for j in range(1, len(mzs) - 1)  # skip the sentinels
            if 0.005 <= ints[j] / feature_i <= 200
        ]

    def rebucket(self, bin_map):
        """
        rebucket the ions to merge adjacent occupied bins.
        :param bin_map:
        :return: None, but updates self.int_dict and self.frag_mass_to_dbsmiles_dict.
        """
        self.buckets = pd.Series(self.bins).apply(bin_map.get).to_numpy()

    def get_bucket_list(
        self,
    ):  # This is a bit of a pain because sometimes there are sentinels and sometimes not.
        return [b for b in self.buckets if 0 < b <= 300000]

    def subset_buckets(self, to_keep):
        """
        Used by MS2 consensus to throw away buckets that aren't interesting (e.g., low population.)
        :param to_keep: list of buckets to keep.
        :return: None; updates Spectrum object.
        """
        inx = (pd.Series(self.buckets).isin(to_keep)).to_numpy()
        self.masses = self.masses[inx]
        self.intensities = self.intensities[inx]
        self.buckets = self.buckets[inx]
        self.bins = self.bins[inx]

    def score_buckets(self, total):
        """
        build a dictionary ranking ions (buckets) by intensity.  Used by MS2Consensus to pick normalizer.
        :param total: The total number of distinct ions in all observed spectra.
        :return: dictionary mapping bin to rank.
        """
        lis = [(i, b) for i, b in zip(list(self.intensities), list(self.buckets))]
        lis = lmap(itemgetter(1), sorted(lis))  # list of buckets ranked by intensity.
        missing = [0] * (total - len(lis))  # add fake buckets for missing buckets
        lis = missing + lis
        rank_dict = dict(zip(lis, range(1, len(lis) + 1)))  # assign position to each bucket.
        if 0 in rank_dict:  # throw out dummies
            del rank_dict[0]
        return rank_dict  # return dictionary

    def normalize(self, normalizer):
        """
        Create a dictionary of intensities relative to the normalizer.
        :param normalizer: basis of relative intensity.  Supplied by MS2Consensus.
        :return: None; creates a dictionary mapping bucket to relative intensity.
        """
        self.norm_int = self.intensities[self.buckets == normalizer][0]
        self.rel_intensities = self.intensities / self.norm_int
        self.norm_dict = dict(zip(self.buckets, self.rel_intensities))
        if max(self.norm_dict.values()) > 1:
            logging.warn("Weirdness normalizing spectrum")
            logging.warn(f"{self.norm_dict}")
