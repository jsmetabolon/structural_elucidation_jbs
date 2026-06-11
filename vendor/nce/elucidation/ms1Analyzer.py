from collections import defaultdict, Counter
from functools import reduce
from operator import add, itemgetter

import numpy as np
import pandas as pd
import regex
from elucidation.binning import bin_for_mass, within_5ppm
from elucidation.buddyFunction import build_delta_df
from elucidation.chemFacts import (
    split_fmla,
    molecular_ion_name,
    neutralize_precursor_mz,
    count_for_ratio,
)
from elucidation.libraryService import LibraryService
from elucidation.ms1Buddy import MS1Buddy
from elucidation.conveniences import is_empty_df, qw, lmap


class MS1Analyzer:
    def __init__(self, feature, spectra):
        self.spectra = spectra
        self.feature = feature
        self.buddy_lib_df = None
        self.isotope_df = None
        self.mz = feature.mz
        self.ion_names = None
        self.inferred_masses = []
        self.bright_features = []
        self.buddy_df = None
        self.iso_count = defaultdict(int)
        self.labeling_df = None

    def as_dict(self):
        """
        Turns instance into a dict for easier (version-free) pickling.
        @return: A dict including all the important results of the MS1 analysis.
        @rtype:
        """
        return {
            attr: self.__getattribute__(attr)
            for attr in qw(
                "buddy_df buddy_lib_df isotope_df ion_names inferred_masses bright_features"
            )
        }

    def extract_buddies_from_one_spectrum(self, sp):
        """
        Given an MS1 spectrum, extracts all its (delta, ratio) tuples.
        @param sp: an MS1 spectrum
        @type sp: Spectrum
        @return: a list of (delta, ratio) tuples.
        @rtype: list of tuples
        """
        return sp.get_dr_tuples(self.mz)

    def set_up_binning_for_labeling(self, ppm):
        """
        Given the labeling df from buddy function, set up windows for binning deltas.
        @param ppm: ppm (usually at precursor mass.)
        @type ppm: float
        @return: edges for delta-binning bins
        @rtype: list of floate
        """

        labeling_df = self.labeling_df
        ppm_window = ppm * 1e-6 * np.maximum(self.mz, self.mz + labeling_df.theor_delta)
        labeling_df["lo"] = labeling_df.theor_delta - ppm_window  # feature mass
        labeling_df["hi"] = labeling_df.theor_delta + ppm_window
        bin_edges = sorted(set(list(labeling_df.lo) + list(labeling_df.hi)))
        labeling_df["bin_lo"] = np.digitize(labeling_df.lo, bins=bin_edges)
        labeling_df["bin_hi"] = np.digitize(labeling_df.hi, bins=bin_edges)
        labeling_df["bin_range"] = [
            set(range(bin_lo, bin_hi))
            for bin_lo, bin_hi in zip(labeling_df.bin_lo, labeling_df.bin_hi)
        ]
        self.labeling_df = labeling_df
        return bin_edges

    def is_good_isotope(self, s, good_features):  # Can't be (e.g.) C13-2 if no C13-1.
        if s.endswith("-2]") and s[:-3] + "-1]" not in good_features:
            return False
        if s.endswith("-3]") and (
            s[:-3] + "-1]" not in good_features or s[:-3] + "-2]" not in good_features
        ):
            return False
        return True

    def extract_buddies(self, ppm=2):
        """
        Takes the feature's MS1 spectra and aggregates them to form a consensus.
        When possible, labels the ions with appropriate ion names (adducts,isotopes, etc.)
        @param ppm: tolerance to be applied to deltas.
        @type ppm: float
        @return: None, but sets self.buddy_df (and some other stuff.)
        @rtype: NoneType
        """
        self.spectra = self.feature.ms1_spectra
        # list of (delta, feature mass, buddy mass, feature intensity, buddy intensity)
        all_dr_pairs = reduce(add, lmap(self.extract_buddies_from_one_spectrum, self.spectra), [])
        # delta, feature mass, buddy mass, feature intensity, buddy intensity
        self.labeling_df = build_delta_df(self.mz, self.feature.lib_id)

        bin_edges = self.set_up_binning_for_labeling(ppm)  #
        lab_df = self.labeling_df

        tuple_bins = np.digitize(lmap(itemgetter(0), all_dr_pairs), bins=bin_edges)
        delta_bins = defaultdict(list)
        for t, bin in zip(all_dr_pairs, tuple_bins):
            delta_bins[bin].append(t)

        # Round 1: Determine feature labels that have m-H/m+H buddy.
        # This is because we are looking for the relationship of the feature to the metabolite, m-H/m+H.
        mol_ion = molecular_ion_name(self.feature.lib_id)
        good_features = [mol_ion]  # This gives base_ion a boost in the counting.
        for t in lab_df[lab_df.buddy.isin({mol_ion, "x"})].itertuples():
            dr_pairs = reduce(add, [delta_bins[b] for b in t.bin_range], [])
            if dr_pairs and len(dr_pairs) >= 0.1 * len(self.spectra) and len(dr_pairs) > 1:
                new_buddy = MS1Buddy(
                    self, t, dr_pairs
                ).get_me()  # What the heck, do this twice, maybe.
                if new_buddy is not None:
                    good_features.append(t.feature)

        good_features = [s for s in good_features if self.is_good_isotope(s, good_features)]

        good_features = Counter(good_features)
        good_features[mol_ion] += good_features["x"]
        good_features = sorted(good_features.keys(), key=itemgetter(1), reverse=True)
        # We are relying on sorting to put, e.g., x[Cl37-1] before x[Cl37-2]

        self.ion_names = [f for f in good_features if f != ["x"]]

        # Round 2: Actually build the buddies.
        good_features.append("x")
        processed_bins = set()
        buddies = []
        for t in lab_df[lab_df.feature.isin(good_features)].itertuples():
            # for t in lab_df.itertuples():
            dr_pairs = reduce(add, [delta_bins[b] for b in t.bin_range], [])
            if dr_pairs and len(dr_pairs) >= 0.1 * len(self.spectra) and len(dr_pairs) > 1:
                new_buddy = MS1Buddy(self, t, dr_pairs).get_me()
                if new_buddy is not None:
                    buddies.append(new_buddy)
                    processed_bins |= t.bin_range

        buddy_df = pd.DataFrame([bud.get_summary() for bud in buddies if bud is not None])
        if len(buddy_df) > 0:
            self.bright_features = list(buddy_df[buddy_df.ratio > 1.001].feature)

        # Collect any deltas that have not been accounted for by labels.  Bin these uniformly.
        extra_deltas = reduce(
            add, [lis for b, lis in delta_bins.items() if b not in processed_bins], []
        )
        buddy_bins = defaultdict(list)
        for t in extra_deltas:
            buddy_bins[bin_for_mass(t[2])].append(t)
        extra_buddies = [
            MS1Buddy(self, None, deltas).get_me()
            for b, deltas in sorted(buddy_bins.items())
            if len(deltas) >= 0.25 * len(self.spectra) and len(deltas) > 1
        ]

        extra_summaries = [bud.get_summary() for bud in extra_buddies if bud is not None]
        extra_buddy_df = pd.DataFrame(extra_summaries)

        if len(extra_buddy_df) > 0:
            extra_buddy_df = extra_buddy_df[extra_buddy_df.c_o_v < 0.25]
            extra_buddy_df = extra_buddy_df.sort_values("tally").tail(10)
        self.buddy_df = pd.concat([buddy_df, extra_buddy_df], ignore_index=True)
        if len(self.buddy_df) > 0:
            self.buddy_df = self.buddy_df[self.buddy_df.delta.abs() > 0.1]
            self.buddy_df = self.buddy_df[self.buddy_df.tally >= 0.1 * len(self.spectra)]
            self.buddy_df = self.buddy_df.sort_values("delta")

    def build_buddy_library_df(self):
        """
        Looks at the MS1 buddy m/zs and tries to find matching ions in the library.
        Creates a dataframe of these for user's report on feature.
        @return: None; sets self.buddy_lib_df.
        @rtype: NoneType
        """
        buddy_df = self.buddy_df
        if is_empty_df(buddy_df):
            return
        inferred_mass = neutralize_precursor_mz(self.mz, self.feature.lib_id)
        im_bins = defaultdict(list)
        for m in buddy_df.inferred_m:
            if isinstance(m, float) and not np.isnan(m):
                im_bins[bin_for_mass(m)].append(m)
        self.inferred_masses = [-inferred_mass] + sorted(lmap(np.median, im_bins.values()))
        lib_df = LibraryService.get_instance().get_library_ions(self.feature.lib_id)
        med_ri = self.feature.peak_df.ri.median()
        lib_df = lib_df[
            (lib_df.ri - lib_df.ri_window <= med_ri) & (med_ri <= lib_df.ri + lib_df.ri_window)
        ].sort_values("neutral_mass")
        dfs = [
            lib_df[within_5ppm(np.fabs(m), lib_df.neutral_mass)]
            for m in self.inferred_masses
            if type(m) is float and not np.isnan(m)
        ]
        if len(dfs) == 0:
            return
        lib_df = pd.concat(dfs, ignore_index=True).drop_duplicates()

        lib_df = lib_df[
            qw("comp_id ri chemical_id neutral_mass mol_formula name")
        ].drop_duplicates()
        self.buddy_lib_df = lib_df.reset_index(drop=True)

    def add_prediction(self, row, iso, r):
        """
        Helps to build self.isotope_df for further analysis and  user's report.
        @param row: This is a dict to which a prediction is to be added (if any made),
        e.g., predicted number of C given C13 ratio.
        @type row: dict
        @param iso: name of an isotope, such as "[C13-1]"
        @type iso: str
        @param r: observed intensity ratio with feature m/z.
        @type r: float
        @return: None; adds to row dict
        @rtype: NoneType
        """
        atom_ranges = {
            "C": (0.5, 80),
            "Br": (0.5, 3),
            "Cl": (0.5, 3),
            "O": (-2, -1),
            "S": (0.5, 3),
            "N": (-2, -1),
        }  # no O, N.
        atom_count = count_for_ratio(iso, r)
        atom = regex.findall(r"[A-Z][a-z]?", iso)[0]
        atom_lo, atom_hi = atom_ranges[atom]
        if atom_lo <= atom_count <= atom_hi:
            row["predict"] = atom_count
            row[atom + "_predict"] = atom_count

    def analyze_isotopes(self):
        """
        Searches for isotopic relationships in self.buddy_df and builds self.isotope_df for further analysis and
        elucidation report.
        @return: None; sets self.isotope_df
        @rtype: NoneType
        """

        if is_empty_df(self.buddy_df):
            return
        # trim down buddy_df just to ease debugging
        iso_df = self.buddy_df[qw("delta feature buddy tally ratio c_o_v")][
            self.buddy_df.feature.str.startswith("x")
        ].copy()
        iso_df.columns = qw("delta feature buddy tally ratio c_o_v")
        if len(iso_df) == 0:
            return
        rows = []

        iso_set = qw("C13 Cl37 Br81 S34 O18 N15")
        for iso in iso_set:
            # First "has":  Does the feature have this isotope as a buddy?
            iso1 = f"[{iso}-1]"
            df = iso_df[(iso_df.buddy == f"x{iso1}") & (iso_df.feature == "x")]
            if len(df) > 0:  # Yes... probably.
                d = dict(list(df.iterrows())[0][1])
                row = {s: d[s] for s in qw("tally delta ratio c_o_v")}
                row["isotope"] = d["buddy"]
                row["has_is"] = "has"
                self.add_prediction(row, iso, row["ratio"])  # add a number (maybe)
                if "predict" in row:
                    rows.append(row)
                    self.iso_count[iso] = row[
                        "predict"
                    ]  # convenient to have them in this form, too.
            # Now "is":  Is this feature an isotope of one of its buddies?
            df = iso_df[(iso_df.buddy == "x") & (iso_df.feature == ("x" + iso1))]
            if len(df) > 0:
                d = dict(list(df.iterrows())[0][1])
                row = {s: d[s] for s in qw("tally delta ratio c_o_v")}
                row["isotope"] = d["feature"]
                row["has_is"] = "is"
                self.add_prediction(row, iso, 1 / row["ratio"])
                if "predict" in row:
                    rows.append(row)
                    self.iso_count[iso] = row[
                        "predict"
                    ]  # convenient to have them in this form, too.
        # Deuterium.
        # Deuterium is always artificially introduced, so we can't make a prediction based on natural abundance.
        # First "has"
        df = iso_df[(iso_df.buddy == "x[D-1]") & (iso_df.feature == "x")]
        if len(df) > 0:
            d = dict(list(df.iterrows())[0][1])
            row = {s: d[s] for s in qw("tally delta ratio c_o_v")}
            row["isotope"] = d["buddy"]
            row["has_is"] = "has"
            rows.append(row)
        # Now "is"
        df = iso_df[(iso_df.buddy == "x") & (iso_df.feature == "x[D-1]")]
        if len(df) > 0:
            d = dict(list(df.iterrows())[0][1])
            row = {s: d[s] for s in qw("tally delta ratio c_o_v")}
            row["isotope"] = d["feature"]
            row["has_is"] = "is"

        if len(rows) > 0:
            df = pd.DataFrame(rows)
            df = df.sort_values("delta")
            self.isotope_df = df
            return
        return

    def is_formula_isotopically_ok(self, fmla):
        """
        Test if a molecular formula meets isotope requirements, e.g., if a Cl37 isotope has been identified,
        does the formula have a Cl?
        @param fmla: a molecular formula like C13H20Cl
        @type fmla: str
        @return: True if formula has necessary isotopes; False if any are missing.
        @rtype: bool
        """
        charge, d = split_fmla(fmla)
        if abs(d["C"] - self.iso_count["C13"]) > 4:
            return False
        if abs(d["Cl"] - self.iso_count["Cl37"]) > 0.5:
            return False
        if abs(d["Br"] - self.iso_count["Br81"]) > 0.5:
            return False
        # S is special because S34 doesn't always show up. If it does, we demand it in the formula.
        if self.iso_count["S34"] >= 0.5:
            if abs(d["S"] - self.iso_count["S34"]) > 0.5:
                return False
        return True

    def library_unknown_extract_buddies(self, comp_id, lib_id, ppm=2):
        """
        This is only used for elucidating X-compounds.
        It's a big hack to use the library to build a buddy_df to fake out the rest of the MS1 analysis.
        It could be greatly simplified if the requirements were looked at more carefully.
        @param comp_id: comp id of X-compound
        @type comp_id: int
        @param lib_id: library of X-compound
        @type lib_id: int
        @param ppm: parts per million for analyzing deltas.
        @type ppm: float
        @return: None, but sets self.buddy_df (principally)
        @rtype: NoneType
        """
        df = LibraryService.get_instance().get_ms1(comp_id, lib_id)
        q_ion = next(df[df.name.isin(qw("m-H m+H"))].itertuples())
        df = df[qw("mass weighting")][~df.name.isin(qw("m-H m+H"))]
        all_dr_pairs = [
            (t.mass - q_ion.mass, t.weighting / 100, t.mass, t.weighting / 100)
            for t in df.itertuples()
        ]

        self.labeling_df = build_delta_df(q_ion.mass, lib_id)
        bin_edges = self.set_up_binning_for_labeling(ppm)  #
        lab_df = self.labeling_df
        tuple_bins = np.digitize(lmap(itemgetter(0), all_dr_pairs), bins=bin_edges)
        delta_bins = defaultdict(list)
        for t, bin in zip(all_dr_pairs, tuple_bins):
            delta_bins[bin].append(t)

        # Round 1: Determine feature labels that have m-H/m+H buddy.
        # This is because we are looking for the relationship of the feature to the metabolite, m-H/m+H.
        mol_ion = molecular_ion_name(self.feature.lib_id)
        good_features = [mol_ion]  # This gives base_ion a boost in the counting.
        for t in lab_df[lab_df.buddy.isin({mol_ion, "x"})].itertuples():
            dr_pairs = reduce(add, [delta_bins[b] for b in t.bin_range], [])
            if dr_pairs:
                new_buddy = MS1Buddy(
                    self, t, dr_pairs
                ).get_me()  # What the heck, do this twice, maybe.
                if new_buddy is not None:
                    good_features.append(t.feature)

        good_features = [s for s in good_features if self.is_good_isotope(s, good_features)]
        good_features = Counter(good_features)
        good_features[mol_ion] += good_features["x"]
        good_features = sorted(good_features.keys(), key=itemgetter(1), reverse=True)
        # We are relying on sorting to put, e.g., x[Cl37-1] before x[Cl37-2]

        self.ion_names = [f for f in good_features if f != ["x"]]

        # Round 2: Actually build the buddies.
        good_features.append("x")
        processed_bins = set()
        buddies = []
        for t in lab_df[lab_df.feature.isin(good_features)].itertuples():
            # for t in lab_df.itertuples():
            dr_pairs = reduce(add, [delta_bins[b] for b in t.bin_range], [])
            if dr_pairs:
                new_buddy = MS1Buddy(self, t, dr_pairs).get_me()
                if new_buddy is not None:
                    buddies.append(new_buddy)
                    processed_bins |= t.bin_range

        buddy_df = pd.DataFrame([bud.get_summary() for bud in buddies if bud is not None])
        if len(buddy_df) > 0:
            self.bright_features = list(buddy_df[buddy_df.ratio > 1.001].feature)

            # Collect any deltas that have not been accounted for by labels.  Bin these uniformly.
        extra_deltas = reduce(
            add, [lis for b, lis in delta_bins.items() if b not in processed_bins], []
        )
        buddy_bins = defaultdict(list)
        for t in extra_deltas:
            buddy_bins[bin_for_mass(t[2])].append(t)
        extra_buddies = [
            MS1Buddy(self, None, deltas).get_me() for b, deltas in sorted(buddy_bins.items())
        ]

        # Take only the 10 best extras.
        extra_summaries = [bud.get_summary() for bud in extra_buddies if bud is not None]
        extra_buddy_df = pd.DataFrame(extra_summaries)
        self.buddy_df = pd.concat([buddy_df, extra_buddy_df], ignore_index=True)
        if len(self.buddy_df) > 0:
            self.buddy_df = self.buddy_df.sort_values("delta")
