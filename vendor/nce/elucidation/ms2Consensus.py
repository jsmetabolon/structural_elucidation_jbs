import logging
from collections import Counter, defaultdict
from functools import reduce
from operator import add

import elucidation.chemFacts as cf
import numpy as np
import pandas as pd
from elucidation.binning import bin_for_mass, bin_base
from elucidation.binning import create_bucket_map
from elucidation.chemFacts import isotope_mass_difference
from elucidation.chemontPredictor import ChemontPredictor
from elucidation.conveniences import lmap, qw, lzip, is_empty_df
from elucidation.formulaFinder import FormulaFinder
from elucidation.modelingAssistant import ModelingAssistant
from elucidation.ms2Ion import MS2Ion
from elucidation.ms2Matcher import Ms2Matcher
from elucidation.postgresSmilesDb import PostgresSmilesDb
from elucidation.predictedSpectrumMatcher import PredictedSpectrumMatcher
from elucidation.rdkitUtils import (
    fmla_for_smiles,
    smiles_mass,
    canonicalize_smiles,
    modify_smiles,
)
from elucidation.riPredictor import RiPredictor
from elucidation.smilesPredictor import SmilesPredictor


def analyze_ms2_isotopes(buckets):
    """
    Finds possible isotopic relationships among sets of (MS2) ions.
    :param buckets: list of buckets conta
    :return:
    """
    isotopes = defaultdict(list)
    for isotope in qw("C13 Cl37 Br81 S34 O18 N15"):
        diff = isotope_mass_difference(isotope)
        iso_buckets = [bin_for_mass(diff + bin_base(bucket)) for bucket in buckets]
        iso_buckets = reduce(add, [[b - 1, b, b + 1] for b in iso_buckets])
        for b in iso_buckets:
            isotopes[b].append(isotope)
    return [",".join(isotopes[b]) for b in buckets]


class MS2Consensus(object):
    def __init__(self, feature_):
        """
        Creates an MS2Consensus object for the consensus MS2 spectrum.
        :param feature_: The feature whose consensus is being constructed.
        """
        # If you add something here, does it need to be exported to the feature's pickle
        self.support = []
        self.normalizer = None
        self.ms2_ions = None
        self.ms2_ion_dict = {}
        self.feature = feature_  # feature to which this spectrum belongs.
        self.consensus_df = None
        self.ms2_library_report = None
        self.min_ms2_intensity = None
        self.ms2_counts = defaultdict(int)  # sample name -> count
        self.ms2_smiles_predictions = None
        self.structure_prediction_report = None
        self.chemont_predictions = None
        self.predicted_spectrum_match_report = None
        self.modified_spectrum_match_report = None
        self.model_versions = {}
        self.rediscovery_of = []
        self.rediscovery_report = None
        self.modeled_ions, self.unmodeled_ions = None, None

    def as_dict(self):
        """
        Returns a dict representing the object for pickling.
        @return: a dict to be pickled.  I just said that.
        @rtype: dict
        """
        for klass, name in [
            (FormulaFinder, "formula finder"),
            (ChemontPredictor, "chemont"),
            (PredictedSpectrumMatcher, "smiles2ms2"),
            (SmilesPredictor, "ms2smiles"),
        ]:
            obj = klass.get_instance(self.feature.lib_id)
            self.model_versions[name] = (
                obj.get_gitbranch(),
                obj.get_timestamp(),
            )

        d = {
            attr: self.__getattribute__(attr)
            for attr in qw(
                "normalizer consensus_df ms2_library_report "
                + "ms2_smiles_predictions chemont_predictions "
                + "structure_prediction_report predicted_spectrum_match_report "
                + "model_versions modified_spectrum_match_report rediscovery_report rediscovery_of"
            )
        }
        d["ms2_ions"] = {k: v.as_dict() for k, v in self.ms2_ion_dict.items()}
        return d

    def build(self, spectra, min_ms2_intensity):
        """
        Constructs the consensus spectrum.
        :param min_ms2_intensity: The minimum relative intensity that will be retained.
        :return: None; updates self.ms2_ions (and other stuff)
        """
        self.support = spectra
        if len(self.support) == 0:
            self.normalizer = 0
            return
        self.min_ms2_intensity = min_ms2_intensity
        all_bins = Counter(reduce(add, [list(sp.bins) for sp in self.support], []))
        bucket_map, _ = create_bucket_map(set(all_bins))
        for sp in self.support:
            sp.rebucket(bucket_map)
        bucket_counts = Counter(
            reduce(add, [list(sp.get_bucket_list()) for sp in self.support], [])
        )
        important_buckets = [k for k, v in bucket_counts.items() if v >= len(self.support) / 2]
        if not important_buckets:
            logging.info("Can't build consensus spectrum -- too little support for ions.")
            logging.info(f"{len(self.support)} spectra in support for {self.feature.serial}")
            return
        rank_totals = Counter()  # Not really using it as a counter.
        for sp in self.support:
            sp.subset_buckets(important_buckets)
            rank_totals.update(sp.score_buckets(len(important_buckets)))
        self.normalizer = max(rank_totals.keys(), key=rank_totals.get)
        self.support = [sp for sp in self.support if self.normalizer in sp.get_bucket_list()]
        for sp in self.support:
            sp.normalize(self.normalizer)
        all_buckets = sorted(set(reduce(add, [sp.get_bucket_list() for sp in self.support], [])))
        all_buckets = [b for b in all_buckets if b > 0]  # TODO: Figure out why this happens.
        self.ms2_ion_dict = {b: MS2Ion(b, self).get_me() for b in all_buckets}
        self.ms2_ion_dict = {k: v for k, v in self.ms2_ion_dict.items() if v}  # eliminate None.
        # Throw out small relative intensities.
        self.ms2_ion_dict = {
            k: v for k, v in self.ms2_ion_dict.items() if v.med_ratio >= self.min_ms2_intensity
        }
        self.ms2_ions = lmap(self.ms2_ion_dict.get, sorted(self.ms2_ion_dict.keys()))

    def create_ion_consistency_report(self, fmlas):  # TODO  What did I want to do...?
        """
        Creates the ion consistency report, a table showing which ions in the consensus are consistent with
        which candidate formulas in the list fmlas.
        :param fmlas: candidate formulas for the precursor.
        :return: A dataframe in which rows are ions and columns are candidate precursor formulas.
        """
        ion_reports = [ms2.ion_consistency_report(fmlas) for ms2 in self.ms2_ions]
        ion_report = pd.DataFrame(ion_reports).sort_values("rel_intens", ascending=False)
        ion_report["isotopy"] = analyze_ms2_isotopes(ion_report.bucket)
        ion_report = ion_report[
            qw("bucket mz mass_range isotopy scan_count rel_intens") + fmlas + ["cand_fmlas"]
        ]
        return ion_report

    def create_consensus_df(self):
        """
        Creates a dataframe summarizing the consensus spectrum.
        :return: None; sets self.ion_df
        """
        if not self.ms2_ions:  # We are sad, because none of the ions had enough support.
            return
        dicts = [ion.consensus_stats() for ion in self.ms2_ions]

        self.consensus_df = pd.DataFrame(dicts).sort_values("mz")
        # change below to use masses (?)  # TODO
        self.consensus_df["isotopy"] = analyze_ms2_isotopes(self.consensus_df.bucket)

    def build_ms2_library_report(self):
        """
        Causes the consensus spectrum to be matched against library spectra and summarizes results.
        :return: None; sets self.ms2_library_report.
        """
        matcher = Ms2Matcher.get_instance(self.feature.lib_id)
        self.ms2_library_report = matcher.get_best_ms2_library_matches(
            self.consensus_df.mz, self.consensus_df.norm_int, self.feature.mz, top=10
        )

    def build_rediscovery_report(self):
        """
        Causes the consensus spectrum to be matched against library spectra and summarizes results.
        :return: None; sets self.ms2_library_report.
        """
        matcher = Ms2Matcher.get_instance(self.feature.lib_id)
        self.rediscovery_report = matcher.get_rediscovery_matches(
            self.consensus_df.mz, self.consensus_df.norm_int, self.feature.mz, self.feature.ri,
        )
        self.rediscovery_of = list(self.rediscovery_report.comp_id)

    def build_predicted_spectrum_match_report(self):
        """
        Build report of matches of MS2 to SMILES strings in the database,
        and to common modifications of SMILES strings in the database.
        @return: None; it sets attributes
        @rtype: NoneType
        """
        self.predicted_spectrum_match_report = self.help_build_spectrum_match_report(
            self.feature.smilesdb_df
        )
        dfs = [self.build_modified_spectrum_match_report(t) for t in self._get_reactions()]
        dfs = [df for df in dfs if not is_empty_df(df)]
        if dfs:
            self.modified_spectrum_match_report = pd.concat(dfs, ignore_index=True)

    def _get_reactions(self):
        """
        returns SMARTS strings, etc., for the common modifications used in searching.
        @return: list of tuples
        @rtype: list of tuples
        """
        return [
            ("sulfate", "[O!H0:2]  >> [O:2]-S(=O)(=O)O", 79.95681486000001),
            ("glucuronide", "[O!H0:2]  >> [O:2]-C1OC(C(=O)O)C(O)C(O)C(O)1", 176.032087976,),
            ("acetyl", "[O!H0:2]  >> [O:2]-C(=O)C", 42.010564684),
            ("glycine", "[C:1](=O)O >> [C:1](=O)-NCC(=O)O", 57.021463716),
            ("glutamine", "[C:1](=O)O >> [C:1](=O)-NC(CCC(=O)N)C(=O)O", 128.058577496),
            ("carnitine", "[C:1](=O)O >> [C:1](=O)-OC(C[N+](C)(C)C)CC(=O)[O-]", 143.094628656,),
            ("taurine", "[C:1](=O)O >> [C:1](=O)-NCCS(=O)(=O)O", 107.0040994),
            (
                "carboxyethyl",
                "[C:2](C(=O)O)[#7;H2:1] >> [C:2](C(=O)O)[#7:1]C(C)C(=O)O",
                72.021129368,
            ),
        ]

    def build_modified_spectrum_match_report(self, modification):
        """
        Builds a match report for SMILES strings in matching database that can be modified with this modification
        into the something with the feature's mass.
        @param modification: a tuple, e.g. ("acetyl", "[O!H0:2]  >> [O:2]-C(=O)C", 42.010564684),
        with name of modification, SMARTS reaction for modification, and mass added by modification.
        @type modification: tuple
        @return: a DataFrame with a report for this modification.
        @rtype: DataFrame
        """
        mod_name, smarts, delta_m = modification
        pgsmiles = PostgresSmilesDb.get_instance()
        # unsulfated_mz = self.feature.mz - smiles_mass("*OS(=O)(=O)O") + smiles_mass("*O")
        unmodified_mz = self.feature.mz - delta_m
        df = pgsmiles.get_records_for_neutral_mass(
            cf.neutralize_precursor_mz(unmodified_mz, self.feature.lib_id),
            ppm_mass=self.feature.mz,
        )
        df = df[qw("source name smiles")]
        modifications = []
        for t in df.itertuples():
            for mod_smiles in modify_smiles(smarts, t.smiles):
                d = dict(
                    smiles=mod_smiles,
                    name=f"{mod_name} + {t.name}",
                    source=f"{mod_name} + {t.source}",
                    formula=fmla_for_smiles(mod_smiles),
                    exact_mass=smiles_mass(mod_smiles),
                )
                modifications.append(d)

        return self.help_build_spectrum_match_report(pd.DataFrame(modifications))

    def help_build_spectrum_match_report(self, smilesdb_df):
        """
        Given a dataframe with a column of SMILES strings, this method:
        - consolidates names and sources of rows with the same SMILES string.
        - uses PredictedSpectrumMatcher to score the strings and set up matching information to
          display the match graphically in the feature report.
        @param smilesdb_df: The dataFrame to be aggregated and scored.
        @type smilesdb_df: DataFrame
        @return: Scored DataFrame
        @rtype: DataFrame
        """

        def group_matches(tab):
            names = [n for n in tab["name"] if type(n) is str]
            t = (
                tab.smiles.max(),
                sorted(set(names)),  # dammit!
                sorted(set(tab.source)),
                tab.formula.max(),
                tab.exact_mass.max(),
            )
            return pd.DataFrame([t], columns=qw("smiles name source formula mass"))

        matcher = PredictedSpectrumMatcher.get_instance(self.feature.lib_id)

        if is_empty_df(smilesdb_df):
            return

        smilesdb_df = smilesdb_df.groupby("smiles").apply(group_matches).reset_index(drop=True)

        spectrum = lzip(self.consensus_df.mz, self.consensus_df.norm_int)
        parent_mz = self.feature.mz
        ma = ModelingAssistant.get_instance(self.feature.lib_id)
        q_spec, q_map = ma.bin_spectrum(spectrum, parent_mz)

        self.modeled_ions, self.unmodeled_ions = matcher.check_modeling(q_spec)
        # This one only computes scores.
        # Scores will be recomputed with matchings in a sec.
        # Classic space/time tradeoff... If there are a gazillion SMILES, we'll filter first.

        score_dict = matcher.score_many_smiles_for_spectrum(
            set(smilesdb_df.smiles), q_spec, self.modeled_ions
        )

        # First of all, report <= 20 results.
        good_scores = sorted(score_dict.values(), reverse=True)[:20]
        if good_scores[0] < 0.1:
            return

        # And second, they shouldn't be horrible compared to best.
        min_allowed = max(0.1, 0.5 * good_scores[0])
        good_scores = [score for score in good_scores if score >= min_allowed]
        score_dict = {
            smiles: score for smiles, score in score_dict.items() if score >= good_scores[-1]
        }
        smiles_df = smilesdb_df[smilesdb_df.smiles.isin(set(score_dict.keys()))].copy()
        smiles_df["score"] = smiles_df.smiles.apply(score_dict.get)

        # This one will score AND record matched ions, adding match column
        df = matcher.match_spectrum_to_smiles_df(smiles_df, parent_mz, spectrum, self.modeled_ions)
        rip = RiPredictor.get_instance(self.feature.lib_id)
        df["ri_pred"] = df.smiles.apply(rip.predict_for_smiles, other_ri=self.feature.ri)

        def remap_to_mz(bin_list):
            return tuple([(q_map[b], q_spec[b]) for b in bin_list])

        successes = list(df.spectrum_match)
        failures = [tuple(self.modeled_ions - success) for success in successes]
        successes = lmap(tuple, successes)
        unmodeled = [tuple(self.unmodeled_ions)] * len(df)

        successes = lmap(remap_to_mz, successes)
        failures = map(remap_to_mz, failures)
        unmodeled = map(remap_to_mz, unmodeled)

        df.spectrum_match = lmap(tuple, zip(successes, failures, unmodeled))
        return df

    def build_ms2_smiles_prediction_report(self):
        """
        Use SmilesPredictor and ML models to predict SMILES substructures based on ions in spectrum.
        @return: no return; sets an attribute to save a DataFrame with report.
        @rtype: NoneType
        """
        spectrum = lzip(self.consensus_df.mz, self.consensus_df.norm_int)
        parent_mass = self.feature.mz
        predictor = SmilesPredictor.get_instance(self.feature.lib_id)
        self.ms2_smiles_predictions = predictor.predict_for_spectrum(parent_mass, spectrum)

    def build_chemont_prediction_report(self):
        """
        Uses ChemontPredictor and its ML models to predict Chemical Ontology categories from ions in MS2.
        @return: None; puts a report in self.chemont_predictions
        @rtype:
        """
        spectrum = lzip(self.consensus_df.mz, self.consensus_df.norm_int)
        parent_mass = self.feature.mz
        predictor = ChemontPredictor.get_instance(self.feature.lib_id)
        self.chemont_predictions = predictor.predict_for_spectrum(parent_mass, spectrum)

    # currently not doing ion-loss combos.
    def build_structure_prediction_report(self):
        all_structure_predictions = reduce(
            add, [ms2_ion.matched_smiles for ms2_ion in self.ms2_ions], []
        )
        logging.debug(f"all_struct pred {len(all_structure_predictions)}")
        if not all_structure_predictions:
            return
        all_smiles, all_probs = zip(*all_structure_predictions)
        all_smiles = lmap(canonicalize_smiles, all_smiles)
        ion_loss_score = defaultdict(float)
        for smiles, p in zip(all_smiles, all_probs):
            ion_loss_score[smiles] += 1 + p
        df = pd.DataFrame(
            [
                dict(
                    smiles=s,
                    ion_loss_score=ion_loss_score[s],
                    formula=fmla_for_smiles(s),
                    name="",
                    source="ion_loss_combo",
                    exact_mass=smiles_mass(s),
                )
                for s in all_smiles
            ]
        )
        df = df.sort_values("ion_loss_score", ascending=False)
        if len(df) > 30:
            df = df.head(30)
        # self.structure_prediction_report = self.help_build_spectrum_match_report(df)

    # I used this for a slide in a presentation.  Run it and look at the results if you want...
    def mass_triangle(self):
        a = np.array(self.consensus_df.mz).round(4)
        arr = [[a[j] - a[i] for i in range(15)] for j in range(15)]
        arr = np.array(arr)
        arr2 = arr.astype(str)
        ff = self.feature.lds.get_fmla_finder()
        for i in range(15):
            for j in range(i):
                arr2[j][i] = ""
                df = ff.fmlas_for_ion_mz(arr[i][j])
                if len(df) > 0:
                    fmlas = [
                        f for f in df.cand_fmla if not (set(f) & set("FSP"))
                    ]  # TODO sometimes we want these.
                    if fmlas:
                        arr2[j][i] = fmlas[0]
                if abs(arr[i][j] - 1.0033) < 0.001:
                    arr2[j][i] = "[13C]"
        cols = list(a.round(4).astype(str))
        df = pd.DataFrame(arr2, columns=cols, index=cols)
        return df
