import logging
from functools import reduce
from operator import add, itemgetter

import pandas as pd
from elucidation.modelingAssistant import ModelingAssistant
from elucidation.rdkitUtils import smiles_mass, fmla_for_smiles, fragment_smiles
from elucidation.conveniences import qw, lmap
from elucidation.tableLoaders import smiles2ms2_models_loader
from elucidation.conveniences import funcTime


class PredictedSpectrumMatcher(object):
    """

    """

    __instances__ = {}  # singleton-ish: one per library.

    @staticmethod
    def get_instance(lib_id):
        """
        Used to get a PredictedSpectrumMatcher object for the library.  Enforces a singleton-per-library discipline.
        @param lib_id: library: 209, etc.
        @type lib_id: int
        @return: PredictedSpectrumMatcher object
        @rtype: PredictedSpectrumMatcher
        """
        if lib_id not in PredictedSpectrumMatcher.__instances__:
            PredictedSpectrumMatcher.__instances__[lib_id] = PredictedSpectrumMatcher(lib_id)
        return PredictedSpectrumMatcher.__instances__[lib_id]

    def __init__(self, lib_id):
        self.ma = ModelingAssistant.get_instance(lib_id)
        if lib_id in PredictedSpectrumMatcher.__instances__:
            raise Exception("try again with PredictedSpectrumMatcher.get_instance(lib_id)")
        self.lib_id = lib_id
        self.models = smiles2ms2_models_loader(lib_id)
        self.timestamp = self.models["timestamp"]
        self.gitbranch = self.models["git_branch"]
        self.models = self.models["models"]
        self.model_dict = {m["y"]: m for m in self.models}
        self.ma = ModelingAssistant.get_instance(lib_id)
        # self.df = ma.get_modeling_df()
        PredictedSpectrumMatcher.__instances__[lib_id] = self

    def match_spectrum_to_smiles_db(self, parent_mz, spectrum, smiles_df, model_bins):
        logging.debug(len(smiles_df))

        def group_matches(tab):
            t = (
                tab.smiles.max(),
                list(tab["name"]),  # dammit!
                list(tab.source),
                tab.formula.max(),
                tab.exact_mass.max(),
            )
            return pd.DataFrame([t], columns=qw("smiles name source formula mass"))

        smiles_df = smiles_df.groupby("smiles").apply(group_matches)

        matched_df = self.match_spectrum_to_smiles_df(
            smiles_df, parent_mz, spectrum, model_bins=model_bins
        )

        max_score = matched_df.score.max()
        logging.debug(f"max score {max_score}")
        matched_df = matched_df[(matched_df.score >= 0.5 * max_score) & (max_score > 0)]
        matched_df = matched_df.sort_values("score", ascending=False).head(50)
        matched_df = matched_df[qw("name source smiles formula score mass hits")]
        # matched_df.columns = qw("name source smiles formula score mass spectrum_match")
        logging.debug(matched_df)
        return matched_df.copy()

    def check_modeling(self, q_spec):
        """
        Looks at the ions in q_spec and determines which are modeled and which not.
        @param q_spec: a bin spectrum: dict(bin1234=0.05234, bin4567=1.0, ...)
        @type q_spec: dict
        @return: two sets: ions with/without models.
        @rtype: tuple
        """
        modeled_ions = set([m["y"] for m in self.models if m["y"] in q_spec.keys()])
        unmodeled_ions = set(q_spec.keys()) - modeled_ions

        return modeled_ions, unmodeled_ions

    def predict_ions_from_a_smiles(self, smiles, model_bins=None):
        """
        Given a SMILES and a list of models for ms2 ions, finds the ions predicted by the fragmentation of the SMILES.
        @param smiles: SMILES
        @type smiles: str
        @param model_bins: list of bins for which models must be checked.  Generally, we will not want to predict
        every possible ion for the SMILES; only the ones in the query spectrum.
        @type model_bins: list or set
        @return: set of bins with positive prediction
        @rtype: set
        """
        if model_bins == set():
            return set()
        if model_bins is None:
            models = self.models
        else:
            models = [self.model_dict[b] for b in model_bins]

        mol_m = smiles_mass(smiles)
        # Start by predicting anything that can result from breakage of a single bond,
        # just like the old days with MS2SMILES.
        fragments = set(fragment_smiles(smiles))
        formulas = set(map(fmla_for_smiles, fragments))
        single_bond_bins = reduce(add, [self.ma.fragment_bin_range(frag) for frag in fragments], [])
        single_bond_bins += reduce(
            add, [self.ma.lost_fragment_bin_range(frag, mol_m) for frag in fragments], []
        )

        # Now get bins matched by ML models
        # Ignore probabilities.
        ml_bins = lmap(itemgetter(0), self.ma.find_matching_models(models, formulas))
        old_ml_bins = lmap(itemgetter("y"), ml_bins)
        ml_bins = [m[0]["y"] for m in self.ma.find_matching_models(models, formulas)]
        assert old_ml_bins == ml_bins
        return set(ml_bins + single_bond_bins)

    def score_a_smiles_for_spectrum(self, smiles, spectrum, model_bins=None):
        """
        Assigns a match score for a SMILES and an MS2 spectrum based on intensities of matched ions.
        @param smiles: SMILES to match
        @type smiles: str
        @param spectrum: bin spectrum: dict(bin1234=0.234, bin3454=1.0, bin43423=0.43, ...)
        @type spectrum: dict
        @param model_bins: set or list of bins.  Could get this from spectrum.
        @type model_bins: set or list
        @return: match score
        @rtype: float
        """
        if model_bins == set():
            return 0
        smiles_ions = self.predict_ions_from_a_smiles(smiles, model_bins=model_bins)
        score = sum([spectrum[b] for b in smiles_ions if b in spectrum])
        return score

    def score_many_smiles_for_spectrum(self, smiles_list, spectrum, model_bins=None):
        """
        Find scores for a bunch of SMILES.  Doesn't return the actual matchings.
        @param smiles_list: list of SMILES strings
        @type smiles_list: list(str)
        @param spectrum: a bin spectrum dict(bin1234=0.234, bin3454=1.0, bin43423=0.43, ...)
        @type spectrum: dict
        @param model_bins: list or set of models applicable to spectrum.
        @type model_bins: list or set
        @return: dict(SMILES1=score1, ...)
        @rtype: dict
        """
        d = {
            s: self.score_a_smiles_for_spectrum(s, spectrum, model_bins=model_bins)
            for s in smiles_list
        }
        return d

    @funcTime
    def match_spectrum_to_smiles_df(self, smiles_df, parent_mz, spectrum, model_bins=None):
        """
        Given a dataframe with a smiles column, score the SMILES against the query spectrum,
        and fills in a column with scores and another with the bins matched.
        @param smiles_df: Dataframe with SMILES column
        @type smiles_df: DataFrame
        @param parent_mz: m/z of the precursor to the MS2 scan.
        @type parent_mz: float
        @param spectrum: raw spectrum for the query.
        @type spectrum: list of pairs
        @param model_bins: bins in the spectrum with models.
        @type model_bins: list/set
        @return: DataFrame; actually, the original dataframe is updated and returned.
        @rtype: DataFrame
        """
        matched_df = smiles_df
        q_spec, q_map = self.ma.bin_spectrum(spectrum, parent_mz)
        predictions = matched_df.smiles.apply(
            lambda smiles: self.predict_ions_from_a_smiles(smiles, model_bins)
        )
        q_keys = set(q_spec.keys())
        matched_df["spectrum_match"] = [q_keys & set(pred) for pred in predictions]
        return matched_df.sort_values("score", ascending=False).copy()

    def get_timestamp(self):
        """
        This will return the timestamp of the model file, for inclusion in reports.
        @return: something like 2021-03-04-1609
        @rtype: str
        """
        return self.timestamp

    def get_gitbranch(self):
        """
        This will return the git branch of the program that built the model file, for inclusion in reports.
        @return: something like master@wcejc23c
        @rtype: str
        """
        return self.gitbranch
