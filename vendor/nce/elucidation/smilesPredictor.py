# coding: utf-8

import pandas as pd
from elucidation.modelingAssistant import ModelingAssistant
from elucidation.rdkitUtils import smiles_mass, fmla_for_smiles
from elucidation.tableLoaders import smiles_predictor_model_loader


class SmilesPredictor:
    """
    Class SmilesPredictor is for prediction SMILES fragments from MS2 spectra -- the original MS2SMILES.
    It relies on decision models created by train_ms2smiles_models.py.
    SmilesPredictors must be obtained with get_instance to enforce a singleton regime.
    """

    __instances__ = dict()

    @staticmethod
    def get_instance(lib_id):
        """
        This is a singleton (sort of) class; one object per library.
        This prevents re-reading table and model files.
        @param lib_id: yeah, a library id, right?
        @type lib_id: int
        @return: the singleton instance for this library; created if needed.
        @rtype: SmilesPredictor
        """
        if lib_id not in SmilesPredictor.__instances__:
            SmilesPredictor.__instances__[lib_id] = SmilesPredictor(lib_id)
        return SmilesPredictor.__instances__[lib_id]

    def __init__(self, lib_id):
        """
        Creates an instance.  Reads ML model files, as well as a file of extended text descriptions of categories.
        @param lib_id:
        @type lib_id: int
        """
        if lib_id in SmilesPredictor.__instances__:
            raise Exception("try again with SmilesPredictor.get_instance(lib_id)")
        self.lib_id = lib_id
        model_d = smiles_predictor_model_loader(lib_id)
        self.models = model_d["models"]
        self.timestamp = model_d["timestamp"]
        self.gitbranch = model_d["git_branch"]
        self.ma = ModelingAssistant.get_instance(lib_id)

    def predict_for_spectrum(self, parent_mass, spectrum):
        """
        Builds a dataframe of SMILES fragment predictions based on an MS2 spectrum.
        @param parent_mass: The mass of the precursor ion.
        @type parent_mass:
        @param spectrum: The "raw" (unbinned) spectrum to be matched" [(mz1,intensity1), (mz2, intensity2), ...]
        @type spectrum dict
        @return: dataframe of predictions
        @rtype: pandas.DataFrame
        """
        masses, intensities = zip(*spectrum)
        tspec = self.ma.range_spectrum(spectrum, parent_mass)
        matches = self.ma.find_matching_models(self.models, tspec)
        if not matches:
            return

        rows = [
            dict(
                subsmiles=model["y"],
                fmla=fmla_for_smiles(model["y"]),
                mass=smiles_mass(model["y"]),
                mass_frac=smiles_mass(model["y"]) / parent_mass,
                score=p * model["precision"],
            )
            for model, p in matches
        ]
        return pd.DataFrame(rows).sort_values("mass_frac", ascending=False)

    def get_timestamp(self):
        """
        This will return the timestamp of the model file, for inclusion in reports.
        @return: something like 2021-03-04-1609
        @rtype: str
        """
        return self.timestamp

    def get_gitbranch(self):
        """
        This will return the timestamp of the model file, for inclusion in reports.
        @return: something like 2021-03-04-1609
        @rtype: str
        """
        return self.gitbranch
