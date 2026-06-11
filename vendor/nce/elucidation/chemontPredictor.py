# coding: utf-8

from operator import itemgetter

import pandas as pd
from elucidation.modelingAssistant import ModelingAssistant
from elucidation.tableLoaders import chemont_definition_loader, chemont_model_loader


class ChemontPredictor:
    """
    ChemontPredictor is a class for predicting Chemical Ontology categories such as "lipid" from MS2 spectra.
    It relies on decision models created by train_classyfile_models.py.
    ChemontPredictors must be obtained with get_instance to enforce a singleton regime.
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
        @rtype: ChemontPredictor
        """
        if lib_id not in ChemontPredictor.__instances__:
            _ = ChemontPredictor(lib_id)
        return ChemontPredictor.__instances__[lib_id]

    def __init__(self, lib_id):
        """
        Creates an instance.  Reads ML model files, as well as a file of extended text descriptions of categories.
        @param lib_id:
        @type lib_id: int
        """
        # Enforce singleton.
        if lib_id in ChemontPredictor.__instances__:
            raise Exception("try again with ChemontPredictor.get_instance(lib_id)")
        self.lib_id = lib_id
        model_d = chemont_model_loader(lib_id)
        self.models = sorted(model_d["models"], key=itemgetter("y"))
        self.timestamp = model_d["timestamp"]
        self.gitbranch = model_d["git_branch"]
        defs = chemont_definition_loader()
        defs.chemont_id = defs.chemont_id.apply(lambda id: "CHEMONTID:%07d" % id)
        self.chemont_def_dict = dict(zip(defs.chemont_id, defs["def"]))
        self.ma = ModelingAssistant.get_instance(lib_id)
        ChemontPredictor.__instances__[lib_id] = self

    def predict_for_spectrum(self, parent_mass, spectrum):
        """
        Builds a dataframe of Chemical Ontology prediction based on an MS2 spectrum.
        @param parent_mass: The mass of the precursor ion.
        @type parent_mass:
        @param spectrum: The "raw" (unbinned) spectrum to be matched.
        @type spectrum dict
        @return: dataframe of predictions
        @rtype: pandas.DataFrame
        """
        tspec = self.ma.range_spectrum(spectrum, parent_mass)
        matches = self.ma.find_matching_models(self.models, tspec)
        if not matches:
            return
        rows = [
            dict(
                chemont_id=model["y"],
                chemont_name=model["name"],
                chemont_def=self.chemont_def_dict[model["y"]],
                score=p * model["precision"],
            )
            for model, p in matches
        ]
        return pd.DataFrame(rows).sort_values("score", ascending=False)

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
