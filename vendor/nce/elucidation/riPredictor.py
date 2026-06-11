from argparse import Namespace
from bisect import bisect

import numpy as np
from elucidation.chemFacts import split_fmla, dbe
from elucidation.rdkitUtils import smiles_mass, fragment_smiles, fmla_mass, fmla_for_smiles, logp
from elucidation.tableLoaders import ri_model_loader

"""
This class is for predicting RI of a compound from its molecular formula or its SMILES string.
The class is initialized for the given library by reading pickled sklearn objects from a file.
The predictors return two values.
The first value is a predicted RI.
The second value is a confidence score of sorts. It is based on absolute residuals collected during the
cross-validation process.  The confidence score reflects the probability that the predicted RI is a likely estimate.
For example, if the prediction is 751 for a given SMILES string, and the observed RI of the feature being elucidate is 
760, this is obviously a very close match.  The score in this case is 1, meaning that fewer than 10% of estimates
collected during cross validation were this good.  If the observed RI is 948, the score (for library 209) is 5, meaning
than fewer than 50% of estimates from cross-validation were this good.
The observed RI (other_ri) is optional.  If absent, the methods predict RI and report a score of -1.
"""


class RiPredictor:

    __instances__ = {}

    @staticmethod
    def get_instance(lib_id):
        """
        enforces singleton-per-library discipline
        @param lib_id: library for which instance is required
        @type lib_id: int
        @return: RiPredictor instance for the library.
        @rtype: RiPredictor
        """
        if lib_id is None:  # No instances have been created, so load tables into class variables.
            lib_id = 0
        if lib_id not in RiPredictor.__instances__:
            RiPredictor.__instances__[lib_id] = RiPredictor(lib_id)

        return RiPredictor.__instances__[lib_id]

    def __init__(self, lib_id, test_dict=None):
        """
        Create an instance.  Following a singleton/library regime.
        @param lib_id: library
        @type lib_id: int
        @param test_dict: If this is provided, models are taken from here instead of loaded from file.
        This is done by the model-building program for testing.
        @type test_dict: dict
        """
        if lib_id in RiPredictor.__instances__:
            raise Exception(f"try again with RiPredictor.get_instance({lib_id})")
        self.lib_id = lib_id
        if test_dict is not None:
            d = test_dict
        else:  # usual case.
            d = ri_model_loader(lib_id)
        ns = Namespace(**d)
        self.smiles_vectorizer = ns.smiles_vectorizer
        self.smiles_regressor = ns.smiles_regressor
        self.smiles_scores = ns.smiles_scores + [1e6]
        self.formula_vectorizer = ns.formula_vectorizer
        self.formula_regressor = ns.formula_regressor
        self.formula_scores = ns.formula_scores + [1e6]
        RiPredictor.__instances__[lib_id] = self

    def predict_for_smiles(self, smiles, other_ri=None):
        """
        makes an RI prediction for a SMILES string
        @param smiles: SMIlES string
        @type smiles: str
        @param other_ri: measured RI to which the estimate is to be compared.
        @type other_ri: float
        @return: pair: (predicted RI, confidence)  see above.
        @rtype: tuple(int, int)
        """
        formula = fmla_for_smiles(smiles)
        spfmla = split_fmla(formula)[1]
        X_fmla = self.formula_vectorizer.transform(spfmla)
        frags = set([fmla_for_smiles(s) for s in fragment_smiles(smiles)])
        X_sfrags = self.smiles_vectorizer.transform([" ".join(frags)]).toarray()
        X_dbe = np.array([[dbe(formula)]])
        X_mass = np.array([[smiles_mass(smiles)]])
        X_logp = np.array([[logp(smiles)]])
        ri_pred = int(
            self.smiles_regressor.predict(
                np.concatenate([X_sfrags, X_fmla, X_dbe, X_mass, X_logp], axis=1)
            )[0]
        )
        if other_ri is not None:
            score = min([10, bisect(self.smiles_scores, abs(ri_pred - other_ri))])
            return ri_pred, score
        else:
            return ri_pred

    def predict_for_formula(self, formula, other_ri=None):
        """
        Makes an RI prediction for a formula.
        @param formula:  a molecular formula, e.g. C16H34
        @type formula: str
        @return: pair: (predicted RI, confidence)  see above.
        @rtype: tuple(int, int)
        """
        spfmla = split_fmla(formula)[1]
        X_fmla = self.formula_vectorizer.transform(spfmla)
        X_dbe = np.array([[dbe(formula)]])
        X_mass = np.array([[fmla_mass(formula)]])
        ri_pred = int(
            self.formula_regressor.predict(np.concatenate([X_fmla, X_dbe, X_mass], axis=1))[0]
        )
        if other_ri is not None:
            score = min([10, bisect(self.formula_scores, abs(ri_pred - other_ri))])
            return ri_pred, score
        else:
            return ri_pred
