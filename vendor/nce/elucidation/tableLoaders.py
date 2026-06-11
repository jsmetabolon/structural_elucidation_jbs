import logging
import os

import dotenv
import pandas as pd
from elucidation.conveniences import pickle_load
from elucidation.config import BASE_DIR, REFERENCE_DIR

dotenv.load_dotenv()
from glob import glob

# Local directories and pickled files
system_base_dir = BASE_DIR
_system_reference_dir = REFERENCE_DIR


def get_config_var(key):
    return os.getenv(key, os.path.join(os.getcwd(), key))


def pickle_loader(s):
    """
    returns a function that will load the latest version of a reference pickle file.
    @param s: a string specifying the format of the filename
    @type s: str
    @return: function that returns contents of relevant pickle file.
    @rtype: func
    """

    def loader(*args):  # args could be things like library.
        needed = s % args
        pat = f"{os.path.join(_system_reference_dir, needed)}*.pkl"
        logging.debug(f"Trying to find {pat}")
        file = max(glob(pat))
        logging.info(
            f"loading {file} for {needed}"
        )  # TODO arrange to record this in elucidation report.
        return pickle_load(file)

    return loader


def csv_loader(s):
    """
    returns a function that will load the latest version of a reference CSV file.
    @param s: a string specifying the format of the filename
    @type s: str
    @return: function that returns contents of relevant CSV file.
    @rtype: func
    """

    def loader(*args):  # args could be things like library.
        needed = s % args
        pat = f"{os.path.join(_system_reference_dir, needed)}*.csv"
        logging.debug(f"Trying to find {pat}")
        file = max(glob(pat))
        logging.info(
            f"loading {file} for {needed}"
        )  # TODO arrange to record this in elucidation report.
        return pd.read_csv(file)

    return loader


fmla_finder_cache_loader = pickle_loader("fmla_finder_cache")
harmony_cache_loader = pickle_loader("harmony")
chemont_model_loader = pickle_loader("chemont_models_%d")
ff_dbe_loader = pickle_loader("ff_dbe")
library_fragments_loader = pickle_loader("library_fragments")
chemont_definition_loader = csv_loader("chemont_defs")
smiles_predictor_model_loader = pickle_loader("ms2smiles_models_%d")
smiles2ms2_models_loader = pickle_loader("smiles2ms2_models_%d")
ri_model_loader = pickle_loader("ri_models_%d")
