#!/usr/bin/env python
# coding: utf-8
"""Train per-library retention-index RandomForest models from library metadata.

Fits formula- and fragment-based RI regressors and writes ri_models_{lib}.pkl
loaded by RiPredictor during feature elucidation.
"""

import argparse
import sys
from collections import Counter
from functools import reduce
from operator import add, itemgetter

import numpy as np
import pandas as pd
from elucidation.chemFacts import dbe, split_fmla
from elucidation.modelingAssistant import ModelingAssistant
from elucidation.rdkitUtils import logp, fmla_for_smiles, smiles_mass
from elucidation.conveniences import pickle_dump, get_git_branch, now_string
from elucidation.riPredictor import RiPredictor
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_extraction import DictVectorizer
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.model_selection import KFold


def get_arguments(argv):
    """Parse CLI options for RI model training."""
    parser = argparse.ArgumentParser(description="""RI prediction training.""")
    parser.add_argument(
        "-l",
        "--library",
        choices=[209, 400, 402, 305],
        help="library, only 209, 305, 400, 402 are allowed",
        type=int,
        required=True,
    )
    return parser.parse_args(argv[1:])


def build_subformula_matrix(lib_df):
    lib_df["frag_fmlas"] = lib_df.fragments.apply(
        lambda L: list(set([fmla_for_smiles(s) for s in L]))
    )
    all_frags = Counter(reduce(add, lib_df.frag_fmlas))
    pct1 = 0.01 * len(lib_df)
    all_frags = set([k for k, v in all_frags.items() if v >= pct1 and v <= len(lib_df) - pct1])
    lib_df.frag_fmlas = lib_df.frag_fmlas.apply(lambda L: [s for s in L if s in all_frags])
    lib_df["frag_string"] = lib_df.frag_fmlas.apply(" ".join)
    cv = CountVectorizer(lowercase=False, token_pattern=r"(?u)\b\w+\b")
    fit = cv.fit_transform(list(lib_df.frag_string))
    X_sfrags = fit.toarray()
    return cv, X_sfrags


def build_formula_matrix(lib_df):
    split_fmlas = lib_df.formula.apply(split_fmla).apply(itemgetter(1))
    dv = DictVectorizer(sparse=False)
    return dv, dv.fit_transform(split_fmlas)


def test_a_model(X_list, ri, splits):
    reg = RandomForestRegressor(n_estimators=100)
    X = np.concatenate(X_list, axis=1)
    dfs = []
    for trainx, testx in splits:
        print(tag, end=" ")
        reg.fit(X[trainx], ri[trainx])
        ri_pred = reg.predict(X[testx])
        abs_resid = np.abs(ri[testx] - ri_pred)
        dfs.append(pd.DataFrame(dict(ri_pred=ri_pred, abs_resid=abs_resid, ri=ri[testx])))
    print()
    return pd.concat(dfs, ignore_index=True)


if __name__ == "__main__":
    lib_id = get_arguments(sys.argv).library
    git_branch = (get_git_branch(),)
    time_str = now_string()
    ma = ModelingAssistant(lib_id)
    lib_df = ma.get_modeling_df()
    lib_df = lib_df[lib_df.smiles != ""]

    # Build training matrix for RI, mass, dbe, and logp
    ri = np.array(lib_df.ri)
    X_dbe = np.array(lib_df.formula.apply(dbe)).reshape((-1, 1))
    X_mass = np.array(lib_df.smiles.apply(smiles_mass)).reshape((-1, 1))
    X_logp = np.array(lib_df.smiles.apply(logp)).reshape((-1, 1))

    # Build subformula training matrix
    smiles_vectorizer, X_sfrags = build_subformula_matrix(lib_df)

    # ## Build formula matrix
    formula_vectorizer, X_fmla = build_formula_matrix(lib_df)

    # Perform cross-validation and return residuals.
    splits = list(KFold(n_splits=5, shuffle=True, random_state=13).split(ri))
    results = {}
    for (tag, X_list) in [
        ("fdm", [X_fmla, X_dbe, X_mass]),
        ("sfdml", [X_sfrags, X_fmla, X_dbe, X_mass, X_logp]),
    ]:
        results[tag] = test_a_model(X_list, ri, splits)

    d = {"ri": list(results.values())[0].ri}
    for tag, df in results.items():
        d[f"abs_resid_{tag}"] = df.abs_resid
    result_df = pd.DataFrame(d)
    scores_df = result_df.describe([0.1 * i for i in range(10)])
    print(scores_df)
    smiles_scores = list(scores_df["abs_resid_sfdml"])[4:]
    formula_scores = list(scores_df["abs_resid_fdm"])[4:]

    # Build final SMILES prediction model using all data.
    smiles_regressor = RandomForestRegressor(n_estimators=100)
    smiles_regressor.fit(np.concatenate([X_sfrags, X_fmla, X_dbe, X_mass, X_logp], axis=1), ri)

    # Build final formula prediction model using all data.
    formula_regressor = RandomForestRegressor(n_estimators=100)
    formula_regressor.fit(np.concatenate([X_fmla, X_dbe, X_mass], axis=1), ri)

    d = dict(
        smiles_vectorizer=smiles_vectorizer,
        smiles_regressor=smiles_regressor,
        smiles_scores=smiles_scores,
        formula_vectorizer=formula_vectorizer,
        formula_regressor=formula_regressor,
        formula_scores=formula_scores,
        time_str=time_str,
        git_branch=git_branch,
    )

    rip = RiPredictor(lib_id, test_dict=d)

    smiles = lib_df.smiles[0]
    fmla = fmla_for_smiles(smiles)
    ri = lib_df.ri[0]
    print(f"RI {ri}")
    print(smiles, rip.predict_for_smiles(smiles))
    print(smiles, rip.predict_for_smiles(smiles, other_ri=ri))
    print(fmla, rip.predict_for_formula(fmla))
    print(fmla, rip.predict_for_formula(fmla, other_ri=ri))

    pkl = f"c:/cdtemp/references/ri_models_{lib_id}_{time_str}.pkl"
    print("writing", pkl)
    pickle_dump(pkl, d)
