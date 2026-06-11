# coding: utf-8

"""
This program trains models to predict Chemical Ontology categories from MS2 spectra.
The results are used by ChemontPredictor.
"""
import argparse
import datetime
import sys
from collections import Counter
from functools import reduce
from operator import add

import numpy as np
import pandas as pd
from elucidation.conveniences import qw, pickle_load, lmap
from elucidation.modelingAssistant import ModelingAssistant
from elucidation.conveniences import resource_string, pickle_dump, now_string, qw, get_git_branch
from mdsutils.easyQuery import EasyQuery
from sklearn import linear_model, metrics

exec("import elucidation.initPackageOptions")  # won't be optimized away by pycharm


def elapsed_time(comment):
    s = resource_string()
    print(f"{s}: {comment}")


def getArguments(argv):
    parser = argparse.ArgumentParser(description="""MS2SMILES training.""")
    parser.add_argument("-l", "--library", type=int, help="library")
    return parser.parse_args(argv[1:])


def build_chemont_dict(lib_id):
    """
    Builds dictionaries of known annotations of compounds in MLIMS.
    @param lib_id:
    @type lib_id:
    @return: two dictionaries.
    The first maps MLIMS comp ids to lists of annotated categories as downloaded by Edi a long time ago.
    The second maps ChemOnt ids to ChemOnt category names.
    @rtype: dict, dict
    """
    d = pickle_load("chemont_prediction_training.pkl")

    classyfire_df = d["classyfire"]
    chemont_name_dict = dict(zip(classyfire_df.chemont_id, classyfire_df.chemont_name))

    cle_df = d["annotated_chro_lib_entry"]
    cle_df = cle_df[cle_df.lib_id == lib_id]
    # Throw out categories that are too rare or too frequent.
    c = Counter(reduce(add, list(cle_df.chemont)))
    good_categories = set([id for id, n in c.items() if (5 <= n <= 1000)])
    cle_df.chemont = cle_df.chemont.apply(lambda lis: list(set(lis) & good_categories))
    # Everything in cle_df at start had some annotation, so we assume these were annotated by the Wishart Lab.
    # The rest of MLIMS we've assumed was not even considered by the Wishart Lab (though maybe this is wrong.)
    # Anything left with no annotations now is still good as a negative example for every category remaining.

    good_categories = sorted(good_categories)
    category_to_column = {cat: i for i, cat in enumerate(good_categories)}
    cle_df["cat_inx"] = cle_df.chemont.apply(lambda lis: lmap(category_to_column.get, lis))
    return cle_df, category_to_column, chemont_name_dict


if __name__ == "__main__":
    today = datetime.date.today()
    today = "%4d-%02d-%02d" % (today.year, today.month, today.day)
    time_str = now_string()
    git_branch = get_git_branch()

    # We will read these from the command line.
    opt = getArguments(sys.argv)
    lib_id = opt.library

    elapsed_time("starting")
    ma = ModelingAssistant.get_instance(lib_id)

    d = pickle_load("chemont_prediction_training.pkl")
    cle_df, category_to_column, chemont_name_dict = build_chemont_dict(lib_id)
    # Get spectra from library, but only for annotated compounds.
    X, X_d, comp_ids = ma.get_matrix_for_bin_spectrum(comp_ids=cle_df.comp_id)
    cle_df = cle_df[
        cle_df.comp_id.isin(comp_ids)
    ]  # maybe some annotated comp_ids didn't have spectra... it happened.
    # Now we must realign X to match cle_df.
    vec_dict = dict(zip(comp_ids, list(X)))
    X = cle_df.comp_id.apply(vec_dict.get)
    X = np.array(lmap(list, list(X)))
    print(X)

    print(X.shape, "training matrix")
    Y = np.zeros((len(cle_df), len(category_to_column)))
    for i, t in enumerate(cle_df.itertuples()):
        for j in t.cat_inx:
            Y[i, j] = 1

    elapsed_time(f"training matrices built; X:{X.shape}; Y:{Y.shape}")
    models = []
    for trace, (chemont_id, inx) in enumerate(category_to_column.items()):
        y = Y[:, inx]
        count = int(y.sum())
        if count < 5:  # Have to screen for this again, apparently.
            continue
        logreg = linear_model.LogisticRegression(C=1, solver="liblinear")
        logreg.fit(X, y)
        logreg_ypred = logreg.predict(X)
        recall = metrics.recall_score(y, logreg_ypred)
        precision = metrics.precision_score(y, logreg_ypred)
        count = int(y.sum())
        # Not restricting this right now.
        if precision > 0.8 and recall > 0.1:
            models.append(
                dict(
                    lib_id=lib_id,
                    y=chemont_id,
                    name=chemont_name_dict[chemont_id],
                    x_cols=X_d,
                    precision=precision,
                    recall=recall,
                    count=count,
                    model=logreg,
                    timestamp=time_str,
                    git_branch=git_branch,
                )
            )
        elapsed_time(
            f"{trace + 1}/{len(category_to_column)} {chemont_id} recall={recall:.3f} "
            + f"precision={precision:.3f} count={count} {chemont_name_dict[chemont_id]}    "
        )

    out_df = pd.DataFrame(models)
    out_df = out_df.drop("model", axis=1)
    out_df.to_csv(f"performance/chemont_performance_{lib_id}_{time_str}.csv", index=True)
    print(f"{len(models)} successful models")

    pkl = f"c:/cdtemp/references/chemont_models_{lib_id}_{time_str}.pkl"
    print("writing", pkl)
    pickle_dump(
        pkl, dict(models=models, timestamp=time_str, git_branch=git_branch),
    )
