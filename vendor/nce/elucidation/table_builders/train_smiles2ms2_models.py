# coding: utf-8
import argparse
import datetime
import sys

import numpy as np
import pandas as pd
from elucidation.modelingAssistant import ModelingAssistant
from elucidation.conveniences import resource_string, pickle_dump, now_string, get_git_branch

exec("import elucidation.initPackageOptions")  # won't be optimized away by pycharm
# from sklearn.model_selection import train_test_split
from sklearn import linear_model, metrics


def elapsed_time(comment):
    s = resource_string()
    print(f"{s}: {comment}")


def getArguments(argv):
    parser = argparse.ArgumentParser(description="""train smile2ms2 to predict MS2 from SMILES.""")
    parser.add_argument("-l", "--library", type=int, help="library")
    return parser.parse_args(argv[1:])


if __name__ == "__main__":
    today = datetime.date.today()
    today = "%4d-%02d-%02d" % (today.year, today.month, today.day)
    time_str = now_string()
    git_branch = get_git_branch()

    # We will read these from the command line.
    opt = getArguments(sys.argv)
    lib_id = opt.library

    elapsed_time("starting")

    # Get spectra,etc., from library and massage into dictionaries.
    ma = ModelingAssistant.get_instance(lib_id)
    X0, X0_d, _ = ma.get_matrix_for_fragments()
    X1, X1_d, _ = ma.get_matrix_for_frag_formulas()
    print(len(list(X0)), len(list(X1)))
    X = np.array([l1 + l2 for l1, l2 in zip(map(list, list(X0)), map(list, list(X1)))])

    print(X0.shape, X1.shape, X.shape)
    X_d = X0_d.copy()
    X_d.update({k: v + len(X0_d) for k, v in X1_d.items()})
    Y, Y_d, _ = ma.get_matrix_for_range_spectrum()

    # Try using only formulas.
    X, X_d = X1, X1_d

    elapsed_time("training matrices built")

    print(X.shape, Y.shape)

    models = []
    for i, c in enumerate(Y_d.keys()):
        y = Y[:, Y_d[c]]
        logreg = linear_model.LogisticRegression(C=1, solver="liblinear")
        logreg.fit(X, y)
        logreg_ypred = logreg.predict(X)
        recall = metrics.recall_score(y, logreg_ypred)
        precision = metrics.precision_score(y, logreg_ypred)
        count = int(((y + 1).sum()) / 2)
        if precision > 0.5 and recall > 0.1:
            models.append(
                dict(
                    lib_id=lib_id,
                    y=c,
                    x_cols=X_d,
                    precision=precision,
                    recall=recall,
                    count=count,
                    model=logreg,
                    git_branch=git_branch,
                    time_stamp=time_str,
                )
            )
        elapsed_time(
            f"{i + 1}/{len(Y_d)} {c} recall={recall:.3f} "
            + f"precision={precision:.3f} count={count}       "
        )

    out_df = pd.DataFrame(models)
    out_df = out_df.drop("model", axis=1)
    out_df.to_csv(f"performance/smiles2ms2_performance_{lib_id}_{time_str}.csv", index=True)
    print(f"{len(models)} successful models")

    pkl = f"c:/cdtemp/references/smiles2ms2_models_{lib_id}_{time_str}.pkl"
    print("writing", pkl)
    pickle_dump(
        pkl, dict(models=models, timestamp=time_str, git_branch=git_branch,),
    )
