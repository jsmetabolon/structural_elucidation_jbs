"""One-off cohort statistics script for study 21676 feature annotations.

Correlates integrated peak areas with clinical metadata and writes
annotations.csv. Retained vendor utility, not part of the spec2mol inference CLI.
"""

from collections import Counter
from glob import glob
import pandas as pd
import numpy as np
from elucidation.conveniences import pickle_load, qw
from scipy.stats import ttest_ind, pearsonr

meta_df = pd.read_csv(r"C:\CDtemp\21676\meta_data_21676.csv")
meta_d = dict()
for c in qw("gender age bmi creatinine cystatin abs_gfr rel_gfr egfr_mdrd"):
    meta_d[c] = dict(zip(meta_df.sample_name, meta_df[c]))

pickles = glob("C:/cdtemp/21676/209/clusters/F*.pkl")
cols = "abs_gfr age bmi creatinine cystatin rel_gfr egfr_mdrd"

rows = []
for pkl in pickles:
    f = pickle_load(pkl)
    d = dict(serial=f["serial"])
    peak_df = f["peak_df"]
    half_min = peak_df.area[peak_df.area > 1].min() / 2
    peak_df.area = np.maximum(peak_df.area, half_min)
    peak_df = peak_df.merge(meta_df, on="sample_name", how="left")

    boys = peak_df.area[peak_df.gender == "M"]
    girls = peak_df.area[peak_df.gender == "F"]
    d["p_gender"] = ttest_ind(boys, girls).pvalue

    for col in qw(cols):
        missing = np.isnan(peak_df[col])
        cor, p = pearsonr(peak_df.area[~missing], peak_df[col][~missing])
        d["p_" + col] = p
        d["R_" + col] = cor
    d["p_min"] = min([v for k, v in d.items() if k.startswith("p_")])
    rows.append(d)

cols = (
    qw("serial p_min p_gender ")
    + qw("R_" + cols.replace(" ", " R_"))
    + qw("p_" + cols.replace(" ", " p_"))
)
df = pd.DataFrame(rows)[cols]
df.to_csv("C:/cdtemp/21676/209/annotations.csv", index=False)
