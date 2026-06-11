import numpy as np
from elucidation.conveniences import qw, lmap
from statsmodels.stats.weightstats import DescrStatsW  # need this for std err of mean
from elucidation.chemFacts import count_for_ratio
import weightedstats as ws  # need this one for median
import regex as re


class MS1Buddy:
    def __init__(self, ms1, label, dr_tuples):
        self.me = self
        self.label = label
        self.tally = len(dr_tuples)
        self.ms1 = ms1
        self.mz_feature = ms1.mz
        # delta, intensity ratio, buddy_mass, buddy intensity
        self.deltas, self.ratios, self.mz_buddies, self.int_buddies = zip(*dr_tuples)
        self.int_features = [
            i / r for i, r in zip(self.int_buddies, self.ratios)
        ]  # :( work backwards.
        weights = self.int_features
        self.delta = ws.weighted_median(self.deltas, weights=weights)
        self.ratio = ws.weighted_median(self.ratios, weights=weights)
        self.mz_buddy = ws.weighted_median(self.mz_buddies, weights=weights)
        self.int_buddy = ws.weighted_median(self.int_buddies, weights=weights)
        ratio_logs = lmap(np.log, self.ratios)
        w_stats = DescrStatsW(ratio_logs, weights=weights, ddof=0)
        w_mean = w_stats.mean
        if w_mean == 0:  # feature == feature
            self.me = None
            return
        self.c_o_v = w_stats.std / abs(w_mean)
        if self.c_o_v > 0.5:
            self.me = None
            return
        if label is not None:
            # print(label)
            self.feature, self.buddy = label.feature, label.buddy
            feature, buddy = label.feature, label.buddy
            self.theor_delta = label.theor_delta
            base_ions = ["m-H", "m+H", "x"]
            ratio = self.ratio
            # The base intensity must be bigger than the derivative's.
            if (buddy in base_ions and ratio < 0.99) or (
                feature in base_ions and ratio > 1.01
            ):  # The base intensity must be bigger than the derivative's.
                self.me = None
                return
            m = re.match(r"x\[([CSBON]\w?\d\d)-1\]$", buddy)  # not deuterium.
            if m and count_for_ratio(m.group(1), ratio) < 0.5:
                self.me = None
                return
            if (
                m
                and re.match(r"x\[(Cl|S|Br)\d\d-1\]", buddy)
                and count_for_ratio(m.group(1), ratio) > 2.5
            ):
                self.me = None
                return

            fmla = label.inferred_m
            # We don't really need sympy for this.
            self.inferred_m = eval(
                fmla.replace("d", str(self.delta)).replace("ff", str(self.mz_feature))
            )
        else:
            self.feature, self.buddy = ("", "")
            self.inferred_m = None
            self.theor_delta = np.nan

    def get_me(self):
        return self.me

    def get_summary(self):
        return {
            s: self.__getattribute__(s)
            for s in qw(
                "mz_buddy mz_feature theor_delta delta buddy feature inferred_m tally ratio c_o_v"
            )
        }
