"""Post-elucidation clustering of features by inferred neutral mass and RI.

Groups related features across a batch summary and emits HTML cluster reports.
Called from elucidator.py after all features in a run are processed.
"""

from functools import reduce
from operator import add, itemgetter

import numpy as np
import regex as re
from discovery.precursorClusterer import PrecursorClusterer
from elucidation.binning import minus_5ppm, plus_5ppm
from elucidation.htmlUtils import FileReporter, boxit
from elucidation.libraryService import LibraryService
from elucidation.reformatter import reformat
from elucidation.conveniences import *
from elucidation.tableLoaders import get_config_var
from elucidation.unionFind import UnionFind
from elucidation.config import BASE_DIR

local_base_dir = BASE_DIR


class Ms1ClusterAnalyzer:
    """Cluster elucidated features that share inferred masses within an RI window."""

    def __init__(self, summary_df, debug=1e9):
        df = summary_df[summary_df["inferred_masses"].notnull()]
        quads = lzip(df.serial, df.ri, df.mz, df.inferred_masses)
        f_to_masses = {extract_feature_name(f): extract_masses(t) for f, ri, mz, t in quads}
        f_to_ris = {extract_feature_name(f): ri for f, ri, mz, t in quads}
        triples = [[(f_to_ris[f], m, f) for m in f_to_masses[f]] for f in f_to_ris.keys()]
        triples = reduce(add, triples)
        # Cluster, and retrieve just the feature names from the clusters.
        pc = PrecursorClusterer(triples, MAX_RI_SPAN=10)
        clusters = [set(map(itemgetter(2), g)) for g in pc.clusters]
        #  Throw out singleton clusters sooner rather than later.
        clusters = [cl for cl in clusters if len(cl) > 1]
        # At this point, we expect features to be in more than one cluster, because clusters are tied to a
        # single inferred mass.  So we have to merge clusters with common features.
        uf = UnionFind()
        for cl in clusters:
            uf.union(cl)
        self.clusters = uf.get_clusters()
        return


def extract_masses(t):
    return set([np.fabs(float(s)) for s in t.split(", ")])


def extract_feature_name(s):
    m = re.search(r"([FGHJ]\d\d\d\d)", s)
    return m.group(1)


def filter_buddy_df(f):
    df = f["buddy_df"]
    df = df[df.feature != ""]
    df["serial"] = f["serial"]
    df = df[qw("mz_buddy mz_feature buddy feature inferred_m serial")]
    return df


def split_group_by_mass(group):  # list of (ri,mz) pairs
    if not group:  # shouldn't happen...
        return []
    group = sorted(group)
    start = 0
    groups = []
    gap = group[0] * 10e-6
    for i in range(1, len(group)):
        if group[i] > group[i - 1] + gap:
            groups.append(group[start:i])
            start = i
            gap = group[i] * 10e-6

    groups.append(group[start:])
    return groups


def assign_masses_to_groups(masses):
    masses = sorted(masses)
    group_num = 0
    groups = [[masses[0]]]
    d = {masses[0]: 0}
    gap = masses[0] * 10e-6
    for i in range(1, len(masses)):
        if masses[i] > masses[i - 1] + gap:
            group_num += 1
            gap = masses[i] * 10e-6
            groups.append([])
        d[masses[i]] = group_num
        groups[-1].append(masses[i])
    medians = lmap(np.median, groups)
    d = {m: medians[g] for m, g in d.items()}
    return d


def find_serial(mz, mz_and_serial):
    lis = [(mz1, serial) for mz1, serial in mz_and_serial if minus_5ppm(mz) <= mz1 <= plus_5ppm(mz)]
    return lis[0][1] if lis else ""


def foo(tab, mz_and_serial):
    d = dict(lzip(tab.mz_buddy, tab.buddy) + lzip(tab.mz_feature, tab.feature))
    groups = split_group_by_mass(d.keys())
    med_inf_m = np.median(tab.inferred_m)
    rows = []
    for g in groups:
        med_mz = np.median(g)
        rows.append(
            {
                "neutral_mass": med_inf_m,
                "ion_mz": med_mz,
                "serial": find_serial(med_mz, mz_and_serial),
                "names": sorted(set([d[m] for m in g])),
            }
        )
    df = pd.DataFrame(rows)[qw("neutral_mass ion_mz serial names")]
    return df


def bar(fd):
    buddy_df = fd["buddy_df"]
    df = pd.DataFrame(
        dict(
            mz=list(buddy_df.mz_buddy) + [fd["mz"]],
            ratio=list(buddy_df.ratio) + [1],
            feature=([fd["serial"]] * (1 + len(buddy_df))),
        )
    )
    return df


class Ms1Cluster:
    def __init__(self, serial, feature_dicts):
        feat0 = feature_dicts[0]
        self.feature_dicts = feature_dicts
        self.serial = serial
        self.ri = np.median([f["ri"] for f in feature_dicts])
        self.lib_id = feat0["lib_id"]
        self.feature_names = sorted([f["serial"] for f in feature_dicts])
        self.internal_standard = feat0["internal_standard"]
        self.me = self
        self.study_id = 21676  # TODO feat0["study_id"]

        buddy_dfs = [filter_buddy_df(f) for f in feature_dicts]
        buddy_df = pd.concat(buddy_dfs, ignore_index=True)
        # Correct all masses to median of similar masses.
        all_masses = list(buddy_df.mz_buddy) + list(buddy_df.mz_feature) + list(buddy_df.inferred_m)
        if not all_masses:
            self.me = None
            return
        d = assign_masses_to_groups(all_masses)
        for c in qw("mz_buddy mz_feature inferred_m"):
            buddy_df[c] = buddy_df[c].apply(d.get).round(4)
        mz_and_serial = [(f["mz"], f["serial"]) for f in feature_dicts]  # should make this go away!
        buddy_df = (buddy_df.groupby("inferred_m").apply(foo, mz_and_serial)).reset_index(drop=True)

        buddy_df.names = buddy_df.names.apply(lambda lis: [s.replace("x", "m-H") for s in lis])
        pvt = buddy_df.pivot(index="ion_mz", columns="neutral_mass", values="names")
        mass_cols = pvt.columns
        pvt = pvt.reset_index()
        pvt["serial"] = pvt.ion_mz.apply(lambda mz: find_serial(mz, mz_and_serial))
        for col in pvt.columns:
            pvt[col] = ["" if type(x) is float and np.isnan(x) else x for x in pvt[col]]
        self.label_df = pvt[["serial", "ion_mz"] + list(mass_cols)].reset_index(drop=True)

        area_df = pd.DataFrame({ft["serial"]: ft["peak_df"].area for ft in self.feature_dicts})[
            self.feature_names
        ]
        self.area_corr_df = area_df.corr()
        self.area_stats_df = area_df.describe().apply(np.int32)

        buddy_dfs = [bar(f) for f in feature_dicts]
        buddy_df = pd.concat(buddy_dfs, ignore_index=True)
        d = assign_masses_to_groups(buddy_df.mz)
        buddy_df.mz = buddy_df.mz.apply(d.get)
        buddy_df = buddy_df.drop_duplicates()
        pvt = buddy_df.pivot_table(index="mz", columns="feature", values="ratio")
        pvt = pvt.reset_index()

        pvt.index.name = None
        self.merged_spectra = pvt
        return

    def get_me(self):
        return self.me


class Ms1ClusterReporter:
    def __init__(self, c):
        self.c = c
        self.rep = None
        self.reminder = None
        self.report_ri = None  # This should actually be in the feature.
        self.report_mass = None  # ditto, should be in feature.
        self.summary = None

    def img(self, _):
        return FileReporter.img(_)

    def write_report(self, opts):
        """
        The top level function for generating a report.
        It opens an html file for output and generates various elements from dataframes stored in the Feature object.
        :param opts: stuff off the command line.
        :return: a dict with summary statistics for the feature; it will eventually be included in the spreadsheet.
        """
        c = self.c
        report_file_name = f"ms1cluster_{opts.tag}_ri{c.ri}_{c.serial}.html"
        report_file = f"{opts.feature_dir}/{report_file_name}"
        self.report_file = report_file
        # rank_dict, breaks_list = get_sample_plot_sort(eq,task_df)
        summary = {
            "ri": c.ri,
            "internal_standard": c.internal_standard,
            "report": report_file,
            "num_features": len(c.feature_names),
            "features": c.feature_names,
            "serial": c.serial,
        }
        c.summary = summary  # attach it to the feature.
        self.summary = summary  # and to this object.
        # Add sample peak areas to summary.

        # Open report file.
        rep = FileReporter(report_file)
        self.rep = rep
        if True:
            # Add title block to top of report.
            rep.htmls(self.build_title())
            rep.start_accordion("Label Table")
            rep.comment("<h3>Column = Neutral Mass</h3>")
            df = c.label_df
            for col in df.columns:
                df[col] = df[col].apply(lambda v: ", ".join(v) if type(v) is list else v)
            rep.htmls(reformat(df))
            rep.end_accordion()
            rep.start_accordion("Area Correlations and Distributions")
            rep.htmls(
                reformat(c.area_corr_df, index=True),
                "&emsp;&emsp;",
                reformat(c.area_stats_df, index=True),
            )
            rep.end_accordion()
            rep.start_accordion("Feature MS1s")
            pvt = c.merged_spectra.copy()
            for col in pvt.columns:
                pvt[col] = ["" if type(x) is float and np.isnan(x) else x for x in pvt[col]]
            rep.htmls(reformat(pvt, index="ion_mz"), rep.img(self.plot_ms1s()))
            rep.end_accordion()

        # except Exception as e:
        #     print(e)
        #     rep.close_report()
        #     raise
        rep.close_report()
        return summary

    def build_title(self):
        c = self.c
        feature_str = ", ".join(sorted(c.feature_names))
        study_title = LibraryService.get_instance().get_study_title(c.study_id)
        template = """
         <table>
         <tr> <td colspan=4> %s </td> </tr>
         <tr> <td> %s </td> <td> %s </td> <td> %s </td> <td> %s </td></tr>
         <tr> <td colspan=3> %s </td> <td> %s </td> </tr>
         </table>
         """ % (
            boxit(study_title, tag="h3"),
            boxit(f"Library<br>{c.lib_id}", tag="h3"),
            boxit(f"MS1 Cluster<br>{c.serial}", tag="h3"),
            boxit(f"RI<br>{c.ri}", tag="h3"),
            boxit(time.ctime() + "<br>Version " + get_git_branch(), tag="h3"),
            boxit(f"Features: {feature_str}", tag="h3"),
            boxit(f"Nearby Internal Standard: {c.internal_standard}"),
        )
        return template

    def plot_ms1s(self):
        rep = self.rep
        df = self.c.merged_spectra
        cols = list(df.columns)[1:]
        n = len(cols)
        x = df.mz
        fig = plt.figure(figsize=(10, n * 3))
        axes = fig.subplots(n, 1, sharex=True)
        for ax, col in zip(axes, cols):
            ax.set_ylabel(col)
            y = list(df[col])
            for m, i in zip(x, y):
                ax.plot([m, m], [0, i], lw=2, color="r")
            ax.plot(x, y, ".b")
