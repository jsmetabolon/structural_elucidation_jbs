import logging
from argparse import Namespace
from functools import reduce
from operator import add

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import regex as re
import seaborn as sns
from elucidation.binning import bin_base
from elucidation.dfExplainer import build_explainer as build_explainer0
from elucidation.feature import Feature
from elucidation.htmlUtils import boxit, FileReporter, wrapit
from elucidation.libraryService import LibraryService
from elucidation.reformatter import reformat as reformat0
from elucidation.conveniences import (
    plt_reset,
    qw,
    is_empty_df,
    leng,
    pickle_load,
)

"""
This module generates a report for a feature.  
This is factored out from Feature like this so that it is easy to generate reports in different formats.
"""

"""
This is a list of columns that could be interesting in a debugging report, but not for elucidator.
"""
columns_to_drop = qw(
    "bucket cand_fmla neutral_mass lib_fraction supporting_pairs total_prob c_o_v predict details"
) + qw("recall precision prob theor_delta mass_range")


def reformat(df_in, **kwargs):
    """
    This calls Reformatter.reformat after removing the non-user columns.
    @param df_in:
    @type df_in:
    @param kwargs:
    @type kwargs:
    @return:
    @rtype:
    """
    if df_in is None:
        return reformat0(df_in, **kwargs)
    # return reformat0(df_in, **kwargs)
    df = df_in.drop(columns_to_drop, axis=1, errors="ignore")
    return reformat0(df, **kwargs)


def build_explainer(df_in, **kwargs):
    """
    This calls its counterpart in dfExplainer.py after removing the non-user columns.
    The "explaining" is the table of column explanations, generally to the right of the table.
    @param df_in:
    @type df_in:
    @param kwargs:
    @type kwargs:
    @return:
    @rtype:
    """
    if df_in is None:
        return build_explainer0(df_in, **kwargs)
    # return reformat0(df_in, **kwargs)
    df = df_in.drop(columns_to_drop, axis=1, errors="ignore")
    return build_explainer0(df, **kwargs)


class UserFeatureReporter:
    def __init__(self, feature):
        """
        UserFeatureReporter generates a freestanding HTML report on a single feature.
        @param feature: We need a dict, but if it's a Feature instance, we can extract the dict.
        @type feature: dict or Feature
        """
        # We want a dict, one way or another.
        if type(feature) is dict:
            pass
        if type(feature) is Feature:
            feature = feature.as_dict()
        if type(feature) is str:
            feature = pickle_load(feature)
        self.f = Namespace(**feature)
        self.rep = None
        self.reminder = None
        self.summary = None

    def img(self, _):
        """
        Puts the current matplotlib plot into the report.  It's kind of funky because the argument is a dummy.
        The plot is really being passed in the global in pyplot.  This should be sorted out someday.  #TODO
        @param _: dummy
        @type _:
        @return: None; causes HTML to be emitted.
        @rtype: Nonetype
        """
        return FileReporter.img(_)

    def write_report(self, opts):
        """
        The top level function for generating a report.
        It opens an html file for output and generates various elements from dataframes stored in the Feature object.
        :param opts: stuff off the command line.
        :return: a dict with summary statistics for the feature; it will eventually be included in the spreadsheet.
        """
        f = self.f
        report_mass = "%.4f" % f.mz
        report_ri = int(f.ri)
        file_serial = re.search(r"^[a-zA-Z0-9-]+", f.serial).group(0)
        report_file_name = f"feature_{opts.tag}_mz{report_mass}_ri{report_ri}_{file_serial}.html"
        report_file = f"{opts.feature_dir}/{report_file_name}"
        # summary will be the value eventually returned for inclusion in the summary spreadsheet.
        # More stuff will be gathered into this later.
        summary = {
            "mz": report_mass,
            "ri": report_ri,
            "internal_standard": f.internal_standard,
            "report": report_file,
            "num_peaks": len(f.peak_df[f.peak_df.area > 1]),
            "num_spectra": sum(f.ms2_counts.values()),
            "num_library_matches": leng(f.library_report),
            "cand_fmlas": leng(f.fmla_consistency_report),
            "max_harmony": f.max_harmony,
            "serial": f.serial,
            # "time": f.elucidation_time,
            "wins_ms1": "*" if len(f.bright_features) == 0 else "",
        }
        self.ion_tracker = f.discoverer_version == "IonTracker"
        ion_tracker = self.ion_tracker

        f.summary = summary  # attach it to the feature.
        self.summary = summary  # and to this object.
        # Add sample peak areas to summary.
        summary["areas"] = dict(zip(f.peak_df.sample_name, f.peak_df.area))
        summary["areas"]["serial"] = f.serial

        # Open report file.
        rep = FileReporter(report_file)
        self.rep = rep

        # Add title block to top of report.
        self.title_block()
        if ion_tracker:
            rep.htmls(wrapit(f"<h3>Abbreviated Report for Comp ID {f.comp_id} / {f.serial}</h3>"))
        if f.comments:
            rep.comment("<h2>Comments</h2>")
            for c in f.comments:
                rep.comment(f"<h3>{c}</h3><br>")
        if opts.annot_df is not None:
            df = opts.annot_df
            if "serial" in set(df.columns):
                df = df[df.serial == f.serial].drop("serial", axis=1)
            elif "fxxxx" in set(df.columns):
                df = df[df.fxxxx == f.serial].drop("fxxxx", axis=1)
            if leng(df) > 0:
                summary.update({c: df[c].min() for c in df.columns})
                rep.comment("<h2>Annotations</h2>")
                rep.htmls("&emsp;&emsp;", reformat(df.T, index=True))

        if f.consensus_df is None:
            rep.htmls("<h1>Failed to build MS2 consensus spectrum.</h1>", "<b>Done</b>")
            logging.warning(f"Failed to build MS2 consensus spectrum for {f.serial}")
            rep.close_report()
            return {}

        # When we are debugging, we would like for the program to crash if there's a problem.
        # In production, we will want to do a try-except to move on gracefully to the next feature.
        try:
            # The reminder is the little block that reappears repeatedly on the left margin.
            self.reminder = (
                f"<br>[{f.serial} RI={f.ri} m/z={f.mz:.4f}]</p>"
                # f'<p style="font-size:18px">&emsp;&emsp;[{f.serial} RI={f.ri} m/z={f.mz:.4f}]</p>'
            )
            reminder = self.reminder

            if not self.ion_tracker:
                self.report_on_rediscovery()
                summary["rediscovery_of"] = f.rediscovery_of

            self.report_on_spectral_predictions()

            if not ion_tracker:
                self.peak_statistics()
                self.integrated_peaks()
                self.peak_plots()
            self.report_on_ms1()
            self.mass_library_report()
            self.candidate_formula_consistency_report()
            self.ms2_spectrum_report()
            if not self.ion_tracker:
                self.ms2_ion_consistency_report()
            self.report_on_smiles_predictions()
            self.report_on_chemont_predictions()
            self.report_on_structure_assignments()
            self.report_on_ms2_spectrum_matching()

        except Exception as e:
            logging.warning(f"{self.f.serial}: {str(e)}")
            rep.close_report()
            rep = None
            raise
        finally:
            rep.close_report()
            rep = None
            return summary

    def title_block(self):
        """
        Creates and emits HTML for the title block at the beginning of the report:
        m/z, RI, number of peaks, etc.
        @return: None
        @rtype: NoneType
        """
        f = self.f
        peak_df = f.peak_df
        n_peaks = len(peak_df[peak_df.area > 1])
        n_spectra = self.summary["num_spectra"]

        versions = self.get_software_versions()
        versions.software_version = [
            s if s.startswith("master") or s == "IonTracker" else f'<b style="color:red">{s}</b>'
            for s in versions.software_version
        ]

        self.rep.htmls(
            """
            <table>
            <tr> <td colspan=7> %s </td> </tr>
            <tr> <td> %s </td> <td> %s </td> <td> %s </td> <td> %s </td> <td> %s </td><td> %s </td></tr>
            <tr> <td></td><td><h3>Software<br>Status</h3></td><td colspan=5>%s</td></tr>
            </table>
            """
            % (
                boxit(LibraryService.get_instance().get_study_title(f.study_id), tag="h2"),
                boxit(f"Library<br>{f.lib_id}", tag="h2"),
                boxit(f"Feature<br>{f.serial}", tag="h2"),
                boxit(f"m/z<br>{f.mz}", tag="h2"),
                boxit(f"RI<br>{f.ri}", tag="h2"),
                boxit(f"{n_peaks}<br>peaks", tag="h2"),
                boxit(f"{n_spectra}<br>MS2 spectra", tag="h2"),
                boxit(versions.to_html(index=False, col_space=200, escape=False), tag="h4"),
            )
        )

    def get_software_versions(self):
        f = self.f
        versions = [
            dict(element="Elucidation", software_version=f.gitbranch, run_date=f.elucidation_time,),
            dict(
                element="Discovery",
                software_version=f.discoverer_version,
                run_date=f.discovery_time,
            ),
        ]
        for k, (v, t) in sorted(f.model_versions.items()):
            versions.append(
                dict(element=k + " models", software_version=v + "&nbsp;&nbsp;", run_date=t,)
            )
        versions = pd.DataFrame(versions)
        return versions

    def report_on_chemont_predictions(self):
        """
        Adds HTML for the Chemical Ontology prediction section of the report
        @return: None
        @rtype: None
        """
        rep, f = self.rep, self.f
        df = f.chemont_predictions
        rep.start_accordion(
            f"Chemical Ontology Predictions ({leng(df)} predictions)" + self.reminder
        )
        f.summary["chemont_preds"] = leng(df)
        if not is_empty_df(df):
            rep.htmls(reformat(df), build_explainer(df))
        else:
            rep.comment("No Chemical Ontology Predictions")
        rep.end_accordion()

    def report_on_smiles_predictions(self):
        """
        Adds HTML for the SMILES fragment prediction section of the report
        @return: None
        @rtype: None
        """
        rep, f = self.rep, self.f
        df = f.ms2_smiles_predictions
        rep.start_accordion(
            f"Substructure SMILES Predictions ({leng(df)} predictions)" + self.reminder
        )
        if not is_empty_df(df):
            rep.htmls(reformat(df), build_explainer(df))
            f.summary["num_smiles_preds"] = leng(df)
            f.summary["largest_smiles_prediction"] = df.head(1).subsmiles.min()
            f.summary["fraction_predicted"] = df.mass_frac.max()
        else:
            rep.comment("No Substructure SMILES Predictions")
        rep.end_accordion()

    def candidate_formula_consistency_report(self):
        """
        Adds HTML for the Candidate Formula Consistency section of the report
        @return: None
        @rtype: None
        """
        rep, f = self.rep, self.f
        df = f.fmla_consistency_report
        rep.start_accordion(
            f"Candidate Formula Consistency Report ({leng(df)} formulas)" + self.reminder
        )

        iso_df = f.isotope_df
        summary = self.summary
        summary["isotopes"] = self.create_isotope_summary_string(iso_df, -1)
        if leng(iso_df) > 0:
            rep.htmls(
                reformat(iso_df, alt="No MS1 isotopes detected"), build_explainer(iso_df),
            )
        else:
            rep.comment("No MS1 isotopes detected")
        df = f.fmla_consistency_report
        if leng(df) > 0:
            rep.htmls(reformat(df, index="cand_fmla"), build_explainer(df))
        else:
            rep.comment("No Candidate Formulas")
        rep.end_accordion()

    def peak_plots(self):
        """
        Plots the peak from every sample with sample label.  Too much ink if many samples.
        :return: None
        """
        rep, f, df = self.rep, self.f, self.f.peak_df

        if not (1 <= leng(df) <= 250):
            return
        rep.start_accordion("Peak Plots" + self.reminder)
        df = df.sort_values("area", ascending=False)
        rep.htmls(self.img(self.plot_every_peak()))
        rep.end_accordion()

    def integrated_peaks(self):
        """
        Emits HTML for the table of integrated peaks.
        Only samples with peaks are displayed.  They are sorted by decreasing area.
        Note that DataFrame.to_html() is used.
        The styling stuff used by Reformatter.reformat is a GIANT performance problem for a big table.
        @return:
        @rtype:
        """
        # Now show the dataframe of peaks for all samples.
        f, rep, df = self.f, self.rep, self.f.peak_df
        if leng(df) <= 1:
            rep.comment("No peaks")
            return
        df = df[df.area > 1].sort_values("area", ascending=False)
        if leng(df) == 0:  # not sure if we need to check again.
            rep.comment("No peaks")
            return
        rep.start_accordion("Integrated Peaks" + self.reminder)
        for col in qw("ri area height start_scan apex_scan end_scan"):
            df[col] = np.int32(df[col])
        df.ms2 = ["" if np.isnan(x) else "***" for x in df.ms2]
        df = df[qw("sample_name task_id ri area mass apex_scan apex_rt ms2 raw_file")]
        df = df[qw("task_id ri mass apex_scan apex_rt sample_name area ms2 raw_file")]

        # rep.htmls(reformat(df, index="sample_name", alt="no peaks"))  ## No, no, no!
        rep.htmls(df.to_html(index=False, col_space=100))
        rep.end_accordion()

    def peak_statistics(self):
        """
        Emit HTML for summary statistics about peaks: mass, ri, apex_rt, area.
        @return:
        @rtype:
        """

        f, rep, df = self.f, self.rep, self.f.peak_df
        if leng(df) <= 1:
            return
        rep.start_accordion("Peak Statistics" + self.reminder)
        try:
            f.peak_stats
        except:
            f.peak_stats = pd.concat(
                [f.peak_df[c].describe().to_frame().T for c in qw("mass ri apex_rt area")],
                ignore_index=True,
            ).T
            f.peak_stats.columns = qw("mass ri apex_rt area")
        peak_stats = f.peak_stats
        peak_stats["stat"] = peak_stats.index
        peak_stats = peak_stats[qw("stat mass ri apex_rt area")]
        rep.htmls(
            reformat(peak_stats, alt="no peaks stats"), self.img(self.plot_peak_histogram()),
        )
        # rep.htmls(img(plot_peak_stats(peak_df)))
        rep.end_accordion()

    def report_on_ms1(self):
        """
        Emits HTML for isotope and buddy tables, as well as library matches to masses inferred from buddies.
        :return: None
        """
        rep = self.rep
        f = self.f
        summary = self.summary
        buddy_df = f.buddy_df
        if leng(buddy_df) == 0:
            return

        rep.start_accordion("MS1 Analysis" + self.reminder)
        if f.ion_names is not None:
            summary["ion_names"] = ", ".join(sorted(f.ion_names))
            summary["inferred_masses"] = ", ".join(["%.4f" % x for x in f.inferred_masses])

        summary["lib_buddies"] = (
            ", ".join(list(set(f.buddy_lib_df.name))) if f.buddy_lib_df is not None else ""
        )

        rep.htmls(self.img(self.plot_ms1_consensus(buddy_df)))
        buddy_df = buddy_df[qw("mz_buddy tally ratio theor_delta inferred_m buddy feature")]
        rep.htmls(
            reformat(buddy_df, alt="No MS1 features detected"), build_explainer(buddy_df),
        )
        lib_df = f.buddy_lib_df
        rep.comment("<h3>Library Matches to Inferred Masses</h3>")
        if is_empty_df(lib_df):
            rep.comment("No library matches to inferred masses")
        else:
            rep.htmls(
                reformat(lib_df, alt="No library matches to inferred masses"),
                build_explainer(lib_df),
            )
        rep.end_accordion()
        return

    def create_single_isotope_string(self, t, pos):
        """
        Creates a string for the isotope summary from a row (itertuple) of isotope_df.
        :param pos:
        :param t: a row of isotope_df
        :return: A string such as 'has[C13-1]/2.9'
        """
        if "H" in t.isotope or "D" in t.isotope:
            return f"{t.has_is} {t.isotope[2:-1]}"
        lis = list(t)
        predict = [x for x in lis if type(x) is float and not np.isnan(x)][pos]
        # count, r2, ratio
        return f"{t.has_is} {t.isotope[2:-1]}/{predict:.1f}"

    def create_isotope_summary_string(self, isotope_df, pos):
        """
        Creates a string summarizing information in isotope_df for inclusion on the spreadsheet.
        :param pos:
        :param isotope_df: The dataframe of isotopes found in MS1 analysis.
        :return: summary string.
        """
        if isotope_df is None:
            return ""
        else:
            return ", ".join(
                [self.create_single_isotope_string(t, pos) for t in isotope_df.itertuples()]
            )

            # Report on stuff in the library with the same mass as the feature.

    def mass_library_report(self):
        """
        Emits html for an according displaying the Mass Library Report; see reports for example.
        @return: None; emits html
        @rtype: NoneType
        """
        rep = self.rep
        df = self.f.library_report
        rep.start_accordion(f"Mass Library Report ({leng(df)} matches)" + self.reminder)
        if is_empty_df(df):
            rep.comment("No matches")
        else:
            rep.htmls(
                reformat(df, index="comp_id", alt="No library matches"), build_explainer(df),
            )
        rep.end_accordion()

    def ms2_spectrum_report(self):
        """
        Emits html for an accordion displaying the MS2 consensus as table and graphic.
        @return: None; emits html
        @rtype: NoneType
        """
        f, rep = self.f, self.rep
        consensus_df = f.consensus_df
        rep.start_accordion(f"Consensus MS2 Spectrum ({leng(consensus_df)} ions)" + self.reminder)
        consensus_df = f.consensus_df
        rep.htmls(self.img(self.plot_consensus(consensus_df)), reformat(consensus_df))
        rep.htmls(build_explainer(consensus_df))
        rep.end_accordion()

    def ms2_ion_consistency_report(self):
        """
        Emits html for an according displaying the Ion Consistency Report; see reports for example.
        @return: None; emits html
        @rtype: NoneType
        """
        f, rep, df = self.f, self.rep, self.f.ion_consistency_report
        if leng(df) == 0:
            rep.comment("No MS2 Ion Consistency Report")
            return
        rep.start_accordion("MS2 Ion Consistency Report" + self.reminder)
        rep.htmls(reformat(df), build_explainer(df))
        rep.end_accordion()

    def report_on_spectral_predictions(self):
        """
        Combines predictions from spectral matching from three sources, ion+loss, database, database+modification.
        Creates a single table for report, and adds it to an accordion
        @return:
        @rtype:
        """

        dfs = [
            self.f.structure_prediction_report,
            self.f.modified_spectrum_match_report,
            self.f.predicted_spectrum_match_report,
        ]
        dfs = [df for df in dfs if not is_empty_df(df)]
        if dfs:
            # Give results a tiny boost depending on their degree of connection to the known world.
            for i, df in enumerate(dfs):
                df.score += i * 1e-4
            df = pd.concat(dfs, ignore_index=True)

            def agglomerate(tab):
                """
                The rows of tab all have the same SMILES.  Agglomerate names and sources into a single row and return
                the one-row table.
                @param tab: DataFrame to agglomerate.
                @type tab: DataFrame
                @return: DataFrame with one row.
                @rtype: DataFrame
                """
                if len(tab) == 1:
                    return tab
                rows = list(tab.itertuples())
                row1 = rows[0]
                return pd.DataFrame(
                    [
                        dict(
                            smiles=row1.smiles,
                            name=reduce(add, [row.name for row in rows]),
                            source=reduce(add, [row.source for row in rows]),
                            formula=row1.formula,
                            mass=row1.mass,
                            score=row1.score,
                            spectrum_match=row1.spectrum_match,
                            ri_pred=row1.ri_pred,
                        )
                    ]
                )

            df = df.groupby("smiles").apply(agglomerate)
            df = df.sort_values("score", ascending=False).reset_index(drop=True)
            df = df[df.score >= df.score.max() * 0.66]
            df = df[df.score >= 0.1]
            df = df[qw("smiles spectrum_match score formula mass ri_pred source name")]

        else:
            df = None

        self.rep.start_accordion(
            f"Database Matches by Spectral Prediction ({leng(df)} matches)" + self.reminder
        )
        self.rep.htmls(reformat(df, index="score", alt="<h3>No Matches</h3>"), build_explainer(df))
        summary = self.summary  # this is also f.summary
        if not is_empty_df(df):
            max_score = df.score.max()
            summary["best_spectrum_prediction"] = df.smiles[df.score == max_score].max()
            summary["max_spectrum_score"] = round(max_score, 3)
            summary["num_best_predictions"] = (df.score == max_score).sum()
        else:
            summary["best_spectrum_prediction"] = ""
            summary["max_spectrum_score"] = 0
            summary["num_best_predictions"] = 0
        self.rep.end_accordion()

    def report_on_structure_assignments(self):
        """
        Adds a subsection on each MS2 ion to the report including SMILES analysis.
        :return: None.
        """
        rep = self.rep
        f = self.f
        consensus_df = f.consensus_df
        ms2_ion_dict = f.ms2_ions
        df = f.structure_prediction_report
        # if not is_empty_df(df):
        #     max_support = df.supporting_pairs.max()
        # else:
        #     max_support = 0
        rep.start_accordion(f"Fragment Structure Assignment" + self.reminder)
        # rep.comment("<h3>Ion+Loss Predictions</h3>")
        # self.summary["ion_loss_predictions"] = leng(df)
        # self.summary["max_ion_loss_score"] = df.score.max()
        # if not is_empty_df(df):
        #     rep.comment(f"{len(df)} paired predictions.")
        #     # if len(df) > 10:
        #     #     rep.comment(f"{len(df)} paired predictions.  Showing only 10.")
        #     rep.htmls(
        #         reformat(df.head(10), alt="<h3>No paired prediction</h3>"), build_explainer(df)
        #     )
        # All of this is just to figure out if an explainer is needed.
        for t in consensus_df.itertuples():
            smiles = ms2_ion_dict[t.bucket]["ms2smiles"]
            if not is_empty_df(smiles):
                rep.htmls(build_explainer(smiles))
                break
        # ;padding:12px;font-size:12px;;width:40%
        for t in consensus_df.itertuples():
            df = consensus_df[consensus_df.bucket == t.bucket]
            df = reformat(df.T, index=True)
            smiles = ms2_ion_dict[t.bucket]["ms2smiles"]
            loss_smiles = ms2_ion_dict[t.bucket]["loss_smiles"]
            rep.start_small_accordion(
                f"""<i style="color:blue;font-size:12px">
                            m/z = {t.mz:.4f}&emsp;&emsp;&emsp; intensity = {t.norm_int:.3f}&emsp;&emsp;&emsp;
                            ({leng(smiles)} ion SMILES;  {leng(loss_smiles)} loss SMILES)</i>"""
            )
            report_smiles = reformat(smiles, alt="<h3>No SMILES<br>fragments<br>match mass</h3>")

            report_loss_smiles = reformat(
                loss_smiles, alt="<h3>No SMILES<br>fragments<br>match loss</h3>"
            )
            rep.htmls(
                df, "&emsp;&emsp;&emsp;", report_smiles, "&emsp;&emsp;&emsp;", report_loss_smiles,
            )
            rep.end_accordion()
        rep.end_accordion()

    # TODO This guy and the others should probably return HTML elements rather than receiving a FileReporter argument.

    def report_on_ms2_spectrum_matching(self):
        """
        Add results of matching MS2 spectrum against library.
        :return: None.
        """
        rep, f = self.rep, self.f
        df = f.ms2_library_report
        if df is not None:
            df = df[df.match_score > 0.1].copy()
        rep.start_accordion(
            f"ChemInventory Matching by Spectral Similarity ({leng(df)} matches)" + self.reminder
        )
        if is_empty_df(df) is None:
            rep.comment("No matches to MS2 library")
        else:
            rep.htmls(reformat(df), build_explainer(df))
        rep.end_accordion()

    def report_on_rediscovery(self):
        """
        Add results of matching MS2 spectrum against X- and  NMD- compound of same mass.
        :return: None.
        """
        rep, f = self.rep, self.f
        df = f.rediscovery_report
        rep.start_accordion(f"Rediscovery Report ({leng(df)} matches)" + self.reminder)
        if is_empty_df(df) is None:
            rep.comment("New Discovery")
        else:
            rep.htmls(reformat(df), build_explainer(df))
        rep.end_accordion()

    def plot_peak_histogram(self):
        """
        Plots a histogram of peak areas.  Leaves the plot pending in matplotlib.  There's surely a better way.
        @return: None
        @rtype: None
        """
        f = self.f
        df = f.peak_df
        sns.distplot(
            a=df.area,
            bins=50,
            hist=True,
            kde=False,
            hist_kws=dict(color="blue", alpha=1),
            rug=True,
            rug_kws=dict(height=0.025, color="red", alpha=0.25),
        )
        plt.title(f"Feature {f.serial}")

    def plot_every_peak(self):
        logging.debug("plot_every_peak")
        f = self.f
        df = f.peak_df
        height = min(
            65535 / 200, 0.5 + 0.2 * len(df)
        )  # matplotlib imposes a max height = 2^16 pixels ~ 27 feet!
        plt.figure(figsize=(10, height))
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        plt.yticks(list(range(1, len(df) + 1)), list(df.sample_name), fontsize=7)
        plt.plot(list(df.area), range(1, len(df) + 1), ".")
        x_mean = df.area.mean()
        plt.axvline(x=x_mean, color="r")
        plt.axvline(x=x_mean + np.std(df.area), color="g")
        plt.axvline(x=x_mean - np.std(df.area), color="g")
        plt.title(f"Feature {f.serial}")
        plt.ylabel("sample name")
        plt.xlabel("area")

    def plot_consensus(self, consensus_df, show=False):
        """
        Plots consensus MS2 spectrum.
        There is a vertical bar at each m/z with height corresponding to relative intensity.
        :param consensus_df: has a row for each MS2 ion.
        :param show: If true, plt.show() is called.  If false, the plot is "left hanging" to be put in an HTML report.
        :return: None.
        """
        logging.debug("plot_consensus")

        masses = consensus_df.bucket.apply(bin_base)
        ratios = consensus_df.norm_int
        plt_reset()
        plt.figure(figsize=(10, 5))
        for m, i in zip(masses, ratios):
            plt.plot([m, m], [0, i], lw=2, color="r")
            y = i + (1 - i) / 2
            plt.plot([m, m], [i, y], lw=0.5, color="lightgray")
            plt.text(m, y, f"{m:.2f}", fontsize=7, rotation=90)
        plt.plot(
            masses, ratios, ".b"
        )  # Put caps on the bars.  Useful when two masses are very close.

        if show:
            plt.show()

    def plot_ms1_consensus(self, ion_df, show=False):
        """
        Plots consensus MS2 spectrum.
        There is a vertical bar at each m/z with height corresponding to relative intensity.
        :param ion_df: has a row for each MS2 ion.
        :param show: If true, plt.show() is called.  If false, the plot is "left hanging" to be put in an HTML report.
        :return: None.
        """
        logging.debug("plot_peak_histogram")
        masses = list(ion_df.mz_buddy) + [self.f.mz]
        ratios = list(ion_df.ratio) + [1]
        top = max(ratios)
        plt_reset()
        plt.figure(figsize=(10, 5))
        for m, i in zip(masses, ratios):
            plt.plot([m, m], [0, i], lw=2, color="r")
            y = i + (top - i) / 2
            plt.plot([m, m], [i, y], lw=0.5, color="lightgray")
            plt.text(m, y, f"{m:.2f}", fontsize=7, rotation=90)
        plt.plot(
            masses, ratios, ".b"
        )  # Put caps on the bars.  Useful when two masses are very close.

        if show:
            plt.show()
