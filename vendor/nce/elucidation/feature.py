import logging
from time import time

import elucidation.chemFacts as cf
import numpy as np
import pandas as pd
from elucidation.binning import plus_5ppm, minus_5ppm
from elucidation.formulaFinder import FormulaFinder
from elucidation.formulaScorer import FormulaScorer
from elucidation.libraryService import LibraryService
from elucidation.ms1Analyzer import MS1Analyzer
from elucidation.ms2Consensus import MS2Consensus
from elucidation.postgresSmilesDb import PostgresSmilesDb
from elucidation.rdkitUtils import canonicalize_smiles, fmla_for_smiles
from elucidation.conveniences import (
    qw,
    is_empty_df,
    reorganize_columns,
    leng,
    now_string,
)
from elucidation.riPredictor import RiPredictor
from elucidation.spectrum import Spectrum


class Feature(object):
    def __init__(self, d, opts):
        """
        Creates a new feature from a list of (RI,m/z) pairs.
        :param d: (pickled) dictionary from which discovered feature will be rebuilt for elucidation.
        """
        # The following 7 lines are just to make pycharm happy. The values are copied from d.
        self.ri = None
        self.mz = None
        self.peak_df = None
        self.lib_id = None
        self.ms2_counts = None
        self.comments = []
        self.discoverer_version = None

        for k, v in d.items():
            self.__setattr__(k, v)

        self.ms1_spectra = [Spectrum(d) for d in self.ms1_spectra]
        self.ms2_spectra = [Spectrum(d) for d in self.ms2_spectra]

        if "raw_file" not in self.peak_df.columns:  # This should be done during discovery.
            LibraryService.get_instance().add_filenames_to_peaks(self.peak_df)

        self.mz = self.peak_df.mass.median()
        self.fmlas = []  # This is superfluous.
        self.fmla_df = None
        self.consensus = None
        self.library_report = None
        self.ion_consistency_report = None
        self.fmla_consistency_report = None
        self.max_harmony = -1
        self.internal_standard = ""
        self.iso_count = None
        self.smilesdb_df = None
        self.ms1 = None
        self.elucidation_time = now_string()
        self.gitbranch = opts.gitbranch
        pgsmiles = PostgresSmilesDb.get_instance()
        self.smilesdb_df = pgsmiles.get_records_for_neutral_mass(
            cf.neutralize_precursor_mz(self.mz, self.lib_id)
        )
        self.smilesdb_df = self.smilesdb_df[qw("smiles name source id formula exact_mass")]

    def as_dict(self):
        """
        Turn the object into a dict for pickling.
        @return: dict to be pickled.
        @rtype: dict
        """
        self.iso_count = dict(self.iso_count)  # no defaultdict to pickle, please!
        d = {
            attr: self.__getattribute__(attr)
            for attr in qw("ri mz peak_df lib_id ms2_counts")
            + qw("fmlas fmla_df library_report fmla_consistency_report")
            + qw("ion_consistency_report serial")
            + qw("max_harmony internal_standard iso_count")
            + qw("discoverer_version study_id discovery_time comments elucidation_time gitbranch")
        }
        d.update(self.ms1.as_dict())
        d.update(self.consensus.as_dict())
        return d

    def get_candidate_fmlas(self):
        """
        Uses formula finder to generate candidate formulas for this mass.
        :return: sets fmla_df and fmlas
        """
        t0 = time()
        fmla_df = FormulaFinder.get_instance(self.lib_id).fmlas_for_precursor_mz(self.mz)
        if is_empty_df(fmla_df):  # no formulas to analyze.
            return
        fmla_df["neutral_mass"] = fmla_df.cand_fmla.apply(cf.fmla_mass)

        # It's important to hack away as many formulas as possible as soon as possible before doing expensive stuff.
        fmla_df["neutral"] = fmla_df.cand_fmla.apply(cf.neutralize_fmla)
        logging.debug(f"get_candidate_fmlas neutralize {time() - t0}, {leng(fmla_df)}")

        fmla_df["dbe"] = fmla_df.neutral.apply(cf.dbe)
        # Want dbe to be positive and whole (or just non-negative?)
        fmla_df = fmla_df[((fmla_df.dbe > 0) & (fmla_df.dbe == np.int32(fmla_df.dbe)))]
        logging.debug(f"get_candidate_fmlas dbe {time() - t0}, {leng(fmla_df)}")

        fs = FormulaScorer.get_instance(self.lib_id)
        fmla_df["harmony"] = np.minimum(10, fmla_df.neutral.apply(fs.fmla_harmony))
        fmla_df = fmla_df[(fmla_df.harmony > 1)]
        logging.debug(f"get_candidate_fmlas harmony {time() - t0}, {leng(fmla_df)}")

        if is_empty_df(fmla_df):  # done filtering.
            return

        self.max_harmony = fmla_df.harmony.max()
        self.fmla_df = fmla_df
        self.fmlas = list(self.fmla_df.cand_fmla)  # ff always returns a dataframe

    def create_library_report(self):
        """
        Create a dataframe of all ions already in the library that match this mass.
        :return: sets library_report
        """
        lo_mass = minus_5ppm(self.mz)
        hi_mass = plus_5ppm(self.mz)
        df = LibraryService.get_instance().get_library_ions(self.lib_id)
        df = reorganize_columns(
            df[(lo_mass <= df.ion_mass) & (df.ion_mass <= hi_mass)],
            qw("chemical_id comp_id neutral_mass ion_mass ri mol_formula smiles ion_name name"),
            drop_leftovers=True,
        )
        df.smiles = df.smiles.apply(canonicalize_smiles)
        self.library_report = df.sort_values("ri").reset_index(drop=True)

    def perform_ms1_analysis(self):
        """
        performs the steps of MS1 analysis.
        but it is here right now to make debugging a little more convenient.
        :return: None
        """
        self.ms1 = MS1Analyzer(self, self.ms1_spectra)
        ms1 = self.ms1
        ms1.extract_buddies(ppm=3)
        ms1.analyze_isotopes()
        ms1.build_buddy_library_df()
        self.iso_count = ms1.iso_count
        self.complete_candidate_fmlas()

    def complete_candidate_fmlas(self):
        fmla_df = self.fmla_df
        smilesdb_df = None if leng(self.smilesdb_df) == 0 else self.smilesdb_df.copy()
        neutral_mass = cf.neutralize_precursor_mz(self.mz, self.lib_id)
        fs = FormulaScorer.get_instance(self.lib_id)
        if leng(fmla_df) == 0 and leng(smilesdb_df) == 0:
            return

        rip = RiPredictor.get_instance(self.lib_id)
        if leng(smilesdb_df) > 0:
            smilesdb_df["ppm"] = np.abs(1e6 * ((smilesdb_df.exact_mass / neutral_mass) - 1))
            smilesdb_df.formula = smilesdb_df.smiles.apply(fmla_for_smiles)
            smilesdb_df["dbe"] = smilesdb_df.formula.apply(cf.dbe)
            smilesdb_df["harmony"] = np.minimum(10, smilesdb_df.formula.apply(fs.fmla_harmony))
            smilesdb_df["ri_pred"] = smilesdb_df.smiles.apply(
                rip.predict_for_smiles, other_ri=self.ri
            )
            smilesdb_df["cand_fmla"] = [""] * len(smilesdb_df)
            smilesdb_df.columns = qw(
                "smiles name source id neutral neutral_mass ppm dbe harmony ri_pred cand_fmla"
            )
            # no names in this table.
            smilesdb_df = smilesdb_df[
                qw("cand_fmla ppm neutral_mass neutral dbe harmony ri_pred source id smiles")
            ]
        if leng(fmla_df) > 0:
            fmla_df["ri_pred"] = fmla_df.neutral.apply(rip.predict_for_formula, other_ri=self.ri)

        dfs = [fmla_df, smilesdb_df]
        fmla_df = pd.concat([df for df in dfs if leng(df) > 0], ignore_index=True)

        # Throw away isotopic nonsense -- but not when analyzing X-xxxxx.
        if self.discoverer_version != "IonTracker":
            fmla_df = fmla_df[fmla_df.neutral.apply(self.ms1.is_formula_isotopically_ok)]

        db = LibraryService.get_instance()
        fmla_df["library_matches"] = fmla_df.neutral.apply(
            lambda fmla: leng(db.get_chemicals_for_formula(fmla, self.lib_id))
        )

        # Sort rows by #C
        fmla_df["foo"] = fmla_df.neutral.apply(cf.get_c_count)
        fmla_df = fmla_df.sort_values(
            qw("ri_pred harmony foo ppm"), ascending=(True, False, True, True)
        ).drop("foo", axis=1)

        self.fmla_df = fmla_df
        self.fmla_consistency_report = fmla_df  # These seem to be the same thing now... #TODO
        self.fmlas = sorted(set(self.fmla_df.cand_fmla) - {""})  # And this seems redundant.

    def perform_ms2_analysis(self, opts=None):
        """
        performs the steps of MS2 analysis.  This, too, should be pushed down into Feature, but it is here right
        now to make debugging a little more convenient.
        :return: None
        """

        self.peak_df["ms2"] = self.peak_df.sample_name.apply(self.ms2_counts.get)
        self.consensus = MS2Consensus(self)
        min_ms2_int = 0.01 if opts is None else opts.min_ms2_intensity

        ms2 = self.consensus
        ms2.build(self.ms2_spectra, min_ms2_int)
        ms2.create_consensus_df()
        if ms2.consensus_df is None:
            return
        self.ion_consistency_report = ms2.create_ion_consistency_report(self.fmlas)

        ms2.build_ms2_library_report()
        ms2.build_rediscovery_report()
        # ms2.build_general_ms2_smiles_match_report()
        ms2.build_ms2_smiles_prediction_report()
        # ms2.build_zapped_ms2_smiles_match_report()
        ms2.build_structure_prediction_report()
        ms2.build_chemont_prediction_report()
        ms2.build_predicted_spectrum_match_report()

    def count_ion_matches(self, fmla):
        """
        Counts the number of MS2 ions consistent with this formula.
        :param fmla: a candidate formula for the feature.
        :return: the number of consistent ions.
        """
        df = self.ion_consistency_report
        return sum([1 for f in df[fmla] if type(f) is list and len(f) > 0])

    def count_weighted_ion_matches(self, fmla):
        """
        Counts the number of MS2 ions consistent with this formula,
                weighted by the number of spectra containing the MS2 ion.
        :param fmla: a candidate formula for the feature.
        :return: the weighted count of consistent MS2 ions.
        """
        df = self.ion_consistency_report
        return sum([sc for f, sc in zip(df[fmla], df.scan_count) if type(f) is list and len(f) > 0])

    def create_fmla_consistency_reports(self):
        """
        Creates a dataframe with a row for each candidate formula.  Columns include number of consistent MS2 ions
        and number of library chemicals with the formula.  Formulas are filtered out if they don't match isotopes
        found/not found in MS1 analysis.
        :return: sets fmla_consistency_report and filters ion_consistency report.
        """
        t0 = time()
        fmla_df = self.fmla_df

        if leng(fmla_df) > 0 and self.discoverer_version != "IonTracker":
            fmla_df["ions"] = fmla_df.cand_fmla.apply(self.count_ion_matches)
            fmla_df["weighted_ions"] = fmla_df.cand_fmla.apply(self.count_weighted_ion_matches)
        logging.info(f"create_fmla_consistency_reports ion counts {time() - t0}, {leng(fmla_df)}")

        # Some formulas have been thrown out, so those columns must be removed
        # from ion_consistency_report.
        if self.discoverer_version != "IonTracker":
            good_cols = qw("bucket mass_range scan_count rel_intens cand_fmlas isotopy") + list(
                self.fmla_consistency_report.cand_fmla
            )
            self.ion_consistency_report = self.ion_consistency_report[good_cols]
        logging.info(f"create_fmla_consistency_reports finished {time() - t0}, {leng(fmla_df)}")
