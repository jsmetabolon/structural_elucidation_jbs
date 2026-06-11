import logging
from collections import defaultdict
from itertools import product as cartesian

import elucidation.chemFacts as cf
import numpy as np
import pandas as pd
import weightedstats as ws
from elucidation.formulaFinder import FormulaFinder
from elucidation.fragmentFinder import FragmentFinder
from elucidation.rdkitUtils import fmla_for_smiles, smiles_mass, glue_smiles_fragments
from elucidation.conveniences import qw, sort_df_by_vector


class MS2Ion(object):
    def __init__(self, bucket_, consensus_):
        """
        Builds one ion for the MS2 consensus.
        :param bucket_: the bucket containing the observed ions to be aggregated.
        :param consensus_: The consensus object of which this ion is a part.
        """
        self.me = self
        self.bucket = bucket_
        self.fmlas = None
        self.matched_smiles = []
        self.ms2smiles = None
        self.loss_smiles = None
        self.loss_m = -1

        self.consensus = consensus_  # This is the consensus spectrum of which this ion is a part.
        self.lib_id = self.consensus.feature.lib_id
        # Spectrum object contains a dictionary of normalized intensities (intensities relative to a "normalizer" ion.)
        # Get the ones that belong to this bucket.
        self.ratios = [
            sp.norm_dict[bucket_] for sp in consensus_.support if bucket_ in sp.norm_dict
        ]
        # Also, get the ABSOLUTE intensities of the normalizers.
        self.norm_intensities = [
            sp.get_bucket_int(consensus_.normalizer)
            for sp in consensus_.support
            if bucket_ in sp.norm_dict
        ]
        if len(self.norm_intensities) < 2:
            self.me = None
            return

        self.tally = len(self.ratios)
        # Find the masses for this ion.
        masses = [sp.get_bucket_mz(bucket_) for sp in consensus_.support if bucket_ in sp.norm_dict]
        masses = [m for m in masses if m is not None]

        # Initialize isc stuff.
        self.mass_range = (min(masses), max(masses))
        self.med_ratio = ws.weighted_median(self.ratios, weights=self.norm_intensities)
        self.mz = np.median(masses)

        # Generate the SMILES report.
        # Get smiles and facts for this ion. They come from this and adjacent bins in the ion_facts table.
        # TODO THis would be a place where we could refine the mass calculation
        self.build_structure_assignment_reports()

    def get_me(self):
        return self.me

    def consensus_stats(self):
        """
        Create a summary row for the ion that the consensus parent can add to its dataframe.
        """
        return {
            "mass_range": "[%.4f,%.4f]" % self.mass_range,
            "mz": self.mz,
            "norm_int": self.med_ratio,
            # "norm_pct_dev": self.med_dev,
            "bucket": self.bucket,
            "num_spectra": self.tally,
        }

    def _consistent_ifmlas(self, pfmla):
        """
        Returns a list of the candidate formulas of this ion that are also consistent with candidate
        precursor formula pfmla.  For ion consistency report.
        :param pfmla: a candidate precursor formuula.
        :return: list of compatible formulas.
        """
        return [ifmla for ifmla in self.fmlas if cf.subset_fmla(ifmla, pfmla)]

    def ion_consistency_report(self, pfmlas):
        """
        builds a dictionary from precursor formulas to consistent formulas for this ion.
        Uses formula finder to find formulas for this ion.
        :param pfmlas: list of candidate precursor formulas.
        :return: dictionary from precursor formaulas to lists of ion formulas.
        """
        # Change this to build a dataframe of ions.
        mass_str = "[%.4f,%.4f]" % self.mass_range
        d = {
            "bucket": self.bucket,
            "mass_range": mass_str,
            "mz": self.mz,
            "scan_count": self.tally,
            "rel_intens": self.med_ratio,
        }
        # Could narrow this search to list of ion masses actually observed rather than bucket range.  # TODO
        ff = FormulaFinder.get_instance(self.lib_id)
        ff_df = ff.fmlas_for_ion_mz(self.mz)
        ff_df = sort_df_by_vector(ff_df, ff_df.cand_fmla.apply(cf.get_c_count))
        if self.consensus.feature.ms1.iso_count["[Cl37-1]"] == 0:  # Cl not supported by isotopes.
            ff_df = ff_df[~ff_df.cand_fmla.str.contains("Cl")]
        if self.consensus.feature.ms1.iso_count["[Br81-1]"] == 0:  # Cl not supported by isotopes.
            ff_df = ff_df[~ff_df.cand_fmla.str.contains("Br")]

        self.fmlas = ff_df.cand_fmla
        d["cand_fmlas"] = len(self.fmlas)  # TODO
        # if len(self.fmlas) > 0:
        if len(pfmlas) > 0:
            d.update({pfmla: self._consistent_ifmlas(pfmla) for pfmla in pfmlas})
        else:
            d.update({k: "" for k in pfmlas})
        return d

    def build_structure_assignment_reports(self):

        mz = self.mz
        ff = FragmentFinder.get_instance(self.lib_id)

        rows = list(ff.get_fragments_for_ion_mz(mz))
        if rows:
            df = pd.DataFrame(rows, columns=qw("mass subsmiles lib_occurrences")).sort_values(
                "lib_occurrences", ascending=False
            )
            df["fmla"] = df.subsmiles.apply(fmla_for_smiles)
            total = df.lib_occurrences.sum()
            df["lib_fraction"] = df.lib_occurrences / total
            df["neutral_mass"] = df.subsmiles.apply(smiles_mass)
            self.ms2smiles = df

        # This one is neutral; we've subtracted the charged mass.
        precursor_mz = self.consensus.feature.mz
        self.loss_m, rows = ff.get_fragments_for_loss(mz, precursor_mz)
        if rows:
            df = pd.DataFrame(rows, columns=qw("mass subsmiles lib_occurrences")).sort_values(
                "lib_occurrences", ascending=False
            )
            df["fmla"] = df.subsmiles.apply(fmla_for_smiles)
            total = df.lib_occurrences.sum()
            df["lib_fraction"] = df.lib_occurrences / total
            df["neutral_mass"] = df.subsmiles.apply(smiles_mass)
            self.loss_smiles = df

        if self.ms2smiles is not None and self.loss_smiles is not None:
            pairs = cartesian(self.ms2smiles.itertuples(), self.loss_smiles.itertuples())
            self.matched_smiles = [
                (
                    glue_smiles_fragments(t1.subsmiles, t2.subsmiles),
                    t1.lib_fraction * t2.lib_fraction,
                )
                for t1, t2 in pairs
            ]
            logging.debug(f"{len(self.matched_smiles)} pairs")

    def as_dict(self):
        return {
            attr: self.__getattribute__(attr)
            for attr in qw("bucket mz loss_m ms2smiles loss_smiles matched_smiles")
        }
