"""Offline MLIMS snapshot loader for elucidation without live Oracle access.

Reads mlims_tables.pkl from REFERENCE_DIR and mirrors MlimsLibraryService
queries from cached DataFrames. Used when spec2mol inference runs air-gapped.
"""

import os

import pandas as pd
from elucidation.conveniences import pickle_load, qw
from dotenv import load_dotenv
from mdsutils.easyQuery import EasyQuery


class CachedMlimsLibraryService:
    """Singleton library service backed by a pickled MLIMS table dump."""

    __instance__ = None
    __name_dict__ = None

    def __init__(self):
        load_dotenv()
        d = pickle_load(f'{os.environ["REFERENCE_DIR"]}/mlims_tables.pkl')
        self.cle_df = d["chro_lib_entry"]
        self.clm_df = d["chro_lib_masses"]
        self.cce_df = d["chro_comp_entry"]
        self.chem_df = d["chemicals"]
        self.syn_df = d["chemical_synonym"]
        self.ss_df = d["sample_sets"]

    @staticmethod
    def get_instance():
        if CachedMlimsLibraryService.__instance__ is None:
            CachedMlimsLibraryService.__instance__ = CachedMlimsLibraryService()
        return CachedMlimsLibraryService.__instance__

    def get_name_dict(self) -> dict:
        if CachedMlimsLibraryService.__name_dict__ is not None:
            return CachedMlimsLibraryService.__name_dict__
        d = dict()

        d.update(zip(self.cce_df.comp_id, self.cce_df.name))
        df = self.cce_df[~self.cce_df.report_name.isnull()]
        d.update(zip(df.comp_id, df.report_name))

        foo = self.cle_df.merge(self.chem_df[~self.chem_df.name.isnull()], on="chemical_id")
        d.update(zip(foo.comp_id, foo.name))
        foo = self.cle_df.merge(self.syn_df[~self.syn_df.name.isnull()], on="chemical_id")
        d.update(zip(foo.comp_id, foo.name))
        d = {k: v.replace("ß", "beta") for k, v in d.items()}  # There's always something.
        CachedMlimsLibraryService.__name_dict__ = d
        return d

    def getMs2(self, lib_id):
        clm_df = self.clm_df[self.clm_df.lib_id == lib_id]
        clm_parent = clm_df[clm_df.name.isin(qw("m-H m+H")) & (clm_df.quant_mass == "Y")][
            qw("lib_id comp_id mass name chro_lib_masses_id")
        ]
        clm_parent.columns = qw("lib_id comp_id parent_mass parent_name id")
        clm_child = clm_df[~clm_df.parent_id.isnull()][qw("parent_id mass mass_ratio")]
        clm_child.columns = qw("id mass mass_ratio")
        clm_df = clm_parent.merge(clm_child, on="id")

        cle_df = self.cle_df
        cle_df = cle_df[(cle_df.lib_id == lib_id) & (cle_df.ri > 0) & (cle_df.confidence >= 50)]

        chem_df = self.chem_df
        chem_df = chem_df[~chem_df.smiles.isnull() & (chem_df.smiles != "")].drop(
            qw("created_on status"), axis=1
        )
        smiles_df = chem_df.merge(cle_df[cle_df.lib_id == lib_id], on="chemical_id")[
            qw("comp_id smiles")
        ]

        comp_ids = smiles_df.comp_id

        clm_df = clm_df[clm_df.comp_id.isin(comp_ids)].copy()
        clm_df = clm_df.merge(smiles_df, on="comp_id")
        clm_df["name"] = clm_df.comp_id.apply(self.get_name_dict().get)
        return clm_df

    def get_study_title(self, study_id):
        s = f"Study {study_id}: "
        try:
            return s + self.ss_df[self.ss_df.sample_set_id == int(study_id)].name.min()
        except:
            return s

    def get_uniq_fmlas(self, lib_id):
        chem_df = self.chem_df
        chem_df = chem_df[~chem_df.mol_formula.isnull()]
        chem_df = chem_df.merge(self.cle_df[self.cle_df.lib_id == lib_id], on="chemical_id")
        return chem_df.mol_formula

    def get_isotope_from_mlims(self, iso, lib_id, comp_ids):
        clm_df = self.clm_df
        clm_df = clm_df[
            (clm_df.name == f"m-H[{iso}-1]")
            & (clm_df.lib_id == lib_id)
            & clm_df.comp_id.isin(comp_ids)
        ]
        return clm_df

    def get_library_ions(self, lib_id: int) -> pd.DataFrame:

        # Should we be reporting ALL libraries?
        cle_df = self.cle_df
        cle_df = cle_df[(cle_df.lib_id == lib_id) & (cle_df.ri > 0) & (cle_df.confidence > 0)]
        # This relates to MM's "defunct" compound, and you'd have to ask him why it's a good idea.
        cle_df = cle_df[
            cle_df.retired_for_comp_id.isnull() | cle_df.retired_for_comp_id.isin({62832})
        ]  # DEFUNCT
        cle_df = cle_df[qw("comp_id chemical_id target_mass ri ri_window rt rt_window")]
        cle_df.columns = qw("comp_id chemical_id neutral_mass ri ri_window rt rt_window")

        clm_df = self.clm_df

        clm_df = clm_df[(clm_df.lib_id == lib_id) & ~clm_df.name.isnull() & (clm_df.weighting >= 1)]
        clm_df = clm_df[qw("lib_id comp_id name mass mass_ratio quant_mass weighting")]
        clm_df.columns = qw("lib_id comp_id ion_name ion_mass mass_ratio quant_mass weighting")

        chem_df = self.chem_df[qw("chemical_id mol_formula smiles")]

        df = clm_df.merge(cle_df, on="comp_id").merge(chem_df, on="chemical_id", how="left")
        df["name"] = df.comp_id.apply(self.get_name_dict().get)
        return df

    def get_chemicals_for_formula(self, fmla: str, lib_id: int) -> pd.DataFrame:
        # Casson: I could not find any examples of these in the chemicals table
        # Rex: that's correct.  fmlaFinder may return a formula with Br' = [81]Br.
        # The correct response is "no chemicals in database".
        # Return empty frame if formula contains a single quote
        if "'" in fmla:  # e.g., Br'
            cols = [
                "comp_id",
                "chemical_id",
                "ri",
                "target_mass",
                "mol_formula",
                "name",
                "ion_mass",
            ]
            return pd.DataFrame(data=None, columns=cols)
        chem_df = self.chem_df[self.chem_df.mol_formula == fmla][qw("chemical_id mol_formula")]
        cle_df = self.cle_df[self.cle_df.lib_id == lib_id][qw("chemical_id ri target_mass comp_id")]
        df = chem_df.merge(cle_df, on="chemical_id")
        df["name"] = df.comp_id.apply(self.get_name_dict().get)
        return df

    def get_known_ms2(self, lib_id: int, mass_ratio_lower_limit: float = 1):
        clm_df = self.clm_df[self.clm_df.lib_id == lib_id]
        clm_parent = clm_df[clm_df.name.isin(qw("m-H m+H")) & (clm_df.quant_mass == "Y")][
            qw("lib_id comp_id mass name chro_lib_masses_id")
        ]
        clm_parent.columns = qw("lib_id comp_id parent_mass parent_name id")
        clm_child = clm_df[
            ~clm_df.parent_id.isnull() & (clm_df.mass_ratio >= mass_ratio_lower_limit)
        ][qw("parent_id mass mass_ratio")]
        clm_child.columns = qw("id mass mass_ratio")
        clm_df = clm_parent.merge(clm_child, on="id")

        cle_df = self.cle_df
        cle_df = cle_df[(cle_df.lib_id == lib_id) & (cle_df.ri > 0) & (cle_df.confidence >= 50)][
            qw("comp_id chemical_id ri")
        ]
        chem_df = self.chem_df
        chem_df = chem_df[~chem_df.smiles.isnull() & (chem_df.smiles != "")]
        chem_df = chem_df[qw("chemical_id mol_formula smiles iupac_name")]

        cle_chem_df = cle_df.merge(chem_df, on="chemical_id")

        df = clm_df.merge(cle_chem_df, on="comp_id")
        df["name"] = df.comp_id.apply(self.get_name_dict().get, "")
        df = df.sort_values("name")
        df = df.reset_index(drop=True)  # ??
        return df

    def getUnknownMs2(self, lib_id):  # Gets SOME unknowns, maybe not all.
        clm_df = self.clm_df[self.clm_df.lib_id == lib_id]
        clm_parent = clm_df[clm_df.name.isin(qw("m-H m+H")) & (clm_df.quant_mass_par == "Y")][
            qw("lib_id comp_id mass name chro_lib_masses_id")
        ]
        clm_child = clm_df[~clm_df.parent_id.isnull()]["parent_id mass mass_ratio"]
        clm_df = clm_parent.merge(
            clm_child,
            left_on="chro_lib_masses_id",
            right_on="parent_id",
            suffixes=qw("_par _child"),
        )

        cle_df = self.cle_df
        cle_df = cle_df[(cle_df.lib_id == lib_id) & (cle_df.ri > 0) & (cle_df.confidence > 0)]
        cle_df = pd.concat(
            [cle_df[cle_df.retired_for_comp_id.isnull()], cle_df[cle_df.retired_for_comp_id == 0]],
            ignore_index=True,
        )

        clm_df = clm_df[clm_df.comp_id.isin(cle_df.comp_ids)].copy()
        clm_df["name"] = clm_df.comp_id.apply(self.get_name_dict().get)
        clm_df = clm_df[qw("lib_id comp_id mass_par name parent_id mass_child mass_ratio_child")]
        clm_df.columns = qw("lib_id comp_id parent_mass parent_name id mass mass_ratio")
        return clm_df

    def get_unknown_ms2(self, lib_id: int, mass_ratio_lower_limit: float = 1) -> pd.DataFrame:

        cle_df = self.cle_df
        cle_df = cle_df[(cle_df.lib_id == lib_id) & (cle_df.ri > 0) & (cle_df.confidence == 100)]
        cle_df = pd.concat(
            [cle_df[cle_df.retired_for_comp_id.isnull()], cle_df[cle_df.retired_for_comp_id == 0]],
            ignore_index=True,
        )
        cle_df["name"] = cle_df.comp_id.apply(self.get_name_dict().get)
        cle_df = cle_df[cle_df.name.str.match(r"^X-\d+$")]
        cle_df = cle_df[qw("comp_id name putative_details")]

        clm_df = self.clm_df[
            (self.clm_df.lib_id == lib_id) & self.clm_df.comp_id.isin(cle_df.comp_id)
        ]
        clm_parent = clm_df[clm_df.name.isin(qw("m-H m+H")) & (clm_df.quant_mass == "Y")][
            qw("lib_id comp_id mass chro_lib_masses_id")
        ]
        clm_child = clm_df[
            ~clm_df.parent_id.isnull() & (clm_df.mass_ratio >= mass_ratio_lower_limit)
        ][qw("parent_id mass mass_window mass_ratio")]
        clm_df = clm_parent.merge(
            clm_child,
            left_on="chro_lib_masses_id",
            right_on="parent_id",
            suffixes=qw("_par _child"),
        )

        clm_df = clm_df.merge(cle_df, on="comp_id")
        clm_df = clm_df[
            qw("lib_id comp_id mass_par name parent_id mass_child mass_ratio putative_details")
        ]
        clm_df.columns = qw("lib_id comp_id parent_mass parent_name id mass mass_ratio details")
        clm_df.details = clm_df.details.apply(clean_details)
        clm_df = clm_df.reset_index(drop=True)
        return clm_df

    def get_chemicals(self):
        return self.chem_df

    def add_filenames_to_peaks(self, peak_df):
        # These filenames should hove been attached in discovery...
        with EasyQuery(db="mlims") as eq:
            df = eq.easy_query(
                "select task_id, filename from worksheet_entry where task_id in ($BATCH)",
                BATCH=peak_df.task_id,
            )
        d = dict(zip(df.task_id, df.filename))  # I guess this could be a merge.
        peak_df["raw_file"] = peak_df.task_id.apply(d.get)


def clean_details(details):
    details = str(details)
    if details == "nan":
        details = ""
    else:
        details = details.replace("\n", " | ").replace("\r", "")
        # if '|' in details: print(details[:-3])
        if "|" in details[-3:]:
            details = details[:-3]
    return details
