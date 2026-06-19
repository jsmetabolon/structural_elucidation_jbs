"""Live Oracle MLIMS queries for library matching during elucidation.

Singleton service fetching named/unknown compounds, MS2 spectra, ion tables,
and rediscovery data from production MLIMS. Default backend for LibraryService.
"""

import os
import pandas as pd
from dotenv import load_dotenv
from mdsutils.easyQuery import EasyQuery


class MlimsLibraryService:
    """Singleton Oracle-backed access to MLIMS chromatography library tables."""

    __get_named_ms2__ = {}
    __get_unnamed_nmd__ = {}
    __get_named_nmd__ = {}
    __get_unknown_ms2__ = {}
    __get_library_ions__ = {}
    __get_uniq_fmlas__ = {}
    __get_unnamed_ms2__ = {}
    __get_ms1_rediscovery__ = {}
    __get_known_ms2__ = {}
    __instance__ = None
    __name_dict__ = None

    def __init__(self):
        load_dotenv()
        self.mlims_credentials = dict(db=os.environ["DB_NAME"], user=os.environ["DB_USER"],)

    @staticmethod
    def get_instance():
        if MlimsLibraryService.__instance__ is None:
            MlimsLibraryService.__instance__ = MlimsLibraryService()
        return MlimsLibraryService.__instance__

    def get_name_dict(self) -> dict:
        """
        Gets a dictionary of compound names from MLIMS.
        @return: dict: comp_id -> name
        @rtype: dict: int->str
        """
        if MlimsLibraryService.__name_dict__ is not None:
            return MlimsLibraryService.__name_dict__
        d = dict()
        with EasyQuery(**self.mlims_credentials) as eq:
            foo = eq.easy_query("select cce.comp_id, name from limsuser.chro_comp_entry cce")
            d.update(zip(foo.comp_id, foo.name))
            foo = eq.easy_query(
                "select cce.comp_id, report_name from limsuser.chro_comp_entry cce where report_name is not null"
            )
            d.update(zip(foo.comp_id, foo.report_name))

            foo = eq.easy_query(
                """
                select cle.lib_id, cle.comp_id, chem.name 
                from limsuser.chro_lib_entry cle join limsuser.chemicals chem on cle.chemical_id = chem.chemical_id
                where name is not null
                """
            )
            d.update(zip(foo.comp_id, foo.name))

            foo = eq.easy_query(
                """
                select cle.comp_id, cs.name 
                from limsuser.chro_lib_entry cle join limsuser.chemical_synonym cs on cle.chemical_id = cs.chemical_id
                where cs.name is not null and cs.report_name = 'T'
                """
            )
            d.update(zip(foo.comp_id, foo.name))

            d = {k: v.replace("ß", "beta") for k, v in d.items()}  # There's always something.
            MlimsLibraryService.__name_dict__ = d
            return d

    def get_named_ms2(self, lib_id):
        """
        Gets a dataframe of MS2 ions for every named compound in the library.
        Excluded are non-NMD compounds with confidence < 50 and compounds with no SMILES string.
        The result is suitable for training models (once augmented by get_modeling_df.)
        @param lib_id: library
        @type lib_id: library
        @return: dataframe of MS2 ions.
        @rtype: dataframe
        """
        if lib_id in MlimsLibraryService.__get_named_ms2__:
            return MlimsLibraryService.__get_named_ms2__[lib_id]

        # Retired NMD compounds have comp IDs in the normal sequence, but their confidence is set to 0 to avoid
        # unintended curation and release.  Since everything else is filtered by confidence, we need to identify these.
        with EasyQuery(**self.mlims_credentials) as eq:
            nmd_comp_ids = eq.easy_query(
                """select retired_for_comp_id from chro_lib_entry 
                where comp_id>=2000000 and retired_for_comp_id is not null"""
            ).retired_for_comp_id

            query = """
            select
                ch.smiles,
                ch.mol_formula formula,
                parent.comp_id,
                parent.lib_id,
                cle.ri,
                parent.mass parent_mass,
                parent.name parent_name,
    
                child.parent_id id,
                child.mass mass,
                child.mass_window mass_window,
                child.mass_ratio mass_ratio
    
            from
                chro_lib_masses child 
                join chro_lib_masses parent on child.parent_id = parent.chro_lib_masses_id
                join chro_lib_entry cle on parent.comp_id = cle.comp_id and parent.lib_id = cle.lib_id
                join chemicals ch on cle.chemical_id = ch.chemical_id
            where  cle.lib_id in ($LIB_ID) 
                and cle.ri > 0
                and (cle.confidence >= 50 or cle.comp_id in ($BATCH))
                and parent.name in ('m-H','m+H')
                and parent.quant_mass in ('Y','y')
                and length(smiles) > 0
    
            order by comp_id
            """
            df = eq.easy_query(query, LIB_ID=lib_id, BATCH=nmd_comp_ids)
            df["name"] = df.comp_id.apply(self.get_name_dict().get)
            MlimsLibraryService.__get_named_ms2__[lib_id] = df
            return df

    def get_unnamed_ms2(self, lib_id):  # Gets SOME unknowns, maybe not all.
        if lib_id in MlimsLibraryService.__get_unnamed_ms2__:
            return MlimsLibraryService.__get_unnamed_ms2__[lib_id]
        with EasyQuery(**self.mlims_credentials) as eq:
            query = """
            select
    
                parent.comp_id,
                parent.lib_id,
                parent.mass parent_mass,
                parent.name parent_name,
                child.parent_id id,
                child.mass mass,
                child.mass_ratio mass_ratio,
                cle.confidence confidence,
                cle.ri ri,
                cle.created_on, 
                cle.comments,
                cle.putative_details, 
                cle.retired_for_comp_id
    
            from
                chro_lib_masses child 
                join chro_lib_masses parent on child.parent_id = parent.chro_lib_masses_id and child.lib_id = parent.lib_id
                join chro_lib_entry cle on parent.comp_id = cle.comp_id and parent.lib_id = cle.lib_id
    
            where 
                child.lib_id = $LIB_ID
                and cle.ri > 0
                and cle.confidence > 0
                and (cle.retired_for_comp_id is null or cle.retired_for_comp_id = 0 or cle.comp_id >= 2000000)
                and parent.name in ('m+H', 'm-H')
                and parent.quant_mass in ('Y','y')
            """
            df = eq.easy_query(query, LIB_ID=lib_id)
            name_dict = self.get_name_dict()
            df["name"] = df.comp_id.apply(name_dict.get)
            comp_names = df.comp_id.apply(name_dict.get)
            df["unknown_name"] = comp_names
            with open("foo.txt", "w") as f:
                print(sorted(set(comp_names)), file=f)
            df = df[comp_names.str.match(r"^X-\d+$") | comp_names.str.match(r"^NMD")]
            MlimsLibraryService.__get_unnamed_ms2__[lib_id] = df
            return df

    def get_ms1_rediscovery(
        self, lib_id
    ):  # Gets MS1 ions to repeat filter done in discovery in case of reintroduction in peak pulling step.
        with EasyQuery(**self.mlims_credentials) as eq:
            query = """
            select 1.000 match_score, clm.lib_id, clm.comp_id, ri, ri_window, clm.name, clm.name comments,
                    clm.mass parent_mass, clm.mass_window, clm.weighting, cle.retired_for_comp_id
                from chro_lib_entry cle
                join chro_lib_masses clm on clm.lib_id=cle.lib_id and cle.comp_id = clm.comp_id
                where clm.lib_id = $LIB_ID and clm.name is not null 
                  and (RETIRED_FOR_COMP_ID is null OR RETIRED_FOR_COMP_ID = 62832)
                  and weighting >= 0.1
                  and rt >= 0.5
                  and cle.CONFIDENCE >0 
            """
            df = eq.easy_query(query, LIB_ID=lib_id)
            return df

    def get_study_title(self, study_id):
        """
        Get study name from MLIMS, if available.
        @param study_id: a sample set id
        @type study_id: int
        @return: study name string
        @rtype: str
        """
        s = f"Study {study_id}: "
        try:
            with EasyQuery(**self.mlims_credentials) as eq:
                return (
                    s
                    + eq.easy_query(
                        f"select * from sample_sets where SAMPLE_SET_ID = {study_id}"
                    ).name.min()
                )
        except:
            return s

    def get_uniq_fmlas(self, lib_id):
        """
        This is used to train the FormulaScorer
        @param lib_id: yes, a library
        @type lib_id: int
        @return: a set of formulas
        @rtype: set of str
        """
        if lib_id in MlimsLibraryService.__get_uniq_fmlas__:
            return MlimsLibraryService.__get_uniq_fmlas__[lib_id]
        with EasyQuery(**self.mlims_credentials) as eq:
            df = set(
                list(
                    eq.easy_query(
                        f"""
                                select mol_formula from chemicals ch 
                                join chro_lib_entry cle on ch.chemical_id=cle.chemical_id
                                where lib_id = {lib_id} and mol_formula is not null"""
                    ).mol_formula  # just 209 should be fine.
                )
            )
        MlimsLibraryService.__get_uniq_fmlas__[lib_id] = df
        return df

    def get_isotope_from_mlims(self, iso, lib_id, comp_ids):
        """
        Retrieves MS1 ion entries of isotopes of comp_id found in library.
        @param iso: specifies the isotope, e.g. "C13", "Br81"
        @type iso: str
        @param lib_id: library
        @type lib_id: int
        @param comp_ids: list of comp_ids for which to retrieve isotope.
        @type comp_ids: iterable of ints
        @return: dataframe of MS1 ions for isotopes
        @rtype: DataFrame
        """
        with EasyQuery(**self.mlims_credentials) as eq:
            foo = eq.easy_query(
                """select * from chro_lib_masses 
                where name in ('m-H[$ISO-1]', 'm+H[$ISO1-1]') 
                and comp_id in ($BATCH) and lib_id=$LIBID""",
                BATCH=comp_ids,
                ISO=iso,
                ISO1=iso,
                LIBID=lib_id,
            )
            return foo

    def get_library_ions(self, lib_id: int) -> pd.DataFrame:
        if lib_id in MlimsLibraryService.__get_library_ions__:
            return MlimsLibraryService.__get_library_ions__[lib_id]
        with EasyQuery(**self.mlims_credentials) as eq:
            lib_df = eq.easy_query(
                """
                    select clm.lib_id, clm.comp_id, rt, rt_window, ri, ri_window, ch.chemical_id, clm.name ion_name, 
                            cle.target_mass neutral_mass, ch.MOL_FORMULA, ch.smiles,
                            clm.mass ion_mass, clm.mass_window, clm.mass_ratio, clm.quant_mass, clm.weighting
                        from chro_lib_entry cle
                        join chro_lib_masses clm on clm.lib_id=cle.lib_id and cle.comp_id = clm.comp_id
                        left join chemicals ch on cle.CHEMICAL_ID = ch.CHEMICAL_ID
                        where clm.lib_id = $LIBID and clm.name is not null 
                          and (RETIRED_FOR_COMP_ID is null OR RETIRED_FOR_COMP_ID = 62832) -- defunct
                          and weighting >= 1
                          and cle.CONFIDENCE >0 
                          and ri > 0
                """,
                LIBID=lib_id,
            )
            lib_df["name"] = lib_df.comp_id.apply(self.get_name_dict().get)
            MlimsLibraryService.__get_library_ions__[lib_id] = lib_df
            return lib_df

    def get_chemicals_for_formula(self, fmla: str, lib_id: int) -> pd.DataFrame:
        """
        Gets all compounds from library with the given formula.
        @param fmla: formula, e.g. C5H10O
        @type fmla: str
        @param lib_id: lib_id
        @type lib_id: int
        @return: Dataframe of matching compounds
        @rtype: dataframe
        """
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
        with EasyQuery(**self.mlims_credentials) as eq:
            df = eq.easy_query(
                """
                select ch.chemical_id, ri, target_mass, comp_id, mol_formula 
                from chemicals ch join chro_lib_entry cle on ch.chemical_id = cle.chemical_id
                where lib_id = $LIBID and mol_formula = '$FMLA'""",
                LIBID=lib_id,
                FMLA=fmla,
            )
            df["name"] = df.comp_id.apply(self.get_name_dict().get)
            return df

    def get_known_ms2(self, lib_id: int, mass_ratio_lower_limit: float = 1):
        if (lib_id, mass_ratio_lower_limit) in MlimsLibraryService.__get_known_ms2__:
            return MlimsLibraryService.__get_known_ms2__[(lib_id, mass_ratio_lower_limit)]
        comp_id_to_name_dict = self.get_name_dict()
        query = """
            select
                cle.lib_id,
                cle.comp_id,
                chems.chemical_id,
                cle.ri,
                chems.mol_formula,
                chems.smiles,
                chems.iupac_name,
                parent.mass parent_mass,
                parent.name parent_name,
                child.parent_id id,
                child.mass mass,
                child.mass_window,
                child.mass_ratio              
            from
                limsuser.chro_lib_masses child,
                limsuser.chro_lib_masses parent,
                limsuser.chro_lib_entry cle,
                limsuser.chemicals chems
            where child.parent_id = parent.chro_lib_masses_id
              and parent.comp_id = cle.comp_id
              and parent.lib_id = cle.lib_id
              and parent.comp_id = child.comp_id
              and cle.chemical_id = chems.chemical_id
              and child.lib_id = $LIBID
              and cle.confidence >= 50
              and cle.ri > 0
              and parent.name in ('m-H','m+H')
              and parent.quant_mass in ('Y','y')
              and length(chems.smiles) > 0
        """
        with EasyQuery(**self.mlims_credentials) as eq:
            df = eq.easy_query(query, LIBID=lib_id)
        df = df[df.mass_ratio >= mass_ratio_lower_limit]
        df["name"] = df.comp_id.apply(comp_id_to_name_dict.get, "")
        df = df.sort_values("name")
        df = df.reset_index(drop=True)
        MlimsLibraryService.__get_known_ms2__[(lib_id, mass_ratio_lower_limit)] = df
        return df

    def get_unknown_ms2(self, lib_id: int, mass_ratio_lower_limit: float = 1) -> pd.DataFrame:
        if (lib_id, mass_ratio_lower_limit) in MlimsLibraryService.__get_unknown_ms2__:
            return MlimsLibraryService.__get_unknown_ms2__[(lib_id, mass_ratio_lower_limit)]
        query = """
            select
                cle.lib_id,
                cle.comp_id,
                parent.mass parent_mass,
                parent.name parent_name,
                child.parent_id id,
                child.mass mass,
                child.mass_window,
                child.mass_ratio,
                cle.confidence,
                cle.putative_details details, 
                cle.ri
            from
                limsuser.chro_lib_entry cle,
                limsuser.chro_lib_masses parent,
                limsuser.chro_lib_masses child
    
            where child.parent_id = parent.chro_lib_masses_id
                and cle.lib_id = $LIB_ID
                and parent.lib_id = $LIB_ID
                and child.lib_id = $LIB_ID
                and cle.comp_id = child.comp_id
                and cle.comp_id = parent.comp_id
                and cle.ri > 0 
                and cle.confidence = 100 -- and cle.confidence > 0  -- 1/29/18
                and cle.rt > 0 -- 1/29/18
                and cle.RETIRED_FOR_COMP_ID is null
                and parent.name in ('m+H', 'm-H')
                and parent.quant_mass in ('Y','y','1','2')
        """
        with EasyQuery(**self.mlims_credentials) as eq:
            df = eq.easy_query(query, LIB_ID=lib_id)
        df = df[df.mass_ratio >= mass_ratio_lower_limit]

        df["name"] = df.comp_id.apply(self.get_name_dict().get, "")
        df = df[~df.name.isnull()]
        df = df[df.name.str.match(r"^X-\d{5}$")]
        df.details = df.details.apply(clean_details)
        df = df.reset_index(drop=True)
        MlimsLibraryService.__get_unknown_ms2__[(lib_id, mass_ratio_lower_limit)] = df
        return df

    def get_ms1(self, comp_id, lib_id):
        """
        Gets MS1 ions for given compound and library.  This is used to elucidate library X-compounds.
        @return: dataframe with one row per ion
        @rtype: dataframe
        """
        with EasyQuery(**self.mlims_credentials) as eq:
            return eq.easy_query(
                f"""
                    select * from chro_lib_masses
                    where comp_id = {comp_id} and lib_id={lib_id} and parent_id is null
                    order by mass"""
            )

    def add_filenames_to_peaks(self, peak_df):
        """
        Adds raw file names to peak_df to help elucidators load up peak in tools.
        These should have been and are now added in discovery.
        @param peak_df: The peak dataframe from discovery.
        @type peak_df: DataFrame
        @return: None.  Column added to dataframe
        @rtype: NoneType
        """
        # These filenames should hove been attached in discovery...
        with EasyQuery(db="mlims") as eq:
            df = eq.easy_query(
                "select task_id, filename from worksheet_entry where task_id in ($BATCH)",
                BATCH=peak_df.task_id,
            )
        d = dict(zip(df.task_id, df.filename))  # I guess this could be a merge.
        peak_df["raw_file"] = peak_df.task_id.apply(d.get)

    def cache_library(self, cache_dir="."):
        """
        Call this to cache stuff for CachedMlimsLibraryService.
        #TODO: I guess these should be dated.
        @param cache_dir: directory to hold
        @type cache_dir: str
        @return: None
        @rtype: NoneType
        """
        for lib_id in (209, 400, 402):
            # TODO sort out all the redundancy of these methods.
            self.get_named_ms2(lib_id).to_csv(f"{cache_dir}/getMS2_{lib_id}.csv")
            self.get_unnamed_ms2(lib_id).to_csv(f"{cache_dir}/get_unnamed_ms2_{lib_id}.csv")
            self.get_known_ms2(lib_id).to_csv(f"{cache_dir}/get_known_ms2_{lib_id}.csv")
            self.get_unknown_ms2(lib_id).to_csv(f"{cache_dir}/get_unknown_ms2_{lib_id}.csv")
            self.get_library_ions(lib_id).to_csv(f"{cache_dir}/get_unknown_ms2_{lib_id}.csv")


def clean_details(details):
    """
    Clean up the "putative details" and "comments" to eliminate line feeds.
    @param details: string containing details
    @type details: str
    @return: cleaned-up string.
    @rtype: str
    """
    details = str(details)
    if details == "nan":
        details = ""
    else:
        details = details.replace("\n", " | ").replace("\r", "")
        # if '|' in details: print(details[:-3])
        if "|" in details[-3:]:
            details = details[:-3]
    return details
