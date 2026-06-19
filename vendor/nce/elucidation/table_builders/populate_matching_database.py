"""Rebuild the Postgres smiles_db compound-matching database from source files.

CLI and function library for dropping and repopulating smiles_db from ZINC,
HMDB, MLIMS, LipidMaps, and other public structure sources. Seeds neutral-mass
lookups used by PostgresSmilesDb during elucidation.
"""


import os
import subprocess
import sys
import urllib.request
from argparse import Namespace

import numpy as np
import pandas as pd
import regex as re
from bs4 import BeautifulSoup
from bs4.element import NavigableString
from elucidation.aminoAcids import (
    build_dipeptide_namer,
    build_tripeptide_namer,
    build_tetrapeptide_namer,
)
from elucidation.library_services.mlimsLibraryService import MlimsLibraryService as LibraryService

# from pprint import pprint
from elucidation.rdkitUtils import (
    canonicalize_smiles as canonicalize,
    smiles_mass,
    fmla_for_smiles,
)
from elucidation.conveniences import get_postgres_password, qw, pickle_load
from mdsutils.easyQuery import EasyQuery
from sqlalchemy import create_engine

sql_string = "postgresql+psycopg2://compounddiscoveryservice:{}@ra-pgsql.metabolon.com:5432/compounddiscoveryservice"


def recreate_database():
    """
    drop the table and recreate empty table.
    :return: None
    :rtype: None
    """
    if input("Really delete the database? Type YES: ") != "YES":
        sys.exit("Aborted by user")
    conn.execute("drop index if exists smiles_db_mass_index")
    conn.execute("drop table if exists smiles_db")
    conn.execute(
        """    
    CREATE TABLE SMILES_DB (
    source  text,
    name  text,
    exact_mass  double precision,
    formula  text,
    smiles  text,
    id text,
    source_record text)"""
    )


def create_mass_index():
    """
    indexes the database on exact mass.  It is faster to run this at the end rather
    than maintaining the index through a million insertions.
    :return: None
    :rtype: NoneType
    """
    conn.execute("create index smiles_db_mass_index on smiles_db (exact_mass)")


def add_df_to_smiles_db(source, df, block_size=1000):
    """
    Adds all the items in df to the database, tagged as coming from the
    named source.
    :param source: source tag, e.g. "HMBD"
    :type source: str
    :param df: Dataframe of items to add.
    Required columns (or at least expected):
    - name, a string, is the compound name that wil be displayed in reports, etc.
    - exact_mass, float, is ... the exact match.  MANDATORY
    - formula, str, is the molecular formula
    - smiles, canonical SMILES representation.  MANDATORY
    - id, str, an ID, such as HMDB1234567
    - source_record, a dict converted to a string.  This (generally) contains ALL
      the fields in the source with their original names. To use this column, you'd
      simply call eval(source_record) to recover the dictionary.
    :type df: DataFrame
    :param block_size: number of records to be inserted by each query.
    :type block_size: int
    :return: The orginal DataFrame with the source added.  (Not typically used.)
    :rtype: DataFrame
    """
    print(source, len(df), "records")
    df["source"] = [source] * len(df)
    dfs = np.array_split(df, 1 + (int(len(df) / block_size)))
    for df in dfs:
        print(".", end="")
        df.to_sql("smiles_db", con=eng, index=False, if_exists="append")
    print("Done")
    return df


def insert_zinc(file="ZINC_DB.csv"):
    """
    All of the insertion routines look about the same, and it looks like it should be
    pretty easy to abstract and eliminate similar code, but each database also has
    its peculiarities, so I didn't try.
    :param file: location of ZINC_DB.csv download spreadsheet.  It seems like it would
    be difficult to download it programmatically. I'm currently using the "world" subset,
    which includes about 5000 drugs approved in any "major jurisdiction".
    :type file: str
    :return: None, but could return the DataFrame, I guess
    :rtype: NoneType
    """
    df = pd.read_csv(file)
    rows = []
    for r in df.iterrows():
        d = dict(r[1])
        ns = Namespace(**d)
        try:
            exact_mass = smiles_mass(ns.smiles)
            assert 50 < exact_mass < 1000
            row = dict(
                source_record=str(d),
                name=ns.preferred_name,
                smiles=canonicalize(ns.smiles),
                id=ns.zinc_id,
                formula=fmla_for_smiles(ns.smiles),
                exact_mass=exact_mass,
            )
            rows.append(row)
        except Exception as e:
            pass
    add_df_to_smiles_db("ZINC", pd.DataFrame(rows))


def insert_peptides():
    for namer in [
        build_dipeptide_namer(pyr=True, ggl=True, cyclo=True),
        build_tripeptide_namer(pyr=True, ggl=True),
        build_tetrapeptide_namer(pyr=True, ggl=True),
    ]:
        df = pd.DataFrame({"name": list(namer.values()), "smiles": list(namer.keys())})
        df["id"] = ["pep"] * len(df)
        df["formula"] = df.smiles.apply(fmla_for_smiles)
        df["exact_mass"] = df.smiles.apply(smiles_mass)
        df["source_record"] = ["{}"] * len(df)
        add_df_to_smiles_db("peptides", df)


def insert_mlims():
    conn.execute("delete from smiles_db where source='MLIMS'")
    ls = LibraryService.get_instance()
    name_dict = ls.get_name_dict()

    def condense_lib_ids(t):
        libs = ",".join(sorted(map(str, t.lib_id)))
        t0 = t.head(1).copy()
        t0["mlims_libraries"] = [libs]
        # print(t0)
        return t0

    with EasyQuery(db="mlims") as eq:
        mlims_df = eq.easy_query(
            """select lib_id ,comp_id
        id, cle.chemical_id mlims_chemical_id, smiles, mol_formula formula,
          (putative_details || comments) thoughts,
          iupac_name, exact_mass, cas, inchi, inchikey
         from chro_lib_entry cle join chemicals ch on cle.chemical_id = ch.chemical_id
         where smiles is not null and lib_id in (209,400,402,305) 
         and confidence>0 and ri>0 
         and retired_for_comp_id is null"""
        )
    mlims_df = mlims_df.groupby("id").apply(condense_lib_ids).reset_index(drop=True)

    mlims_df["name"] = mlims_df.id.apply(name_dict.get)
    mlims_df = mlims_df[qw("name exact_mass smiles id")]
    mlims_df.smiles = mlims_df.smiles.apply(canonicalize)
    mlims_df["formula"] = mlims_df.smiles.apply(fmla_for_smiles)
    mlims_df["source_record"] = ["{}"] * len(mlims_df)
    add_df_to_smiles_db("MLIMS", mlims_df)
    return mlims_df
    pass


def run_biotransformer_cyp450(infile="ZINC_DB.csv"):
    """
    This guy takes all the compounds from the file and runs biotransformer on the them
    with reactions for CYP450 enzymes, which are important metabolyzers of drugs.
    This was pulled out of a notebook, and might need some fine tuning.  It wasn't rerun
    because it takes several hours, and I really don't imagine it will be necessary to
    rerun it in the future.  biotransformer must be installed, and the input file must be
    as biotransformer expects; for this, see the biotransformer docs.
    We were given permission to use biotransformer by its creator
    Yannick Djoumbou Feunang <djoumbou@ualberta.ca> in May 2020.  However, we do NOT
    have permission to use the EnviPath transformations, which are owned by someone else.
    :param infile:  Contains the compounds to be transformed.
    :type infile: CSV file
    :return: A dataframe containing all the in silico metabolites.
    :rtype: DataFrame
    """
    zinc_df = pd.read_csv(infile)
    zinc_df["exact_mass"] = zinc_df.smiles.apply(smiles_mass)
    zinc_df = zinc_df[(50 <= zinc_df.exact_mass) & (zinc_df.exact_mass < 1000)]
    save_cwd = os.getcwd()
    os.chdir("biotransformerjar/")

    def transform_smiles(smiles):
        try:
            os.remove("foo.csv")
        except Exception:
            pass
        try:
            cmd = f'java -jar biotransformer-1.1.5.jar -k pred -b cyp450 -ismi "{smiles}" -s 1 -a -ocsv foo.csv'
            out = subprocess.run(qw(cmd), capture_output=True)
        except Exception as e:
            print(f"Couldn't transform {smiles}: {e}")
            return None
        try:
            df = pd.read_csv("foo.csv")
            # print(f"{len(df)} results for {smiles}")
            return df
        except Exception as e:
            # print()
            # print(f"no results for {smiles}")
            return None

    results = []
    for i, t in enumerate(zinc_df.itertuples()):
        result = transform_smiles(t.smiles)
        if result is not None:
            results.append((t, result))
            print(f"{i+1}/{len(zinc_df)}: {len(result)}", end="")

    dfs = []
    for t, df in results:
        df["precursor_preferred_name"] = [t.preferred_name] * len(df)
        df["precursor_zinc_id"] = [t.zinc_id] * len(df)
        dfs.append(df)
    big_df = pd.concat(dfs, ignore_index=True)
    os.chdir(save_cwd)
    big_df.to_csv("btf_zinc.csv", index=False)
    return big_df


def insert_cyp450_zinc(file="btf_zinc.csv"):
    df = pd.read_csv(file)
    rows = []
    for r in df.iterrows():
        d = dict(r[1])
        ns = Namespace(**d)
        smiles = ns.SMILES
        exact_mass = smiles_mass(smiles)
        name = f"{d['Reaction ID']}({ns.precursor_preferred_name})"
        id = f"{d['Reaction ID']}({ns.precursor_zinc_id})"
        try:
            assert 50 < exact_mass < 1000
            row = dict(
                source_record=str(d),
                name=name,
                smiles=canonicalize(smiles),
                id=id,
                formula=fmla_for_smiles(smiles),
                exact_mass=smiles_mass(smiles),
            )
            rows.append(row)
        except Exception as e:
            pass
    df = pd.DataFrame(rows)
    add_df_to_smiles_db("CYP450(ZINC)", df)
    return df


def insert_hmdb(file="chemont_prediction_training.pkl"):
    d = pickle_load(file)
    df = d["hmdb_summary_df"]
    rows = []
    for r in df.iterrows():
        d = dict(r[1])
        ns = Namespace(**d)
        smiles = ns.smiles
        try:
            assert 50 < ns.exact_mass < 1000
            row = dict(
                source_record=str(d),
                name=ns.name,
                smiles=canonicalize(smiles),
                id=ns.hmdb_id,
                formula=fmla_for_smiles(smiles),
                exact_mass=smiles_mass(smiles),
            )
            rows.append(row)
        except Exception as e:
            pass
    add_df_to_smiles_db("HMDB", pd.DataFrame(rows))
    pass


def insert_blood_exposeome(file="BloodExposome_v_1.0.xlsx"):
    df = pd.read_excel(file, header=2)
    df["Compound Name"] = df["Compound Name"].apply(str)
    df.columns = [re.sub(r"[ -]+", "_", k) for k in df.columns]
    rows = []
    for r in df.iterrows():
        d = dict(r[1])
        ns = Namespace(**d)
        smiles = ns.CanonicalSMILES
        try:
            assert 50 < ns.ExactMass < 1000
            assert ns.Charge == 0
            row = dict(
                source_record=str(d),
                name=ns.Compound_Name,
                smiles=canonicalize(smiles),
                id=f"PubChem:{ns.PubChem_CID}",
                formula=fmla_for_smiles(smiles),
                exact_mass=smiles_mass(smiles),
            )
            rows.append(row)
        except Exception as e:
            pass
    add_df_to_smiles_db("BloodExposome", pd.DataFrame(rows))


def scrape_pmdb():
    """
    Scrapes the PMDB website using BeautifulSoup.
    After querying for a list of compound, the sit has to be hit again for every compound.
    It's only about 1000 compounds, so it doesn't take a super-long time.
    It doesn't look like there is any development going on on this site.
    :return:
    :rtype:
    """
    # Get all the PMDB ids for the entire database.
    html_doc = urllib.request.urlopen("http://scbt.sastra.edu/pmdb/ebrowse.php")
    soup = BeautifulSoup(html_doc, "html.parser")
    lines = str(soup).replace("\r", "").split("\n")
    pmdb_ids = [
        re.search(r"(pmdb\d{5})", s).group(1)
        for s in str(soup).split("\n")
        if re.search(r"(pmdb\d{5})", s)
    ]

    # Query ID by ID and build a list of beautiful soups.
    soups = []
    for id in pmdb_ids:
        html_doc = urllib.request.urlopen(f"http://scbt.sastra.edu/pmdb/edetail.php?is={id}")
        soups.append(BeautifulSoup(html_doc, "html.parser"))

    # Fields of the soup we care about.
    fields = set(
        (
            "PMDB ID:|COMMON NAME:|CAS REGISTRY NUMBER:|MOLECULAR FORMULA:|MOLECULAR WEIGHT:|IUPAC NAME:|"
            + "SMILES(Isomeric):|SMILES(Canonical):|KEGG ID:|CHEBI ID:|PUBCHEM COMPOUND ID:|PUBCHEM SUBSTANCE ID"
        ).split("|")
    )

    # Go through all the soups and build a list of dictionaries.
    rows = []
    for soup in soups:
        ps = [s.contents for s in soup.find_all(id="tablebig")]
        d0 = {}
        key = None
        for p in ps:
            if p[0] == "MOL FILE:":
                break
            if p[0] in fields:
                key = p[0]
                if key[-1] == ":":
                    key = key[:-1]  # drop colon
            else:
                s = "".join([t if type(t) is NavigableString else str(t.contents[0]) for t in p])
                val = re.sub(r"[\r\n\t]", "", s)
                d0[key] = val
                key = None
        del d0[None]
        if "SMILES(Canonical)" not in d0:
            continue
        smiles = canonicalize(d0["SMILES(Canonical)"].strip())
        d = dict(
            id=d0["PMDB ID"].strip(),
            name=d0["COMMON NAME"].strip(r"\r\t "),
            smiles=smiles,
            formula=fmla_for_smiles(smiles),
            exact_mass=smiles_mass(smiles),
            source_record=str(d0),
        )
        rows.append(d)
    # Turn it all into a dataframe and save.
    df = pd.DataFrame(rows)
    df.to_csv("pmdb_scraping_results.csv", index=False)
    return df


def insert_pmdb(file="pmdb_scraping_results.csv"):
    add_df_to_smiles_db("PMDB", pd.read_csv(file))


def scrape_plantcyc():
    """
    Scrapes the PlantCyc website using BeautifulSoup.
    After querying for a list of species, the site has to be hit again for every species.
    There is obviously a great deal of overlap among species.
    I coallesced the species for each compound into a list
    It's <5000 compounds, so it takes only moderately long time.
    This site appears to get on-going but infrequent updates.
    :return: None; dumps result to a file.
    :rtype: NoneType
    """
    html_doc = urllib.request.urlopen(
        "ftp://ftp.plantcyc.org/pmn/Pathways/Data_dumps/PMN14_January2020/compounds"
    )
    s = BeautifulSoup(html_doc, "html.parser")
    files = re.findall("\w+compounds[.]\d{8}", str(s))

    dfs = []

    os.makedirs("plantcyc", mode=0o777, exist_ok=True)
    for f in files:
        url = f"ftp://ftp.plantcyc.org/pmn/Pathways/Data_dumps/PMN14_January2020/compounds/{f}"
        species = re.match(r"([A-Z]\w+)_compounds", f).group(1)
        print(species)
        html_doc = urllib.request.urlopen(url)
        s = BeautifulSoup(html_doc, "html.parser")
        csv = f"plantcyc/{species}.csv"
        with open(csv, "w", encoding="utf-8") as fp:
            print(s, file=fp)
        df = pd.read_csv(csv, sep="\t")
        df["Species"] = [species] * len(df)
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True)

    def coallesce_species(tab):
        d = {col: tab[col].apply(str).max() for col in tab.columns}
        d["Species"] = sorted(set(tab.Species))
        return d

    df = pd.DataFrame(list(df.groupby("Compound_id").apply(coallesce_species)))
    df.to_csv("plantcyc_scrape.csv", index=False)


def insert_plantcyc(file="plantcyc_scrape.csv"):
    df = pd.read_csv(file, low_memory=False)
    print(df.columns)
    # PlantCyc has a lot of bad SMILES.
    rows = []
    for r in df.iterrows():
        d = dict(r[1])
        ns = Namespace(**d)
        smiles = ns.Smiles
        try:
            mass = smiles_mass(smiles)
            assert 50 < mass < 1000
            row = dict(
                source_record=str(d),
                name=ns.Compound_common_name,
                smiles=canonicalize(smiles),
                id=ns.Compound_id,
                formula=fmla_for_smiles(smiles),
                exact_mass=smiles_mass(smiles),
            )
            rows.append(row)
        except Exception as e:
            print(e)
            pass
    add_df_to_smiles_db("PlantCyc", pd.DataFrame(rows))


def insert_foodb(file="FooDB_Compound.csv"):
    df = pd.read_csv(file, low_memory=False)
    # The column headings are totally jumbled.
    df.columns = qw(
        """id public_id name state annotation_quality description cas_number moldb_smiles
    moldb_inchi moldb_mono_mass moldb_inchikey moldb_iupac kingdom superklass klass subklass"""
    )
    rows = []
    for r in df.iterrows():
        d = dict(r[1])
        ns = Namespace(**d)
        smiles = ns.moldb_smiles
        try:
            assert 50 < ns.moldb_mono_mass < 1000
            row = dict(
                source_record=str(d),
                name=ns.name,
                smiles=canonicalize(smiles),
                id=ns.public_id,
                formula=fmla_for_smiles(smiles),
                exact_mass=smiles_mass(smiles),
            )
            rows.append(row)
        except Exception as e:
            pass
    df = pd.DataFrame(rows)
    add_df_to_smiles_db("FooDB", df)
    return df


def insert_lipidmaps(file="structures.sdf"):
    with open(file, encoding="utf-8") as f:
        records = f.read().split("$$$$")[:-1]
    rows = []
    for record in records:
        try:
            d = dict([line.strip().split("\n") for line in record.split("> ")[1:]])
            d = {k[1:-1]: v for k, v in d.items()}
            ns = Namespace(**d)
            smiles = ns.SMILES
            assert 50 < float(ns.EXACT_MASS) < 1000
            row = dict(
                source_record=str(d),
                name=ns.NAME,
                smiles=canonicalize(smiles),
                id=ns.LM_ID,
                formula=fmla_for_smiles(smiles),
                exact_mass=smiles_mass(smiles),
            )
            rows.append(row)
        except Exception as e:
            pass
    df = pd.DataFrame(rows)
    add_df_to_smiles_db("LipidMaps", df)
    return df


def insert_binding_db(file="binding_db_good.csv"):
    # Looking for the code that turns BindingDB_All.csv in binding_db_good.csv
    df = pd.read_csv(file)
    rows = []
    for i, r in enumerate(df.iterrows()):
        if i % 1000 == 0:
            print(i, len(rows))
        d = dict(r[1])
        ns = Namespace(**d)
        smiles = ns.canonical_smiles
        try:
            assert 50 < ns.exact_mass < 1000
            row = dict(
                source_record=str(d),
                name=ns.name,
                smiles=canonicalize(smiles),
                id=f"BNDDB_{ns.bindingdb_id}",
                formula=fmla_for_smiles(smiles),
                exact_mass=smiles_mass(smiles),
            )
            rows.append(row)
        except Exception as e:
            print(e)
            pass
    df = pd.DataFrame(rows)
    add_df_to_smiles_db("BindingDB", df)
    return df


def insert_metacyc(file="metacyc_compounds.txt"):
    with open("metacyc_compounds.txt") as fp:
        s = fp.read()
    # Throw away comment lines.
    lines = s.split("\n")
    lines = [line for line in lines if not (len(line) == 0 or line.startswith("#"))]
    # Rejoin into long string, then divide into records for single compounds
    records = re.split("//", "\n".join(lines))
    # Process records to build database records.
    rows = []
    for record in records:
        try:
            d = dict(re.findall(r"(.+?) - (.+?)\n", record))
            smiles = d["SMILES"]
            exact_mass = smiles_mass(smiles)
            assert 50 < exact_mass < 1000
            row = dict(
                source_record=str(d),
                name=d["COMMON-NAME"],
                smiles=canonicalize(smiles),
                id=d["UNIQUE-ID"],
                formula=fmla_for_smiles(smiles),
                exact_mass=exact_mass,
            )
            rows.append(row)
        except Exception as e:
            pass
    df = pd.DataFrame(rows)
    add_df_to_smiles_db("MetaCyc", df)
    return df
    pass


def insert_chebi():  # TODO
    pass


def insert_refmet(file="refmet_with_smiles.csv"):
    df = pd.read_csv(file)
    m = df.smiles.apply(smiles_mass)
    df = df[(50 < m) & (m < 1000)]

    rows = []
    for i, r in enumerate(df.iterrows()):
        if i % 1000 == 0:
            print(i, len(rows))
        d = dict(r[1])
        ns = Namespace(**d)
        try:
            exact_mass = smiles_mass(ns.smiles)
            row = dict(
                source_record=str(d),
                name=ns.refmet_name,
                smiles=canonicalize(ns.smiles),
                id=f"PubChem:{ns.pubchem_cid}",
                formula=ns.formula,
                exact_mass=exact_mass,
            )
            rows.append(row)
        except Exception as e:
            print("insert_refmet", str(e))
            print(ns)
            pass
    df = pd.DataFrame(rows)
    add_df_to_smiles_db("RefMet", df)
    return df


def insert_cip():  # TODO
    pass


if __name__ == "__main__":

    eng = create_engine(sql_string.format(get_postgres_password()))
    conn = eng.connect()

    # recreate_database()
    # insert_zinc()
    # insert_cyp450_zinc()
    # insert_peptides()
    # insert_plantcyc()
    # insert_mlims()
    # insert_hmdb(file="chemont_prediction_training.pkl")
    # insert_blood_exposeome(file="BloodExposome_v_1.0.xlsx")
    # insert_pmdb()
    # insert_foodb()
    # insert_lipidmaps()
    #
    # insert_binding_db()
    # insert_metacyc()
    # insert_chebi()  # TODO
    # insert_refmet()
    # insert_cip()  #TODO; just a couple of hundred
    # create_mass_index()

    for row in conn.execute("select source, count(*) n from smiles_db group by source"):
        print(row)
    print(len(list(conn.execute("select smiles from smiles_db"))), "canonical SMILES")
    print(len(set(conn.execute("select smiles from smiles_db"))), "distinct canonical SMILES")
    # print(len(set(conn.execute("select smiles from smiles_db where source != 'RefMet'"))), "distinct canonical SMILES")
