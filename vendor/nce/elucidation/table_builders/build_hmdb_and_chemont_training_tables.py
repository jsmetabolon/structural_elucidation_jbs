"""Build Chemont training pickle from HMDB and ClassyFire downloads.

Maps HMDB metabolites to MLIMS chemical ids and ClassyFire ontology labels for
use by train_chemont_models.py and ChemontPredictor at inference time.
"""

from collections import defaultdict

import pandas as pd
from defusedxml.ElementTree import parse
from elucidation.rdkitUtils import canonicalize_smiles, smiles_mass
from elucidation.conveniences import pickle_dump
from mdsutils.easyQuery import EasyQuery

xml_file = "hmdb_metabolites.xml"
spreadsheet_file = "hmdb_spreadsheet.csv"
classyfire_annotation_file = (
    "c:/CompoundDiscoveryExploratory/ClassyFire/HMDB_36_classyfire_21_annotations.csv"
)
output_pkl_file = "chemont_prediction_training.pkl"

if __name__ == "__main__":
    if False:  # This will build the hmdb spreadsheet from the download from www.hmdb.ca
        root = parse(xml_file).getroot()
        rows = []
        for hm in root:
            row = {}
            for elt in hm.findall("{http://www.hmdb.ca}accession"):
                row["hmdb_id"] = elt.text

            for elt in hm.findall("{http://www.hmdb.ca}secondary_accessions"):
                row["secondary_ids"] = ";".join(
                    [elt1.text for elt1 in elt.findall("{http://www.hmdb.ca}accession")]
                )

            for elt in hm.findall("{http://www.hmdb.ca}smiles"):
                row["smiles"] = elt.text
            for elt in hm.findall("{http://www.hmdb.ca}name"):
                row["name"] = elt.text
            rows.append(row)
        df = pd.DataFrame(rows)
        df["canonical_smiles"] = df.smiles.apply(canonicalize_smiles)
        df["exact_mass"] = df.canonical_smiles.apply(smiles_mass)
        df.to_csv(spreadsheet_file, index=False)

    df = pd.read_csv(spreadsheet_file)
    df.secondary_ids = df.secondary_ids.apply(str)

    hmdb_to_smiles = defaultdict(str)
    for i, t in enumerate(df.itertuples()):
        for id in [t.hmdb_id] + t.secondary_ids.split(";"):
            hmdb_to_smiles[id] = t.canonical_smiles

    annotations = pd.read_csv(classyfire_annotation_file)
    annotations.columns = "hmdb_id chemont_id chemont_name parent_name".split()

    smiles_to_annotations = defaultdict(list)
    for t in annotations.itertuples():
        smiles_to_annotations[hmdb_to_smiles[t.hmdb_id]].append(t.chemont_id)

    with EasyQuery(db="mlims") as eq:
        chem_df = eq.easy_query(
            "select chemical_id, smiles from chemicals where smiles is not null"
        )
    chem_df["canonical_smiles"] = chem_df.smiles.apply(canonicalize_smiles)
    chem_df["chemont"] = chem_df.canonical_smiles.apply(smiles_to_annotations.get)
    chem_df = chem_df[~chem_df.chemont.isnull()]

    with EasyQuery(db="mlims") as eq:
        cle_df = eq.easy_query(
            """select cle.chemical_id, comp_id, lib_id from chro_lib_entry cle
        join chemicals ch on ch.chemical_id=cle.chemical_id where cle.chemical_id in ($BATCH)
        and lib_id in (209, 305, 402, 400)""",
            BATCH=list(chem_df.chemical_id),
        )
    cle_df = cle_df.merge(chem_df, on="chemical_id")

    d = dict(
        hmdb_download="2021-03-15",
        hmdb_summary_df=df,
        hmdb_to_smiles=hmdb_to_smiles,
        classyfire_version=classyfire_annotation_file,
        classyfire=annotations,
        smiles_to_annotations=smiles_to_annotations,
        annotated_chemical_df=chem_df,
        annotated_chro_lib_entry=cle_df,
    )

    pickle_dump(output_pkl_file, d)
