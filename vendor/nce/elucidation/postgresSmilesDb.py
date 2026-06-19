"""Postgres smiles_db access for neutral-mass compound lookup during elucidation.

Singleton SQLAlchemy client querying the compound-discovery smiles database for
structures near a feature's inferred neutral mass. Populated by table_builders
populate_matching_database.py.
"""

import pandas as pd
from elucidation.binning import plus_5ppm, minus_5ppm
from elucidation.rdkitUtils import canonicalize_smiles, smiles_mass
from elucidation.conveniences import get_postgres_password, leng
from sqlalchemy import create_engine


class PostgresSmilesDb:
    """Singleton reader for the Postgres compound-discovery smiles_db table."""

    __instance__ = None

    @staticmethod
    def get_instance():
        if PostgresSmilesDb.__instance__ is None:
            PostgresSmilesDb.__instance__ = PostgresSmilesDb()
        return PostgresSmilesDb.__instance__

    def __init__(self):
        if PostgresSmilesDb.__instance__ is not None:
            raise Exception("use PostgresSmilesDb.get_instance()")
        self.sql_string = (
            "postgresql+psycopg2://compounddiscoveryservice:"
            + get_postgres_password()
            + "@ra-pgsql.metabolon.com:5432/compounddiscoveryservice"
        )
        self.connection = create_engine(self.sql_string).connect()

    def get_records_for_neutral_mass(self, mass, ppm_mass=None):
        """
        Gets records from postgres at given mass +/- 5ppm.  If ppm_mass is given, ppm is found from
        that mass.  This is needed in the case of modifications: If we think we have X glucuronide and we are searching
        for X, the mass we have measured is X glucuronide, and that determines our accuracy.
        @param mass: Neutral mass at which compounds are to be found.
        @type mass: float
        @param ppm_mass: Mass for computing 5ppm.
        @type ppm_mass: float
        @return: DataFrame of compounds found.
        @rtype: DataFrame
        """
        if ppm_mass is None:
            ppm_mass = mass
        lo = minus_5ppm(ppm_mass) - ppm_mass + mass
        hi = plus_5ppm(ppm_mass) - ppm_mass + mass
        query = f"""select * from smiles_db where {lo} <= exact_mass and exact_mass <= {hi}"""
        return pd.read_sql_query(query, self.connection)

    def get_records_for_smiles(self, smiles):
        """
        Gets compounds with mass near mass of SMILES.
        @param smiles: SMILES to be "matched"
        @type smiles: str
        @return: DataFrame of compounds
        @rtype: DataFrame
        """
        df = self.get_records_for_neutral_mass(smiles_mass(smiles))  # indexed on mass
        if leng(df) == 0:
            return None
        can_smiles = canonicalize_smiles(smiles)
        return df[df.smiles == can_smiles]

    def get_a_smiles_name(self, smiles):
        """
        Gets some name to attach to a SMILES string for reports.
        Needs to be refined.  Removed from use pending fix.  #TODO
        @param smiles: SMILES string of compound
        @type smiles: str
        @return: some name for the compound
        @rtype: str
        """
        smiles = canonicalize_smiles(smiles.replace("*", "[H]"))
        # Since table is not indexed on smiles, we use mass.
        mass = smiles_mass(smiles)
        lo = minus_5ppm(mass)
        hi = plus_5ppm(mass)
        query = f"""select smiles, name from smiles_db where smiles = '{smiles}'
                    and {lo} <= exact_mass and exact_mass <= {hi}"""
        print("get_a_smiles_name", query)
        df = pd.read_sql_query(query, self.connection)
        print("get_a_smiles_name", df)
        if len(df) > 0:
            return min(list(df["name"]), key=len)
        else:
            return ""
