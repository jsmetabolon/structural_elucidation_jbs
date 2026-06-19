"""Resolve InChIKeys to InChI, RDKit molecules, and SMILES via ChemBL UniChem.

Fallback structure lookup when a compound identifier cannot be resolved from
the local Postgres smiles_db. Used sparingly during elucidation reporting
when only an InChIKey is available.
"""

from chembl_webresource_client.unichem import unichem_client as unichem
from rdkit.Chem import (
    MolToSmiles,
    MolFromInchi,
)


def get_inchi_for_inchikey(inchikey):  # TODO: exception handling
    """
    Hits ChemBL API to get an InChI for an InChIKey
    @param inchikey: An InChIKey, e.g. ADFCQWZHKCXPAJ-UHFFFAOYSA-N  (equol)
    @type inchikey: str
    @return: an InChI, e.g.:
      InChI=1S/C15H14O3/c16-13-4-1-10(2-5-13)12-7-11-3-6-14(17)8-15(11)18-9-12/h1-6,8,12,16-17H,7,9H2
    @rtype: str
    """
    return unichem.inchiFromKey(inchikey)[0]["standardinchi"]


def get_mol_for_inchikey(inchikey, lenient=True):  # TODO: exception handling
    """
    Given an InChIKey, hits chembl to get an inchi, then creates
    @param inchikey: An InChIKey, e.g. ADFCQWZHKCXPAJ-UHFFFAOYSA-N  (equol)
    @type inchikey: str
    @param lenient: How to handle exceptions.
    @type lenient: bool
    @return: molecule as represented by rdkit.
    @rtype: Molecule datatype from rdkit
    """
    return MolFromInchi(get_inchi_for_inchikey(inchikey))


def get_smiles_for_inchikey(inchikey, lenient=True):  # TODO: exception handling
    """
    Given an InChIKey, returns the corresponding SMILES string.
    @param inchikey: An InChIKey, e.g. ADFCQWZHKCXPAJ-UHFFFAOYSA-N  (equol)
    @type inchikey: str
    @param lenient: How to handle exceptions.
    @type lenient: bool
    @return: a SMILES string, e.g.: OC1=CC=C2C(OCC(C3=CC=C(O)C=C3)C2)=C1  (equol again)
    @rtype: str
    """
    return MolToSmiles(get_mol_for_inchikey(inchikey))
