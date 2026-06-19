"""RDKit helpers for SMILES manipulation, depiction, and Excel export.

Canonicalization, fragmentation, formula/mass derivation, and molecule images
for elucidation reports and table-building scripts in spec2mol v21 inference.
"""

import base64
import logging
from functools import reduce
from io import BytesIO
from operator import add

import numpy as np
import pandas as pd
import rdkit
import regex as re
from elucidation.chemFacts import fmla_mass
from rdkit.Chem import (
    MolToSmiles,
    MolFromSmiles,
    FragmentOnSomeBonds,
    CanonSmiles,
)
from rdkit.Chem.AllChem import ReactionFromSmarts
from rdkit.Chem.Crippen import MolMR, MolLogP
from rdkit.Chem.Draw import MolToImage, rdMolDraw2D
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt

graphics_cache = dict()  # Caching up to 5000 images, then tossing them all.


def to_excel_with_molecules(df, xlsx_name, mol_cols, sheet_name="Sheet1", size=(300, 300)):
    """
    Writes a dataframe into an excel spreadsheet, adding graphic of molecule to specified columns.
    @param df: dataframe to write
    @type df: dataframe
    @param xlsx_name: name of output excel file.  Must end ".xlsx"
    @type xlsx_name: str
    @param mol_cols: a column or list of columns with smiles strings to be expanded.
    If any element of the columns can't be parsed as as SMILES string, it will be left alone.
    @type mol_cols: str or list of str (actually, any iterable of str)
    @param sheet_name: The name of the sheet in the workbook.
    @type sheet_name: str
    @param size: (width,height) of images
    @type size: tuple (could be list)
    @return: None
    @rtype: NoneType
    """
    df = df.copy()

    if type(mol_cols) is str:
        mol_cols = [mol_cols]
    mol_cols = list(mol_cols)

    writer = pd.ExcelWriter(xlsx_name, engine="xlsxwriter")
    df.to_excel(writer, sheet_name=sheet_name, index=False)
    workbook = writer.book
    worksheet = writer.sheets[sheet_name]

    mol_cols = [i for i, col in enumerate(df.columns) if col in mol_cols]
    df.columns = list(range(len(df.columns)))

    for c in mol_cols:
        worksheet.set_column(c, c, width=size[0] / 6)
        for rr, smiles in enumerate(df[c]):
            try:
                r = rr + 1  # row 0 of df will be in row 1 of spreadsheet.
                image_data = BytesIO()
                img = MolToImage(MolFromSmiles(smiles), size=size)
                img.save(image_data, format="PNG")
                worksheet.set_row(r, height=size[1])
                worksheet.insert_image(r, c, "f", {"image_data": image_data})
            except:  # Some are expected to fail, e.g., None when there is no matching molecule.
                pass

    workbook.close()
    return


def modify_smiles(reaction_smarts, smiles):
    """
    Takes the molecule represented by the SMILES and applies the modification specified by the SMARTS whereever
    possible to get a list of SMILES of modified molecules.
    @param smiles: SMILES string to be modified.
    @type smiles: str
    @param reaction_smarts: SMARTS string for the modification. The attachment atom should be specified with [X:1],
    e.g. the SMARTS for sulfation is:  "[S:1](=O)(=O)O"
    @type reaction_smarts: str
    @return: list of SMILES for modified molecules.
    @rtype: list of strings.
    """
    rxn = ReactionFromSmarts(reaction_smarts)
    try:
        new_mols = rxn.RunReactants((MolFromSmiles(smiles),))
        new_smiles = {MolToSmiles(m[0]) for m in new_mols}
        return new_smiles
    except ValueError as e:
        logging.warning(str(e))
        return set()


def fragment_smiles(s, lenient=True):
    """
    Given a SMILES string, breaks it at all single bonds not in cycles, and returns a list of fragments.
    @param s: SMILES string
    @type s: str
    @return: list of SMILES strings for fragments
    @rtype: list
    """
    try:
        m = MolFromSmiles(s)
        idx = [
            b.GetIdx() for b in m.GetBonds() if b.GetBondType() == rdkit.Chem.rdchem.BondType.SINGLE
        ]
        frags = [FragmentOnSomeBonds(m, (i,)) for i in idx]
        smiles_fragmentation = [MolToSmiles(f[0]) for f in frags]
        smiles_frags = reduce(add, [sf.split(".") for sf in smiles_fragmentation if "." in sf], [])
        smiles_frags = [re.sub(r"\[\d+\*\]", "*", sf) for sf in smiles_frags]
        smiles_frags = sorted(set(smiles_frags), key=len)
        return smiles_frags
    except Exception as e:
        logging.warning(f"Couldn't fragment {s}")
        logging.warning(str(e))
        if lenient:
            return []
        raise


def glue_smiles_fragments(s1, s2):
    """
    Given two SMILES fragment, glues them together to form a SMILES string for a complete molecule.
    @param s1: SMILES string for a fragment.  Should have exactly 1 '*'  (connection point)
    @type s1: str
    @param s2: SMILES string for a fragment.  Should have exactly 1 '*'  (connection point)
    @type s2: str
    @return: SMILES string for new molecule
    @rtype: str
    """
    # There's no doubt a way to do this without using Astatine,
    # "the rarest naturally occurring element in the Earth's crust", but I couldn't figure out the reaction SMARTS.
    rxn = ReactionFromSmarts("[*:1]-[At:3] . [*:2]-[At:4] >> [*:1]-[*:2] . [At:3][At:4]")
    mol1 = MolFromSmiles(s1.replace("*", "[At]"))
    mol2 = MolFromSmiles(s2.replace("*", "[At]"))
    new_mol = rxn.RunReactants((mol1, mol2))[0][0]
    return MolToSmiles(new_mol)


def get_mol_for_smiles(smiles, lenient=True):
    # If we want, we can pass an rdkit Molecule instead of a string.
    # The assumption is that it's either a SMILES or a Molecule.
    if type(smiles) is not str:
        return smiles
    smiles = smiles.strip(" \n\r")
    mol = MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"rdkit couldn't parse {smiles}")
    return mol


def smiles_mass(smiles, lenient=True):  # returns the NEUTRAL mass in case of a fragment.
    # This was probably a mistake... a fragment has an electron to share... I think.
    withstars = smiles
    smiles = smiles.replace("*", "[H]")
    try:
        mol = get_mol_for_smiles(smiles, lenient=False)
        mass = CalcExactMolWt(mol)
        n_stars = (len(smiles) - len(withstars)) / 2
        result = mass - n_stars * fmla_mass("H")
        return result
    except Exception as e:
        if lenient:
            logging.warning(str(e))
            return 0
        raise e


def fmla_for_smiles(smiles, lenient=True):
    try:
        # add a soupçon of France.
        smiles = smiles.replace("*", "[Fr]")
        mol = get_mol_for_smiles(smiles, lenient=lenient)
        fmla = CalcMolFormula(mol)
        fmla = fmla.replace("Fr", "")
        return fmla
    except Exception as e:
        logging.warning(str(e))
        if lenient:
            return ""
        raise e


def get_substructure_matches(smiles, smarts, lenient=True):
    try:
        smiles = get_mol_for_smiles(smiles, lenient=lenient)
        smarts = get_mol_for_smiles(smarts, lenient=lenient)
        return smiles.GetSubstructureMatches(smarts)
    except Exception as e:
        if lenient:
            return []
        raise e


def draw_smiles_image_for_html(m, aspect=1, area_per_atom=50 * 50):  # golden ratio
    global graphics_cache
    if len(graphics_cache) > 5000:
        logging.info("**************************************** purging graphics cache")
        graphics_cache = dict()
    if type(m) is str:
        smiles = m
        if (smiles, aspect, area_per_atom) in graphics_cache:
            return graphics_cache[(smiles, aspect, area_per_atom)]
        mol = get_mol_for_smiles(m)
    else:
        mol = m
        smiles = MolToSmiles(m, isomericSmiles=False)

    # scale drawing to ~size of molecule.  Br and Cl get a little extra! :)
    atoms = len(re.sub(r"[^A-Za-z*]", "", smiles))
    area = area_per_atom * atoms
    scale = np.sqrt(area)
    width = int(scale * aspect)
    height = int(scale * 1 / aspect)

    try:
        d = rdMolDraw2D.MolDraw2DCairo(width, height)  # or MolDraw2DCairo to get PNGs
        rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
        txt = d.GetDrawingText()
        encoding = base64.b64encode(txt).decode("utf-8").replace("\n", "")
        graphics_cache[(smiles, aspect, area_per_atom)] = encoding
        return encoding
    except Exception as e:
        logging.warning(str(e))
        graphics_cache[(smiles, aspect, area_per_atom)] = ""
        return ""


def canonicalize_smiles(smiles, lenient=True, chiral=False):
    """
    Returns canonical form of a SMILES string. Use rdkit
    """
    try:
        canonical = CanonSmiles(smiles, useChiral=chiral)
        return canonical
    except Exception as e:
        if lenient:
            return smiles
        raise e


def logp(smiles, clip=True):
    """
    use rdkit to compute logp of molecule specified by SMILES.
    :param smiles: SMILES string of molecule.
    :param clip: If True, clip large and small outliers.
    :return: logp value, a float.
    """
    m = get_mol_for_smiles(smiles)
    if m is None:
        return np.nan

    x = MolLogP(m)
    if not clip:
        return x
    if x < -7.5:
        return -7.5
    if x > 15:
        return 15
    return x


def mr(smiles):
    """
    MR ("molecular reflectivity"?)
    :param smiles: SMILES string
    :return: a mysterious number we are not using right now
    """
    m = get_mol_for_smiles(smiles)
    if m is None:
        return np.nan
    return MolMR(m)
