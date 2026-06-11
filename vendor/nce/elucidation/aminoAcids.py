"""
This is not currently used during elucidation per se, but it is used to populate the matching database.
"""
import io
from collections import defaultdict
from functools import reduce
from itertools import product as cartesian, combinations_with_replacement
from operator import add, itemgetter

import pandas as pd
import regex as re

# from pprint import pprint
from elucidation.binning import bin_for_mass
from elucidation.chemFacts import split_fmla, fmla_mass, unsplit_fmla
from elucidation.rdkitUtils import canonicalize_smiles as canonicalize
from elucidation.conveniences import lmap, qw
from elucidation.riPredictor import RiPredictor

# I don't have a need for these right now... but they were hard to dig out.
mlims_amino_acids = """"CHEMICAL_ID","NAME","MOL_FORMULA","EXACT_MASS","SMILES"
"231","Arg","C6H14N4O2","174.111676","N=C(N)NCCCC(N)C(O)=O"
"234","Asp","C4H7NO4","133.037509","NC(C(O)=O)CC(O)=O"
"340","Gly","C2H5NO2","75.032029","NCC(O)=O"
"355","His","C6H9N3O2","155.069477","NC(C(O)=O)CC1=CN=CN1"  
"376","Ile","C6H13NO2","131.094629","NC(C(O)=O)C(C)CC"
"397","Leu","C6H13NO2","131.0946286667","O=C([C@H](CC(C)C)N)O"
"407","Lys","C6H14N2O2","146.105528","NC(C(O)=O)CCCCN"
"415","Met","C5H11NO2S","149.051051","CSCCC(N)C(O)=O"
"460","Phe","C9H11NO2","165.078979","NC(C(O)=O)CC1=CC=CC=C1"
"480","Pro","C5H9NO2","115.063329","OC(C1CCCN1)=O"
"503","Ser","C3H7NO3","105.042594","OC([C@@H](N)CO)=O"
"561","Glu","C5H9NO4","147.053159","O=C(O)[C@@H](N)CCC(O)=O"
"563","Gln","C5H10N2O3","146.0691421977","OC([C@@H](N)CCC(N)=O)=O"
"564","Thr","C4H9NO3","119.058244","N[C@H](C(O)=O)[C@H](O)C"
"565","Trp","C11H12N2O2","204.089878","O=C(O)[C@@H](N)CC1=CNC2=C1C=CC=C2"
"566","Val","C5H11NO2","117.078979","N[C@H](C(O)=O)C(C)C"
"800","Cys","C3H7NO2S","121.019751","SC[C@H](N)C(O)=O"
"811","Ala","C3H7NO2","89.047679","O=C(O)[C@@H](N)C"
"815","Tyr","C9H11NO3","181.073894","N[C@@H](CC1=CC=C(O)C=C1)C(O)=O"
"917","Asn","C4H8N2O3","132.053493","NC(C[C@H](N)C(O)=O)=O"
"""
mlims_aa_table = pd.read_csv(io.StringIO(mlims_amino_acids))

s = """Alanine C3H5NO	Ala	A	71.037114	71.0779	Ala
Arginine C6H12N4O	Arg	R	156.101111	156.1857	Arg
Asparagine C4H6N2O2	Asn	N	114.042927	114.1026	Asn
Aspartic_acid C4H5NO3	Asp	D	115.026943	115.0874	Asp
Cysteine C3H5NOS	Cys	C	103.009185	103.1429	Cys
Glutamic_acid C5H7NO3	Glu	E	129.042593	129.114	Glu
Glutamine C5H8N2O2	Gln	Q	128.058578	128.1292	Gln
Glycine C2H3NO	Gly	G	57.021464	57.0513	Gly
Histidine C6H7N3O	His	H	137.058912	137.1393	His
Isoleucine C6H11NO	Ile	I	113.084064	113.1576	Ile
Leucine C6H11NO	Leu	L	113.084064	113.1576	Leu
Lysine C6H12N2O	Lys	K	128.094963	128.1723	Lys
Methionine C5H9NOS	Met	M	131.040485	131.1961	Met
Phenylalanine C9H9NO	Phe	F	147.068414	147.1739	Phe
Proline C5H7NO	Pro	P	97.052764	97.1152	Pro
Serine C3H5NO2	Ser	S	87.032028	87.0773	Ser
Threonine C4H7NO2	Thr	T	101.047679	101.1039	Thr
Selenocysteine C3H5NOSe	Sec	U	150.95363	150.0379	SeC
Tryptophan C11H10N2O	Trp	W	186.079313	186.2099	Trp
Tyrosine C9H9NO2	Tyr	Y	163.06332	163.1733	Tyr
Valine C5H9NO	Val	V	99.068414	99.1311	Val
Pyroglutamic_acid C5H5NO2 Glp x 111.03202879428 -1 Pyr"""

"""
This is the string above in a dataframe.
"""
amino_acid_table = pd.DataFrame(
    [tuple(re.split(r"\s+", t)[:-1]) for t in s.split("\n")],
    columns=qw("name fmla three one mono_mass avg_mass"),
)
# We don't guess isoleucine in formulas because we can't distinguish it from leucine by mass.
amino_acid_table = amino_acid_table[amino_acid_table.name != "Isoleucine"]


def build_aa_smiles_dict():
    # We do NOT want them in canonical form.  We want to be able to concatenate them.
    side_chains = dict(
        Ala="(C)",
        Arg="(CCCNC(=N)N)",
        Asn="(CC(N)=O)",
        Asp="(CC(=O)O)",
        Cys="(CS)",
        Gln="(CCC(N)=O)",
        Glu="(CCC(=O)O)",
        Gly="",  # implicit H
        His="(Cc1c[nH]cn1)",
        Ile="(C(C)CC)",
        Lys="(CCCCN)",
        Leu="(CC(C)C)",
        Met="(CCSC)",
        Phe="(Cc1ccccc1)",
        Ser="(CO)",
        Thr="(C(O)C)",
        Trp="(Cc1c[nH]c2ccccc12)",
        Tyr="(Cc1ccc(O)cc1)",
        Val="(C(C)C)",
        BADHIS="(Cc1cnc[nH]1)",  # This is the incorrect histidine from PubChem.  It's here to help identify errors.
        BADARG="(CCCN=C(N)N)",
    )
    d = {k: f"NC{v}C(=O)O" for k, v in side_chains.items()}
    d["Pro"] = "N1CCCC1C(=O)O"
    d["Pyr"] = "N1C(=O)CCC1C(=O)O"
    d["gGl"] = "NC(C(=O)O)CCC(=O)O"
    return d


aa_smiles_dict = build_aa_smiles_dict()
twenty_aas = tuple(
    qw("Ala Arg Asn Asp Cys Gln Glu Gly His Ile Lys Leu Met Phe Pro Ser Thr Trp Tyr Val")
)


def build_monopeptide_namer(aas=twenty_aas, pyr=False, ggl=False, bad_his=False, bad_arg=False):
    aas, first = add_misc_aas(aas, pyr=pyr, ggl=ggl, bad_his=bad_his, bad_arg=bad_arg)
    d = {canonicalize(build_oligo_smiles([p])): p for p in first}
    return d


def build_dipeptide_namer(
    aas=twenty_aas, pyr=False, ggl=False, cyclo=True, bad_his=False, bad_arg=False
):
    aas, first = add_misc_aas(aas, pyr=pyr, ggl=ggl, bad_his=bad_his, bad_arg=bad_arg)
    d = {canonicalize(build_oligo_smiles(p)): "-".join(p) for p in cartesian(first, aas)}
    if cyclo:
        d.update(
            {
                canonicalize(build_cyclic_oligo_smiles(p)): f"cyclo({p[0]}-{p[1]})"
                for p in combinations_with_replacement(aas, 2)
            }
        )
        print(len(d))

    return d


def build_tripeptide_namer(aas=twenty_aas, pyr=False, ggl=False, bad_his=False, bad_arg=False):
    aas, first = add_misc_aas(aas, pyr=pyr, ggl=ggl, bad_his=bad_his, bad_arg=bad_arg)
    d = {canonicalize(build_oligo_smiles(p)): "-".join(p) for p in cartesian(first, aas, aas)}
    return d


def add_misc_aas(aas, pyr=False, ggl=False, bad_his=False, bad_arg=False):
    aas = list(aas)
    if bad_his:
        aas += ["BADHIS"]
    if bad_arg:
        aas += ["BADARG"]
    first = aas.copy()
    if pyr:
        first.append("Pyr")
    if ggl:
        first.append("gGl")
        return aas, first


def build_tetrapeptide_namer(aas=twenty_aas, pyr=False, ggl=False, bad_his=False, bad_arg=False):
    aas, first = add_misc_aas(aas, pyr=pyr, ggl=ggl, bad_his=bad_his, bad_arg=bad_arg)
    print(len(list(cartesian(first, aas, aas, aas))))
    d = {canonicalize(build_oligo_smiles(p)): "-".join(p) for p in cartesian(first, aas, aas, aas)}
    return d


def build_polypeptide_dict():
    """
    This function is called once to build a table that will be used in subsequent queries by mass.
    :return: Using amino acid table above, enumerates all possible polypeptide formulas with up to 4 AAs.
                Note that Pyr can only appear once, at N-terminus.
    """
    aa_tuples = sorted(  # 3-letter abbreviation, exact mass, and molecular formula.
        [(t.three, float(t.mono_mass), t.fmla) for t in amino_acid_table.itertuples()]
    )

    def multiples(t):
        """
        :param t: A tuple from aa_tuples.
        :return:  A list of AA multiples: AA, factor, mass, list of fmla repeated factor times (kludge!)
        """
        aa, amu, fmla = t
        if aa == "Pyr":  # only 1 of these.
            return [(aa, 1, amu, [fmla])]
        return [(aa, i, i * amu, [fmla] * i) for i in range(1, 4)]

    all_options = [multiples(t) for t in aa_tuples]

    def f(options):
        """
        This is a recursive function that builds up a list of lists.  Each inner list represents a polypeptide.
        :param options: A list of lists: Each inner list is a list of multiples for an AA.
            We will choose at most one multiple from each list.
        :return: A list of lists, each inner list representing one polypeptide.
        """
        # recursion bottoms out if options=[]; return the list of a single "empty" polypeptide for technical reasons
        if len(options) == 0:
            return [[]]

        # What can we construct ignoring the first option?
        results = f(options[1:])
        new_results = results.copy()

        # Cool; now add in all possibilities for the first option.
        for multiple in options[0]:
            aa, i, amu, fmlas = multiple
            new_results.extend(
                [
                    r + [multiple]
                    for r in results
                    if aa_count(r) + i <= 6 and aa_weight(r) + amu < 1000
                ]
            )
        return new_results

    h2o = fmla_mass("H2O")

    def aa_weight(options):
        """
        Find mass of a polypeptide given in the form of a list of multiples.
        :param options: a list of multiples.  The masses will be extracted from them and summed.
        :return:
        """
        return sum(
            lmap(itemgetter(2), options), h2o
        )  # There is a loss of water when C and N termini bond.

    def aa_count(options):
        return sum(lmap(itemgetter(1), options), 0)

    def aa_fmla(p):
        """
        Creates an AA formula and a molecular formula for a polypeptide expressed as a list of multiples.
        :param p: A proposed polypeptide in the form of a list of multiples.
        :return: A polypeptide formula and a molecular formula for the multiples combined.
        """
        s = "+".join([f"{i}*{aa}" for aa, i, mass, fmla in p[::-1]])
        fmlas = reduce(add, [fmla for aa, i, mass, fmla in p])
        big_fmla = split_fmla("H2O")[1]
        for fmla in fmlas:
            for k, v in split_fmla(fmla)[1].items():
                big_fmla[k] += v
        return unsplit_fmla(0, big_fmla), re.sub(r"1[*]", "", s)

    # Construct a list representing all polypeptides.  Each polypeptide is represented by a list of multiples.
    polypeptides_as_multiples = f(all_options)
    polypeptide_dict = defaultdict(list)
    for p in polypeptides_as_multiples:
        if p:  # Remember that "empty" polypeptide?  Time to ignore it.
            polypeptide_dict[bin_for_mass(aa_weight(p))].append(aa_fmla(p))

    return polypeptide_dict


polypeptide_dict = None


def guess_polypeptide(mass):
    """
    Gather a dataframe of polypeptides for a neutral mass.
    :param mass: neutral mass for which polypeptide formula is sought.
    :return: A dataframe of candidate polypeptides including molecular formula and AA formula.
    """
    global polypeptide_dict
    if polypeptide_dict is None:
        polypeptide_dict = build_polypeptide_dict()
    b = bin_for_mass(mass)
    df = pd.DataFrame(
        sorted(reduce(add, [polypeptide_dict[b + i] for i in (-1, 0, 1)])),
        columns=qw("fmla polypeptide"),
    )

    df["ri_pred_209"] = df.fmla.apply(RiPredictor.predict_fmla_ri_209)
    df["ri_pred_400"] = df.fmla.apply(RiPredictor.predict_fmla_ri_400)

    # TODO construct SMILES strings and add more accurate RI prediction from SMILES.
    return df


def build_oligo_smiles(aa_list):
    return "".join([aa_smiles_dict[aa][:-1] for aa in aa_list]) + "O"


def build_cyclic_oligo_smiles(aa_list):
    straight = build_oligo_smiles(aa_list)
    # form the cyclic.  '9' is safe, because almost no molecules of any sort have 9 cycles.
    return "N9" + straight[1:-5] + "9(=O)"


if __name__ == "__main__":

    # print(guess_polypeptide(396.212118))  # neutral
    # print(guess_polypeptide(475.1834 + proton_mass))  # m-H
    # print(polypeptide_dict)
    df = mlims_aa_table.sort_values("NAME")
    df["my_smiles"] = df.NAME.apply(aa_smiles_dict.get)
    df["c1"] = df.SMILES.apply(canonicalize)
    df["c2"] = df.my_smiles.apply(canonicalize)
    df["same"] = df.c1 == df.c2
    for t in df.itertuples():
        if not t.same:
            print(f"{t.NAME} {t.SMILES} {t.c1}\n     {t.my_smiles} {t.c2}\n")

    print(df)
    print(len(df))
