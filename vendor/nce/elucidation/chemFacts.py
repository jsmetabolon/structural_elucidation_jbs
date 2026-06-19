"""Core chemistry constants, periodic-table data, and formula utilities.

Loads isotope masses and provides ion-mode-aware neutralization, formula
parsing, DBE checks, and exact-mass calculations used throughout NCE
elucidation and spec2mol v21 inference.
"""
# TODO turn these dictionaries into functions, and perhaps turn it all into a singleton class.
import logging
from collections import defaultdict

import numpy as np
import pandas as pd
import regex as re
from elucidation.conveniences import qw

POS_MODE_LIBRARIES = (400, 402)
NEG_MODE_LIBRARIES = (209, 305)


POS = "pos"
NEG = "neg"

library_mode = {209: NEG, 305: NEG, 400: POS, 402: POS}

proton_mass = 1.007276466
electron_mass = 0.00054858
neutron_mass = 1.008664

# This is a list of all isotopes of all elements stolen from some unrecorded source on the web.  (There must be many.)
# This is probably built into scipy.
# We have: name, isotope symbol mass in amu, natural abundance.
# True!  We are not interested in most of this!


_periodic_table = """Aluminum	Al(27)	26.981541	100.00
Antimony	Sb(121)	120.903824	57.30
Antimony	Sb(123)	122.904222	42.70
Argon	Ar(36)	35.967546	0.34
Argon	Ar(38)	37.962732	0.063
Argon	Ar(40)	39.962383	99.60
Arsenic	As(75)	74.921596	100.00
Barium	Ba(130)	129.906277	0.11
Barium	Ba(132)	131.905042	0.10
Barium	Ba(134)	133.904490	2.42
Barium	Ba(135)	134.905668	6.59
Barium	Ba(136)	135.904556	7.85
Barium	Ba(137)	136.905816	11.23
Barium	Ba(138)	137.905236	71.70
Beryllium	Be(9)	9.012183	100.00
Bismuth	Bi(209)	208.980388	100.00
Boron	B(10)	10.012938	19.80
Boron	B(11)	11.009305	80.20
Bromine	Br(79)	78.918336	50.69
Bromine	Br(81)	80.916290	49.31
Calcium	Ca(40)	39.962591	96.95
Calcium	Ca(42)	41.958622	0.65
Calcium	Ca(43)	42.958770	0.14
Calcium	Ca(44)	43.955485	2.086
Calcium	Ca(46)	45.953689	0.004
Calcium	Ca(48)	47.952532	0.19
Cadmium	Cd(106)	105.906461	1.25
Cadmium	Cd(108)	107.904186	0.89
Cadmium	Cd(110)	109.903007	12.49
Cadmium	Cd(111)	110.904182	12.80
Cadmium	Cd(112)	111.902761	24.13
Cadmium	Cd(113)	112.904401	12.22
Cadmium	Cd(114)	113.903361	28.73
Cadmium	Cd(116)	115.904758	7.49
Carbon	C(12)	12.000000	98.90
Carbon	C(13)	13.003355	1.10
Cerium	Ce(136)	135.907140	0.19
Cerium	Ce(138)	137.905996	0.25
Cerium	Ce(140)	139.905442	88.48
Cerium	Ce(142)	141.909249	11.08
Cesium	Cs(133)	132.905433	100.00
Chlorine	Cl(35)	34.968853	75.77
Chlorine	Cl(37)	36.965903	24.23
Chromium	Cr(50)	49.946046	4.35
Chromium	Cr(52)	51.940510	83.79
Chromium	Cr(53)	52.940651	9.50
Chromium	Cr(54)	53.938882	2.36
Cobalt	Co(59)	58.933198	100.00
Copper	Cu(63)	62.929599	69.17
Copper	Cu(65)	64.927792	30.83
Dysprosium	Dy(156)	155.924287	0.060
Dysprosium	Dy(158)	157.924412	0.10
Dysprosium	Dy(160)	159.925203	2.34
Dysprosium	Dy(161)	160.926939	18.90
Dysprosium	Dy(162)	161.926805	25.50
Dysprosium	Dy(163)	162.928737	24.90
Dysprosium	Dy(164)	163.929183	28.20
Erbium	Er(162)	161.928787	0.14
Erbium	Er(164)	163.929211	1.61
Erbium	Er(166)	165.930305	33.60
Erbium	Er(167)	166.932061	22.95
Erbium	Er(168)	167.932383	26.80
Erbium	Er(170)	169.935476	14.90
Europium	Eu(151)	150.919860	47.80
Europium	Eu(153)	152.921243	52.20
Fluorine	F(19)	18.998403	100.00
Gallium	Ga(69)	68.925581	60.10
Gallium	Ga(71)	70.924701	39.90
Gadolinium	Gd(152)	151.919803	0.200
Gadolinium	Gd(154)	153.920876	2.18
Gadolinium	Gd(155)	154.822629	14.80
Gadolinium	Gd(156)	155.922130	20.47
Gadolinium	Gd(157)	156.923967	15.65
Gadolinium	Gd(158)	157.924111	24.84
Gadolinium	Gd(160)	159.927061	21.86
Germanium	Ge(70)	69.924250	20.50
Germanium	Ge(72)	71.922080	27.40
Germanium	Ge(73)	72.923464	7.80
Germanium	Ge(74)	73.921179	36.50
Germanium	Ge(76)	75.921403	7.80
Gold	Au(197)	196.966560	100.00
Hafnium	Hf(174)	173.940065	0.16
Hafnium	Hf(176)	175.941420	5.20
Hafnium	Hf(177)	176.943233	18.60
Hafnium	Hf(178)	177.943710	27.10
Hafnium	Hf(179)	178.945827	13.74
Hafnium	Hf(180)	179.946561	35.20
Helium	He(3)	3.016029	.0001
Helium	He(4)	4.002603	100.00
Holmium	Ho(165)	164.930332	100.00
Hydrogen	H(1)	1.007825	99.99
Hydrogen	H(2)	2.014102	0.015
Indium	In(113)	112.904056	4.30
Indium	In(115)	114.903875	95.70
Iodine	I(127)	126.904477	100.00
Iridium	Ir(191)	190.960603	37.30
Iridium	Ir(193)	192.962942	62.70
Iron	Fe(54)	53.939612	5.80
Iron	Fe(56)	55.934939	91.72
Iron	Fe(57)	56.935396	2.20
Iron	Fe(58)	57.933278	0.28
Krypton	Kr(78)	77.920397	0.35
Krypton	Kr(80)	79.916375	2.25
Krypton	Kr(82)	81.913483	11.60
Krypton	Kr(83)	82.914134	11.50
Krypton	Kr(84)	83.911506	57.00
Krypton	Kr(86)	85.910614	17.30
Lanthanum	La(138)	137.907114	0.09
Lanthanum	La(139)	138.906355	99.91
Lead	Pb(204)	203.973037	1.40
Lead	Pb(206)	205.974455	24.10
Lead	Pb(207)	206.975885	22.10
Lead	Pb(208)	207.976641	52.40
Lithium	Li(6)	6.015123	7.42
Lithium	Li(7)	7.016005	92.58
Lutetium	Lu(175)	174.940785	97.40
Lutetium	Lu(176)	175.942694	2.60
Magnesium	Mg(24)	23.985045	78.90
Magnesium	Mg(25)	24.985839	10.00
Magnesium	Mg(26)	25.982595	11.10
Manganese	Mn(55)	54.938046	100.00
Mercury	Hg(196)	195.965812	0.15
Mercury	Hg(198)	197.966760	10.10
Mercury	Hg(199)	198.968269	17.00
Mercury	Hg(200)	199.968316	23.10
Mercury	Hg(201)	200.970293	13.20
Mercury	Hg(202)	201.970632	29.65
Mercury	Hg(204)	203.973481	6.80
Molybdenum	Mo(92)	91.906809	14.84
Molybdenum	Mo(94)	93.905086	9.25
Molybdenum	Mo(95)	94.905838	15.92
Molybdenum	Mo(96)	95.904676	16.68
Molybdenum	Mo(97)	96.906018	9.55
Molybdenum	Mo(98)	97.905405	24.13
Molybdenum	Mo(100)	99.907473	9.63
Neodymium	Nd(142)	141.907731	27.13
Neodymium	Nd(143)	142.909823	12.18
Neodymium	Nd(144)	143.910096	23.80
Neodymium	Nd(145)	144.912582	8.30
Neodymium	Nd(146)	145.913126	17.19
Neodymium	Nd(148)	147.916901	5.76
Neodymium	Nd(150)	149.920900	5.64
Neon	Ne(20)	19.992439	90.60
Neon	Ne(21)	20.993845	0.26
Neon	Ne(22)	21.991384	9.20
Nickel	Ni(58)	57.935347	68.27
Nickel	Ni(60)	59.930789	26.10
Nickel	Ni(61)	60.931059	1.13
Nickel	Ni(62)	61.928346	3.59
Nickel	Ni(64)	63.927968	0.91
Niobium	Nb(93)	92.906378	100.00
Nitrogen	N(14)	14.003074	99.63
Nitrogen	N(15)	15.000109	0.37
Osmium	Os(184)	183.952514	0.02
Osmium	Os(186)	185.953852	1.58
Osmium	Os(187)	186.955762	1.60
Osmium	Os(188)	187.955850	13.30
Osmium	Os(189)	188.958156	16.10
Osmium	Os(190)	189.958455	26.40
Osmium	Os(192)	191.961487	41.00
Oxygen	O(16)	15.994915	99.76
Oxygen	O(17)	16.999131	0.038
Oxygen	O(18)	17.999159	0.20
Palladium	Pd(102)	101.905609	1.02
Palladium	Pd(104)	103.904026	11.14
Palladium	Pd(105)	104.905075	22.33
Palladium	Pd(106)	105.903475	27.33
Palladium	Pd(108)	107.903894	26.46
Palladium	Pd(110)	109.905169	11.72
Phosphorus	P(31)	30.973763	100.00
Platinum	Pt(190)	189.959937	0.010
Platinum	Pt(192)	191.961049	0.79
Platinum	Pt(194)	193.962679	32.90
Platinum	Pt(195)	194.964785	33.80
Platinum	Pt(196)	195.964947	25.30
Platinum	Pt(198)	197.967879	7.20
Potassium	K(39)	38.963708	93.20
Potassium	K(40)	39.963999	0.012
Potassium	K(41)	40.961825	6.73
Praseodymium	Pr(141)	140.907657	100.00
Rhenium	Re(185)	184.952977	37.40
Rhenium	Re(187)	186.955765	62.60
Rhodium	Rh(103)	102.905503	100.00
Rubidium	Rb(85)	84.911800	72.17
Rubidium	Rb(87)	86.909184	27.84
Ruthenium	Ru(96)	95.907596	5.52
Ruthenium	Ru(98)	97.905287	1.88
Ruthenium	Ru(99)	98.905937	12.70
Ruthenium	Ru(100)	99.904218	12.60
Ruthenium	Ru(101)	100.905581	17.00
Ruthenium	Ru(102)	101.904348	31.60
Ruthenium	Ru(104)	103.905422	18.70
Samarium	Sm(144)	143.912009	3.10
Samarium	Sm(147)	146.914907	15.00
Samarium	Sm(148)	147.914832	11.30
Samarium	Sm(149)	148.917193	13.80
Samarium	Sm(150)	149.917285	7.40
Samarium	Sm(152)	151.919741	26.70
Samarium	Sm(154)	153.922218	22.70
Scandium	Sc(45)	44.955914	100.00
Selenium	Se(74)	73.922477	0.90
Selenium	Se(76)	75.919207	9.00
Selenium	Se(77)	76.919908	7.60
Selenium	Se(78)	77.917304	23.50
Selenium	Se(80)	79.916521	49.60
Selenium	Se(82)	81.916709	9.40
Silicon	Si(28)	27.976928	92.23
Silicon	Si(29)	28.976496	4.67
Silicon	Si(30)	29.973772	3.10
Silver	Ag(107)	106.905095	51.84
Silver	Ag(109)	108.904754	48.16
Sodium	Na(23)	22.989770	100.00
Strontium	Sr(84)	83.913428	0.56
Strontium	Sr(86)	85.909273	9.86
Strontium	Sr(87)	86.908902	7.00
Strontium	Sr(88)	87.905625	82.58
Sulfur	S(32)	31.972072	95.02
Sulfur	S(33)	32.971459	0.75
Sulfur	S(34)	33.967868	4.21
Sulfur	S(36)	35.967079	0.020
Tantalum	Ta(180)	179.947489	0.012
Tantalum	Ta(181)	180.948014	99.99
Tellurium	Te(120)	119.904021	0.096
Tellurium	Te(122)	121.903055	2.60
Tellurium	Te(123)	122.904278	0.91
Tellurium	Te(124)	123.902825	4.82
Tellurium	Te(125)	124.904435	7.14
Tellurium	Te(126)	125.903310	18.95
Tellurium	Te(128)	127.904464	31.69
Tellurium	Te(130)	129.906229	33.80
Terbium	Tb(159)	158.925350	100.00
Thallium	Tl(203)	202.972336	29.52
Thallium	Tl(205)	204.974410	70.48
Thorium	Th(232)	232.038054	100.00
Thulium	Tm(169)	168.934225	100.00
Tin	Sn(112)	111.904826	0.97
Tin	Sn(114)	113.902784	0.65
Tin	Sn(115)	114.903348	0.36
Tin	Sn(116)	115.901744	14.70
Tin	Sn(117)	116.902954	7.70
Tin	Sn(118)	117.901607	24.30
Tin	Sn(119)	118.903310	8.60
Tin	Sn(120)	119.902199	32.40
Tin	Sn(122)	121.903440	4.60
Tin	Sn(124)	123.905271	5.60
Titanium	Ti(46)	45.952633	8.00
Titanium	Ti(47)	46.951765	7.30
Titanium	Ti(48)	47.947947	73.80
Titanium	Ti(49)	48.947871	5.50
Titanium	Ti(50)	49.944786	5.40
Tungsten	W(180)	179.946727	0.13
Tungsten	W(182)	181.948225	26.30
Tungsten	W(183)	182.950245	14.30
Tungsten	W(184)	183.950953	30.67
Tungsten	W(186)	185.954377	28.60
Uranium	U(234)	234.040947	0.006
Uranium	U(235)	235.043925	0.72
Uranium	U(238)	238.050786	99.27
Vanadium	V(50)	49.947161	0.25
Vanadium	V(51)	50.943963	99.75
Xenon	Xe(124)	123.905894	0.10
Xenon	X(126)	125.904281	0.09
Xenon	Xe(128)	127.903531	1.91
Xenon	Xe(129)	128.904780	26.40
Xenon	Xe(130)	129.903510	4.10
Xenon	Xe(131)	130.905076	21.20
Xenon	Xe(132)	131.904148	26.90
Xenon	Xe(134)	133.905395	10.40
Xenon	Xe(136)	135.907219	8.90
Ytterbium	Yb(168)	167.933908	0.13
Ytterbium	Yb(170)	169.934774	3.05
Ytterbium	Yb(171)	170.936338	14.30
Ytterbium	Yb(172)	171.936393	21.90
Ytterbium	Yb(173)	172.938222	16.12
Ytterbium	Yb(174)	173.938873	31.80
Ytterbium	Yb(176)	175.942576	12.70
Yttrium	Y(89)	88.905856	100.00
Zinc	Zn(64)	63.929145	48.60
Zinc	Zn(66)	65.926035	27.90
Zinc	Zn(67)	66.927129	4.10
Zinc	Zn(68)	67.924846	18.80
Zinc	Zn(70)	69.925325	0.60
Zirconium	Zr(90)	89.904708	51.45
Zirconium	Zr(91)	90.905644	11.27
Zirconium	Zr(92)	91.905039	17.17
Zirconium	Zr(94)	93.906319	17.33
Zirconium	Zr(96)	95.908272	2.78"""


def _rephrase_isotope(s):  # C(13) -> [13C]
    e, n, _ = re.split(r"[()]", s)
    return f"[{n}{e}]"


def _build_isotope_df(periodic_table):
    # Turns the string above into a DataFrame.
    isotope_df = pd.DataFrame(
        [s.split("\t") for s in periodic_table.split("\n")],
        columns=["element", "isotope", "mass", "abundance"],
    )
    # Add a symbol column (element without mass) and convert strings to floats.
    isotope_df["symbol"] = isotope_df.isotope.apply(lambda s: re.sub(r"\(\d+\)", "", s))
    isotope_df.mass = np.float32(isotope_df.mass)
    isotope_df.abundance = np.float32(isotope_df.abundance)
    isotope_df1 = isotope_df.copy()
    isotope_df1.isotope = isotope_df1.isotope.apply(_rephrase_isotope)
    isotope_df = pd.concat([isotope_df, isotope_df1], ignore_index=True)
    return isotope_df


def _build_mass_dict(isotope_df):
    mass_dict = dict(
        isotope_df.groupby("symbol").apply(lambda t: t.mass[t.abundance == t.abundance.max()].max())
    )
    mass_dict.update(dict(zip(isotope_df.isotope, isotope_df.mass)))
    mass_dict["D"] = mass_dict["H(2)"]
    return mass_dict


isotope_df = _build_isotope_df(_periodic_table)
_mass_dict = _build_mass_dict(isotope_df)
_abundance_dict = dict(zip(isotope_df.isotope, isotope_df.abundance))

_abundance_ratio_dict = defaultdict(lambda: 1)
_abundance_ratio_dict.update(
    dict(
        C13=_abundance_dict["[13C]"] / _abundance_dict["[12C]"],
        Cl37=_abundance_dict["[37Cl]"] / _abundance_dict["[35Cl]"],
        Br81=_abundance_dict["[81Br]"] / _abundance_dict["[79Br]"],
        S34=_abundance_dict["[34S]"] / _abundance_dict["[32S]"],
        O18=_abundance_dict["[18O]"] / _abundance_dict["[16O]"],
        N15=_abundance_dict["[15N]"] / _abundance_dict["[14N]"],
    )
)

_iso_mass_diff_dict = defaultdict(lambda s: 0)
_iso_mass_diff_dict.update(
    dict(
        C13=_mass_dict["[13C]"] - _mass_dict["[12C]"],
        Cl37=_mass_dict["[37Cl]"] - _mass_dict["[35Cl]"],
        Br81=_mass_dict["[81Br]"] - _mass_dict["[79Br]"],
        S34=_mass_dict["[34S]"] - _mass_dict["[32S]"],
        O18=_mass_dict["[18O]"] - _mass_dict["[16O]"],
        N15=_mass_dict["[15N]"] - _mass_dict["[14N]"],
    )
)


def element_mass(e):
    return _mass_dict[e]


def isotope_abundance(iso):
    return _abundance_dict[iso]


def count_for_ratio(iso, ratio):
    return ratio / _abundance_ratio_dict[iso]


def isotope_mass_difference(iso):
    return _iso_mass_diff_dict[iso]


def molecular_ion_name(lib_id):
    if lib_id in POS_MODE_LIBRARIES:
        return "m+H"
    if lib_id in NEG_MODE_LIBRARIES:
        return "m-H"
    raise Exception(f"molecular_ion_name: can't establish mode for library {lib_id}")


def neutralize_ion_mz(mz, lib_id):
    if lib_id == 0:  # not sure why we'd do this.
        return mz
    if lib_id in POS_MODE_LIBRARIES:
        return mz + electron_mass
    if lib_id in NEG_MODE_LIBRARIES:
        return mz - electron_mass
    raise Exception(f"neutralize_ion_mass: can't establish mode for library {lib_id}")


def neutralize_precursor_mz(mz, lib_id):
    if lib_id == 0:
        return mz
    if lib_id in POS_MODE_LIBRARIES:
        return mz - proton_mass
    if lib_id in NEG_MODE_LIBRARIES:
        return mz + proton_mass
    raise Exception(f"neutralize_precursor_mass: can't establish mode for library {lib_id}")


def fmla_mass(fmla):
    try:
        charge, d = split_fmla(fmla)
        electrons = -charge * electron_mass
        return sum([v * element_mass(k) for k, v in d.items()]) + electrons
    except Exception as e:
        logging.warning(f"Couldn't find mass for formula {fmla}; returning -1")
        return -1


def split_fmla(s):
    """
    Split a formula into a dictionary of counts.
    :param s: the formula
    :return: a charge and a dict mapping elements to counts.
    """
    c = defaultdict(int)
    charge = 0
    # stuff + ( + digit,maybe + +/- + digit,maybe + )
    m = re.search(r"(.*?)[(](\d?)([+-])(\d?)[)]$", s)
    # New formats coming from matching database.
    # stuff + +/- + digit, maybe
    m1 = re.search(r"(.*)([+-])(\d?)$", s)
    if m:  # There's a charge; work it out.
        # stuff + ( + digit,maybe + +/- + digit,maybe + )
        s, n1, sign, n2 = m.group(1, 2, 3, 4)
        charge = {"+": 1, "-": -1}[sign]
        if len(max(n1, n2)) > 0:  # Could have: (2+) or (+2)...
            charge *= int(max(n1, n2))
    elif m1:
        # stuff + +/- + digit,maybe
        s, sign, n = m1.groups()
        n = 1 if n == "" else int(n)
        charge = n * {"+": 1, "-": -1}[sign]

    lis = re.split(r"([A-Z][a-z]?[']?)", s)[1:]
    for e, n in zip(lis[0::2], lis[1::2]):
        c[e] += 1 if n == "" else int(n)
    return charge, c


# These priorities were extracted from the library.  C will come first in a formula and Pt will come last.
_priorities = defaultdict(int)
_priorities.update(
    zip(qw("Pt Se Si S P O K Mo Ni N I Cs Na Hf Co Zr Mg Fe No F D Cl Br H C"), range(100))
)


def unsplit_fmla(charge, counts):
    """
    Turn a formula in dict form such as returned by split_fmla into a string formula.
    :param charge: The charge as an integer
    :param counts: Dict containing counts of elements.
    :return: Formula in string form.
    """
    charge = "" if charge == 0 else f"({charge})" if charge < 0 else f"(+{charge})"
    pairs = sorted(counts.items(), key=lambda p: _priorities[p[0]], reverse=True)
    s = "".join([f"{e}{k}" if k > 1 else e for e, k in pairs]) + charge
    return s


def canonicalize_fmla(s):
    """
    Organizes the elements in a molecular formula into a canonical order so that it can easily be compared to another.
    :param s: a molecular formula
    :return: a formula with elements in canonical order.
    """
    return unsplit_fmla(*split_fmla(s))


def subset_fmla(f1, f2):
    """
    Determines if f1 can be a fragment of f2, i.e. if f1 is C7, f2 must be C>=7, etc.
    :param f1: candidate fragment formula
    :param f2: precursor formula
    :return: True iff f1 could be a fragment of f2.
    """
    c1, d1 = split_fmla(f1)
    c2, d2 = split_fmla(f2)
    return all([c1 <= d2[e1] for e1, c1 in d1.items()])


def neutralize_fmla(s):
    """
    Given a formula for an ion, i.e., positively or negatively charged, return corresponding neutral formula.
    :param s: a (charged) ion formula.
    :return: the corresponding neutral formula.  (If original formula was neutral, returns it.)
    """  # OK, this is really, really correct.
    try:
        charge, d = split_fmla(s)
        d["H"] -= charge
        return unsplit_fmla(0, d)
    except:
        logging.warning(f"Couldn't neutralize {s}; returning {s}")
        return s


# Computes the Double Bond Equivalent of a formula.
# This number is based on the combinatorics of valence:
# If the formula is legit, it will be a whole number >= 0.
# If the formula is a hydrocarbon (all C and H), this will be the number of double bonds.
# It's an extension of the (Euler's?) observation that the number of vertices of odd degree in a groph must be even.
def dbe(fmla):
    """
    Computes the Double Bond Equivalence of a formula.  If the formula is legit, it will be a whole number >= 0.
    :param fmla: a formula
    :return: a number
    """
    charge, d = split_fmla(neutralize_fmla(fmla))
    # No D from our predictions, but we might evaluate something from MLIMS.
    # To avoid filtering out good formulas, we need to use MAXIMUM valence, i.e., P=5, S=6.
    # In particular, this applies to sulfates and phosphates.
    # N can be 5, but it might generate too many false positives.
    valence = dict(C=4, H=1, F=1, Br=1, Cl=1, I=1, N=3, P=5, S=6, O=2, D=1)
    dbe = sum([d[e] * (valence[e] - 2) for e in qw("C H F Cl Br I N P S D")]) / 2 + 1
    return dbe


def is_good_dbe(fmla):
    """
    Reports if the DBE number of a formula is "good": a whole number >= 0.
    :param fmla: a formula
    :return: True iff the DBE number is a whole number >= 0.
    """
    k = dbe(fmla)  # dbe neutralizes formula
    return k >= 0 and k == int(k)


def get_c_count(fmla):
    """
    Given a formula string, how many carbons?
    :param fmla: formula
    :return: number of carbons.
    """
    charge, d = split_fmla(fmla)
    return d["C"]
