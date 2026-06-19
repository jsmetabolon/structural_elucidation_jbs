"""Offline code generator for buddyFunction.py adduct delta tables.

Uses SymPy to symbolically derive MS1 buddy mass relationships and emits the
build_delta_df source consumed by MS1Analyzer. Run manually when ion chemistry
rules change.
"""

import argparse
import itertools

import pandas as pd
import regex as re
from functools import reduce
from operator import add

from elucidation.chemFacts import (
    electron_mass,
    element_mass,
    fmla_mass,
    proton_mass,
)

# I seem to recall we can't import numpy as np because sympy also defines an np.
# I'm leaving this dead code ^ because it was hard to sort out.
# from numpy import fabs

from regex import compile  # Something similar to above.
from elucidation.conveniences import qw
from sympy import Eq, solve, symbols, sympify

# also hard to sort out:
# from sympy import Float, srepr, Integer

from elucidation.tableLoaders import get_config_var

local_base_dir = get_config_var("LOCAL_BASE_DIR")

# These are groups of ion names that will be treated in the same way.
ion_name_lists = {
    402: dict(
        half=qw(
            "(m+2H)/2 (m+2H[S34-1])/2 (m+2H[C13-1])/2 (m+2H[C13-1][Cl37-1])/2 (m+2H[Cl37-1])/2 (m+2H[C13-2])/2"
        ),
        third=qw(
            "(m+3H)/3 (m+3H[C13-1])/3 (m+3H[C13-1][Cl37-1])/3 (m+3H[Cl37-1])/3 (m+3H[C13-2])/3 (m+3H[S34-1])/3"
        ),
        twice=qw(
            "2m+H 2m+H[C13-1] 2m+H[C13-1][Cl37-1] 2m+H[Cl37-1] 2m+H[C13-2] 2m+H[C13-3] 2m+H[S34-1] 2m+Na 2m+Na[C13-1]"
        ),
        thrice=qw("3m+H 3m+H[C13-1] 3m+H[C13-1][Cl37-1] 3m+H[Cl37-1] 3m+H[S34-1]"),
        once=qw(
            "m+H  m+H[C13-1][Cl37-1] m+H[C13-1] m+H[C13-2] "
            + "m+H[C13-3] m+H[S34-1] m+H[S34-2] m+H[Cl37-1] m+H[Cl37-2] "
            + "m+H-H2O m+H-H2O[C13-1] "
            + "m+H-H4O2 m+H-H4O2[C13-1] "
            + "m+H-CH2O2 m+H-CH2O2[C13-1] "
            + "m+H+NH3 m+H+NH3[C13-1] "
            + "m+Na m+Na[C13-1] "
            + "m+K m+K[C13-1]"
            # m+H[Br81-1] m+H[C13-1][Br81-1]"
        ),
        isotopes=qw("m+H[Br81-1] m+H m+H[O18-1] m+H[D-1]"),
    ),
    209: dict(
        half=qw(
            "(m-2H)/2 (m-2H[S34-1])/2 (m-2H[C13-1])/2 (m-2H[C13-1][Cl37-1])/2 (m-2H[Cl37-1])/2 (m-2H[C13-2])/2"
        ),
        third=qw(
            "(m-3H)/3 (m-3H[C13-1])/3 (m-3H[C13-1][Cl37-1])/3 (m-3H[Cl37-1])/3 (m-3H[C13-2])/3 (m-3H[S34-1])/3"
        ),
        # This was just to test out a hypothesis made about one particular feature.
        # ninth=qw(
        # "(m-9H)/9 (m-9H[C13-1])/9 (m-9H[C13-1][Cl37-1])/9 (m-9H[Cl37-1])/9 (m-9H[C13-2])/9 (m-9H[C13-3])/9 (m-9H[C13-4])/9"
        # (,
        twice=qw(
            "2m-H 2m-H[C13-1] 2m-H[C13-1][Cl37-1] 2m-H[Cl37-1] 2m-H[C13-2] 2m-H[C13-3] 2m+Na-2H 2m-H[S34-1]"
        ),
        thrice=qw("3m-H 3m-H[C13-1] 3m-H[C13-1][Cl37-1] 3m-H[Cl37-1] 3m-H[S34-1]"),
        once=qw(
            "m-H m-2H+Na m-H-CO2 m-H-H2O m-H-H2O[C13-1] m-H[C13-1][Cl37-1] m-H[C13-1] m-H[C13-2] "
            + "m-H[C13-3] m-H[S34-1] m-H[S34-2] m-H[Cl37-1] m-H[Cl37-2] m-H+CH2O3 m-H+CH2O3[C13-1] "
            + "m+Cl35 m+Cl35[C13-1] m+Cl35[C13-1][Cl37-1] m+Cl35[Cl37-1]"  # m-H[Br81-1] m-H[C13-1][Br81-1]"
        ),
        isotopes=qw("m-H[Br81-1] m-H m-H[O18-1] m-H[D-1] m-H[N15-1]"),
    ),
}

# These are the mass substitutions for constants in formulas.
substitutions = [
    (sympify("H"), proton_mass),
    (sympify("D"), element_mass("H(2)") - element_mass("H(1)")),
    (sympify("C13"), element_mass("C(13)") - 12),
    (sympify("O18"), element_mass("O(18)") - element_mass("O(16)")),
    (sympify("Na"), element_mass("Na") - electron_mass),
    (sympify("K"), element_mass("K") - electron_mass),
    (sympify("N15"), element_mass("N(15)") - element_mass("N(14)")),
    (sympify("Cl35"), element_mass("Cl(35)") + electron_mass),
    (sympify("Cl37"), element_mass("Cl(37)") - element_mass("Cl(35)")),
    (sympify("Br81"), element_mass("Br(81)") - element_mass("Br(79)")),
    (sympify("S34"), element_mass("S(34)") - element_mass("S(32)")),
] + [(sympify(s), fmla_mass(s)) for s in qw("H2O CO2 CH2O3 CH2O2 NH3 H4O2")]


# Take an ion name from the library (such as (m-2H[C13-1][Cl37-1])/2) and turns it into a sympy formula.
def process_ion_name(s):
    # [D-2] => (2*D)
    s = compile(r"\[D-(\d)\]").sub(r"+(\1*D)", s)
    # [Cl37-1] => (1*Cl37)
    s = compile(r"\[([A-Z][a-z]?\d+)-(\d)\]").sub(r"+(\2*\1)", s)
    # 3m => (3*m);    3H => (3*H)
    s = compile(r"(\d)([mH])").sub(r"(\1*\2)", s)
    # Should be a good normal algebraic formula now, so simplify and fill in constants.
    s = sympify(s)
    s = s.subs(substitutions)
    return s


# Want: delta in terms of f
# f label, b label
# m in terms of delta  (inferred mass) - this will be a formula to pass back into buddyFunctions with delta:
#                                      - it's just a substitution, but we want encapsulate sympy.
# delta = b-f, m=func(f)


def create_row(f_str, b_str, formula_dict):
    a, d, f, m, ff = symbols("a d f m ff")

    f = formula_dict[f_str]
    b = formula_dict[b_str]
    if (b_str.startswith("m+H") and f_str.startswith("m+H")) or (
        b_str.startswith("m-H") and f_str.startswith("m-H")
    ):
        b_str = "x" + b_str[3:]
        f_str = "x" + f_str[3:]
        if len(b_str) > 1 and len(f_str) > 1:
            return
    if ("Cl" in b_str and "Br" in f_str) or ("Br" in b_str and "Cl" in f_str):
        return

    # a formula for inferred mass (m) in terms of feature mass. (ff)
    m_fmla = solve(Eq(f, ff), m)[0]
    d_fmla = (b - f).subs(m, m_fmla)  # a formula for delta in terms of mass m
    if d_fmla.is_real and -0.001 < d_fmla < 0.001:  # delta is too small.
        return
    if not d_fmla.is_real:
        # print(d, d_fmla, type(d_fmla), f,b)
        # In this case, we can more accurately estimate feature mass and therefore inferred mass from delta.
        ff_fmla = solve(Eq(d_fmla, d), ff)[0]  # a formula for feature mass in terms of delta.
        m_fmla = m_fmla.subs(ff, ff_fmla)  # a formula for inferred mass in terms of delta.
    return {
        "delta": round_in_string(str(d_fmla)),
        "f_str": f_str,
        "f": f,
        "b_str": b_str,
        "b": b,
        "length": len(f_str) + len(b_str),
        "m": f"'{round_in_string(str(m_fmla))}'",
        # this goes in quotes because it will have to be returned to MS1Analyzer
        # and then passed back with an observed value for delta.
        # "b_fmla": b_fmla,
    }


def construct_unified_list(ion_name_groups):
    all_ion_names = set(reduce(add, ion_name_groups.values()))
    all_ion_pairs = list(itertools.permutations(all_ion_names, 2))

    # Saves processing the same ion name repeatedly:
    formula_dict = {ion_name: process_ion_name(ion_name) for ion_name in all_ion_names}
    results = [create_row(f_str, b_str, formula_dict) for f_str, b_str in all_ion_pairs]
    results = [r for r in results if r is not None]
    df = pd.DataFrame(results)
    print(len(df))
    df = df.groupby("delta").apply(clean_tab)
    df = df.drop(["length"], axis=1)
    print(len(df))
    return df


def clean_tab(tab):
    base_ions = qw("x m-H m+H")
    tab1 = tab[tab.f_str.isin(base_ions) | tab.b_str.isin(base_ions)]
    if len(tab1) > 0:
        return tab1
    if len(tab) > 1:
        tab = tab[tab.length == tab.length.min()]
        tab = tab[tab.f_str.str.len() == tab.f_str.str.len().min()]
        tab = tab[tab.f_str == tab.f_str.min()]
        if len(tab) > 1:
            print("Still not unique!")
            print(tab)
            print("=======================================================================")
    return tab


def get_arguments(argv):
    parser = argparse.ArgumentParser(
        description="Generates python code for module for MS1 analysis using sympy symbolic algebra package."
    )
    parser.add_argument(
        "--py",
        type=str,
        default=f"{local_base_dir}/buddyFunctions.py",
        help="The python file that will code the output code.",
    )
    return parser.parse_args(argv[1:])


# These are for debugging and testing.
# little_ion_names = dict(
#     half=qw("(m-2H)/2 (m-2H[S34-1])/2"),
#     once=qw("m-H m-2H+Na m-H-CO2 m-H-H2O"),
#     isotopes=qw("m-H[Br81-1] m-H[C13-1]"),
# )


# These are bits of the code produced.
src_start = """import pandas as pd
from elucidation.conveniences import qw

def build_delta_df(ff, lib_id):
"""

table_start = """
    if lib_id in %s:
        df = pd.DataFrame([
"""

src_end = """], 
        columns = qw("feature buddy theor_delta inferred_m"))
        df = df.sort_values('theor_delta')  # Not sure how this helps.
        df = df.reset_index(drop=True) 
        return df
"""


def round_in_string(s):
    t = re.sub(r"(\.\d{5})\d+", r"\1", s)
    return t


def debug():
    with open("new_buddyFunction.py", "w") as src:
        print(src_start, file=src)
        for lib_ids in ((402, 400), (209, 305)):
            print(table_start % str(lib_ids), file=src)
            df = construct_unified_list(ion_name_lists[lib_ids[0]])
            for t in df.itertuples():
                m_s = round_in_string(t.m)
                delta_s = round_in_string(t.delta)
                print(f"('{t.f_str}', '{t.b_str}', {delta_s}, {m_s}),", file=src)
            print(src_end, file=src)


debug()
