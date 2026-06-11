"""
This program builds a cache of arithmetically possible formulas by mass.
The resulting file is read by FormulaFinder.
If it is run with the -j flag, it will explain itself...
"""

import sys
from collections import Counter
from functools import reduce
from itertools import product
from operator import itemgetter, add
from pprint import pprint

from elucidation.chemFacts import element_mass, fmla_mass
from elucidation.conveniences import pickle_dump, qw, lmap, lzip, now_string, get_git_branch
from elucidation.formulaFinder import (
    _build_large_mass_tables,
    _mass_from_tuple,
    _count_multiplicity,
    pack_fmla_tuple,
    FormulaFinder,
)
from mdsutils.easyQuery import EasyQuery

fmla_finder_cache_file = "fmla_finder_cache_file.pkl"
import numpy as np
import regex as re


def quick_tests():
    return
    # Below doesn't work because the file has not been moved to the reference directory.
    # It uses the old file.
    ff = FormulaFinder.get_instance(0)
    for f in qw("C10H20P2 H2SO4 C H2O C3ClF3 CHBr"):
        df = ff.fmlas_for_neutral_mass(fmla_mass(f))
        if leng(df[df.cand_fmla == f]) > 0:
            print("test of ", f, "passed")
            print("test of ", f, "failed")


if __name__ == "__main__":

    print(
        f"""Building Formula Finder Reference File in {fmla_finder_cache_file}:
    This will take 10-20 minutes.
    """
    )
    justify = sys.argv[-1] in ("-j", "--justify")
    if justify:
        print(
            """First we will repeat the analysis justifying the choices.
        This will add a couple of minutes.
        """
        )
    else:
        print("""If you want the design explanation, rerun with the '-j' option.""")
    mlims_fmlas = None  # Pleases pycharm
    if justify:
        with EasyQuery(db="mlims") as eq:
            mlims_fmlas = eq.get_chemical_formulas()
        print(f"{len(mlims_fmlas)} total formulas; {len(set(mlims_fmlas))} unique formulas")

        # Get all the Br, Cl, I, F, P from formulas
        lis = sorted(
            Counter(
                ["_".join(sorted(re.findall(r"(Br|Cl|P|I|F)", s))) for s in mlims_fmlas]
            ).items(),
            key=itemgetter(1),
            reverse=True,
        )
        # Take all of above that have more than one.
        co_occur = [(s, n) for s, n in lis if "_" in s]
        print(
            "We can group all these together (Br|Cl|P|I|F)\n because there are few co-occurences:",
            co_occur,
        )

        lis = _count_multiplicity("F", mlims_fmlas)
        p, _, r = lis[2]
        print("example: ")
        print("    F", lis)
        print(
            f"    If we limit F to 0 <= F < {p}, it excludes formulas for {r} chemicals in the CHEMICALS table."
        )
        print()

        for elt in qw("Br Cl P I F C H N O S Br"):
            print(elt, _count_multiplicity(elt, mlims_fmlas), "\n")

    print(
        """our preliminary choices: C 0-55,  H 0-94, N 0-11, O 0-22, S 0-2, Br 0-1, Cl 0-2, P 0-3, I 0-1, F 0-3"""
    )

    def build_elt_list(elt, n):
        e_mass = element_mass(elt)
        return [("", 0), (elt, e_mass)] + [(f"{elt}{i}", e_mass * i) for i in range(2, n + 1)]

    opt_lists = {
        elt: build_elt_list(elt, n)
        for elt, n in [
            ("C", 55),
            ("H", 94),
            ("N", 11),
            ("O", 22),
            ("S", 2),
            ("Cl", 2),
            ("P", 2),
            ("I", 1),
            ("F", 3),
            ("Br", 1),
        ]
    }

    x_list = reduce(
        add, lmap(opt_lists.get, qw("Cl P I F Br"))
    )  # Combine into one group.  So we can't have e.g. CLBr in same formula
    # I think this may help with drugs:
    x_list.append(("ClF3", element_mass_dict["Cl"] + 3 * element_mass_dict["F"]))

    opt_lists["x"] = sorted(set(x_list))  # remove extra ("",0)
    options = lmap(opt_lists.get, qw("C H N O S x"))
    pprint(options)

    if justify:
        print(
            f"""our final choices: C 0-55,  H 0-94, N 0-11, O 0-22, S 0-2
        PLUS one option from this list:
        {lmap(itemgetter(0), opt_lists['x'])}"""
        )

        logs = np.ceil(np.log(lmap(len, options)) / np.log(2))
        print(
            f"""giving {' * '.join([str(len(lis)) for lis in options])} = {np.prod(lmap(len, options))} possible formulas,
        which we will represent by {' + '.join(lmap(str, logs))} = {np.sum(logs)} bits in a 32-bit integer."""
        )

    # Form Cartesian product of options lists.  Each element is a 6-tuple of indices into the options lists.
    print("Packing formulas into ints and pairing with masses.  This takes several (8?) minutes.")
    cart_prod = product(*[range(len(lis)) for lis in options])
    all_pairs = sorted([(_mass_from_tuple(t, options), pack_fmla_tuple(t)) for t in cart_prod])

    # eliminate formulas that are too heavy.  Keep the light ones because we might want to look up losses.
    small_pairs = [pair for pair in all_pairs if pair[0] <= 50]
    large_pairs = [pair for pair in all_pairs if 50 < pair[0] < 1000]
    print(
        f"""{len(large_pairs)} formulas 50 < m < 1000 and {len(small_pairs)} formulas 0 < m <= 50
    These will be searched separately.  
    For the moment, we will not handle queries that cross the boundary at 50 very gracefully"""
    )

    # Build query structures.
    # Masses <50 are handled specially. There aren't very many.
    large_ints, large_bins = _build_large_mass_tables(large_pairs)
    small_masses, small_ints = lzip(*small_pairs)

    time_str = now_string()
    pkl = f"{fmla_finder_cache_file[:-4]}_{time_str}.pkl"
    print("writing", pkl)

    pickle_dump(
        pkl,
        dict(
            large_ints=large_ints,
            large_bins=large_bins,
            small_ints=small_ints,
            small_masses=small_masses,
            options=options,
            timestamp=time_str,
            git_branch=get_git_branch(),
        ),
    )
    print("Dump finished.")

    quick_tests()

    if not justify:
        exit()

    print(
        """Continuing with justification.
          What MLIMS formulas are missing from the tables?"""
    )

    def sort_fmla(fmla):
        elts = re.findall(r"[A-Z][a-z]?\d*", fmla)
        return "".join(sorted(elts))

    our_fmlas = [sort_fmla(ff.formula_from_int(i)) for i in large_ints[:-1]] + [
        ff.formula_from_int(i) for i in small_ints
    ]
    mlims_fmlas = lmap(sort_fmla, mlims_fmlas)

    print(f"{len(our_fmlas)} formulas in our search tables.")
    print(f"{len(mlims_fmlas)} formulas in MLIMS whereof {len(set(mlims_fmlas))} are unique")

    our_fmlas = set(our_fmlas)
    missing = [f for f in mlims_fmlas if f not in our_fmlas]
    print(
        f"{len(missing)} MLIMS formulas are missing from our tables ({len(set(missing))} unique formulas."
    )
    print("Missing formulas:")
    pprint(sorted(set(missing)))
    print("Done")
