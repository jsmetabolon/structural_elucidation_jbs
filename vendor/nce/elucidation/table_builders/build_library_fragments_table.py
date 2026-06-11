from collections import Counter
from functools import reduce
from operator import add

import numpy as np
from elucidation.conveniences import pickle_dump, now_string, lmap
from elucidation.rdkitUtils import smiles_mass, fragment_smiles

from elucidation.rdkitUtils import canonicalize_smiles as canonicalize
from elucidation.rdkitUtils import smiles_mass
from mdsutils.easyQuery import EasyQuery

with EasyQuery(db="mlims") as eq:
    all_smiles = list(
        eq.easy_query(
            """select max(smiles) smiles from chemicals ch
        join chro_lib_entry cle on cle.chemical_id = ch.chemical_id
        where smiles is not null
        and lib_id in (209,305,400,402)
        group by comp_id"""
        ).smiles
    )

can_smiles = set(lmap(canonicalize, all_smiles))

all_frags = reduce(add, map(fragment_smiles, can_smiles))

frag_counts = Counter(all_frags)

all_pairs = sorted([(smiles_mass(f), c, f) for f, c in frag_counts.items()])
all_pairs = [p for p in all_pairs if p[0] > 0 and p[0] < 1000]

all_masses, all_counts, all_frags = lmap(np.array, zip(*all_pairs))

time_str = now_string()
pkl = f"c:/cdtemp/library_fragments_{time_str}.pkl"
print("writing", pkl)
pickle_dump(
    pkl, dict(masses=all_masses, fragments=all_frags, counts=all_counts, timestamp=time_str)
)
