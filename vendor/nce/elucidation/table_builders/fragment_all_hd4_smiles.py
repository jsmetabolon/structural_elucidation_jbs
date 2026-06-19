"""Extract unique fragment neutral masses from all HD4 library SMILES.

One-off reference builder that pickles distinct fragment masses in the 50–1000 Da
range for modeling and matching tables.
"""

from functools import reduce
from operator import add

from elucidation.conveniences import pickle_dump, lmap, now_string
from elucidation.chemFacts import proton_mass
from elucidation.rdkitUtils import smiles_mass, canonicalize, fragment_smiles
from mdsutils.easyQuery import EasyQuery
import numpy as np

eq = EasyQuery(db="mlims")

df = eq.easy_query(
    """
select smiles from chemicals join chro_lib_entry cle on chemicals.chemical_id = cle.chemical_id
where smiles is not null and lib_id in (209,305,400,402)"""
)

all_smiles = list(df.smiles)

all_fragments = reduce(add, map(fragment_smiles, all_smiles))
print(len(all_fragments))
all_fragments = set(all_fragments)
print(len(all_fragments))

all_fragment_masses = set(np.array(lmap(smiles_mass, all_fragments)).round(4))
all_fragment_masses = [m for m in all_fragment_masses if 50 < m < 1000]
print(len(all_fragment_masses))

time_str = now_string()
pkl = "C:/cdtemp/references/all_HD4_fragment_masses_{time_str}.pkl"
print("writing", pkl)
pickle_dump(pkl, dict(all_fragment_masses=all_fragment_masses, timestamp=time_str))
