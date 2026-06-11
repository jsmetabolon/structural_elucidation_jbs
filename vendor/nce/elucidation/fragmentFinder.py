import numpy as np
from elucidation.binning import plus_5ppm, minus_5ppm
from elucidation.conveniences import qw
from elucidation.chemFacts import electron_mass, fmla_mass, proton_mass, neutralize_precursor_mz
from elucidation.chemFacts import neutralize_ion_mz
from elucidation.rdkitUtils import smiles_mass
from elucidation.tableLoaders import library_fragments_loader


class FragmentFinder:
    __instances__ = {}

    @staticmethod
    def get_instance(lib_id):
        if lib_id is None:  # No instances have been created, so load tables into class variables.
            lib_id = 0
        if lib_id not in FragmentFinder.__instances__:
            FragmentFinder.__instances__[lib_id] = FragmentFinder(lib_id)

        return FragmentFinder.__instances__[lib_id]

    def __init__(self, lib_id):
        if lib_id in FragmentFinder.__instances__:
            raise Exception("try again with FragmentFinder.get_instance(lib_id)")
        self.lib_id = lib_id
        _d = library_fragments_loader()
        self.masses, self.fragments, self.counts = [_d[k] for k in qw("masses fragments counts")]

    def get_fragments_for_mass(self, m):
        # return []
        m_lo = minus_5ppm(m) if m > 100 else m - 0.0005
        m_hi = plus_5ppm(m) if m > 100 else m + 0.0005
        masses = self.masses
        fragments = self.fragments
        counts = self.counts
        lo = np.searchsorted(masses, m_lo)
        results = []
        for i in range(lo, len(masses)):  # I *think* this is faster than while
            if masses[i] > m_hi:
                break
            results.append((masses[i], fragments[i], counts[i]))
        return results

    def get_fragments_for_ion_mz(self, mz):
        return self.get_fragments_for_mass(neutralize_ion_mz(mz, self.lib_id))

    # This is sort of a mess because of adding and subtracting the offset repeatedly.
    def get_fragments_for_loss(self, ion_mz, precursor_mz):
        precursor = neutralize_precursor_mz(precursor_mz, self.lib_id)
        other_fragment = neutralize_ion_mz(ion_mz, self.lib_id)
        neutral_loss = precursor - other_fragment
        return neutral_loss, FragmentFinder.get_instance(0).get_fragments_for_mass(neutral_loss)


if __name__ == "__main__":
    answer = [
        (99.11737549502459, "*C(C)CCCCC", 2),
        (99.11737549502459, "*CC(C)CCCC", 2),
        (99.11737549502459, "*CC(CC)C(C)C", 2),
        (99.11737549502459, "*CCC(C)C(C)C", 2),
        (99.11737549502459, "*C(CC)CCCC", 3),
        (99.11737549502459, "*C(CCC)CCC", 3),
        (99.11737549502459, "*CCCC(C)CC", 6),
        (99.11737549502459, "*CCCCC(C)C", 22),
        (99.11737549502459, "*CCCCCCC", 701),
    ]

    fs = FragmentFinder(0)
    print(fs.get_fragments_for_mass(smiles_mass("CCCCCCC") - fmla_mass("H")) == answer)

    fs = FragmentFinder(209)
    print(fs.get_fragments_for_mass(smiles_mass("CCCCCCC") - proton_mass) == answer)

    fs = FragmentFinder(400)
    print(
        fs.get_fragments_for_mass(smiles_mass("CCCCCCC") - fmla_mass("H") - electron_mass) == answer
    )
