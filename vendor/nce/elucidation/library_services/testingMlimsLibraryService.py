import pandas as pd
from library_services.cachedMlimsLibraryService import CachedMlimsLibraryService
from library_services.mlimsLibraryService import MlimsLibraryService


class TestingMlimsLibraryService:
    __instance__ = None
    __name_dict__ = None

    def __init__(self):
        print("Creating TestingOracleDataProvider")
        self.odp = MlimsLibraryService.get_instance()
        self.codp = CachedMlimsLibraryService.get_instance()

    @staticmethod
    def get_instance():
        if TestingMlimsLibraryService.__instance__ is None:
            TestingMlimsLibraryService.__instance__ = TestingMlimsLibraryService()
        return TestingMlimsLibraryService.__instance__

    def get_name_dict(self) -> dict:
        print("get_name_dict")
        if TestingMlimsLibraryService.__name_dict__ is not None:
            return TestingMlimsLibraryService.__name_dict__
        d1 = self.odp.get_name_dict()
        d2 = self.codp.get_name_dict()
        if len(d1) != len(d2):
            print(len(d1), len(d2))
            raise Exception("get_name_dict different")
        TestingMlimsLibraryService.__name_dict__ = d1
        print("---> OK")
        return d2

    def getMs2(self, lib_id):
        print("getMS2")
        df1 = self.odp.getMs2(lib_id)
        df2 = self.codp.getMs2(lib_id)
        if len(df1) != len(df2):
            print(f"getMs2 different {len(df1)} {len(df2)}")
            print(df1.columns)
            print(df2.columns)
            raise Exception("")
        print("---> OK")
        return df2

    def get_study_title(self, study_id):
        return self.codp.get_study_title(study_id)
        print("get_study_title")
        s1 = self.odp.get_study_title(study_id)
        s2 = self.codp.get_study_title(study_id)
        if s1 != s2:
            print(f"get_study_title different")
            print(s1)
            print(s2)
            raise Exception("")
        print("---> OK")
        return s2

    def get_uniq_fmlas(self, lib_id):
        print("get_uniq_fmlas")
        df1 = self.odp.get_uniq_fmlas(lib_id)
        df2 = self.codp.get_uniq_fmlas(lib_id)
        if len(df1) != len(df2):
            print(f"get_uniq_fmlas different {len(df1)} {len(df2)}")
            raise Exception("")
        print("---> OK")
        return df2

    def get_isotope_from_mlims(self, iso, lib_id, comp_ids):
        print("get_isotope_from_mlims")
        df1 = self.odp.get_isotope_from_mlims(lib_id, comp_ids)
        df2 = self.codp.get_isotope_from_mlims(lib_id, comp_ids)
        if len(df1) != len(df2):
            print(f"get_isotope_from_mlims different {len(df1)} {len(df2)}")
            raise Exception("")
        print("---> OK")
        return df2

    def get_library_ions(self, lib_id: int) -> pd.DataFrame:
        return self.codp.get_library_ions(lib_id)
        print("get_library_ions")
        df1 = self.odp.get_library_ions(lib_id)
        df2 = self.codp.get_library_ions(lib_id)
        if len(df1) != len(df2):
            print(f"get_library_ions different {len(df1)} {len(df2)}")
            raise Exception("")
        print("---> OK")
        return df2

    def get_chemicals_for_formula(self, fmla: str, lib_id: int) -> pd.DataFrame:
        print("get_chemicals_for_formula")
        df1 = self.odp.get_chemicals_for_formula(fmla, lib_id)
        df2 = self.codp.get_chemicals_for_formula(fmla, lib_id)
        if len(df1) != len(df2):
            print(f"get_chemicals_for_formula different {len(df1)} {len(df2)}")
            raise Exception("")
        print("---> OK")
        return df2

    def get_known_ms2(self, lib_id: int, mass_ratio_lower_limit: float = 1):
        print("get_known_ms2")
        df2 = self.codp.get_known_ms2(lib_id, mass_ratio_lower_limit=mass_ratio_lower_limit)
        df1 = self.odp.get_known_ms2(lib_id, mass_ratio_lower_limit=mass_ratio_lower_limit)
        if len(df1) != len(df2):
            print(f"get_known_ms2 different {len(df1)} {len(df2)}")
            raise Exception("")
        print("---> OK")
        return df2

    def getUnknownMs2(self, lib_id):
        print("getUnknownMs2")
        df1 = self.odp.getUnknownMs2(lib_id)
        df2 = self.codp.getUnknownMs2(lib_id)
        if len(df1) != len(df2):
            print(f"getUnknownMs2 different {len(df1)} {len(df2)}")
            raise Exception("")
        print("---> OK")
        return df2

    def get_unknown_ms2(self, lib_id: int, mass_ratio_lower_limit: float = 1) -> pd.DataFrame:
        print("get_unknown_ms2")
        df1 = self.odp.get_unknown_ms2(lib_id, mass_ratio_lower_limit=mass_ratio_lower_limit)
        df2 = self.codp.get_unknown_ms2(lib_id, mass_ratio_lower_limit=mass_ratio_lower_limit)
        if len(df1) != len(df2):
            print(f"get_unknown_ms2 different {len(df1)} {len(df2)}")
            raise Exception("")
        print("---> OK")
        return df2
