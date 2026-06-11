import os
import pandas as pd
from elucidation.conveniences import qw, pickle_dump, now_string
from dotenv import load_dotenv
from mdsutils.easyQuery import EasyQuery

if __name__ == "__main__":
    load_dotenv()
    creds = dict(db=os.environ["DB_NAME"], user=os.environ["DB_USER"])

    d = dict()

    with EasyQuery(**creds) as eq:
        for s in qw("chro_lib_entry chro_lib_masses"):
            d[s] = eq.easy_query(
                f"select * from {s} where lib_id in (209,305,400,402,9209,9400,9305,9402)"
            )
        d["chemical_synonym"] = eq.easy_query(
            "select * from chemical_synonym where report_name='T'"
        )
        d["sample_sets"] = eq.easy_query(
            "select * from sample_sets where created_on > date '2017-01-01'"
        )
        d["chro_comp_entry"] = eq.easy_query("select * from chro_comp_entry")
        d["chemicals"] = eq.easy_query("select * from chemicals")
        time_str = now_string()
        d["timestamp"] = time_str
        pkl = f'{os.environ["REFERENCE_DIR"]}/mlims_tables_{time_str}.pkl'
        print("writing", pkl)
        pickle_dump(pkl, d)
