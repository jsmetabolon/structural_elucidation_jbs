"""Integration harness that runs miner discovery through API providers.

Wires miner MIO providers against live service endpoints for manual or
dispatcher-driven discovery smoke tests.
"""

import os
import logging
from api import settings
from api.settings import discovery_config

from miner.mio.api.task_provider_api import create_task_info_provider
from miner.mio.api.internal_standards_provider_api import create_is_matching_provider
from miner.mio.api.library_provider_api import create_api_library_provider
from miner.mio.api.peak_provider_api import create_cpd_discovery_peak_provider
from flask import Flask, current_app

# noinspection PyUnresolvedReferences
import miner.initPackageOptions

from dotenv import load_dotenv

load_dotenv()


def discover(csid, lib_id, discovery_config, api_config):

    logging.getLogger("flask_caching").setLevel(logging.DEBUG)

    logging.info(f"Run miner_on_api: discover: csid: {csid}, lib_id: {lib_id}")

    library_provider = create_api_library_provider(api_config)
    get_name_dict = library_provider.get_name_dict
    get_library_ions = library_provider.get_library_ions
    get_chemical_formulas = library_provider.get_chemical_formulas
    get_known_ms2 = library_provider.get_known_ms2
    get_unknown_ms2 = library_provider.get_unknown_ms2
    get_chemicals_for_formula = library_provider.get_chemicals_for_formula
    get_smiles = library_provider.get_smiles

    task_info_provider = create_task_info_provider(api_config)
    get_task_info = task_info_provider.get_task_info

    is_matching_provider = create_is_matching_provider(api_config)
    get_internal_standards_for_csid = is_matching_provider.get_internal_standards_for_csid
    get_internal_standards = is_matching_provider.get_internal_standards

    peak_provider = create_cpd_discovery_peak_provider(api_config)
    get_peaks = peak_provider.get_peaks

    dct = get_name_dict(lib_id)

    chems = get_chemicals_for_formula("C5H10O3", lib_id)
    chems = get_chemicals_for_formula("C2H6O2", lib_id)

    # Get group of internal standards
    # stds = get_internal_standards_for_csid(csid, lib_id)
    #
    # tasks = [
    #     "04788579-9e9c-4ff3-8a0c-9bfaf108c7ae",
    #     "fea58338-4c60-49fb-9b0f-320aa0abe32d",
    #     "627b9ccc-8edb-413f-a3f6-35d2236beef0",
    # ]
    # # Now get stds for individual tasks
    # for i in range(3):
    #     iss = get_internal_standards(tasks[i])
    #
    # # get_library_ions
    # ions = get_library_ions(lib_id)
    #
    # # # get_known_ms2
    # ms2 = get_known_ms2(lib_id)
    #
    # # # get_unknown_ms2
    # ms2 = get_unknown_ms2(lib_id)
    #
    # Tasks
    tasks = get_task_info(csid, [lib_id])

    # chemical_formulas
    ans = get_chemical_formulas([209])
    print(ans)

    peaks = get_peaks(tasks.iloc[0:5, :], lo_mass=100, hi_mass=300, lo_ri=2000, hi_ri=3000)

    ans = get_smiles()

    return ans


# Still not working
if __name__ == "__main__":
    # Define the values in these lines

    import pprint

    pp = pprint.PrettyPrinter()

    print("1", os.environ["LIBRARY_SERVICE_ENDPOINT"])
    dct = dict()
    settings.add_to_dict(settings, dct)
    print(dct)
    pp.pprint(dct)

    csid = "532e0da5-c125-4f16-bb56-363271912235"  # Sample Set ID to curate
    lib_id = 400
    discovery_config["study"] = csid
    discovery_config["lib_id"] = lib_id

    pp.pprint(f"settings\n {settings}")
    app = Flask(__name__)
    with app.app_context():
        print("settings", settings.LIBRARY_SERVICE_ENDPOINT)
        settings.add_to_dict(settings, current_app.config)
        pp.pprint(current_app.config)
        print("2", current_app.config["LIBRARY_SERVICE_ENDPOINT"])
        result = discover(
            csid, lib_id, discovery_config=discovery_config, api_config=current_app.config
        )

    print(f"Done with Work, result is {type(result)}")
