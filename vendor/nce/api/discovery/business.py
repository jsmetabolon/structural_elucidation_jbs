from flask import current_app
from api.testing_miner_on_api import discover
from api import settings


def discovery_for_csid(csid: str, lib_id: int) -> tuple:
    discovery_config = settings.discovery_config
    discovery_config["lib_id"] = lib_id
    discovery_config["ssid"] = csid
    return discover(csid, lib_id, discovery_config, current_app.config)
