"""File Management Service persisters for async discovery result payloads."""

import logging
import os
import socket
import tempfile
from abc import ABC, abstractmethod
from typing import Dict

import requests
import ujson

from api import settings
from api.dispatcher.discovery_results import DiscoveryResults


class Persister(ABC):
    def __init__(self) -> None:
        super().__init__()

    @abstractmethod
    def persist(self, token: str, data: DiscoveryResults):
        """Upload discovery results and return an FMS locator string."""
        pass


class FmsPersister(Persister):
    def __init__(self) -> None:
        super().__init__()
        if not settings.FMS_ENDPOINT:
            raise ValueError("FMS_ENDPOINT configuration setting is undefined")
        self._fms_endpoint = settings.FMS_ENDPOINT
        self._logger = logging.getLogger("persister")

    def persist(self, token: str, data: DiscoveryResults) -> str:
        self._logger.info(f"Persisting discovery results for {token}")

        params = self.create_fms_params(token, data)
        filename = params["fileName"]
        with tempfile.TemporaryDirectory() as temp_dir:
            result_file = os.path.join(temp_dir, filename)
            self._logger.debug(f"Writing results to output file {result_file}")
            with open(result_file, "w") as f:
                ujson.dump(data.results, f)

            self._logger.debug("Sending results to FMS")
            with open(result_file, "r") as f:
                files_dict = {"file": (filename, f, "text/plain")}

                url = f"{self._fms_endpoint}/api/v1/File/PutFile"
                response = requests.post(url, params=params, files=files_dict)
                self._logger.debug(f"Response code: {response.status_code}. Text: {response.text}")
                if response.status_code == requests.codes.ok:
                    locator = ujson.loads(response.content)
                    self._logger.info(
                        f"Persisted curation results for {token} - locator is {locator}"
                    )
                    return locator
                else:
                    response.raise_for_status()
                    raise RuntimeError(f"Got non-ok response from FMS {response.content}")

    def create_fms_params(self, token: str, data: DiscoveryResults) -> Dict[str, str]:
        metadata = {"csid": data.csid, "lib_id": str(data.lib_id), "token": token}
        metadata_str = ujson.dumps(metadata)
        filename = f"discovery-{token}.json"
        params = {
            "source": "Discovery",
            #'sourceVersion': data.version,
            "sourceMachine": socket.gethostname(),
            "fileName": filename,
            "fileType": "Discovery Results",
            "metadataJsonStr": metadata_str,
        }
        return params
