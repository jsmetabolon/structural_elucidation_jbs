"""Container for completed discovery results keyed by curation set and library."""


class DiscoveryResults:
    def __init__(self, csid: str, lib_id: int, results: dict):
        self._csid = csid
        self._lib_id = lib_id
        self._results = results

    @property
    def csid(self) -> str:
        return self._csid

    @property
    def lib_id(self) -> int:
        return self._lib_id

    @property
    def results(self) -> dict:
        return self._results
