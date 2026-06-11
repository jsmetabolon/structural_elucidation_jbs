from typing import Any


class UnitOfWork(object):
    """Information about a single request to be processed."""

    def __init__(self, csid: str, lib_id: int, config: dict, app_config: Any, token: str):
        """
        :param csid: Curation set id, cannot be None.
        :param token: Token to identify this UOW, cannot be None.
        """
        if not csid:
            raise ValueError("csid cannot be None")
        if not lib_id:
            raise ValueError("lib_id cannot be None")
        if not config:
            raise ValueError("config cannot be None or empty")
        if not token:
            raise ValueError("token cannot be None")
        self._csid = csid
        self._lib_id = lib_id
        self._config = config
        self._app_config = app_config
        self._token = token

    @property
    def csid(self) -> str:
        """Curation set id."""
        return self._csid

    def lib_id(self) -> str:
        """Curation set id."""
        return self._lib_id

    @property
    def config(self) -> dict:
        return self._config

    @property
    def app_config(self) -> Any:
        return self._app_config

    @property
    def token(self) -> str:
        """Token returned to the creator to identify this UOW."""
        return self._token
