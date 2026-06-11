import logging
import os
from queue import Queue
from threading import Thread

import ujson

from api.testing_miner_on_api import discover
from cache import cache_util
from .discovery_results import DiscoveryResults
from .unit_of_work import UnitOfWork
from api.async_discovery.notifier.notifier import Notifier
from api.async_discovery.persister.persister import Persister
import api.settings as settings


class Dispatcher(object):
    """Simple consumer that drains UnitOfWork objects from a work queue. Consumption starts once start() is called."""

    def __init__(self, notifier: Notifier, persister: Persister):
        self._logger = logging.getLogger("async_dispatcher")
        self._notifier = notifier
        self._persister = persister
        self._work_queue = Queue()
        self._work_thread = Thread(target=self.dispatch, daemon=True)

    def start(self) -> None:
        """Starts processing."""
        self._logger.info("Starting async discovery dispatcher thread")
        self._work_thread.start()

    def dispatch(self) -> None:
        self._logger.info("Starting async discovery dispatcher loop")
        while True:
            # noinspection PyBroadException
            try:
                unit: UnitOfWork = self.work_queue.get(True)
                self._logger.info(
                    f"Starting unit of work - CSID {unit.csid}, Lib_id: {unit.lib_id}, token {unit.token}"
                )
                self.__process(unit)
                self._logger.info(
                    f"Completed unit of work - CSID {unit.csid}, Lib_id: {unit.lib_id}, token {unit.token}"
                )
            except:
                self._logger.exception("Error processing unit of work:")

    def __process(self, uow: UnitOfWork) -> None:
        """Performs autocuration on the given unit of work."""
        try:
            key = f"results-{uow.csid}-{uow.lib_id}-"

            def discover_uow() -> DiscoveryResults:
                result_dict = discover(uow.csid, uow.lib_id, uow.config, uow.app_config)
                # result_dict = curation_df_to_result_dict(result)
                return DiscoveryResults(uow.csid, uow.lib_id, result_dict)

            discovery_results: DiscoveryResults = cache_util.read_cacheable(key, discover_uow)
            self.__cache_local_results(uow, discovery_results)
            locator = self._persister.persist(uow.token, discovery_results)
            self._notifier.notify_success(uow.csid, uow.lib_id, uow.token, locator)
        except BaseException as ex:
            self._logger.exception(f"Error processing {uow.token}:")
            self._notifier.notify_failure(uow.csid, uow.lib_id, uow.token, str(ex))
            raise

    def __cache_local_results(self, uow: UnitOfWork, discovery_results: DiscoveryResults) -> None:
        # noinspection PyBroadException
        try:
            if settings.LOCAL_RESULT_CACHE:
                local_output_file = os.path.join(
                    settings.LOCAL_RESULT_CACHE, f"{uow.token}-results.json"
                )
                with open(local_output_file, "w") as out_file:
                    self._logger.debug(f"Writing output for {uow.token} to {local_output_file}")
                    ujson.dump(discovery_results.results, out_file)
                    self._logger.info(f"Wrote output for {uow.token} to {local_output_file}")
        except:
            # An error here shouldn't interfere with anything, just log and move on
            self._logger.exception(f"Error writing results to local fs for {uow.token}:")

    @property
    def work_queue(self) -> Queue:
        return self._work_queue
