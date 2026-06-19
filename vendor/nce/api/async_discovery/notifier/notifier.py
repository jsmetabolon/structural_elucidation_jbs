"""SNS notification adapters for async discovery completion and failure."""

import logging
from abc import ABC, abstractmethod
from typing import Optional, Dict

import boto3
import ujson

from api import settings


class Notifier(ABC):
    @abstractmethod
    def notify_success(self, csid: str, lib_id: int, token: str, locator: str):
        """Sends a notification that discovery succeeded."""
        pass

    @abstractmethod
    def notify_failure(self, csid: str, lib_id: int, token: str, message: str) -> None:
        """Sends a notification that discovery failed."""
        pass


class SnsNotifier(Notifier):
    def __init__(self) -> None:
        super().__init__()
        self._logger = logging.getLogger("notifier")
        self._enabled = settings.CORE_AWS_ENVIRONMENT is not None
        if self._enabled:
            region = settings.CORE_AWS_REGION
            if not region:
                raise Exception("CORE_AWS_REGION must be specified.")

            self._logger.info(f"CORE_AWS_REGION is {region}")
            self._client = boto3.client("sns", region_name=region)
            self._complete_arn = self.get_arn("AutocurationComplete")
            self._failure_arn = self.get_arn("AutocurationFailed")

            self._logger.info(f"Publishing curation complete messages to {self._complete_arn}")
            self._logger.info(f"Publishing curation failure messages to {self._failure_arn}")
        else:
            self._logger.warning("CORE_AWS_ENVIRONMENT is not set, SNS notifications disabled.")

    def notify_success(self, csid: str, lib_id: int, token: str, locator: str) -> None:
        if self._enabled:
            self._logger.info(
                f"Sending discovery complete to {self._complete_arn} for CSID {csid}, lib_id {lib_id}, and token {token}. "
                f"FMS locator is {locator}"
            )
            msg: Dict[str, str] = self.to_message(csid, lib_id, token, locator)
            msg_str: str = ujson.dumps(msg)
            response = self._client.publish(TopicArn=self._complete_arn, Message=msg_str)
            message_id = response["MessageId"]
            self._logger.info(f"Discovery complete message for {token} is {message_id}")
        else:
            self._logger.debug(
                f"Not sending complete notification to SNS for csid {csid}/ lib_id {lib_id}/token {token}"
            )

    def notify_failure(self, csid: str, lib_id: int, token: str, message: str) -> None:
        if self._enabled:
            self._logger.warning(
                f"Sending discovery failure to {self._failure_arn} for CSID {csid}, lib_id {lib_id}, and token {token}. "
                f"Message is {message}"
            )
            msg: Dict[str, str] = self.to_message(csid, lib_id, token, message=message)
            msg_str: str = ujson.dumps(msg)
            response = self._client.publish(TopicArn=self._failure_arn, Message=msg_str)
            message_id = response["MessageId"]
            self._logger.info(f"Discovery failure message for {token} is {message_id}")
        else:
            self._logger.debug(
                f"Not sending failure notification to SNS for csid {csid}/lib_id {lib_id}/token {token}"
            )

    def to_message(
        self,
        csid: str,
        lib_id: int,
        token: str,
        locator: Optional[str] = None,
        message: Optional[str] = None,
    ) -> Dict[str, str]:
        msg = {"RequestToken": token, "CurationSetId": csid, "LibraryId": lib_id}
        if locator:
            msg["ResultsLocator"] = locator
        if message:
            msg["Message"] = message
        return msg

    def get_arn(self, topic_base: str) -> str:
        core_env = settings.CORE_AWS_ENVIRONMENT
        topic_name = f"{core_env}_{topic_base}"
        complete_response: Dict[str, str] = self._client.create_topic(Name=topic_name)
        return complete_response["TopicArn"]
