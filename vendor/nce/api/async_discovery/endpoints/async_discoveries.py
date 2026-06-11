import logging
import uuid

from flask import current_app
from flask_restplus import Resource, reqparse
from flask_restplus.reqparse import ParseResult
from api import settings
from api import app_container
from api.dispatcher.unit_of_work import UnitOfWork
from api.restplus import api

ns = api.namespace(
    "async_discovery", description="Asynchronously requests discovery on a discovery set."
)

parser = reqparse.RequestParser()
parser.add_argument("csid", required=True, type=str, help="Curation set id")
parser.add_argument("lib_id", required=True, type=int, help="Library id")  # action='append',


@ns.route("/by-csid-async/")
@api.doc(parser=parser)
@api.response(200, "Discovery request queued.")
class AsyncDiscoveries(Resource):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._logger = logging.getLogger("api_service")

    def get(self) -> str:
        """Requests that the given discovery set be asynchronously curated."""
        args: ParseResult = parser.parse_args()
        csid: str = args["csid"]
        lib_id: int = args["lib_id"]

        config: dict = settings.discovery_config
        config["csid"] = args["csid"]
        config["lib_id"] = lib_id

        token: str = uuid.uuid4().hex
        unit_of_work = UnitOfWork(csid, lib_id, config, current_app.config, token)
        self._logger.info(f"Queuing discovery of {csid} with token {token}")
        app_container.dispatcher.work_queue.put(unit_of_work)
        return token
