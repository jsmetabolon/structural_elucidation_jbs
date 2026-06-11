import logging

from flask_restplus import Resource

from api.discovery.business import discovery_for_csid
from api.discovery.serializers import discover
from api.restplus import api

log = logging.getLogger(__name__)
ns = api.namespace(
    "discovery", description="Discovery - opinions on which peaks might be new compounds"
)


@ns.route("/by-csid/<csid>/<lib_id>")
@api.doc(params={"csid": "Curation set id", "lib_id": "Library id"})
@api.response(404, "Discovery not found.")
class Discoveries(Resource):
    @api.marshal_with(discover)
    def get(self, csid: str, lib_id: int) -> dict:
        """Returns a curated discovery set."""
        logging.info(f"logging.info in Discoveries {csid}, {lib_id}")
        results = discovery_for_csid(csid, lib_id)
        return {"results": results}
