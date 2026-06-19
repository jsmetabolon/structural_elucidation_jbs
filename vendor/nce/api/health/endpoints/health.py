"""HTTP resources implementing shallow and deep health pings.

Returns JSON health payloads used by load balancers and orchestration probes.
"""

import logging
from flask_restplus import Resource
from ..business import perform_shallow_ping, DiscoveryDeepPing
from api.restplus import api

log = logging.getLogger(__name__)
ns = api.namespace("health", description="Shallow and deep pings")

shallow_ping = api.schema_model(
    "shallow_ping",
    {"type": "object", "additionalProperties": {"type": "object", "additionalProperties": {}}},
)

deep_ping = api.schema_model(
    "deep_ping",
    {"type": "object", "additionalProperties": {"type": "object", "additionalProperties": {}}},
)


@ns.route("/_ping")
class ShallowPing(Resource):
    @api.response(200, "Shallow ping complete", shallow_ping)
    @api.doc("ShallowPing")
    def get(self):
        """ Performing shallow ping """
        return perform_shallow_ping()


@ns.route("/_deepping", defaults={"chain": ""})
@ns.route("/_deepping/<string:chain>")
class DeepPing(Resource):
    @api.response(200, "Deep ping complete", deep_ping)
    @api.doc("DeepPing")
    def get(self, chain):
        """ Performing deep ping """
        deepping = DiscoveryDeepPing().ping(chain)
        return deepping
