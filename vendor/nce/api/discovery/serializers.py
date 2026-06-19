"""Flask-RESTPlus response models for discovery API payloads."""

from flask_restplus import fields
from api.restplus import api


discover = api.model(
    "Discovery",
    {
        #'version': fields.Raw(readOnly=True, description='Sally Version'),
        #'curation_level': fields.Raw(readOnly=True, description='Level of discovery Sally performed'),
        "results": fields.Raw(readOnly=True, description="Discovery Results as a dictionary.")
    },
)

pagination = api.model(
    "A page of results",
    {
        "page": fields.Integer(description="Number of this page of results"),
        "pages": fields.Integer(description="Total number of pages of results"),
        "per_page": fields.Integer(description="Number of items per page of results"),
        "total": fields.Integer(description="Total number of results"),
    },
)
