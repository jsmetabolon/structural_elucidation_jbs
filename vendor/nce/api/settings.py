"""Environment-driven configuration for the Compound Discovery API.

Loads Flask, downstream service URIs (metadata, curation, is-matching,
library, FMS), discovery defaults, and logging dict used at app startup.
"""

import os
from os import getenv


# Flask settings
# FLASK_SERVER_NAME = getenv('SERVER_NAME', 'localhost:8083')
FLASK_HOST = getenv("FLASK_HOST", "localhost")
FLASK_PORT = getenv("FLASK_PORT", 8083)
FLAS_ENV = "development"
DEBUG = getenv("DEBUG", "true").lower() == "true"  # Do not use debug mode in production
VERSION = "1.0"

# Flask-Restplus settings
SWAGGER_UI_DOC_EXPANSION = "list"
RESTPLUS_VALIDATE = True
RESTPLUS_MASK_SWAGGER = False
ERROR_404_HELP = False

DATADIR = getenv("DATADIR", os.getcwd())

CACHE_DIR = getenv("CACHE_DIR", DATADIR)
CACHE_EXPIRE = int(getenv("CACHE_EXPIRE", 86400))
CACHE_SIZE_MB = int(getenv("CACHE_SIZE_MB", 0))

discovery_config = {
    "study": None,
    "library": 209,
    "ri_lo": 0,
    "ri_hi": 10000,
    "min_ri_window": 50,
    "ri_window_fraction": 0.75,
    "p_value": 0.05,
    "min_ms2_intensity": 0.01,
    "out_dir": DATADIR,
}
# If set, the sally API will write results to this location.
LOCAL_RESULT_CACHE = getenv("LOCAL_RESULT_CACHE")

# metadata service settings
METADATA_SERVICE_HOST = getenv("METADATA_SERVICE_HOST", "")
METADATA_SERVICE_PORT = getenv("METADATA_SERVICE_PORT", 80)
METADATA_SERVICE_USE_HTTPS = getenv("METADATA_SERVICE_USE_HTTPS", "false").lower() == "true"
METADATA_SERVICE_URI = f'{("https://" if METADATA_SERVICE_USE_HTTPS else "http://")}{METADATA_SERVICE_HOST}:{METADATA_SERVICE_PORT}'
METADATA_SERVICE_TIMEOUT_SECONDS = int(getenv("METADATA_SERVICE_TIMEOUT_SECONDS", 10))

# discovery service settings
CURATION_SERVICE_ENDPOINT = getenv("CURATION_SERVICE_ENDPOINT", "")
CURATION_SERVICE_HOST = getenv("CURATION_SERVICE_HOST", "")
CURATION_SERVICE_PORT = getenv("CURATION_SERVICE_PORT")
CURATION_SERVICE_USE_HTTPS = getenv("CURATION_SERVICE_USE_HTTPS", "false").lower() == "true"
CURATION_SERVICE_URI = f'{("https://" if CURATION_SERVICE_USE_HTTPS else "http://")}{CURATION_SERVICE_HOST}:{CURATION_SERVICE_PORT}'
CURATION_SERVICE_TIMEOUT_SECONDS = int(getenv("CURATION_SERVICE_TIMEOUT_SECONDS", 20))

# IsMatching service settings
ISMATCHING_SERVICE_HOST = getenv("ISMATCHING_SERVICE_HOST", "")
ISMATCHING_SERVICE_PORT = getenv("ISMATCHING_SERVICE_PORT", 80)
ISMATCHING_SERVICE_USE_HTTPS = getenv("ISMATCHING_SERVICE_USE_HTTPS", "false").lower() == "true"
ISMATCHING_SERVICE_TIMEOUT_SECONDS = int(getenv("ISMATCHING_SERVICE_TIMEOUT_SECONDS", "600"))
ISMATCHING_SERVICE_URI = f'{("https://" if ISMATCHING_SERVICE_USE_HTTPS else "http://")}{ISMATCHING_SERVICE_HOST}:{ISMATCHING_SERVICE_PORT}/api/v{VERSION}'

# Library Service settings
LIBRARY_SERVICE_ENDPOINT = getenv("LIBRARY_SERVICE_ENDPOINT", "")

AUTH0_LIB_ENDPOINT = getenv("AUTH0_LIB_ENDPOINT")
AUTH0_LIB_AUDIENCE = getenv("AUTH0_LIB_AUDIENCE")
AUTH0_LIB_CLIENT_ID = getenv("AUTH0_LIB_CLIENT_ID")
AUTH0_LIB_CLIENT_SECRET = getenv("AUTH0_LIB_CLIENT_SECRET")

# Fit score service
FIT_SCORE_SERVICE_HOST = getenv("FIT_SCORE_SERVICE_HOST", "")
FIT_SCORE_SERVICE_PORT = getenv("FIT_SCORE_SERVICE_PORT")
FIT_SCORE_SERVICE_USE_HTTPS = getenv("FIT_SCORE_SERVICE_USE_HTTPS", "false").lower() == "true"
FIT_SCORE_SERVICE_BATCH_SIZE = int(getenv("FIT_SCORE_SERVICE_BATCH_SIZE", "1000"))
FIT_SCORE_SERVICE_TIMEOUT_SECONDS = int(getenv("FIT_SCORE_SERVICE_TIMEOUT_SECONDS", "1800"))

# File management service
# Real Value = http://filemanagementservice-alpha.apps.dev.metabolon.com
FMS_ENDPOINT = getenv("FMS_ENDPOINT", "http://localhost")

# Core integration settings
CORE_AWS_ENVIRONMENT = getenv("CORE_AWS_ENVIRONMENT")
CORE_AWS_REGION = getenv("CORE_AWS_REGION")

api_logging_config_dict = {
    "version": 1,
    "formatters": {"default": {"format": "[%(asctime)s] %(levelname)s in %(module)s: %(message)s"}},
    "handlers": {
        "wsgi": {
            "class": "logging.StreamHandler",
            "stream": "ext://flask.logging.wsgi_errors_stream",
            "formatter": "default",
        }
    },
    "root": {"level": "INFO", "handlers": ["wsgi"]},
}


def add_to_dict(settings_module, config_dict):
    """Copy uppercase module-level settings into a Flask ``config`` dict."""
    for setting in dir(settings_module):
        # you can write your filter here
        if setting == "LIBRARY_SERVICE_ENDPOINT":
            print("setting LSE:", LIBRARY_SERVICE_ENDPOINT)
            print("environ LSE:", os.environ["LIBRARY_SERVICE_ENDPOINT"])
        if setting.isupper():
            config_dict[setting] = getattr(settings_module, setting)
    print("final LSE:", config_dict["LIBRARY_SERVICE_ENDPOINT"])
