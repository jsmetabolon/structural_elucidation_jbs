"""Shared Flask app instance and request-level response caching hooks.

Holds the module-level ``app`` and ``dispatcher`` singletons imported
by API endpoints and the async discovery worker thread.
"""

import logging
from typing import Any, Optional

from flask import Flask, request, Response

from api.dispatcher.dispatcher import Dispatcher
from cache import cache_util

app: Flask = Flask(__name__)
dispatcher: Optional[Dispatcher] = None
cache_logger = logging.getLogger("flask_caching")


def is_cacheable(full_path: str) -> bool:
    cacheable = full_path.startswith("/api/discovery/") or full_path.startswith("/api/discovery/")
    cache_logger.info(f"full path is {full_path}, cacheable = {cacheable}")
    ######################
    return False  # cacheable


@app.before_request
def before_request() -> Any:
    key = request.full_path
    if is_cacheable(key):
        cache_logger.info(f"Request is cacheable, attempting to read cached value for key {key}")
        return cache_util.read_cacheable(key, lambda: None)
    else:
        return None


@app.after_request
def after_request(response: Response) -> Response:
    key = request.full_path
    if is_cacheable(key) and cache_util.cache_enabled() and response.status_code == 200:
        logging.info(f"Request successful, storing cache value for key {key}")
        cache_util.cache_object(key, response.data)
    return response
