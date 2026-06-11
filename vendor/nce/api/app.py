import logging
import os
from logging.config import dictConfig
from flask import Blueprint

from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())

from api import app_container, settings
from api.settings import api_logging_config_dict
from api.async_discovery.endpoints.async_discoveries import ns as async_discovery_namespace
from api.async_discovery.notifier.notifier import SnsNotifier
from api.async_discovery.persister.persister import FmsPersister
from api.discovery.endpoints.discoveries import ns as discovery_namespace
from api.dispatcher.dispatcher import Dispatcher
from api.health.endpoints.health import ns as health_namespace
from api.restplus import api


app = app_container.app

# Set up logging
dictConfig(api_logging_config_dict)


def initialize_app(flask_app):
    notifier = SnsNotifier()
    persister = FmsPersister()
    app_container.dispatcher = Dispatcher(notifier, persister)

    settings.add_to_dict(settings, flask_app.config)
    blueprint = Blueprint("api", __name__, url_prefix="/api")
    api.init_app(blueprint)
    api.add_namespace(async_discovery_namespace)
    api.add_namespace(discovery_namespace)
    api.add_namespace(health_namespace)
    api.config = dict([(key, value) for (key, value) in flask_app.config.items()])
    flask_app.register_blueprint(blueprint)
    app_container.dispatcher.start()


def main():
    initialize_app(app)
    # noinspection PyPep8
    logging.info(
        f'>>>>> Starting development server at http://{app.config["FLASK_HOST"]}:{app.config["FLASK_PORT"]}/api/ <<<<<'
    )
    print("Cache ", os.getenv("CACHE_DIR"))
    print("Cache Size", os.getenv("CACHE_SIZE_MB"))
    app.run(
        debug=app.config["DEBUG"],
        port=app.config["FLASK_PORT"],
        host=app.config["FLASK_HOST"],
        threaded=True,
    )


if __name__ == "__main__":
    main()
