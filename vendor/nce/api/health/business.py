from healthserver.shallow.shallow_ping import ShallowPing
from healthserver.deep.deep_ping_base import DeepPingBase
from healthserver.deep.deep_ping import IDeepPing

from curationservice.client.v1.metabolon_curation_service_api import MetabolonCurationServiceAPI
from metadataservice.client.v1.metabolon_metadata_service_api import MetabolonMetadataServiceAPI
from ismatching.client.is_matching_service import IsMatchingService
from libraryservice.client.metabolon_library_service_api import MetabolonLibraryServiceAPI

from api import settings


def perform_shallow_ping():
    shallow_ping = ShallowPing(settings.VERSION)
    return shallow_ping.ping()


class DiscoveryDeepPing(DeepPingBase):
    def __init__(self):
        super(DiscoveryDeepPing, self).__init__(settings.VERSION)

    @property
    def dependencies(self):
        return [  # CurationServiceDeepPingDependency(),
            MetadataServiceDeepPingDependency(),
            IsMatchingServiceDeepPingDependency(),
            LibraryServiceDeepPingDependency(),
        ]

    @property
    def service_name(self):
        return "Sally Service"


# class CurationServiceDeepPingDependency(IDeepPing):
#     def ping(self, chain):
#         service = MetabolonCurationServiceAPI(settings.CURATION_SERVICE_ENDPOINT)
#         return service.deep_ping(chain)
#
#     @property
#     def service_name(self):
#         return "Curation Service"


class MetadataServiceDeepPingDependency(IDeepPing):
    def ping(self, chain):
        service = MetabolonMetadataServiceAPI(settings.METADATA_SERVICE_URI)
        return service.deep_ping(chain)

    @property
    def service_name(self):
        return "Metadata Service"


class IsMatchingServiceDeepPingDependency(IDeepPing):
    def ping(self, chain):
        service = IsMatchingService(settings.ISMATCHING_SERVICE_URI)
        response = service.deep_ping1(chain)
        return response

    @property
    def service_name(self):
        return "IS Matching Service"


class LibraryServiceDeepPingDependency(IDeepPing):
    def ping(self, chain):
        service = MetabolonLibraryServiceAPI(settings.LIBRARY_SERVICE_ENDPOINT)
        return service.deep_ping(settings.VERSION, chain)

    @property
    def service_name(self):
        return "Library Service"
