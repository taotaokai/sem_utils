import collections

from obspy.clients.fdsn.client import get_bulk_string
from obspy.clients.fdsn.routing.federator_routing_client import FederatorRoutingClient


class FederatorCatalog(FederatorRoutingClient):
    def __init__(self, *args, **kwargs):
        super(FederatorCatalog, self).__init__(*args, **kwargs)

    def query_fedcat(self, bulk=[], **kwargs):
        """
        Get stations from multiple data centers.

        It will pass on most parameters to the federated routing service.
        They will also be passed on to the individual FDSNWS implementations
        if a service supports them.

        This can route on a number of different parameters, please see the
        web site of the `EarthScope (former IRIS) Federator
        <https://service.earthscope.org/irisws/fedcatalog/1/>`_ for details.
        """
        if not bulk:
            selection = [
                kwargs.pop(key, "*")
                for key in (
                    "network",
                    "station",
                    "location",
                    "channel",
                    "starttime",
                    "endtime",
                )
            ]
            bulk = [selection]

        bulk = [["*" if x is None else x for x in req ] for req in bulk]

        params = collections.OrderedDict()
        for k in self.kwargs_of_interest:
            if k in kwargs:
                params[k] = str(kwargs[k])
        params["format"] = "request"

        bulk_str = get_bulk_string(bulk, params)
        r = self._download(
            self._url + "/query", data=bulk_str, content_type="text/plain"
        )
        data = r.content.decode() if hasattr(r.content, "decode") else r.content
        split = self._split_routing_response(data)
        return split

    @staticmethod
    def _split_routing_response(data):
        """
        Splits the routing responses per data center for the federator output.

        Returns a dictionary with the keys being the datacenter name, and dictionary
        of service url and requests as the values.

        :param data: The return value from the earthscope Federator routing service.
        """
        split = collections.defaultdict(dict)
        current_key = None
        for line in data.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith("DATACENTER"):
                datacenter_info = line.split("=")[1]
                current_key, url = datacenter_info.split(",")
                split[current_key] = {}
                continue
            # Anything before the first data center can be ignored.
            if current_key is None:
                continue
            if "http://" in line or "https://" in line:
                service, url = line.split("=")
                split[current_key][service] = url
                continue
            if "requests" not in split[current_key]:
                split[current_key]["requests"] = []
            split[current_key]["requests"].append(line.split())

        return split
