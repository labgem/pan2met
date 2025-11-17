"""
SPARQL Endpoint from BioPAX data
"""
import os

from SPARQLWrapper import SPARQLWrapper, JSON

class SPARQL:

    def __init__(self):
        self.ensure_endpoint_config()
        self.wrapper = SPARQLWrapper(
            os.environ["SPARQL_ENDPOINT"]
        )
        self.wrapper.setReturnFormat(JSON)

    def ensure_endpoint_config(self):
        for variable in ["SPARQL_ENDPOINT", "SPARQL_USER", "SPARQL_PASSWORD"]:       
            assert (variable in os.environ
                    and os.environ[variable] != ""), f"`{variable}' environment variable should be set"


