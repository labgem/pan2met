"""
SPARQL Endpoint from BioPAX data
"""
import os

from SPARQLWrapper import SPARQLWrapper, JSON
import ouisparql

from .sql import PathwayKnowledge

SPARQL_QUERY_FILE_PATH = os.path.join(os.path.dirname(__file__), '../../sparql/queries.rq')

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
          

class SPARQLPathwayKnowledge(PathwayKnowledge):

    def __init__(self):

        self.queries = ouisparql.from_path(SPARQL_QUERY_FILE_PATH)
