"""
A knowledge base on Pathway data based on a SPARQL Endpoint from BioPAX data
"""

import os

from SPARQLWrapper import SPARQLWrapper, JSON
import ouisparql


from ..utils import logger
from ..knowledge_base import KnowledgeBase

SPARQL_QUERY_FILE_PATH = os.path.join(
    os.path.dirname(__file__), "../../sparql/queries.rq"
)


class SPARQL:
    def __init__(self):
        self.wrapper = SPARQLWrapper(os.environ["SPARQL_ENDPOINT"])
        self.wrapper.setReturnFormat(JSON)
        logger.info(
            f"Using SPARQL endpoint: {os.environ['SPARQL_ENDPOINT']} as a source of knowledge."
        )


class SPARQLBackendKnowledgeBase(KnowledgeBase):
    def __init__(self):
        self.queries = ouisparql.from_path(SPARQL_QUERY_FILE_PATH, "sparql_wrapper")
        self.sparql_wrapper = SPARQL().wrapper

    def pathways(self):
        return self.queries.list_pathways(self.sparql_wrapper)

    def reactions_of_pathway(self, pathway_id: str) -> list[str]:
        pass

    def reaction_enzymes(self, reaction_id: str) -> list[str]:
        pass

    def reactions_by_ec_number(self, ec_number: str) -> list[str]:
        ec_number_literal = ouisparql.utils.to_string_literal(
            ec_number
        )  # convert 1.1.1.1 to "1.1.1.1"
        result = self.queries.reactions_by_ec_number(
            self.sparql_wrapper, ec_number=ec_number_literal
        )
        return list(map(lambda item: item["reactionId"]["value"], result))
