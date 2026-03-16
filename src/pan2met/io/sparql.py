"""
A knowledge base on Pathway data based on a SPARQL Endpoint from BioPAX data
"""

import pan2met

import os
import importlib.resources

from SPARQLWrapper import SPARQLWrapper, JSON
import ouisparql


from ..utils import logger
from .knowledge_base import KnowledgeBase


class SPARQL:
    def __init__(self):
        self.wrapper = SPARQLWrapper(os.environ["SPARQL_ENDPOINT"])
        self.wrapper.setReturnFormat(JSON)
        logger.info(
            f"Using SPARQL endpoint: {os.environ['SPARQL_ENDPOINT']} as a source of knowledge."
        )


class SPARQLBackendKnowledgeBase(KnowledgeBase):
    def __init__(self):
        queries_str: str = importlib.resources.read_text(pan2met, "sparql/queries.rq")
        self.queries = ouisparql.from_str(queries_str, "sparql_wrapper")
        self.sparql_wrapper = SPARQL().wrapper

    def pathways(self):
        return self.queries.list_pathways(self.sparql_wrapper)

    def reactions_of_pathway(self, pathway_id: str) -> list[str]:
        raise NotImplementedError()

    def reaction_enzymes(self, reaction_id: str) -> list[str]:
        raise NotImplementedError()

    def reactions_by_ec_number(self, ec_number: str) -> list[str]:
        ec_number_literal = ouisparql.utils.to_string_literal(
            ec_number
        )  # convert 1.1.1.1 to "1.1.1.1"
        result = self.queries.reactions_by_ec_number(
            self.sparql_wrapper, ec_number=ec_number_literal
        )
        return list(map(lambda item: item["reactionId"]["value"], result))
