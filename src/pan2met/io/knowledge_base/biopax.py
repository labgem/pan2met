"""
A SPARQL backend using the 'BioPAX' model.
"""

import importlib.resources
import logging
from collections.abc import Iterable

import ouisparql
from SPARQLWrapper import GET, JSON, SPARQLWrapper

import pan2met
import pan2met.sparql.biopax

from . import KnowledgeBase

logger = logging.getLogger()


def ouisparql_list(response, key) -> list[str]:
    return [item[key]["value"] for item in response]


def format_sparql_uri(uri: str) -> str:
    return f"<{uri}>"


def format_sparql_string_literal(string: str) -> str:
    return f'"{string}"'


ENTITIES = [
    "bp:PhysicalEntity",
    "bp:Pathway",
    "bp:Catalysis",
    "bp:Complex",
    "bp:SmallMolecule",
    "bp:Gene",
    "bp:Pathway",
    "bp:Protein",
    "bp:BiochemicalReaction",
]

DATABASES = ["MetaCyc", "HumanCyc", "EcoCyc", "ChEBI"]  # partial.


class BioPAXSparqlKnowledgeBase(KnowledgeBase):
    def __init__(self, config):
        self.ENDPOINT = config["reference"]["biopax"]["SPARQL_ENDPOINT"]
        self.queries = ouisparql.from_str(
            importlib.resources.read_text(pan2met.sparql.biopax, "queries.rq"),
            "sparql_wrapper",
        )
        self.sparql_wrapper = SPARQLWrapper(
            self.ENDPOINT,
        )
        self.sparql_wrapper.setMethod(GET)
        self.sparql_wrapper.setReturnFormat(JSON)

    def pathways(self) -> list[str]:
        response = self.queries.get_all_pathways(self.sparql_wrapper)
        return ouisparql_list(response, "pathway")

    def get_pathway_standard_name(self, pathway_id: str) -> str:
        response = self.queries.get_standard_name_of_pathway(
            self.sparql_wrapper, pathway_id=format_sparql_uri(pathway_id)
        )
        value = next(response)
        return value["pathway_name"]["value"]

    def get_pathway(self, pathway_id: str) -> dict:
        response = self.queries.get_pathway_summary(
            self.sparql_wrapper, pathway_id=format_sparql_uri(pathway_id)
        )
        return next(response)

    def match_pathway_name(self, pattern: str) -> Iterable[dict[str, str]]:
        response = self.queries.get_pathway_by_partial_text_match(
            self.sparql_wrapper, pattern=f'"{pattern}"'
        )
        return response

    def match_compound_name(self, pattern: str) -> Iterable[dict[str, str]]:
        response = self.queries.get_compound_by_partial_text_match(
            self.sparql_wrapper, pattern=f'"{pattern}"'
        )
        return response

    def match_reaction_name(self, pattern: str) -> Iterable[dict[str, str]]:
        response = self.queries.get_reaction_by_partial_text_match(
            self.sparql_wrapper, pattern=f'"{pattern}"'
        )
        return response

    def search_pathway_by_compounds(
        self, compounds: list[str]
    ) -> Iterable[dict[str, str]]:
        def generate_rdf_triplets(compounds: list[str]) -> Iterable[str]:
            for compound_index, compound_uri in enumerate(compounds):
                yield f"""
                ?pathway bp:pathwayComponent ?compound_{compound_index}_reaction .
                ?compound_{compound_index}_reaction rdf:type bp:BiochemicalReaction ;
                          bp:left|bp:right <{compound_uri}> .
                """

        compound_to_pathway_triplets = "\n".join(generate_rdf_triplets(compounds))
        print(compound_to_pathway_triplets)
        response = self.queries.get_pathway_by_partial_compound_match(
            self.sparql_wrapper,
            compound_to_pathway_triplets=compound_to_pathway_triplets,
        )
        return response

    def count_entity_type(self, entity_type: str) -> int:
        response = self.queries.count_entities(
            self.sparql_wrapper, entity_type=entity_type
        )
        result = next(response)
        return result["count"]["value"]

    def get_reaction_standard_name(self, reaction_id: str) -> str:
        response = self.queries.get_standard_name_of_reaction(
            self.sparql_wrapper, reaction_id=format_sparql_uri(reaction_id)
        )
        value = next(response)
        return value["reaction_name"]["value"]

    def get_compound_description(self, compound_id: str) -> str:
        response = self.queries.get_compound_description(
            self.sparql_wrapper, compound_id=format_sparql_uri(compound_id)
        )
        return next(response)

    def get_entity_by_reference_id(
        self, entity_type: str, reference_id: str, reference_database: str
    ) -> list[str]:
        response = self.queries.get_entity_by_reference_id(
            self.sparql_wrapper,
            entity_type=entity_type,
            reference_id=format_sparql_string_literal(reference_id),
            reference_database=format_sparql_string_literal(reference_database),
        )
        return list(response)

    def get_pathway_by_reference_id(self, pathway_id: str, database: str) -> list[str]:
        return self.get_entity_by_reference_id("bp:Pathway", pathway_id, database)

    def database_summary_counts(self) -> dict[str, int]:
        """
        Count how many entities of some types are found in the triplestore.
        """

        summary_counts = {
            entity_type: self.count_entity_type(entity_type) for entity_type in ENTITIES
        }
        return summary_counts
