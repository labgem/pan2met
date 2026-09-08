"""
A SQL backend using 'metabiantes' model for data imported from MetaCyc.
"""

import importlib.resources
import sqlite3
from typing import Iterable, List, Optional, Tuple

import aiosql
import networkx as nx

import pan2met
import pan2met.sql.metabiantes

from . import KnowledgeBase


class MetabiantesKnowledgeBase(KnowledgeBase):
    def __init__(self, config):
        self.queries = aiosql.from_str(
            importlib.resources.read_text(pan2met.sql.metabiantes, "queries.sql"),
            "sqlite3",
        )
        self.connection = sqlite3.connect(config["metabiantes"]["database"])

    def _aiosql_to_list(self, iterator: Iterable[tuple]):
        if iterator is not None:
            return [item[0] for item in iterator]

    def monomers(self) -> List[str]:
        """
        List monomer polypeptides
        """
        return self._aiosql_to_list(self.queries.get_monomers(self.connection))

    def pathways(self) -> List[str]:
        """
        List the pathways referenced by the knowledge base
        """
        return self._aiosql_to_list(self.queries.get_pathways(self.connection))

    def pathway_name(self, pathway_id: str) -> str:
        """Get the name of a pathway"""
        return self.queries.get_pathway_name(self.connection, pathway_id=pathway_id)

    def reactions_of_pathway(self, pathway_id: str) -> List[str]:
        """
        List reactions of a pathway
        """
        return self._aiosql_to_list(
            self.queries.get_reactions_of_pathway(
                self.connection, pathway_id=pathway_id
            )
        )

    def reactions(self) -> List[str]:
        """
        List all reactions referenced in the knowledge base
        """
        return self._aiosql_to_list(self.queries.get_reactions(self.connection))

    def orphan_reactions(self) -> List[str]:
        """
        List orphan reactions
        """
        return self._aiosql_to_list(self.queries.get_orphan_reactions(self.connection))

    def spontaneous_reactions(self) -> List[str]:
        """
        List spontaneous reactions
        """
        return self._aiosql_to_list(
            self.queries.get_spontaneous_reactions(self.connection)
        )

    def non_spontaneous_reactions_of_pathway(self, pathway_id: str) -> List[str]:
        """
        Non-spontaneous reactions of a pathway
        """
        return self._aiosql_to_list(
            self.queries.get_non_spontaneous_reactions_of_pathway(
                self.connection, pathway_id=pathway_id
            )
        )

    def non_orphan_non_spontaneous_reactions_of_pathway(
        self, pathway_id: str
    ) -> List[str]:
        """
        Non-orphan and non-spontaneous reactions of a pathway
        """
        return self._aiosql_to_list(
            self.queries.get_non_orphan_non_spontaneous_reactions_of_pathway(
                self.connection, pathway_id=pathway_id
            )
        )

    def enzymes_of_reaction(self, reaction_id: str) -> List[str]:
        """
        List enzymes catalyzing a reaction
        """
        return self._aiosql_to_list(
            self.queries.get_enzymes_of_reaction(
                self.connection, reaction_id=reaction_id
            )
        )

    def reactions_by_ec_number(self, ec_number: str) -> List[str]:
        """
        List reactions annotated with given EC-number
        """
        return self._aiosql_to_list(
            self.queries.get_reactions_by_ec_number(
                self.connection, ec_number=f"EC-{ec_number}"
            )
        )

    def pathway_taxonomic_range(self, pathway_id: str) -> Optional[int]:
        """
        Return a NCBI-Taxonomy Taxonomy Identifier number
        """
        tax_id = self.queries.get_pathway_taxonomic_range(
            self.connection, pathway_id=pathway_id
        )
        if tax_id is not None:
            return int(tax_id[0])

    def reaction_is_key(self, pathway_id: str, reaction_id: str) -> bool:
        """
        Return True if the reaction is a key reaction of the pathway
        """
        return self.queries.reaction_is_key(
            self.connection, pathway_id=pathway_id, reaction_id=reaction_id
        )

    def key_reactions_of_pathway(self, pathway_id: str) -> List[str]:
        """
        List all key reactions of a pathway.
        """
        return self._aiosql_to_list(
            self.queries.get_key_reactions_of_pathway(
                self.connection, pathway_id=pathway_id
            )
        )

    def pathways_with_reaction(self, reaction_id) -> List[str]:
        """
        List all pathways with the given reaction identifier.
        """
        return self._aiosql_to_list(
            self.queries.get_pathways_with_reaction(
                self.connection, reaction_id=reaction_id
            )
        )

    def count_pathways_with_reaction(self, reaction_id) -> int:
        """
        List all pathways with the given reaction identifier.
        """
        return self.queries.count_pathways_with_reaction(
            self.connection, reaction_id=reaction_id
        )

    def variants_of_pathway(self, pathway_id: str) -> List[str]:
        """
        List the variants of a pathway.
        """
        return self._aiosql_to_list(
            self.queries.get_variants_of_pathway(self.connection, pathway_id=pathway_id)
        )

    def ontology_parent_class_of_pathway(self, pathway_id: str) -> List[str]:
        """
        List the parent class of a pathway in the ontology of pathway tools
        """
        return self._aiosql_to_list(
            self.queries.get_ontology_parent_class_of_pathway(
                self.connection, pathway_id=pathway_id
            )
        )

    def pathway_reaction_order(self, pathway_id: str) -> List[tuple[str, str]]:
        """
        Return a list of pairs of reactions in the pathway, where the first reaction is a predecessor of the second reaction in the pathway.
        """
        return self.queries.get_pathway_reaction_order(
            self.connection, pathway_id=pathway_id
        )

    def species_evidence_of_pathway(self, pathway_id: str) -> List[int]:
        """
        List NCBI Taxonomy identifiers of species where the pathway presence evidence was found in the literature, according to the knowledge base.
        """
        return list(
            map(
                int,
                self._aiosql_to_list(
                    self.queries.get_pathway_species(
                        self.connection, pathway_id=pathway_id
                    )
                ),
            )
        )

    def reaction_graph_topological_order(self, pathway) -> List[str]:
        """
        Return the topological ordering of a pathway reaction graph.

        Assume that the pathway reaction graph is a directed acyclic graph.
        """
        reaction_order: List[Tuple[str, str]] = self.pathway_reaction_order(pathway)
        graph = nx.DiGraph(reaction_order)
        sort = list(next(nx.all_topological_sorts(graph)))
        return sort

    def complex(self) -> List[str]:
        """
        List the protein complex
        """
        return self._aiosql_to_list(self.queries.get_complex(self.connection))

    def components_of_complex(self, complex_id: str) -> List[str]:
        """
        List the components of a complex
        """
        return self._aiosql_to_list(
            self.queries.get_polypeptide_complex_components(
                self.connection, complex_id=complex_id
            )
        )
