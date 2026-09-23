"""
A SQL backend using 'metabiantes' model for data imported from MetaCyc.
"""

import importlib.resources
import sqlite3
from collections.abc import Iterable

import aiosql
import graph_tool as gt
import graph_tool.topology

import pan2met
import pan2met.sql.metabiantes

from . import KnowledgeBase


class MetabiantesKnowledgeBase(KnowledgeBase):
    def __init__(self, config):
        self.queries = aiosql.from_str(
            importlib.resources.read_text(pan2met.sql.metabiantes, "queries.sql"),
            "sqlite3",
        )
        self.connection = sqlite3.connect(
            config["reference"]["metabiantes"]["database"]
        )

    def _aiosql_to_list(self, iterator: Iterable[tuple]):
        if iterator is not None:
            return [item[0] for item in iterator]

    def monomers(self) -> list[str]:
        """
        List monomer polypeptides
        """
        return self._aiosql_to_list(self.queries.get_monomers(self.connection))

    def pathways(self) -> list[str]:
        """
        List the pathways referenced by the knowledge base
        """
        return self._aiosql_to_list(self.queries.get_pathways(self.connection))

    def pathway_name(self, pathway_id: str) -> str:
        """Get the name of a pathway"""
        return self.queries.get_pathway_name(self.connection, pathway_id=pathway_id)

    def reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List reactions of a pathway
        """
        return self._aiosql_to_list(
            self.queries.get_reactions_of_pathway(
                self.connection, pathway_id=pathway_id
            )
        )

    def reactions(self) -> list[str]:
        """
        List all reactions referenced in the knowledge base
        """
        return self._aiosql_to_list(self.queries.get_reactions(self.connection))

    def orphan_reactions(self) -> list[str]:
        """
        List orphan reactions
        """
        return self._aiosql_to_list(self.queries.get_orphan_reactions(self.connection))

    def spontaneous_reactions(self) -> list[str]:
        """
        List spontaneous reactions
        """
        return self._aiosql_to_list(
            self.queries.get_spontaneous_reactions(self.connection)
        )

    def non_spontaneous_reactions(self) -> list[str]:
        """
        List spontaneous reactions
        """
        return self._aiosql_to_list(
            self.queries.get_non_spontaneous_reactions(self.connection)
        )

    def non_spontaneous_reactions_of_pathway(self, pathway_id: str) -> list[str]:
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
    ) -> list[str]:
        """
        Non-orphan and non-spontaneous reactions of a pathway
        """
        return self._aiosql_to_list(
            self.queries.get_non_orphan_non_spontaneous_reactions_of_pathway(
                self.connection, pathway_id=pathway_id
            )
        )

    def enzymes_of_reaction(self, reaction_id: str) -> list[str]:
        """
        List enzymes catalyzing a reaction
        """
        return self._aiosql_to_list(
            self.queries.get_enzymes_of_reaction(
                self.connection, reaction_id=reaction_id
            )
        )

    def reactions_by_ec_number(self, ec_number: str) -> list[str]:
        """
        List reactions annotated with given EC-number
        """
        return self._aiosql_to_list(
            self.queries.get_reactions_by_ec_number(
                self.connection, ec_number=f"EC-{ec_number}"
            )
        )

    def pathway_taxonomic_range(self, pathway_id: str) -> int | None:
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
        res = self.queries.reaction_is_key(
            self.connection, pathway_id=pathway_id, reaction_id=reaction_id
        )
        return res is not None

    def key_reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List all key reactions of a pathway.
        """
        return self._aiosql_to_list(
            self.queries.get_key_reactions_of_pathway(
                self.connection, pathway_id=pathway_id
            )
        )

    def pathways_with_reaction(self, reaction_id) -> list[str]:
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

    def variants_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List the variants of a pathway.
        """
        return self._aiosql_to_list(
            self.queries.get_variants_of_pathway(self.connection, pathway_id=pathway_id)
        )

    def ontology_parent_class_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List the parent class of a pathway in the ontology of pathway tools
        """
        return self._aiosql_to_list(
            self.queries.get_ontology_parent_class_of_pathway(
                self.connection, pathway_id=pathway_id
            )
        )

    def pathway_reaction_order(self, pathway_id: str) -> list[tuple[str, str]]:
        """
        Return a list of pairs of reactions in the pathway, where the first reaction is a predecessor of the second reaction in the pathway.
        """
        return self.queries.get_pathway_reaction_order(
            self.connection, pathway_id=pathway_id
        )

    def species_evidence_of_pathway(self, pathway_id: str) -> list[int]:
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

    def reaction_graph_topological_order(self, pathway) -> list[str] | None:
        """
        Return the topological ordering of a pathway reaction graph.

        Assume that the pathway reaction graph is a directed acyclic graph.
        """
        reaction_order: list[tuple[str, str]] = self.pathway_reaction_order(pathway)
        graph = gt.Graph(directed=True)
        vmap = {}
        for reaction_1, reaction_2 in reaction_order:
            if reaction_1 not in vmap:
                vmap[reaction_1] = graph.add_vertex()
            if reaction_2 not in vmap:
                vmap[reaction_2] = graph.add_vertex()
            graph.add_edge(vmap[reaction_1], vmap[reaction_2])
        try:
            sort = gt.topology.topological_sort(graph)
        except ValueError:
            # topological sort fails if the graph is not a DAG.
            # in this case we return None
            return None
        reverse_vmap = {value: key for key, value in vmap.items()}
        sort_reactions = [reverse_vmap[vertex] for vertex in sort]
        return sort_reactions

    def complex(self) -> list[str]:
        """
        List the protein complex
        """
        return self._aiosql_to_list(self.queries.get_complex(self.connection))

    def components_of_complex(self, complex_id: str) -> list[str]:
        """
        List the components of a complex
        """
        return self._aiosql_to_list(
            self.queries.get_polypeptide_complex_components(
                self.connection, complex_id=complex_id
            )
        )
