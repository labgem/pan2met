"""
A SQLite backend using 'metabiantes' model for data imported from MetaCyc.
"""

from typing import Iterable
import importlib.resources

import aiosql
import psycopg

import pangenome2panmetabolome

from ...knowledge_base import KnowledgeBase
from ...config import config


class MetabiantesKnowledgeBase(KnowledgeBase):
    def __init__(self):
        self.queries = aiosql.from_str(
            importlib.resources.read_text(
                pangenome2panmetabolome, "sql/metabiantes/queries.sql"
            ),
            "psycopg2",
        )
        self.connection = psycopg.connect(
            f"dbname={config['reference']['postgresql_database']}"
        )

    def aiosql_to_list(self, iterator: Iterable[tuple]):
        if iterator is not None:
            return [item[0] for item in iterator]

    def pathways(self) -> list[str]:
        """
        List the pathways referenced by the knowledge base
        """
        return self.aiosql_to_list(self.queries.get_pathways(self.connection))

    def reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List reactions of a pathway
        """
        return self.aiosql_to_list(
            self.queries.get_reactions_of_pathway(
                self.connection, pathway_id=pathway_id
            )
        )

    def non_spontaneous_reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        Non-spontaneous reactions of a pathway
        """
        return self.aiosql_to_list(
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
        return self.aiosql_to_list(
            self.queries.get_non_orphan_non_spontaneous_reactions_of_pathway(
                self.connection, pathway_id=pathway_id
            )
        )

    def enzymes_of_reaction(self, reaction_id: str) -> list[str]:
        """
        List enzymes catalyzing a reaction
        """
        return self.aiosql_to_list(
            self.queries.get_enzymes_of_reaction(
                self.connection, reaction_id=reaction_id
            )
        )

    def reactions_by_ec_number(self, ec_number: str) -> list[str]:
        """
        List reactions annotated with given EC-number
        """
        return self.aiosql_to_list(
            self.queries.get_reactions_by_ec_number(
                self.connection, ec_number=ec_number
            )
        )

    def pathway_taxonomic_range(self, pathway_id: str) -> int:
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

    def key_reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List all key reactions of a pathway.
        """
        return self.aiosql_to_list(
            self.queries.get_key_reactions_of_pathway(
                self.connection, pathway_id=pathway_id
            )
        )

    def pathways_with_reaction(self, reaction_id) -> list[str]:
        """
        List all pathways with the given reaction identifier.
        """
        return self.aiosql_to_list(
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
        return self.aiosql_to_list(
            self.queries.get_variants_of_pathway(self.connection, pathway_id=pathway_id)
        )
