import pythoncyc

from ...utils import logger
from ..knowledge_base import KnowledgeBase


class PythonCycKnowledgeBase(KnowledgeBase):
    """
    PythonCyc - based source of knowledge on pathways.
    """

    def __init__(self):
        try:
            self.pgdb = pythoncyc.select_organism("meta")
        except Exception as e:
            logger.critical(
                "Make sure you launch PathwayTools python API with "
                "`pathway-tools -lisp -python-local-only-non-strict` before running this command"
            )
            raise e

    def pathways(self) -> list[str]:
        """
        List the pathways referenced by the knowledge base
        """
        raise NotImplementedError()

    def reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List reactions of a pathway
        """
        raise NotImplementedError()

    def non_spontaneous_reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        Non-spontaneous reactions of a pathway
        """
        raise NotImplementedError()

    def non_orphan_non_spontaneous_reactions_of_pathway(
        self, pathway_id: str
    ) -> list[str]:
        """
        Non-orphan and non-spontaneous reactions of a pathway
        """
        raise NotImplementedError()

    def reaction_enzymes(self, reaction_id: str) -> list[str]:
        """
        List enzymes catalyzing a reaction
        """
        raise NotImplementedError()

    def reactions_by_ec_number(self, ec_number: str) -> list[str]:
        """
        List reactions annotated with given EC-number
        """
        raise NotImplementedError()

    def pathway_taxonomic_range(self, pathway_id: str) -> int:
        """
        Return a NCBI-Taxonomy Taxonomy Identifier number
        """
        raise NotImplementedError()

    def reaction_is_key(self, pathway_id: str, reaction_id: str) -> bool:
        """
        Return True if the reaction is a key reaction of the pathway
        """
        raise NotImplementedError()

    def key_reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List all key reactions of a pathway.
        """
        raise NotImplementedError()

    def pathways_with_reaction(self, reaction_id: str) -> list[str]:
        """
        List all pathways with the given reaction identifier.
        """
        raise NotImplementedError()

    def variants_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List the variants of a pathway.
        """
        raise NotImplementedError()
