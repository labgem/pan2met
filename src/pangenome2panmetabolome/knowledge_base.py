"""
Define an asbtract interface to the knowledge base.
Knowledge base backends should implement all the methods of the KnowledgeBase abstract class.
"""

from abc import ABC


class KnowledgeBase(ABC):
    """
    Define the required information a template pathway knowledge base
    is expected to contain to be able to use the inference rules.
    """

    def __init__(self):
        pass

    def pathways(self) -> list[str]:
        """
        List the pathways referenced by the knowledge base
        """
        pass

    def reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List reactions of a pathway
        """
        pass

    def non_spontaneous_reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        Non-spontaneous reactions of a pathway
        """
        pass

    def non_orphan_non_spontaneous_reactions_of_pathway(
        self, pathway_id: str
    ) -> list[str]:
        """
        Non-orphan and non-spontaneous reactions of a pathway
        """
        pass

    def reaction_enzymes(self, reaction_id: str) -> list[str]:
        """
        List enzymes catalyzing a reaction
        """
        pass

    def reactions_by_ec_number(self, ec_number: str) -> list[str]:
        """
        List reactions annotated with given EC-number
        """
        pass

    def pathway_taxonomic_range(self, pathway_id: str) -> int:
        """
        Return a NCBI-Taxonomy Taxonomy Identifier number
        """
        pass

    def reaction_is_key(self, pathway_id: str, reaction_id: str) -> bool:
        """
        Return True if the reaction is a key reaction of the pathway
        """
        pass

    def key_reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List all key reactions of a pathway.
        """

    def pathways_with_reaction(self, reaction_id: str) -> list[str]:
        """
        List all pathways with the given reaction identifier.
        """
        pass

    def variants_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List the variants of a pathway.
        """
        pass
