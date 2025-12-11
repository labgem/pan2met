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

    def reaction_enzymes(self, reaction_id: str) -> list[str]:
        """
        List enzymes catalyzing a reaction
        """
        pass

    def reactions_by_ec_number(self, ec_number: str) -> list[str]:
        """
        List reactions annotated with given EC-number
        """
