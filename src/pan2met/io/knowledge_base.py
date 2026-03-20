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

    def monomers(self) -> list[str]:
        """
        List monomer polypeptides
        """
        raise NotImplementedError()

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

    def spontaneous_reactions(self):
        """
        List spontaneous reactions
        """
        raise NotImplementedError()

    def orphan_reactions(self):
        """
        List orphan reactions
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


def select_kb(choice: str) -> KnowledgeBase:
    match choice:
        case "metabiantes":
            from .metabiantes_kb.kb import MetabiantesKnowledgeBase

            return MetabiantesKnowledgeBase()
        case "pythoncyc":
            from .pythoncyc_kb.kb import PythonCycKnowledgeBase

            return PythonCycKnowledgeBase()
        case "padmet":
            from .padmet_kb.kb import PADMetKnowledgeBase

            return PADMetKnowledgeBase()
        case _:
            raise ValueError(
                f"Cannot load kb for choice {choice}. Not in {set(['metabiantes', 'pythoncyc', 'padmet'])}"
            )
