"""
Define an asbtract interface to the knowledge base.
Knowledge base backends should implement all the methods of the KnowledgeBase abstract class.
"""

from typing import List, Optional
from abc import ABC

class KnowledgeBase(ABC):
    """
    An abstract class to implement an interface to a metabolic pathway knowledgebase.
    """

    def __init__(self):
        pass

    def monomers(self) -> List[str]:
        """
        List monomer polypeptides
        """
        raise NotImplementedError()

    def pathways(self) -> List[str]:
        """
        List the pathways referenced by the knowledge base
        """
        raise NotImplementedError()

    def reactions_of_pathway(self, pathway_id: str) -> List[str]:
        """
        List reactions of a pathway
        """
        raise NotImplementedError()

    def name_of_pathway(self, pathway_id: str) -> str:
        """
        Get the name of a pathway
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

    def non_spontaneous_reactions_of_pathway(self, pathway_id: str) -> List[str]:
        """
        Non-spontaneous reactions of a pathway
        """
        raise NotImplementedError()

    def non_orphan_non_spontaneous_reactions_of_pathway(
        self, pathway_id: str
    ) -> List[str]:
        """
        Non-orphan and non-spontaneous reactions of a pathway
        """
        raise NotImplementedError()

    def reaction_enzymes(self, reaction_id: str) -> List[str]:
        """
        List enzymes catalyzing a reaction
        """
        raise NotImplementedError()

    def reactions_by_ec_number(self, ec_number: str) -> List[str]:
        """
        List reactions annotated with given EC-number
        """
        raise NotImplementedError()

    def pathway_taxonomic_range(self, pathway_id: str) -> Optional[int]:
        """
        Return a NCBI-Taxonomy Taxonomy Identifier number
        """
        raise NotImplementedError()

    def reaction_is_key(self, pathway_id: str, reaction_id: str) -> bool:
        """
        Return True if the reaction is a key reaction of the pathway
        """
        raise NotImplementedError()

    def key_reactions_of_pathway(self, pathway_id: str) -> List[str]:
        """
        List all key reactions of a pathway.
        """
        raise NotImplementedError()

    def pathway_reaction_order(self, pathway_id: str) -> List[str]:
        """
        Return an ordered list of reactions of a pathway.
        """
        raise NotImplementedError()

    def pathways_with_reaction(self, reaction_id: str) -> List[str]:
        """
        List all pathways with the given reaction identifier.
        """
        raise NotImplementedError()

    def variants_of_pathway(self, pathway_id: str) -> List[str]:
        """
        List the variants of a pathway.
        """
        raise NotImplementedError()

    def species_evidence_of_pathway(self, pathway_id: str) -> List[int]:
        """
        List species NCBI Taxonomy identifiers where pathway presence evidence were found.
        """
        raise NotImplementedError()

    def ontology_parent_class_of_pathway(self, pathway_id: str) -> List[str]:
        """
        List the pathway ontology parents of a pathway
        """
        raise NotImplementedError()


def select_kb(config) -> KnowledgeBase:
    match config["reference"]["source"]:
        case "metabiantes":
            from .metabiantes_kb.kb import MetabiantesKnowledgeBase

            return MetabiantesKnowledgeBase(config)
        case "pythoncyc":
            from .pythoncyc_kb.kb import PythonCycKnowledgeBase

            return PythonCycKnowledgeBase()
        case "padmet":
            from .padmet_kb.kb import PADMetKnowledgeBase

            return PADMetKnowledgeBase(config)
        case _:
            raise ValueError(
                f"Cannot load kb for choice {config["reference"]["source"]}. Not in {set(['metabiantes', 'pythoncyc', 'padmet'])}"
            )
