"""
Define an asbtract interface to the knowledge base.
Knowledge base backends should implement all the methods of the KnowledgeBase abstract class.
"""

from abc import ABC


def unsupported(function):
    """
    A decorator to make the function raise NotImplementedError and add the attribute __unsupported.

    __unsupported is True when the function would raise NotImplementedError
    """

    def magic(self):
        raise NotImplementedError(function)
        # We do not call the function it self

    setattr(magic, "__unsupported", True)
    return magic


class KnowledgeBase(ABC):
    """
    An abstract class to implement an interface to a metabolic pathway knowledgebase.
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

    def name_of_pathway(self, pathway_id: str) -> str:
        """
        Get the name of a pathway
        """
        raise NotImplementedError()

    def spontaneous_reactions(self) -> list[str]:
        """
        List spontaneous reactions
        """
        raise NotImplementedError()

    def orphan_reactions(self) -> list[str]:
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

    def enzymes_of_reaction(self, reaction_id: str) -> list[str]:
        """
        List enzymes catalyzing a reaction
        """
        raise NotImplementedError()

    def reactions_by_ec_number(self, ec_number: str) -> list[str]:
        """
        List reactions annotated with given EC-number
        """
        raise NotImplementedError()

    def pathway_taxonomic_range(self, pathway_id: str) -> int | None:
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

    def pathway_reaction_order(self, pathway_id: str) -> list[tuple[str, str]]:
        """
        Return a list of (predecessor, successor) reactions of a pathway.
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

    def species_evidence_of_pathway(self, pathway_id: str) -> list[int]:
        """
        List species NCBI Taxonomy identifiers where pathway presence evidence were found.
        """
        raise NotImplementedError()

    def ontology_parent_class_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List the pathway ontology parents of a pathway
        """
        raise NotImplementedError()

    def complex(self) -> list[str]:
        """
        List the protein complex
        """
        raise NotImplementedError()

    def components_of_complex(self, complex_id: str) -> list[str]:
        """
        List the components of a complex
        """
        raise NotImplementedError()


def select_kb(config) -> KnowledgeBase:
    match config["reference"]["knowledge_base"]:
        case "metabiantes":
            from .metabiantes import MetabiantesKnowledgeBase

            return MetabiantesKnowledgeBase(config)
        case "padmet":
            from .padmet import PADMetKnowledgeBase

            return PADMetKnowledgeBase(config)
        case "biopax_folder":
            from .biopax_folder import BioPAXFolderKnowledgeBase

            return BioPAXFolderKnowledgeBase(config)
        case _:
            raise ValueError(
                f"Cannot load kb for choice {config['reference']['source']}. Not in { {'metabiantes', 'pythoncyc', 'padmet'} }"
            )
