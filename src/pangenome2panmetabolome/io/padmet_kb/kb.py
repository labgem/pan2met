"""
A driver wrapper on PADMet API to serve as a source of knowledge on metabolites

The PADMet file required for this tools must contain only non-orphan reactions.
It seems that PADMet does not store information on the spontaneity of a reaction.

To generate this PADMet file metacyc.padmet for MetaCyc,
first, export the MetaCyc database as flat files:
$ pathway-tools -lisp
EC(0): (select-organism :org-id 'meta)
EC(1): (create-flat-files-for-current-kb)

Then, run padmet's pgdb_to_padmet:
$ padmet pgdb_to_padmet --pgdb=~/.local/share/pathway-tools/aic-export/pgdbs/biocyc/metacyc/29.5/data --output=metacyc.padmet --no-orphan --extract-gene
"""

from queue import Queue

from padmet.classes import PadmetSpec


from ..knowledge_base import KnowledgeBase
from ...config import config


class PADMetKnowledgeBase(KnowledgeBase):
    def __init__(self):
        padmet_filename = config["reference"]["padmet_file"]
        self.padmet_object = PadmetSpec(padmet_filename)
        self.pathways_to_reactions_dict = self.padmet_object.getPathwaysReactions()

    def monomers(self) -> list[str]:
        """
        List monomer polypeptides
        """
        raise NotImplementedError("monomer listing not implemented for PADMet.")

    def pathways(self) -> list[str]:
        """
        List the pathways referenced by the knowledge base
        """
        return self.pathways_to_reactions_dict.keys()

    def reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List reactions of a pathway
        """
        return self.pathways_to_reactions_dict[pathway_id]

    def reactions(self) -> list[str]:
        """
        List all reactions referenced in the knowledge base
        """
        return self.padmet_object.getReactions()

    def non_spontaneous_reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        Non-spontaneous reactions of a pathway
        """
        return NotImplementedError(
            "PADMet does not keep the information on reaction spontaneity, a priori."
        )

    def non_orphan_non_spontaneous_reactions_of_pathway(
        self, pathway_id: str
    ) -> list[str]:
        """
        Non-orphan and non-spontaneous reactions of a pathway
        """
        # TODO raise NotImplementedError("Listing non-orphan non-spontaneous is not implemented ")
        reactions = [
            reaction
            for reaction in self.pathways_to_reactions_dict[pathway_id]
            if (
                "SPONTANEOUS" not in self.padmet_object.dicOfNode[reaction].misc
                or self.padmet_object.dicOfNode[reaction].misc["SPONTANEOUS"][0] == "T"
            )
        ]
        return reactions

    def enzymes_of_reaction(self, reaction_id: str) -> list[str]:
        """
        List enzymes catalyzing a reaction
        """
        raise NotImplementedError(
            "no enzymes_of_reaction() method implemented for PADMet format yet"
        )

    def ec_number_of_reaction(self, reaction_id: str) -> str:
        """
        Get the EC-number of a reaction
        """
        return self.dicOfNode[reaction_id].misc["EC-NUMBER"]

    def reactions_by_ec_number(self, ec_number: str) -> list[str]:
        """
        List reactions annotated with given EC-number

        """
        reactions = [
            reaction.id
            for reaction in self.padmet_object.getReactions(get_id=False)
            if "EC-NUMBER" in reaction.misc and ec_number in reaction.misc["EC-NUMBER"]
        ]
        return reactions

    def pathway_taxonomic_range(self, pathway_id: str) -> int:
        """
        Return a NCBI-Taxonomy taxonomy identifier number
        """
        if "TAXONOMIC-RANGE" in self.padmet_object.dicOfNode[pathway_id].misc:
            tax_id = self.padmet_object.dicOfNode[pathway_id].misc["TAXONOMIC-RANGE"][0]
            tax_id_number = int(tax_id.split("-")[-1])
            return tax_id_number
        else:
            return None

    def key_reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List all key reactions of a pathway.
        """
        if "KEY-REACTIONS" in self.padmet_object.dicOfNode[pathway_id].misc:
            return self.padmet_object.dicOfNode[pathway_id].misc["KEY-REACTIONS"]
        else:
            return []

    def pathways_with_reaction(self, reaction_id: str) -> list[str]:
        """
        List all pathways with the given reaction identifier.
        """
        return [
            pathway
            for pathway, reactions in self.pathways_to_reactions_dict.items()
            if reaction_id in reactions
        ]

    def count_pathways_with_reaction(self, reaction_id: str) -> int:
        """
        Count the pathways having the given reaction identifier.
        """
        return len(self.pathways_with_reaction(reaction_id))

    def variants_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List the variants of a pathway.
        """
        return []
        # raise NotImplementedError("For PADMet: variants of pathway is not implemented.")

    def ontology_parent_class_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List the parent class of a pathway in the ontology of pathway tools
        """
        # Trace back the ontology tree
        parent_classes: list[str] = []
        current_class = pathway_id
        queue = Queue()
        queue.put(current_class)
        visited = set()
        root = "FRAMES"
        while not queue.empty():
            current_class = queue.get()
            parents = [
                class_relation.id_out
                for class_relation in self.padmet_object.dicOfRelationIn[current_class]
                if class_relation.type == "is_a_class"
            ]
            parent_classes = parent_classes + parents
            for parent in parents:
                if parent not in visited and parent != root:
                    queue.put(parent)
                    visited.add(parent)
        return parent_classes
