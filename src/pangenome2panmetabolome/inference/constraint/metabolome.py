"""
Clingo Answer-Set Programming constraint-based metabolome inference from reactome
"""

import logging
import argparse
import importlib

import clyngor

import pangenome2panmetabolome
from ...utils import read_list, write_output
from ...config import config
from ...io.knowledge_base import KnowledgeBase, select_kb


logger = logging.getLogger("pangenome2panmetabolome:inference:constraint")


class ASPPathwayInference:
    """
    Constraint based metabolome inference
    """

    def __init__(self, template: KnowledgeBase, reactome: set[str]):
        self.template = template
        self.reactome = reactome

    def dump_pathway_knowledge_to_asp(self):
        """
        Dump a Clingo ASP file with the following atoms, with information taken from the knowledge based:
        - pathway/1.
        - reaction/1.
        - pathway_reaction/2. (pathway(Pathway, Reaction)).
        """
        pathways: list[str] = self.template.pathways()
        pathway_list_asp = self.dump_pathway_list_to_asp(pathways)

        reactions: list[str] = self.template.reactions()
        reaction_list_asp = self.dump_reaction_list_to_asp(reactions)

        pathway_reactions: list[str] = [
            (pathway, self.template.reactions_of_pathway(pathway))
            for pathway in pathways
        ]

        pathway_reaction_list_asp = self.dump_pathway_reaction_list_to_asp(
            pathway_reactions
        )

        return "\n".join(
            [pathway_list_asp, reaction_list_asp, pathway_reaction_list_asp]
        )

    def dump_pathway_list_to_asp(self, pathways: list[str]) -> str:
        return "pathway(" + ";".join(f'"{pathway}"' for pathway in pathways) + ")."

    def dump_reaction_list_to_asp(self, reactions: list[str]) -> str:
        return "reaction(" + ";".join(f'"{reaction}"' for reaction in reactions) + ")."

    def dump_monomer_list_to_asp(self, monomers: list[str]) -> str:
        return "\n".join(f'monomer("{monomer}")' for monomer in monomers)

    def dump_pathway_reaction_list_to_asp(
        self, pathway_reactions: list[tuple[str, list[str]]]
    ) -> str:
        return "\n".join(
            f'pathway_reaction("{pathway}", "{reaction}").'
            for pathway, reactions in pathway_reactions
            for reaction in reactions
        )

    def reactome_to_asp(self) -> str:
        return "\n".join(f'reactome("{reaction}").' for reaction in self.reactome)

    def inferred_pathways(self, kb_asp_path: str) -> set[str]:
        """
        List all inferred pathways using Answer Set Programming/ASP inference rules.

        Arguments
        ---------

            kb_asp_path: str
                filename of the ASP representation of the knowledge base

        Returns
        -------
            the set of inferred metabolic pathways
        """

        with importlib.resources.path(
            pangenome2panmetabolome, "../asp/minimal_covering_pathway.lp"
        ) as path:
            minimal_covering_pathway_rule_path = path

        inline_asp = self.reactome_to_asp()

        answers = clyngor.solve(
            [minimal_covering_pathway_rule_path, kb_asp_path], inline=inline_asp
        )
        logger.info("Done clingo prediction")
        answer = next(answers)
        inferred: set[str] = set()
        for predicate, value in answer:
            if predicate == "infer_present_pathway":
                identifier = value[0]
                identifier = identifier.replace('"', "")
                inferred.add(identifier)

        return inferred


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reactome", help="List of reaction identifiers")
    parser.add_argument(
        "--dump",
        action=argparse.BooleanOptionalAction,
        help="Whether to refresh the knowledge base dump as ASP.",
    )
    parser.add_argument(
        "--kb-asp",
        help="Path to the Knowledge base as ASP, written when --dump, only used when --no-dump",
        required=True,
    )
    parser.add_argument("-o", "--output", help="Output list of pathways")
    parser.add_argument("--reactome-asp", help="Export reactome as ASP atoms")
    return parser.parse_args()


def main():
    logging.basicConfig(level=logging.DEBUG)
    args = parse_arguments()
    if args.dump:
        kb_asp_dump = ASPPathwayInference.dump_pathway_knowledge_to_asp()
        with open(args.kb_asp, "w") as kb_asp_file:
            kb_asp_file.write(kb_asp_dump)
    else:
        with open(args.kb_asp, "r") as kb_asp_file:
            assert len(kb_asp_file.read()) > 1, (
                "Please provide a non-empty knowledge base ASP input, or use --dump option."
            )
    if args.reactome:
        reactome = set(read_list(args.reactome))
        kb = select_kb(config["reference"]["source"])
        inference = ASPPathwayInference(kb, reactome)
        # Write reactome as ASP atoms if --reactome-asp <path> is set
        if args.reactome_asp is not None:
            with open(args.reactome_asp, "w") as reactome_asp_file:
                reactome_asp_file.write(inference.reactome_to_asp())
        # Infer the metabolome and write to output if -o/--output <path> is set
        if args.output is not None:
            inferred_pathways: list[str] = inference.inferred_pathways(args.kb_asp)
            write_output(args.output, "\n".join(inferred_pathways))


if __name__ == "__main__":
    main()
