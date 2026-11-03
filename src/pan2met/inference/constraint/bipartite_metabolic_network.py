"""
Clingo Answer-Set Programming constraint-based metabolism inference from reactome
"""

import argparse
import importlib
import logging

import pan2met

from ...config import default_config, override_config
from ...io.knowledge_base import KnowledgeBase, select_kb
from ...utils import read_list, write_output


class ASPMetabolicBipartiteNetworkInference:
    """
    Constraint based metabolism inference
    """

    def __init__(self, template: KnowledgeBase):
        self.kb = template

    def dump_metabolic_knowledge_to_asp(self) -> str:
        """
        Dump a Clingo ASP file with the following atoms, with information taken from the knowledge based:
        - reaction(Reaction).
        - spontaneous(Reaction).
        - reaction_reactant(Reaction, Metabolite).
        - reaction_product(Reaction, Metabolite).
        - metabolite(Metabolite).
        """
        spontaneous_reaction_asp = self.dump_asp_atoms_semicolon(
            self.kb.spontaneous_reactions(), "spontaneous"
        )
        reactions_asp = self.dump_asp_atoms_semicolon(self.kb.reactions(), "reaction")

        reaction_sides_asp = ""
        metabolites: set[str] = set()
        for reaction in self.kb.reactions():
            products, reactants = self.kb.sides_of_reaction(reaction)
            for product in products:
                reaction_sides_asp += "\n"
                reaction_sides_asp += f'reaction_product("{reaction}", "{product}").'
            for reactant in reactants:
                reaction_sides_asp += "\n"
                reaction_sides_asp += f'reaction_reactant("{reaction}", "{reactant}").'
            for sets in [products, reactants]:
                for item in sets:
                    metabolites.add(item)
        metabolite_asp = self.dump_asp_atoms_semicolon(metabolites, "metabolite")
        return f"{spontaneous_reaction_asp}\n{reactions_asp}\n{reaction_sides_asp}\n{metabolite_asp}"

    def dump_asp_atoms_semicolon(self, atoms: list[str], type: str) -> str:
        return f"{type}(" + ";".join(f'"{atom}"' for atom in atoms) + ")."

    def inferred_reactions(self, catalysis: list[str]) -> set[str]:
        """
        List inferred reactions of the built metabolic graph using Answer Set Programming/ASP inference rules.
        """

        # We keep clyngor as an optional dependency
        import clyngor

        with importlib.resources.path(
            pan2met, "asp/rules/denovo_bipartite_metabolic_graph.lp"
        ) as denovo_bipartite_metabolic_graph_rule_path:
            inline_asp = self.dump_metabolic_knowledge_to_asp()
            inline_asp += self.dump_asp_atoms_semicolon(catalysis, "catalysis")

            # with open('/tmp/atoms.lp', "w") as f:
            #     f.write(inline_asp)

            answers = clyngor.solve(
                [denovo_bipartite_metabolic_graph_rule_path],
                inline=inline_asp,
                nb_model=1,
            )  # FIXME: we might be interested in more than one model.
            answer = next(answers)
            inferred: set[str] = set()
            for predicate, value in answer:
                if predicate == "selected_reaction":
                    identifier = value[0]
                    identifier = identifier.replace('"', "")
                    inferred.add(identifier)
            return inferred


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="config override")
    parser.add_argument(
        "--catalysis", help="list of catalyzed reaction identifiers", required=True
    )
    parser.add_argument("-o", "--output", help="output list of pathways", required=True)
    return parser.parse_args()


def main():
    logging.basicConfig(level=logging.DEBUG)
    args = parse_arguments()

    if args.config:
        config = override_config(args.config, default_config)
    else:
        config = default_config
    if args.catalysis:
        catalysis = set(read_list(args.catalysis))
        kb = select_kb(config)
        inference = ASPMetabolicBipartiteNetworkInference(kb)

        inferred_reactions: set[str] = inference.inferred_reactions(catalysis)
        write_output(args.output, inferred_reactions)


if __name__ == "__main__":
    main()
