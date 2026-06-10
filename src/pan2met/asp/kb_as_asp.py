"""
Output a metabolic knowledge base in ASP format for use with clingo.

The Answer Set Programming atoms we use are:

- pathway(Pathway)
- is_in_pathway(Reaction, Pathway).
- is_a(Pathway, PathwayClass). / is_a(PathwayClass, ParentPathwayClass).
- in_taxonomic_range(Pathway). (depending on the taxonomic range of the pathway and the taxonomy of the target organism.)
- reactome(Reaction). (depending on the reactome of the target organism, inferred from its genome, by homology.)
- orphan(Reaction). (when the reaction has no known enzyme in the reference knowledge base.)
- spontaneous(Reaction). (when the reaction is spontaneous.)
- pathway_key_reaction(Pathway, Reaction). (when the knowledge base specify that the reaction is key for the pathway.)
- pathway_reaction_order(Pathway, Predecessor, Successor). (indicates the topological ordering of the reactions in the pathway.)
"""

import argparse


from ..io.knowledge_base import KnowledgeBase
from ..io.knowledge_base import select_kb
from ..config import default_config, override_config


def pathway_asp_atom(pathway: str) -> str:
    return f'pathway("{pathway}").'


def is_in_pathway_asp_atom(reaction: str, pathway: str) -> str:
    return f'is_in_pathway("{reaction}", "{pathway}").'


def orphan_asp_atom(reaction: str) -> str:
    return f'orphan("{reaction}").'


def spontaneous_asp_atom(reaction: str) -> str:
    return f'spontaneous("{reaction}").'


def in_taxonomic_range_asp_atom(pathway: str) -> str:
    return f'in_taxonomic_range("{pathway}").'


def reactome_asp_atom(reaction: str) -> str:
    return f'reactome("{reaction}").'


def pathway_reaction_order(pathway: str, predecessor: str, successor: str) -> str:
    return f'pathway_reaction_order("{pathway}", "{predecessor}", "{successor}").'


def pathway_key_reaction_asp_atom(pathway: str, reaction: str) -> str:
    return f'pathway_key_reaction("{pathway}", "{reaction}").'


def pathway_class_asp_atom(pathway: str, pathway_class: str) -> str:
    return f'is_a("{pathway}", "{pathway_class}").'


def write_kb_as_asp(kb: KnowledgeBase, output_file: str):
    with open(output_file, "w") as f:
        # Write pathway/1 atoms
        pathways = kb.pathways()
        for pathway in pathways:
            f.write(pathway_asp_atom(pathway) + "\n")

        # Write is_in_pathway/2 atoms
        for pathway in pathways:
            reactions = kb.reactions_of_pathway(pathway)
            for reaction in reactions:
                f.write(is_in_pathway_asp_atom(reaction, pathway) + "\n")

        # Write is_a/2 atoms for pathway ontology
        for pathway in pathways:
            for pathway_class in kb.ontology_parent_class_of_pathway(pathway):
                f.write(pathway_class_asp_atom(pathway, pathway_class) + "\n")

        # Write orphan/1 atoms
        for reaction in kb.orphan_reactions():
            f.write(orphan_asp_atom(reaction) + "\n")

        # Write spontaneous/1 atoms
        for reaction in kb.spontaneous_reactions():
            f.write(spontaneous_asp_atom(reaction) + "\n")

        # Write pathway reaction order atoms
        for pathway in pathways:
            reaction_order = kb.pathway_reaction_order(pathway)
            for predecessor, successor in reaction_order:
                f.write(pathway_reaction_order(pathway, predecessor, successor) + "\n")

        # Write pathway key reactions atoms
        for pathway in pathways:
            for reaction in kb.key_reactions_of_pathway(pathway):
                f.write(pathway_key_reaction_asp_atom(pathway, reaction) + "\n")


def target_organism_reactions_as_asp(reactome: set[str]) -> str:
    return "\n".join(reactome_asp_atom(reaction) for reaction in reactome)


def pathways_in_taxonomic_range_as_asp(
    pathway_in_taxonomic_range: dict[str, bool],
) -> str:
    return "\n".join(
        in_taxonomic_range_asp_atom(pathway)
        for pathway in pathway_in_taxonomic_range
        if pathway_in_taxonomic_range[pathway]
    )


def main():
    parser = argparse.ArgumentParser(
        description="Output a metabolic knowledge base to ASP atoms"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="The output file to write the ASP atoms to",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="Config override",
        required=False,
        default=None
    )
    args = parser.parse_args()
    if args.config:
        config = override_config()
    else:
        config = default_config
    kb = select_kb(config)
    write_kb_as_asp(kb, args.output)


if __name__ == "__main__":
    main()
