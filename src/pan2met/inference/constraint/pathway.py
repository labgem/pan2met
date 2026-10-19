"""
Identify candidate de novo metabolic pathways.

What do we call "metabolic pathway" in this context?
We call a candidate metabolic pathway
a connected subgraph in the reaction graph,
that has at least one catalysis realizes by genes that are adjacent in the pangenome graph.

Currently, we use only the reaction graph edges identified in the metabolic pathway knowledgebase, so they reflect directly the already identified relevant reaction succession.
"""

import argparse
import csv
import importlib.resources

import clyngor
import graph_tool as gt

import pan2met
from pan2met.config import default_config, override_config
from pan2met.io.knowledge_base import KnowledgeBase, select_kb
from pan2met.utils import read_mapping


def dump_reaction_graph_to_asp(kb: KnowledgeBase) -> str:
    """
    Dump a Clingo ASP string with the following atoms, from information taken from the knowledge base:
    - reaction/1.
    - reaction_graph_edge/1.

    :return: the asp atoms as a string
    """
    edges: set[tuple[str, str]] = set()
    reactions: set[str] = set()
    for pathway in kb.pathways():
        for reaction_from, reaction_to in kb.pathway_reaction_order(pathway):
            edge = (reaction_from, reaction_to)
            edges.add(edge)
            reactions.add(reaction_from)
            reactions.add(reaction_to)
    reaction_asp = "\n".join(f'reaction("{reaction}").' for reaction in reactions)
    edge_asp = "\n".join(f'reaction_graph_edge("{u}", "{v}").' for u, v in edges)
    "\n".join(f'spontaneous("{reaction}").' for reaction in kb.spontaneous_reactions())
    return f"{reaction_asp}\n\n{edge_asp}"


def dump_pangenome_graph_to_asp(
    pangenome_graph: gt.Graph, families_to_reactions: dict[str, set[str]]
) -> str:
    """
    Dump a Clingo ASP string with the following atoms from a (pan)genome graph with the following vertex property: nid, a gene family identifiers.

    - catalysis(Reaction, GeneFamilyVertex).
    - pangenome_graph_edge(GeneFamilyVertex1, GeneFamilyVertex2).
    """
    pangenome_graph_edge_asp = ""
    catalysis_asp = ""

    for edge in pangenome_graph.edges():
        u = edge.source()
        v = edge.target()
        pangenome_graph_edge_asp += "\n"
        pangenome_graph_edge_asp += f'pangenome_graph_edge("{u}", "{v}").'

        u_family = pangenome_graph.vp["nid"][u]
        v_family = pangenome_graph.vp["nid"][v]

        if u_family in families_to_reactions:
            for reaction in families_to_reactions[u_family]:
                catalysis_asp += "\n"
                catalysis_asp += f'catalysis("{reaction}", "{u_family}").'
        if v_family in families_to_reactions:
            for reaction in families_to_reactions[v_family]:
                catalysis_asp += "\n"
                catalysis_asp += f'catalysis("{reaction}", "{v_family}").'
    return f"{pangenome_graph_edge_asp}\n\n{catalysis_asp}"


def denovo_pathway_identification(
    pangenome_graph: gt.Graph,
    families_to_reactions: dict[str, set[str]],
    kb: KnowledgeBase,
    limit=100,
) -> list[set[str]]:
    """
    Use clingo to identify candidate metabolic pathway from the pangenome graph.
    """
    pangenome_catalysis_asp = dump_pangenome_graph_to_asp(
        pangenome_graph, families_to_reactions
    )
    reaction_asp = dump_reaction_graph_to_asp(kb)

    with importlib.resources.path(
        pan2met, "asp/rules/pangenome_pathway_extraction.lp"
    ) as pangenome_pathway_extraction_path:
        inline_asp = (
            pangenome_catalysis_asp + "\n\n" + reaction_asp + f"#const n={limit}."
        )
        answers = clyngor.solve(
            [pangenome_pathway_extraction_path], inline=inline_asp, nb_model=1
        )
        answer = next(answers)
        pathway_to_reactions: dict[int, set[str]] = {}
        for predicate, value in answer:
            if predicate == "reaction_is_in_pathway":
                reaction = value[1]
                reaction = reaction.replace('"', "")
                pathway = value[0]
                pathway = int(pathway)
                if pathway in pathway_to_reactions:
                    pathway_to_reactions[pathway].add(reaction)
                else:
                    pathway_to_reactions[pathway] = {reaction}
    return list(pathway_to_reactions.values())


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", help="path to a configuration yaml file", required=False
    )
    parser.add_argument(
        "--gene-reaction", help="path to gene-family - reaction mapping", required=True
    )
    parser.add_argument(
        "-g",
        "--pangenome-graph",
        help="path to the (pan)-genome graph in graph-tool binary gt format with gene family (nid)",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="output filename, a two columns TSV file with column 1 a pathway numeric identifier and column two a reaction of belonging to this pathway, according to the de novo pathway construction.",
        required=True,
    )
    parser.add_argument(
        "-l",
        "--limit",
        help="limit the number of pathway to search, a lower value will often lead to faster solving, but could miss more pathway candidates (Default 100).",
        type=int,
        default=100,
    )
    return parser


def main():
    parser = argument_parser()
    args = parser.parse_args()

    if args.config:
        config = override_config(args.config, default_config)
    else:
        config = default_config

    pangenome_graph: gt.Graph = gt.load_graph(args.pangenome_graph)
    families_to_reaction: dict[str, set[str]] = read_mapping(args.gene_reaction)

    kb: KnowledgeBase = select_kb(config)

    pathway_candidates = denovo_pathway_identification(
        pangenome_graph, families_to_reaction, kb, limit=args.limit
    )

    with open(args.output, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        for pathway_number, pathway_candidate_reaction_set in enumerate(
            pathway_candidates
        ):
            for reaction in pathway_candidate_reaction_set:
                writer.writerow([pathway_number, reaction])


if __name__ == "__main__":
    main()
