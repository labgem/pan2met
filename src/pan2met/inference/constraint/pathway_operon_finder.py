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
import logging
from collections.abc import Iterable

import clyngor
import graph_tool as gt
import graph_tool.topology
import numpy as np

import pan2met
from pan2met.config import default_config, override_config
from pan2met.io.knowledge_base import KnowledgeBase, select_kb
from pan2met.utils import logger, read_mapping

SAVE_ASP: str | None = "/tmp/asp_pathway_instance.lp"


logging.basicConfig(level=logging.DEBUG)
logger.setLevel(logging.DEBUG)


def dump_reaction_graph_to_asp(kb: KnowledgeBase) -> str:
    """
    Dump a Clingo ASP string with the following atoms, from information taken from the knowledge base:
    - reaction/1.
    - spontaneous/1.
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
    spontaneous_asp = "\n".join(
        f'spontaneous("{reaction}").' for reaction in kb.spontaneous_reactions()
    )
    reaction_graph = reaction_edges_to_reaction_graph(edges)
    connected_components_asp = reaction_connected_components_asp(reaction_graph)
    return f"{reaction_asp}\n\n{spontaneous_asp}\n\n{connected_components_asp}"


def reaction_edges_to_reaction_graph(
    reaction_graph_edges: Iterable[tuple[str, str]],
) -> gt.Graph:
    """
    Build a graph-tool directed reaction graph from a list of edges.

    :param reaction_graph_edges: the edges
    :return: a directed graph
    """
    # Build the graph_tool graph
    reaction_graph = gt.Graph(directed=True)
    reaction_to_vertex: dict[str, gt.Vertex] = {}
    vertex_to_reaction: gt.VertexPropertyMap = reaction_graph.new_vertex_property(
        "string"
    )
    for reaction_from, reaction_to in reaction_graph_edges:
        if reaction_from not in reaction_to_vertex:
            u = reaction_graph.add_vertex()
            reaction_to_vertex[reaction_from] = u
            vertex_to_reaction[u] = reaction_from
        else:
            u = reaction_to_vertex[reaction_from]
        if reaction_to not in reaction_to_vertex:
            v = reaction_graph.add_vertex()
            reaction_to_vertex[reaction_to] = v
            vertex_to_reaction[v] = reaction_to
        else:
            v = reaction_to_vertex[reaction_to]
        reaction_graph.add_edge(u, v)
    reaction_graph.vp["reaction"] = vertex_to_reaction
    return reaction_graph


def reaction_connected_components(reaction_graph: gt.Graph) -> gt.VertexPropertyMap:
    """
    Get the reaction graph strongly connected components as sets of reactions
    """
    connected_components_vertex_map, _hist = gt.topology.label_components(
        reaction_graph, directed=False
    )
    return connected_components_vertex_map


def reaction_connected_components_asp(reaction_graph: gt.Graph) -> str:
    """
    Generate a set of Answer Set Programmings atoms reaction_graph_strongly_connected_components/2:
    reaction_graph_connected_components(Component, Reaction). for every reaction Reaction in strongly connected component Component.
    """
    asp = ""
    connected_components_vertex_map = reaction_connected_components(reaction_graph)
    for vertex in reaction_graph.vertices():
        component = connected_components_vertex_map[vertex]
        reaction = reaction_graph.vp["reaction"][vertex]
        asp += "\n"
        asp += (
            f'reaction_graph_strongly_connected_components({component}, "{reaction}").'
        )
    return asp


def reaction_reachability(reaction_graph: gt.Graph):
    # Use Johnson's algorithm to compute the pairwise shortest distance
    # Johnson's algorithm should be faster than Floyd Warshall algorithm for pairwise distance
    # when the graph is sparse, which is expected to be the case for reaction graph
    dist: np.ndarray = gt.topology.shortest_distance(reaction_graph, dense=False)

    # Check what are the pair of vertices that are connected
    for u in reaction_graph.vertices():
        for v in reaction_graph.vertices():
            if dist[u][v] != np.inf:
                yield reaction_graph.vp["reaction"][u], reaction_graph.vp["reaction"][v]


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
        # As the pangenome graph is awfully large, and to make the answer set program solvable in a reasonable time range,
        # we limit the edges to the set of edges with both extremities linked with at least one reaction catalysis
        if u_family in families_to_reactions and v_family in families_to_reactions:
            pangenome_graph_edge_asp += "\n"
            pangenome_graph_edge_asp += f'pangenome_graph_edge("{u}", "{v}").'
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
            pangenome_catalysis_asp
            + "\n\n"
            + reaction_asp
            + f"#const n={limit}."
            + "#show reaction_is_in_pathway/2."
        )
        if SAVE_ASP:
            with open(SAVE_ASP, "w") as f:
                f.write(inline_asp)
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

    logger.info(f"Loading the pangenome graph from {args.pangenome_graph}")
    pangenome_graph: gt.Graph = gt.load_graph(args.pangenome_graph)
    logger.info("Pangenome graph loaded.")
    logger.info("Loading the gene families to reaction catalysis mapping")
    families_to_reaction: dict[str, set[str]] = read_mapping(args.gene_reaction)
    logger.info("Gene families to reaction catalysis mapping loaded.")

    kb: KnowledgeBase = select_kb(config)

    logger.info("Launching the de novo pathway finder.")
    pathway_candidates = denovo_pathway_identification(
        pangenome_graph, families_to_reaction, kb, limit=args.limit
    )
    logger.info(f"Writing the results to the output file {args.output}")

    with open(args.output, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        for pathway_number, pathway_candidate_reaction_set in enumerate(
            pathway_candidates
        ):
            for reaction in pathway_candidate_reaction_set:
                writer.writerow([pathway_number, reaction])
    logger.info("Done.")


if __name__ == "__main__":
    main()
