"""
Find metabolic pathways 'de novo'.

We consider the following type of metabolic pathways:

- connected components in the reaction graph
- each non spontaneous reactions of the pathway is inside a connected component subgraph of the pangenome.
"""

import argparse
import csv
import itertools
from collections.abc import Iterable

import graph_tool as gt

from pan2met.config import default_config, override_config
from pan2met.io.knowledge_base import KnowledgeBase, select_kb
from pan2met.utils import logger, read_mapping


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-g",
        "--pangenome-graph",
        help="path to pangenome graph in graph-tool binary gt format.",
        required=True,
    )
    parser.add_argument(
        "--gene-reaction",
        help="gene reaction TSV table (column 1: gene family, column 2: reaction identifier)",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="output path to the output file, a TSV table (column 1: pathway numeric identifier, column 2: reaction identifier, column 3: gene family identifier of the transitive closure catalyzing the reaction)",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--config",
        help="path to configuration file to override the default config",
        required=False,
    )
    return parser


def build_reaction_graph(
    reaction_graph_edges: Iterable[tuple[str, str]],
) -> tuple[gt.Graph, dict[str, gt.Vertex]]:
    """
    Build a graph tool directed graph from a set of edges

    :return: a graph_tool directed graph, with vertex property map 'reaction' with the reaction identifier, and a dictionnary mapping reaction to the corresponding vertex
    """
    # Build the graph_tool graph
    reaction_graph = gt.Graph(directed=True)
    reaction_to_vertex: dict[str, gt.Vertex] = {}
    vertex_to_reaction: dict = reaction_graph.new_vertex_property("string")
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
    return (reaction_graph, reaction_to_vertex)


def check_reaction_is_connected(
    anchor_reactions: set[str],
    reaction: str,
    reaction_graph: gt.Graph,
    reaction_to_vertex: dict[str, set[gt.Vertex]],
) -> bool:
    if reaction not in reaction_to_vertex:
        return False
    reaction_vertex = reaction_to_vertex[reaction]
    for anchor_reaction in anchor_reactions:
        if anchor_reaction in reaction_to_vertex:
            anchor_reaction_vertex = reaction_to_vertex[anchor_reaction]
            if reaction_graph.edge(
                anchor_reaction_vertex, reaction_vertex
            ) or reaction_graph.edge(reaction_vertex, anchor_reaction_vertex):
                return True
    return False


# def find_candidate_pathways(
#     pangenome_graph: gt.Graph,
#     reaction_graph: gt.Graph,
#     reaction_to_vertex: dict[str, gt.Vertex],
#     families_to_reaction: dict[str, set[str]],
#     max_depth: int = 3,
# ) -> Iterable[set[str]]:
#     """
#     Find candidate operonic pathways.
#     We call an "operonic pathway" a pathway whose reactions
#     are either spontaneous or catalyzed by enzymes found in the same transitive closure in the graph.
#     We retain all pangenome graph edge.

#     :param pangenome_graph: a pangenome graph, gene family nodes with vertex property map 'nid' for gene family identifier.
#     :param reaction_graph: a reaction graph with vertex property map 'reaction' for the reaction identifier.
#     :param reaction_to_vertex: a dictionnary mapping the reaction name to the corresponding node in the reaction graph.
#     :param max_depth: maximum depth to extend the pangenome graph search space, to extend the pathway operon search space.

#     :yield: sets of reactions belonging to a pathway whose reactions forms a connected components, of spontaneous reactions or reactions catalyzed in an operon.
#     """


#     visited_pangenome_nodes: set[gt.Vertex] = set()
#     visited_reactions: set[str] = set()
#     frontier: queue.PriorityQueue[tuple[int, gt.Vertex]] = queue.PriorityQueue()

#     candidate_node_reactions: set[tuple[gt.Vertex, set[str]]] = set() # candidate nodes with reactions that could be introduced in the building pathway

#     for pangenome_node in pangenome_graph.vertices():
#         if pangenome_node not in visited_pangenome_nodes:
#             gene_family: str = pangenome_graph.vp["nid"].get(pangenome_node, None)
#             if gene_family is not None and gene_family in families_to_reaction:
#                 candidate_reaction_anchors = families_to_reaction[gene_family]
#                 # We use this pangenome node as an anchor if it catalyses a reaction.
#                 # We perform the exploration starting from each reaction.
#                 # If the starting pangenome node catalyses several reactions, we keep them in the pathway too.
#                 for reaction_anchor in candidate_reaction_anchors:

#                     pathway_reaction_anchor =

#                     depth: int = 0
#                     frontier.put((pangenome_node, depth))
#                     while not frontier.empty():
#                         depth, current_node= frontier.get()
#                         if current_node not in visited_pangenome_nodes:
#                             visited_pangenome_nodes.add(current_node)
#                         if depth <= max_depth:
#                             next_depth = depth + 1
#                             for successor_vertex in current_node.all_neighbors():
#                                 frontier.put((next_depth, successor_vertex))


def find_neighbor_reactions(
    pangenome_graph: gt.Graph,
    reaction_graph: gt.Graph,
    reaction_to_vertex: dict[str, gt.Vertex],
    families_to_reaction: dict[str, set[str]],
):
    """
    We start with a simpler problem.
    Are there pairs of gene family nodes that catalyse successive reactions in a pathway?
    """
    for starting_pangenome_node in pangenome_graph.vertices():
        starting_pangenome_node_gene_family = pangenome_graph.vp["nid"][
            starting_pangenome_node
        ]
        starting_pangenome_node_reactions = families_to_reaction.get(
            starting_pangenome_node_gene_family, None
        )
        if starting_pangenome_node_reactions is not None:
            for neighbor_pangenome_node in starting_pangenome_node.all_neighbors():
                neighbor_pangenome_node_gene_family = pangenome_graph.vp["nid"][
                    neighbor_pangenome_node
                ]
                neighbor_pangenome_node_reactions = families_to_reaction.get(
                    neighbor_pangenome_node_gene_family, None
                )
                if neighbor_pangenome_node_reactions is not None:
                    for from_reaction, to_reaction in itertools.product(
                        starting_pangenome_node_reactions,
                        neighbor_pangenome_node_reactions,
                    ):
                        from_reaction_vertex = reaction_to_vertex.get(
                            from_reaction, None
                        )
                        to_reaction_vertex = reaction_to_vertex.get(to_reaction, None)
                        if (
                            from_reaction_vertex is not None
                            and to_reaction_vertex is not None
                            and reaction_graph.edge(
                                from_reaction_vertex, to_reaction_vertex
                            )
                            is not None
                        ):
                            yield (
                                starting_pangenome_node_gene_family,
                                neighbor_pangenome_node_gene_family,
                                from_reaction,
                                to_reaction,
                            )


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

    logger.info("Building the reaction graph")
    kb: KnowledgeBase = select_kb(config)
    edges: set[tuple[str, str]] = set()
    for pathway in kb.pathways():
        for reaction_from, reaction_to in kb.pathway_reaction_order(pathway):
            edge = (reaction_from, reaction_to)
            edges.add(edge)
    reaction_graph, reaction_to_vertex = build_reaction_graph(edges)
    logger.info("")

    with open(args.output, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        header = ["gene family 1", "gene family 2", "reaction 1", "reaction 2"]
        writer.writerow(header)
        for record in find_neighbor_reactions(
            pangenome_graph, reaction_graph, reaction_to_vertex, families_to_reaction
        ):
            writer.writerow(record)
    logger.info("Done.")


if __name__ == "__main__":
    main()
