"""
Reaction graph
"""

from collections.abc import Iterable

import graph_tool as gt

from .knowledge_base import KnowledgeBase


def reaction_graph_from_knowledge_base(
    kb: KnowledgeBase,
) -> tuple[gt.Graph, dict[str, gt.Vertex]]:
    """
    Build a reaction graph, based on the content of the knowledgebase
    """
    edges: set[tuple[str, str]] = set()
    for pathway in kb.pathways():
        for reaction_from, reaction_to in kb.pathway_reaction_order(pathway):
            edge = (reaction_from, reaction_to)
            edges.add(edge)
    reaction_graph, reaction_to_vertex = build_reaction_graph(edges)
    return reaction_graph, reaction_to_vertex


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
