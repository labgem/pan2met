"""
Find all linear sub-pathways from a reaction graph.

We build the reaction graph from the set of edges
in the reaction graphs of the metabolic pathways in the knowledgebase.
In effect, this will make possibly possibly unwanted gaps in the metabolic network,
when reactions are in the boundary of metabolic pathways.

"""

import argparse
import csv
from collections.abc import Iterable

import graph_tool as gt

from ..config import default_config, override_config
from ..io.knowledge_base import KnowledgeBase, select_kb
from ..io.reaction_graph import reaction_graph_from_knowledge_base
from ..utils import logger

REMOVE_EXTREMITIES = True


def find_all_non_branching_path(G: gt.Graph) -> Iterable[list[int]]:
    """
    Enumerate all maximal non-branching path in a graph.

    A non branching path is a path between two edge whose inner nodes have in-degree and out-degree equal to one.

    :param G: a graph_tool directed graph
    :yield: linear path, as an ordered list of vertices
    """
    visited: set[gt.Vertex] = set()
    for v in G.get_vertices():
        if v in visited:
            continue
        else:
            visited.add(v)
        # Depth first search while vertex is not branching in both direction.
        in_edges = G.get_in_edges(v)
        out_edges = G.get_out_edges(v)

        path: list[gt.Vertex] = [v]
        # Extend backward
        while len(in_edges) == 1:
            current = in_edges[0][0]
            if current in visited:
                break
            else:
                visited.add(current)
            in_edges = G.get_in_edges(current)
            path = [current] + path
        # Extend forward
        while len(out_edges) == 1:
            current = out_edges[0][0]
            if current in visited:
                break
            else:
                visited.add(current)
            visited.add(current)
            out_edges = G.get_in_edges(current)
            path.append(current)
        if len(path) > 1:
            yield path


def reactions_of_linear_path(
    linear_path: list[int], reaction_graph: gt.Graph
) -> Iterable[str]:
    """
    Enumerate the reactions of a linear path in a reaction graph

    :param linear_path: vertex list, path in the reaction graph
    :param metabolic_graph: graph_tool undirected graph, with vertex property attribute `.vp` with key "type" taking value in {"metabolite", "reaction"}
        and key "id" taking value in the the name of the reaction, or compound in the reference database / sbml
    :yield: the name of reaction in the linear path.
    """
    if REMOVE_EXTREMITIES:
        linear_path = linear_path[1:-1]
    for item in linear_path:
        yield reaction_graph.vp["reaction"][item]


def associate_reaction_set_to_pathway(
    kb: KnowledgeBase, reaction_set: set[str]
) -> Iterable[str]:
    """
    Associate a linear path to a (set of) pathway(s)

    :param linear_path: a vertex list with non-branching path in the metabolic graph
    :param padmet_object: a PADMet format object with information on pathway and reactions
    :param metabolic_graph: a graph_tool object with vertex property dictionnary in attribute

    :yield: name of metabolic pathways where a reaction of the linear path is involved in
        # TODO: filter those pathway to retain only the pathway where the whole linear sub-pathway is present.
    """
    for reaction_set_item in reaction_set:
        for pathway in kb.pathways():
            for pathway_reaction in kb.reactions_of_pathway(pathway):
                if reaction_set_item == pathway_reaction:
                    yield pathway
                    break


def find_linear_subpathways(
    reaction_graph: gt.Graph, kb: KnowledgeBase
) -> Iterable[tuple[list[gt]]]:
    """
    Find linear sub-pathways.

    A linear sub-pathway is the longest non-branching simple path in a pathway that appears only in one pathway of the metabolome.

    It is defined in Unipathway article (doi: 10.1093/nar/gkr1023) as:

    .. quote::

        A 'Linear Subpathway (ULS)', also simply called
        'subpathway', is a chemical transformation from a
        multi-set of initial compounds [substrate(s)] to a multi-
        set of final compounds [product(s)] that does not contain
        any branching reaction or cycle (Figure 2c). More pre-
        cisely, if we deﬁne the 'reaction graph' as a graph where
        vertices are reactions and two vertices are linked by an
        A set of subpathways (ULS) and their connecting
        compounds edge where the product of one reaction is the substrate
        of the next, then ULS are simply paths in this graph.
        Technically, an ULS is therefore an ordered sequence of
        UERs (unipathway reaction). UERs within a ULS are linked via their primary
        metabolites, which are deﬁned according to the context of
        the pathway [...].

    :param reaction_graph: a directed reaction graph
    :param kb: a metabolism knowledgebase
    """
    for linear_path in find_all_non_branching_path(reaction_graph):
        reaction_set: set[str] = set(
            reactions_of_linear_path(linear_path, reaction_graph)
        )
        pathways: Iterable[str] = associate_reaction_set_to_pathway(kb, reaction_set)
        pathways: list[str] = list(pathways)
        if len(pathways) == 1:
            yield reaction_set, pathways[0]
        elif len(pathways) > 1:
            logger.warning(
                f"Found a linear path matching several pathways: {reaction_set}, in {len(pathways)} pathway"
            )


def argument_parser():
    parser = argparse.ArgumentParser(
        description="Identify the sets of linear reaction path that are specific to a pathway."
    )
    parser.add_argument(
        "-c", "--config", help="override default config", required=False
    )
    parser.add_argument(
        "-o",
        "--output",
        help="output filename (defaults to stdout)",
        required=False,
        default="/dev/stdout",
    )
    return parser


def main():
    parser = argument_parser()
    args = parser.parse_args()
    if args.config:
        config = override_config(args.config, default_config)
    else:
        config = default_config

    kb: KnowledgeBase = select_kb(config)
    reaction_graph, _reaction_to_vertex = reaction_graph_from_knowledge_base(kb)

    with open(args.output, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        header = ["pathway", "linear path reactions"]
        writer.writerow(header)
        for linear_path_reactions, pathway in find_linear_subpathways(
            reaction_graph, kb
        ):
            record = [pathway, ",".join(linear_path_reactions)]
            writer.writerow(record)


if __name__ == "__main__":
    main()
