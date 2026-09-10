"""
Identify a (pan)-genome genomic context of a pathway
"""

import logging
import queue
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Set, Tuple

import graph_tool as gt

logger = logging.getLogger()


def get_pangenome_graph_nid_to_vertex_mapping(
    pangenome_graph: gt.Graph,
) -> Dict[str, List[gt.Vertex]]:
    """
    Get a mapping from a gene family to the list of pangenome graph vertex associated with this gene family.
    """
    mapping = defaultdict(list)
    for v in pangenome_graph.vertices():
        family_id = pangenome_graph.vp["nid"][v]
        mapping[family_id].append(v)
    return dict(mapping)


def local_edge_jaccard(
    pangenome_graph: gt.Graph, edge: gt.Edge, genomes: Set[str]
) -> float:
    """
    A Jaccard index is computed on each edge.
    For two families $u$ and $v$ connected by an edge $e_{u, v}$,
    the index $J_{u, v}$ is defined as the ratio of the number of genomes
    in which the edge $e_{u, v}$ exists over the number of genomes in which
    at least one of the two gene families is present.

    :param pangenome_graph: a graph_tool graph datastructure
        with vertex property "strains" listing the organisms having the gene family vertex
        and the edge property "strains" listing the organisms having the edge.
    :param edge: the edge between vertex u and v in `pangenome_graph` datastructure.
    :param genomes: the set of genomes to consider, i.e., the set of genomes where at least one gene families coding
                    for a reaction of the pathway of interest is found.
    :return: the Jaccard index
    """
    u = edge.source()
    v = edge.target()
    genomes_u: List[str] = pangenome_graph.vp["strains"][u]
    genomes_v: List[str] = pangenome_graph.vp["strains"][v]

    edge_genomes: List[str] = pangenome_graph.ep["strains"][edge]
    edge_genomes_pathway: Set[str] = set(edge_genomes).intersection(genomes)

    union_vertices_genomes = set(genomes_u).union(set(genomes_v))
    union_vertices_genomes_pathway = union_vertices_genomes.intersection(genomes)

    if len(edge_genomes_pathway) == 0:
        return 0
    if len(union_vertices_genomes_pathway) == 0:
        return 0
    jaccard_index = len(edge_genomes_pathway) / len(union_vertices_genomes_pathway)
    # logger.debug(
    #     f"Local Edge Jaccard index J({edge}) = {len(edge_genomes_pathway)} / {len(union_vertices_genomes)} = {jaccard_index}"
    # )
    return jaccard_index


def list_genomes_linked_to_the_pathway(
    pathway_reactions: List[str],
    reaction_to_gene_families: Dict[str, Set[str]],
    gene_family_to_pangenome_graph_nodes: Dict[str, List[gt.Vertex]],
    pangenome_graph: gt.Graph,
) -> Set[str]:
    """
    List the genomes where at least one reaction of a pathway have at least one enzyme that catalyzes it.
    :param pathway_id: a pathway identifier
    :param kb: a metabolic pathway knowledge base
    :param reaction_to_gene_families: a mapping of reaction identifier to pangenome gene families identifier
    :param gene_family_to_pangenome_nodes: a mapping of gene family identifier to pangenome graph nodes
    :return: a set of such genomes
    """
    genomes: Set[str] = set()
    genome_reactions: Dict[str, Set[str]] = {}
    for reaction in pathway_reactions:
        genome_reactions[reaction] = set()
        if reaction in reaction_to_gene_families:
            for gene_family in reaction_to_gene_families[reaction]:
                if gene_family in gene_family_to_pangenome_graph_nodes:
                    for node in gene_family_to_pangenome_graph_nodes[gene_family]:
                        for strain in pangenome_graph.vp["strains"][node]:
                            genome_reactions[reaction].add(strain)
        genomes = genomes.union(genome_reactions[reaction])
    # for reaction in pathway_reactions:
    #     if len(genome_reactions[reaction]) > 0:
    #         genome_reactions_reaction_only = genome_reactions[reaction].copy()
    #         for other_reaction in pathway_reactions:
    #             genome_reactions_reaction_only.difference(genome_reactions[other_reaction])

    #         genome_reaction_contribution_percentage = len(genome_reactions_reaction_only) / len(genomes) * 100
    # logger.info(f"Reaction {reaction} contributed {genome_reaction_contribution_percentage}% of the total listed genomes for the pathway, among {len(pathway_reactions)} reactions in the pathway.")
    return genomes


def pathway_transitive_closure_connected_components(
    pathway: str,
    pangenome_graph: gt.Graph,
    pathway_reactions: List[str],
    reaction_to_gene_families: Dict[str, Set[str]],
    gene_family_to_pangenome_graph_nodes: Dict[str, List[gt.Vertex]],
    distance: int,
    local_edge_jaccard_threshold: float,
) -> Iterable[List[Tuple[bool, List[gt.Vertex]]]]:
    """

    :param pangenome_graph:
    :param pathway_reactions: a list of reactions of the pathway
    :param distance: maximum number of gaps between two nodes with a gene family encoding an enzyme catalyzing a reaction of the pathway (e.g, 3)
    :param local_edge_jaccard_threshold: the minimum threshold value for the Jaccard index of a edge we will consider
    :return: an iterable of lists of tuple with a Boolean flag stating whether the pangenome graph gene family vertex belongs to the pathway known enzymes, listing the pangenome graph vertices belonging to each of the transitive closures.
    """
    pathway_genomes = list_genomes_linked_to_the_pathway(
        pathway_reactions,
        reaction_to_gene_families,
        gene_family_to_pangenome_graph_nodes,
        pangenome_graph,
    )
    logger.info(
        f"We consider only {len(pathway_genomes)} genomes for the edge local Jaccard index for this pathway"
    )

    def filter_edge_test(edge: gt.Edge) -> bool:
        """
        Check whether a edge pass the local Jaccard index criterion.
        """
        return (
            local_edge_jaccard(pangenome_graph, edge, pathway_genomes)
            > local_edge_jaccard_threshold
        )

    def reset_depth(successor_node: gt.Vertex) -> bool:
        """
        Check whether a vertex is associated with a gene family coding an enzyme catalyzing a reaction of the pathway.
        """
        family_id: str = pangenome_graph.vp["nid"][successor_node]
        for reaction in pathway_reactions:
            if reaction in reaction_to_gene_families:
                gene_families = reaction_to_gene_families[reaction]
                if family_id in gene_families:
                    return True
        return False

    def transitive_closure(
        seed_node: gt.Vertex, visited: Set[gt.Vertex]
    ) -> Optional[List[gt.Vertex]]:
        if seed_node in visited:
            return None
        closure = [(True, seed_node)]
        q: queue.PriorityQueue[Tuple[int, gt.Vertex]] = queue.PriorityQueue()
        depth = 0  # We start with a gene family node belonging to the set of enzyme catalyzing a reaction of the pathway
        q.put((depth, seed_node))
        while not q.empty():
            depth, node = q.get()
            if depth <= distance:
                for edge in node.out_edges():
                    if filter_edge_test(edge):
                        successor_depth = depth + 1
                        successor_node = edge.target()
                        flag = False
                        if reset_depth(successor_node):
                            successor_depth = 0
                            flag = True
                        if successor_node not in visited:
                            if successor_depth <= distance:
                                closure.append((flag, successor_node))
                                q.put((successor_depth, successor_node))
            visited.add(node)
        return closure

    visited: Set[gt.Vertex] = set()
    for reaction in pathway_reactions:
        if reaction in reaction_to_gene_families:
            reaction_enzyme_gene_families = reaction_to_gene_families[reaction]
            for gene_family in reaction_enzyme_gene_families:
                if gene_family in gene_family_to_pangenome_graph_nodes:
                    for seed_node in gene_family_to_pangenome_graph_nodes[gene_family]:
                        if seed_node not in visited:
                            closure = transitive_closure(seed_node, visited)
                            if closure:
                                yield closure
