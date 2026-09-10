"""
We use the transitive closures on the (pan)-genome graph
on gene family nodes and a list of candidate gene families with uncertain catalysis.
When a candidate protein catalyzing a reaction of the pathway is within a closure of the pathway,
we get more confidence on this candidate catalyzis.
"""

import logging
from collections import namedtuple
from typing import Dict, Iterable, List, Set, Tuple

import graph_tool as gt

import pan2met.inference.genomic_context
from pan2met.io.knowledge_base import KnowledgeBase
from pan2met.utils import reverse_mapping

GeneFamilyClosureAssignation = namedtuple(
    "GeneFamilyClosureAssignation", ["pathway", "vertex", "closure"]
)


logger = logging.getLogger()


def assign_candidate_enzymes_catalysis(
    pathway: str,
    candidate_catalyzis: Dict[str, Set[str]],
    gene_family_to_pangenome_graph_nodes: Dict[str, List[gt.Vertex]],
    pathway_closures: Dict[str, List[List[Tuple[bool, gt.Vertex]]]],
    kb: KnowledgeBase,
) -> Dict[str, Dict[str, List[GeneFamilyClosureAssignation]]]:
    """
    Assign a candidate gene family with uncertain catalysis
    to a pangenome graph closure on gene family nodes near gene
    family coding enzymes catalyzing reactions of the pathway.


    :param closures: a list of closures, a closure is represented as a list of tuples, representing nodes composing the closure where the first member of the tuple
                    is a flag that evaluates to True when the reaction catalyzed by the gene of the node is a member of the pathway
    :param pangenome_graph: the pangenome graph
    """
    result = {}
    for enzyme, reactions in candidate_catalyzis.items():
        for reaction in reactions:
            result[enzyme] = {}
            assignations = list(
                assign_candidate_enzyme_to_pathway_reactions(
                    pathway,
                    enzyme,
                    reaction,
                    gene_family_to_pangenome_graph_nodes,
                    pathway_closures,
                    kb,
                )
            )
            if assignations != []:
                result[enzyme][reaction] = assignations
    return result


def is_within_a_closure(
    vertex: gt.Vertex, closure: List[Tuple[bool, gt.Vertex]]
) -> bool:
    """
    Check whether a vertex belongs to a listed closure
    """
    for item in closure:
        node = item[1]
        if node == vertex:
            return True
    return False


def assign_candidate_enzyme_to_pathway_reactions(
    pathway: str,
    enzyme_gene_family: str,
    enzyme_reaction: str,
    gene_family_to_pangenome_graph_nodes: Dict[str, List[gt.Vertex]],
    pathway_closures: Dict[str, List[List[Tuple[bool, gt.Vertex]]]],
    kb: KnowledgeBase,
) -> Iterable[GeneFamilyClosureAssignation]:
    logger.info(
        f"Attempt to assign a closure node to {enzyme_reaction} in pathway {pathway}"
    )
    if enzyme_gene_family in gene_family_to_pangenome_graph_nodes:
        for gene_family_vertex in gene_family_to_pangenome_graph_nodes[
            enzyme_gene_family
        ]:
            candidate_pathways = kb.pathways_with_reaction(enzyme_reaction)
            for pathway in candidate_pathways:
                if pathway in pathway_closures:
                    closures = pathway_closures[pathway]
                    for closure in closures:
                        if is_within_a_closure(gene_family_vertex, closure):
                            yield GeneFamilyClosureAssignation(
                                pathway, gene_family_vertex, closure=closure
                            )


def pathway_operon_filler(
    pangenome_graph: gt.Graph,
    kb: KnowledgeBase,
    confident_enzyme_catalyzis: Dict[str, Set[str]],
    less_confident_enzyme_catalyzis: Dict[str, Set[str]],
    distance: int,
    local_edge_jaccard_threshold: float,
) -> dict[str, dict[str, list[GeneFamilyClosureAssignation]]]:
    """

    :return: a dictionnary mapping gene family identifier to a dictionnary with reaction identifier as key and GeneFamilyClosureAssignation as value
    """
    pathways = kb.pathways()
    # Select a subset of the pathways that have at least one listed reaction
    reaction_pathways: set[str] = set()
    for _enzyme, reactions in confident_enzyme_catalyzis.items():
        for reaction in reactions:
            for pathway in kb.pathways_with_reaction(reaction):
                reaction_pathways.add(pathway)
    selected_pathways = set(pathways).intersection(reaction_pathways)

    # Get the reversed mapping of enzyme / reactions
    reaction_to_gene_families = reverse_mapping(confident_enzyme_catalyzis)

    # Get the mapping from gene family identifier to graph tool vertex
    gene_family_to_pangenome_graph_nodes = (
        pan2met.inference.genomic_context.get_pangenome_graph_nid_to_vertex_mapping(
            pangenome_graph
        )
    )

    # Compute the transitive closures for these pathways
    closures = {}
    for pathway in selected_pathways:
        pathway_reactions = kb.reactions_of_pathway(pathway)
        pathway_closures = pan2met.inference.genomic_context.pathway_transitive_closure_connected_components(
            pathway,
            pangenome_graph,
            pathway_reactions,
            reaction_to_gene_families,
            gene_family_to_pangenome_graph_nodes,
            distance,
            local_edge_jaccard_threshold,
        )
        closures[pathway] = list(pathway_closures)
    # Assign the less confident reaction annotations to potential pathway closure
    result = assign_candidate_enzymes_catalysis(
        pathway,
        less_confident_enzyme_catalyzis,
        gene_family_to_pangenome_graph_nodes,
        closures,
        kb,
    )
    return result
