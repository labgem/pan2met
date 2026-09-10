from collections import defaultdict
from typing import List, Tuple, Iterable

import pytest
import graph_tool as gt

import pan2met.inference.genomic_context


@pytest.fixture
def pangenome_graph():
    graph = gt.Graph(directed=False)
    graph.vp["strains"] = graph.new_vertex_property("vector<string>")
    graph.ep["strains"] = graph.new_edge_property("vector<string>")
    graph.vp["nid"] = graph.new_vertex_property("string")

    # Create a dummy graph
    node1 = graph.add_vertex()
    node2 = graph.add_vertex()
    node3 = graph.add_vertex()
    node4 = graph.add_vertex()
    node5 = graph.add_vertex()

    gene1 = "gene1"
    gene2 = "gene2"
    gene3 = "gene3"
    gene4 = "gene4"
    gene5 = "gene5"

    edge1 = graph.add_edge(node1, node2)
    edge2 = graph.add_edge(node2, node3)
    edge3 = graph.add_edge(node3, node4)
    edge4 = graph.add_edge(node4, node5)

    genomes = ["strain1", "strain2", "strain3"]
    for vertex, gene_family in zip(
        [node1, node2, node3, node4, node5], [gene1, gene2, gene3, gene4, gene5]
    ):
        graph.vp["nid"][vertex] = gene_family
        graph.vp["strains"][vertex] = genomes
    for edge in [edge1, edge2, edge3, edge4]:
        graph.ep["strains"][edge] = genomes
    return graph


@pytest.fixture()
def gene_family_to_pangenome_graph_nodes(pangenome_graph):
    mapping = defaultdict(list)
    for node in pangenome_graph.vertices():
        family = pangenome_graph.vp["nid"][node]
        mapping[family].append(node)
    return dict(mapping)


def test_genomic_context(pangenome_graph, gene_family_to_pangenome_graph_nodes):
    reaction1 = "rxn-1"
    reaction2 = "rxn-2"
    reaction3 = "rxn-3"

    pathway_reactions = [reaction1, reaction2, reaction3]

    reaction_to_gene_families = {
        "rxn-1": ["gene1"],
        "rxn-2": ["gene2"],
        "rxn-3": ["gene3"],
    }

    closures = pan2met.inference.genomic_context.pathway_transitive_closure_connected_components(
        pangenome_graph,
        pathway_reactions,
        reaction_to_gene_families,
        gene_family_to_pangenome_graph_nodes,
        distance=3,
        local_edge_jaccard_threshold=0.8,
    )

    assert len(list(closures)) == 1


def closure_to_gene_list(
    pangenome_graph, closure: List[Tuple[bool, str]]
) -> Iterable[str]:
    for _flag, node in closure:
        yield pangenome_graph.vp["nid"][node]


def test_unreachable_reaction_genomic_context(
    pangenome_graph, gene_family_to_pangenome_graph_nodes
):
    reaction1 = "rxn-1"
    reaction2 = "rxn-2"
    reaction3 = "rxn-3"

    pathway_reactions = [reaction1, reaction2, reaction3]

    reaction_to_gene_families = {
        "rxn-1": ["gene1"],
        "rxn-2": ["gene2"],
        "rxn-3": ["gene5"],
    }

    closures = pan2met.inference.genomic_context.pathway_transitive_closure_connected_components(
        pangenome_graph,
        pathway_reactions,
        reaction_to_gene_families,
        gene_family_to_pangenome_graph_nodes,
        distance=1,
        local_edge_jaccard_threshold=0.8,
    )
    closures = list(closures)

    assert len(closures) == 2

    closures_genes = list(
        map(
            lambda closure: list(closure_to_gene_list(pangenome_graph, closure)),
            closures,
        )
    )

    assert closures_genes == [["gene1", "gene2", "gene3"], ["gene5", "gene4"]]


def test_transitively_reachable_reaction_genomic_context(
    pangenome_graph, gene_family_to_pangenome_graph_nodes
):
    reaction1 = "rxn-1"
    reaction2 = "rxn-2"
    reaction3 = "rxn-3"

    pathway_reactions = [reaction1, reaction2, reaction3]

    reaction_to_gene_families = {
        "rxn-1": ["gene1"],
        "rxn-2": ["gene2"],
        "rxn-3": ["gene5"],
    }

    closures = pan2met.inference.genomic_context.pathway_transitive_closure_connected_components(
        pangenome_graph,
        pathway_reactions,
        reaction_to_gene_families,
        gene_family_to_pangenome_graph_nodes,
        distance=2,
        local_edge_jaccard_threshold=0.8,
    )
    closures = list(closures)

    assert len(closures) == 1

    closures_genes = list(
        map(
            lambda closure: list(closure_to_gene_list(pangenome_graph, closure)),
            closures,
        )
    )

    assert closures_genes == [["gene1", "gene2", "gene3", "gene4", "gene5"]]
