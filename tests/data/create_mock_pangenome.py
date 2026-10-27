import graph_tool as gt


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


graph = pangenome_graph()

graph.save("mock_pangenome_graph.gt")
