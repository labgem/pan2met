
reaction("RXN-1").
reaction("RXN-2").
reaction("RXN-3").
reaction("RXN-4").

spontaneous("RXN-2").

reaction_graph_edge("RXN-1", "RXN-2").
reaction_graph_edge("RXN-2", "RXN-3").
reaction_graph_edge("RXN-3", "RXN-4").

catalysis("RXN-1", "GENE-1").
catalysis("RXN-1", "GENE-1-prime").
catalysis("RXN-2", "GENE-2").
catalysis("RXN-3", "GENE-3").
catalysis("RXN-4", "GENE-4").

pangenome_graph_edge("GENE-1", "GENE-2").
pangenome_graph_edge("GENE-2", "GENE-3").
pangenome_graph_edge("GENE-3", "GENE-4").

% #show pangenome_connected_reactions/3.
#show pathway/1.
#show reaction_is_in_pathway/2.
