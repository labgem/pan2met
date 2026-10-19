% Test find several distinct connected components pathways
#const n=3.

reaction("RXN-1").
reaction("RXN-2").
reaction("RXN-3").
reaction("RXN-4").
reaction("RXN-5").
reaction("RXN-6").

spontaneous("RXN-2").

reaction_graph_edge("RXN-1", "RXN-2").
reaction_graph_edge("RXN-2", "RXN-3").
reaction_graph_edge("RXN-3", "RXN-4").
reaction_graph_edge("RXN-4", "RXN-5").
reaction_graph_edge("RXN-5", "RXN-6").

catalysis("RXN-1", "GENE-1").
catalysis("RXN-2", "GENE-2").
catalysis("RXN-3", "GENE-3").
catalysis("RXN-4", "GENE-4").
catalysis("RXN-5", "GENE-5").
catalysis("RXN-6", "GENE-6").

pangenome_graph_edge("GENE-1", "GENE-2").
pangenome_graph_edge("GENE-2", "GENE-3").
pangenome_graph_edge("GENE-4", "GENE-5").
pangenome_graph_edge("GENE-5", "GENE-6").

% #show pangenome_connected_reactions/3.
#show reaction_is_in_pathway/2.
