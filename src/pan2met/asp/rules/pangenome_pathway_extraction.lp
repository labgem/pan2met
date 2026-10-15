% pangenome metabolic pathway extraction:
% identifying candidate pathway modules from pangenomes
%
% in this model, we consider a pathway module to be
% a connected component in a reaction graph
% where reactions are either spontaneous, or catalyzed by a known enzyme.
%
% We ignore orphan reactions.
% As a consequence, no metabolic pathway with such a reaction could appear
% in any answer set of this model.

% We use the following input data atoms:
% - reaction(Reaction). : Reaction is a reaction.
% - spontaneous(Reaction). : Reaction is a spontaneous reaction
% - catalysis(Reaction, GeneFamily). : Reaction is catalyzed by an enzyme belonging to the GeneFamily
% - reaction_graph_edge(Reaction1, Reaction2) : directed, Reaction1 produces a metabolite consumed by Reaction2.
% - pangenome_graph_edge(GeneFamily1, GeneFamily2) : undirected, GeneFamily1 and GeneFamily2 are linked in the pangenome gene graph
%
% The output is encoded as follows:
% - pathway(Pathway). : The answer set contains a pathway named Pathway (Pathway is a numeric identifier with no particular meaning)
% - reaction_is_in_pathway(Pathway, Reaction). : The reaction Reaction belongs to the proposed pathway identified by Pathway

#const n=100.

pangenome_graph_edge(GeneFamily1, GeneFamily2) :- pangenome_graph_edge(GeneFamily2, GeneFamily1).

pathway(1).

% a pathway identifier is choosen from numbers in range 1 to n

%
% What makes a pathway a "pathway"
%

2 { reaction_is_in_pathway(Pathway, Reaction): reaction(Reaction) } :- pathway(Pathway).
% a pathway is composed of at least two reactions

:- pathway(Pathway)
                                                        , reaction_is_in_pathway(Pathway, Reaction1)
                                                        , reaction_is_in_pathway(Pathway, Reaction2)
                                                        , not reactions_are_connected(Reaction1, Reaction2).
% a pathway is not a connected component if it has a two not connected reactions

reactions_are_connected(Reaction1, Reaction2) :- reaction_graph_path(Reaction1, Reaction2).
reactions_are_connected(Reaction1, Reaction2) :- reaction_graph_path(Reaction2, Reaction1).
% to check if two reactions are connected, we check if there is a directed edge in either direction,

reaction_graph_path(Reaction1, Reaction2) :- Reaction1 = Reaction2 , reaction(Reaction1) , reaction(Reaction2).
reaction_graph_path(Reaction1, Reaction2) :- reaction_graph_edge(Reaction1, Reaction2).
reaction_graph_path(Reaction1, Reaction3) :- reaction_graph_path(Reaction1, Reaction2) , reaction_graph_path(Reaction2, Reaction3).
% a directed path between Reaction1 and Reaction2 follows the directed edges, edge to edge.

% we reject pathway that are not connected components

%
% What makes a pathway an 'operonic' pathway in the pangenome
%

% an operonic pathway is a pathway with reactions found contiguous in the pangenome graph

pangenome_connected_reactions(Pathway, Reaction1, Reaction3) :- pangenome_connected_reactions(Pathway, Reaction1, Reaction2)
                                                     , pangenome_connected_reactions(Pathway, Reaction2, Reaction3).
pangenome_connected_reactions(Pathway, Reaction, Reaction) :- reaction_is_in_pathway(Pathway, Reaction).
pangenome_connected_reactions(Pathway, Reaction1, Reaction2) :- catalysis(Reaction1, Gene1)
                                                     , catalysis(Reaction2, Gene2)
                                                     , pangenome_graph_edge(Gene1, Gene2)
                                                     , reaction_is_in_pathway(Pathway, Reaction1)
                                                     , reaction_is_in_pathway(Pathway, Reaction2).

% As a starting point, we consider only the transitive closure with no gaps.
% Two reactions catalyzed by adjacent enzymes belongs to the same transitive closure.

:- not pangenome_connected_reactions(Pathway, Reaction1, Reaction2)
   , reaction(Reaction1)
   , reaction(Reaction2)
   , reaction_is_in_pathway(Pathway, Reaction1)
   , reaction_is_in_pathway(Pathway, Reaction2)
   , not spontaneous(Reaction1).

:- not pangenome_connected_reactions(Pathway, Reaction1, Reaction2)
   , reaction(Reaction1)
   , reaction(Reaction2)
   , reaction_is_in_pathway(Pathway, Reaction1)
   , reaction_is_in_pathway(Pathway, Reaction2)
   , not spontaneous(Reaction2).

%
% What makes a pathway a 'good' candidate
%

pathway_size(N, Pathway) :- N = #count { Reaction: reaction_is_in_pathway(Pathway, Reaction) } , pathway(Pathway).

#maximize{N: pathway_size(N, P)}.
% % we want to put the maximum number of reactions in the pathway module

% N > 1 :- pathway_size(N, Pathway) , pathway(Pathway).
% % a pathway has at least to reactions

