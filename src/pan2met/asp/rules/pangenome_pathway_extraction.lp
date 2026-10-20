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
% % - reaction_graph_edge(Reaction1, Reaction2). : directed, Reaction1 produces a metabolite consumed by Reaction2.
% - reaction_graph_strongly_connected_component(Component, Reaction) : if two reactions belongs to the same strongly connected component, there exists a path from the first reaction to the second, or in the other direction.
% - reaction_are_connected(Reaction1, Reaction2).
% - pangenome_graph_edge(GeneFamily1, GeneFamily2) : undirected, GeneFamily1 and GeneFamily2 are linked in the pangenome gene graph
%
% The output is encoded as follows:
% - pathway(Pathway). : The answer set contains a pathway named Pathway (Pathway is a numeric identifier with no particular meaning)
% - reaction_is_in_pathway(Pathway, Reaction). : The reaction Reaction belongs to the proposed pathway identified by Pathway

pangenome_graph_edge(GeneFamily1, GeneFamily2) :- pangenome_graph_edge(GeneFamily2, GeneFamily1).

{ pathway(1..n) }.

% a pathway identifier is choosen from numbers in range 1 to n

%
% What makes a pathway a "pathway"
%

2 { reaction_is_in_pathway(Pathway, Reaction): reaction(Reaction) } :- pathway(Pathway).
% a pathway is composed of at least two reactions

:- pathway(Pathway)
   , #count { Component :
              reaction_is_in_pathway(Pathway, Reaction),
              reaction_graph_strongly_connected_component(Component, Reaction)
            } > 1.
% a pathway cannot contain two reactions that do not belong to the same connected component


%
% What makes a pathway an 'operonic' pathway in the pangenome
%

% an operonic pathway is a pathway with reactions found contiguous in the pangenome graph

% Two reactions are said to be "pangenome connected" if they are in the same pathway
% and that they are on the same operon

% Two reacions

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

% We do not admit duplicates
% Two pathways are different if there exists a reaction in one but not the other
has_a_different_reaction(Pathway1, Pathway2) :-
    pathway(Pathway1),
    pathway(Pathway2),
    reaction(R),
    reaction_is_in_pathway(Pathway1, R),
    not reaction_is_in_pathway(Pathway2, R).

has_a_different_reaction(Pathway1, Pathway2) :-
    pathway(Pathway1),
    pathway(Pathway2),
    reaction(R),
    reaction_is_in_pathway(Pathway2, R),
    not reaction_is_in_pathway(Pathway1, R).

% Same set of reactions = not different
same_set_of_reactions(Pathway1, Pathway2) :-
    pathway(Pathway1),
    pathway(Pathway2),
    not has_a_different_reaction(Pathway1, Pathway2).

% Eliminate duplicates (keep only one representative)
Pathway1 = Pathway2 :-
    same_set_of_reactions(Pathway1, Pathway2).

% number of reactions in the pathway
pathway_size(N, Pathway) :- N = #count { Reaction: reaction_is_in_pathway(Pathway, Reaction) } , pathway(Pathway).

% number of distinct pathways
pathway_count(N) :- N = #count { Pathway: pathway(Pathway) }.

N >= 4 :- pathway(Pathway) , pathway_size(N, Pathway).
% Instead of optimizing for the size of the pathway, simply search for pathway with size at least 4.

% #maximize{1@2, N: pathway_count(N)}.
% #maximize{1@1, N: pathway_size(N, P), pathway(P)}.
% % we want to put the maximum number of reactions in the pathway module

% N > 1 :- pathway_size(N, Pathway) , pathway(Pathway).
% % a pathway has at least to reactions

