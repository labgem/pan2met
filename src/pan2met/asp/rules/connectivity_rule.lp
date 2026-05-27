% Maximize connectivity
%
% Connectivity maximisation is dedfine as the minimisation of gaps within a pathway
%
% Required data atoms:
% - pathway(Pathway). % List all pathway identifiers
% - is_in_pathway(Reaction, Pathway).
% - catalysis(Reaction). % Reaction is catalyzed
% - spontaneous(Reaction). % Reaction is spontaneous
% - pathway_reaction_order(Pathway, PredecessorReaction, SuccessorReaction). % In the directed reaction graph, the edge (PredecessorReaction, SuccessorReaction) exists.

% We define a gap as a reaction, or a set of reactions, that lies between two other reactions of the reaction graph.
% When these gap reactions are lost, the path between their wrapping reactions is broken.
% gap_reaction(Reaction, Pathway) means Reaction is a gap in the Pathway

reactome(Reaction) :- catalysis(Reaction).
reactome(Reaction) :- spontaneous(Reaction).

gap_reaction(Reaction, Pathway) :- pathway(Pathway)
                                , not reactome(Reaction)
                                , is_in_pathway(Reaction, Pathway)
                                , has_predecessor(Reaction, Pathway)
                                , has_successor(Reaction, Pathway).
% Propagate predecessor and successors
expanded_predecessor(Pathway, Reaction, PredecessorReaction) :-
                              pathway_reaction_order(Pathway, PredecessorReaction, Reaction).
expanded_predecessor(Pathway, Reaction, AnotherPredecessorReaction) :-
                              pathway_reaction_order(Pathway, AnotherPredecessorReaction, PredecessorReaction).

expanded_successor(Pathway, Reaction, SuccessorReaction) :-
                              pathway_reaction_order(Pathway, Reaction, SuccessorReaction).
expanded_successor(Pathway, Reaction, AnotherSuccessorReaction) :-
                              pathway_reaction_order(Pathway, PredecessorReaction, AnotherSuccessorReaction).

% Identify dead ends nodes, according to the 'reactome'
source_reaction(Pathway, Reaction) :- reactome(Reaction)
                                  , not pathway_reaction_order(Pathway, _, Reaction).
sourced_reaction(Pathway, Reaction) :- source_reaction(Pathway, Reaction).
sourced_reaction(Pathway, Reaction) :- source_reaction(Pathway, SourceReaction) , expanded_successor(Pathway, SourceReaction, Reaction).

source_reaction(Pathway, Reaction) :- reactome(Reaction)
                                  , not sourced_reaction(Pathway, PredecessorReaction)
                                  , pathway_reaction_order(Pathway, PredecessorReaction, Reaction).

sink_reaction(Pathway, Reaction) :- reactome(Reaction)
                                  , not pathway_reaction_order(Pathway, Reaction, _).
sinked_reaction(Pathway, Reaction) :- sink_reaction(Pathway, Reaction).
sinked_reaction(Pathway, Reaction) :- sink_reaction(Pathway, SinkReaction) , expanded_predecessor(Pathway, Reaction, SinkReaction).

sink_reaction(Pathway, Reaction) :- reactome(Reaction)
                                  , not sinked_reaction(SuccessorReaction)
                                  , pathway_reaction_order(Pathway, Reaction, SuccessorReaction).

% We identified source and sink reaction, wrapping a subgraph of thze reaction graph with some of the reactions known to be catalyzed.
% FIXME: orphan reaction ?

% Now we can identify gap following out definition.
% A gap reaction is a reaction in a path between a source and a sink node not being catalyzed.
gap_reaction(Pathway, Reaction) :- not reactome(Reaction)
                                , source_reaction(Pathway, SourceReaction)
                                , sink_reaction(Pathway, SinkReaction)
                                , expanded_successor(Pathway, Reaction, SinkReaction)
                                , expanded_predecessor(Pathway, Reaction, SourceReaction).


% We want to minimize the number of gaps in the pathways, maximizing the "connectivity"$
% We have two options: either counting the tuples Pathway, GapReaction, or the unique Reaction found.

pathway_reaction_gap_count(N) :- N = #count { Pathway, Reaction : gap_reaction(Pathway, Reaction) }.
% We counted the number of tuple 〈 Pathway,GapReactio 〉

reaction_gap_count(N) :- N = #count { Reaction : gap_reaction(Pathway, Reaction) }.
% We counted the number of reaction identifiers involved in a gap of at least one pathway (maybe many.)

% Option 1.
#minimize { N: pathway_reaction_gap_count(N) }.
% Option 2.
% #minimize { N: reaction_gap_count(N) }.
