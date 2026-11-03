% Build a bipartite metabolite graph 'de novo'
%
% atoms:
% - metabolite(Molecule).
% - reaction(Reaction).
% - reaction_reactant(Reaction, Metabolite).
% - reaction_product(Reaction, Metabolite).
% - catalysis(Reaction).
% - spontaneous(Reaction).

% The goal is to build a connected subgraph of the whole possible reaction graph
% by selecting a set of connected reactions (the one that are catalyzed + the spontaneous reactions + a selection of "holes" that are required to connect danggling the reactions together).
%
% The output of this logic program is the set of predicate selected_reaction(Reaction),
% where Reaction is taken from reaction/1. and are either a catalyzed reaction,
% a spontaneous reaction, or a pathway hole.


% The set of selected reaction is a selection among the set of reactions.
{ selected_reaction(Reaction): reaction(Reaction) }.

% A reaction cannot be selected if it is neither spontaneous, catalyzed nor considered a hole
:- selected_reaction(Reaction) , not spontaneous(Reaction) , not catalysis(Reaction) , not reaction_hole(Reaction).


% We want to select the metabolic network with the largest number of reactions
count_reaction(N) :- N = #count { Reaction: selected_reaction(Reaction) }.
#maximize { 2@N: count_reaction(N) }.


% We want to minimize the number of reaction holes
count_reaction_holes(N) :- N = #count { Reaction: reaction_hole(Reaction) }.
#minimize { 1@N: count_reaction_holes(N) }.

% A reaction hole is a reaction that is neither spontaneous, nor catalyzed, but selected
reaction_hole(Reaction) :- selected_reaction(Reaction)
                          , not spontaneous(Reaction)
                          , not catalysis(Reaction).

% The set of selected reactions forms a connected graph of reactions
reaction_path(Reaction1, Reaction2) :- selected_reaction(Reaction1)
                                      , selected_reaction(Reaction2)
                                      , Reaction1 = Reaction2.
reaction_path(Reaction1, Reaction2) :- selected_reaction(Reaction1)
                                      , selected_reaction(Reaction2)
                                      , reaction_product(Reaction1, Metabolite)
                                      , reaction_reactant(Reaction2, Metabolite).
reaction_path(Reaction1, Reaction3) :- selected_reaction(Reaction1)
                                      , selected_reaction(Reaction2)
                                      , selected_reaction(Reaction3)
                                      , reaction_path(Reaction1, Reaction2)
                                      , reaction_path(Reaction2, Reaction3).

bidirection_reaction_path(Reaction1, Reaction2) :- reaction_path(Reaction1, Reaction2).
bidirection_reaction_path(Reaction2, Reaction1) :- reaction_path(Reaction1, Reaction2).

:- selected_reaction(Reaction1) , selected_reaction(Reaction2) , not bidirection_reaction_path(Reaction1, Reaction2).
