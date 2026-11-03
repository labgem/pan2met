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

% We want to select the metabolic network with the largest number of reactions
count_reaction(N) :- N = #count { Reaction: selected_reaction(Reaction) }.
#maximize { 2@N: count_reaction(N) }.

% We want to minimize the number of reaction holes
count_reaction_holes(N) :- N = #count { Reaction: reaction_hole(Reaction) }.

% A reaction hole is a reaction that is neither spontaneous, nor catalyzed, but selected
reaction_hole(Reaction) :- selected_reaction(Reaction)
                          , not spontaneous(Reaction)
                          , not catalysis(Reaction).


% define a seed reaction with seed(Reaction).

reachable(Reaction) :- seed(Reaction) , selected_reaction(Reaction).
reachable(Reaction1) :- selected_reaction(Reaction1)
                     , selected_reaction(Reaction2)
                     , reaction_product(Reaction1, Metabolite)
                     , reaction_reactant(Reaction2, Metabolite)
                     , reachable(Reaction2).

reachable(Reaction1) :- selected_reaction(Reaction1)
                     , selected_reaction(Reaction2)
                     , reaction_reactant(Reaction1, Metabolite)
                     , reaction_product(Reaction2, Metabolite)
                     , reachable(Reaction2).

% ensure all selected reaction are in the same component
:- selected_reaction(Reaction)
   , not reachable(Reaction).

:- selected_reaction(Reaction)
   , potential_hole(Reaction).

%
potential_hole(Reaction) :- reaction(Reaction)
                          , not catalysis(Reaction)
                          , not spontaneous(Reaction).

% Only allow selection of catalyzed, spontaneous, or potential holes
{ selected_reaction(Reaction) : catalysis(Reaction) }.
{ selected_reaction(Reaction) : spontaneous(Reaction) }.
{ selected_reaction(Reaction) : potential_hole(Reaction) }.

#show selected_reaction/1.
