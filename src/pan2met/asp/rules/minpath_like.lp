% A pathway inferrence inspired by MinPath by Ye and Doak 2009
% A Parsimony Approach to Biological Pathway Reconstruction/Inference for Genomes and Metagenomes. DOI: 10.1371/journal.pcbi.1000465
%
% In Ye and Doak, 2009, they use Mixed Integer Linear Programming to select the pathways
% that contains a reaction identified in the reactome
% minimizing the set of pathway identified.
%
% This Answer Set Programming implementation tries to reproduce
% this idea, using logic programming instead of mixed integer linear programming.
%
% We use the following atoms:
%
% catalysis(Reaction). % For reactions expected to be catalyzed
% is_in_pathway(Reaction, Pathway). % When the pathway Pathway contains the reaction Reaction.
% pathway(Pathway). % Pathway is an identifier of a pathway.
%
% That's all we need in the input knowledge base.

% We identify the set of pathways
% that have at least one reaction in the reactome.
candidate_pathway(Pathway) :- pathway(Pathway)
                           , catalysis(Reaction)
                           , is_in_pathway(Reaction, Pathway).

% The inferred set of pathway is a subset of candidate_pathway/1.
{ inferred_pathway(Pathway) : candidate_pathway(Pathway) }.

% MinPath idea is to identify the set of pathways that is sufficient
% to explain the presence of a reaction in the reactome,
% by an involvement in a pathway.
% So, we minimize the number of inferred_pathway.
#minimize { N : pathway_count(N) }.
pathway_count(N) :- N=#count { Pathway : inferred_pathway(Pathway) }.

% However, we still want to explain all reactions by an involvement in a pathway.
% Hopefully the pathway knowledge base cover is large enough to contain all reaction we can annotate.
covered_reaction(Reaction) :- inferred_pathway(Pathway) , is_in_pathway(Reaction, Pathway).
:- reactome(Reaction) , not covered_reaction(Reaction).

% pathway("PWY1").
% pathway("PWY2").
% pathway("PWY3").
% pathway("PWY4").
% is_in_pathway("RXN1", "PWY1").
% is_in_pathway("RXN2", "PWY2").
% is_in_pathway("RXN3", "PWY1").
% is_in_pathway("RXN4", "PWY3").
% is_in_pathway("RXN4", "PWY1").
% reactome("RXN1"; "RXN2"; "RXN3").
% #show inferred_pathway/1.
% % Expect PWY4 not to be selected in the inferred pathway set.
