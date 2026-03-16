% A PathoLogic Like implementation of the decision rules in answer set programming (Clingo dialect)
% We take as reference the Pathway Tools v29.5 User Guide PDF p. 188, section "Metabolic Pathway Prediction"

% What do we need?
% Here is a listing of what we require as input data atoms for the prediction.
% - is_in_pathway(Reaction, Pathway). --
% - reaction(Reaction). -- Reaction is a reaction
% - reactome(Reaction). -- We know the reaction is in the organism reactome, or at least we have evidence of a putative enzyme catalyzing it in the proteome.
% - spontaneous(Reaction). -- We assume the reaction does not need a catalyzis to be present
% - orphan(Reaction). -- We do not know any enzyme catalyzing the reaction
% - pathway_key_reaction(Pathway, KeyReaction). -- KeyReaction is considered key for the pathway
% - pathway_in_taxonomic_range(Pathway). -- The target taxon is in the expected taxonomic range of the pathway
% - pathway_reaction_order(Pathway, Reaction, Order). -- If Reaction is the first in a topological ordering of the reaction graph, Order is 1, Order is maximum if the reaction is last in the topological ordering.
%   Used for a special rule for catabolysis and biosynthesis pathways.
% - pathway_passing_score(Pathway).


% In the PathoLogic algorithm, the decision rules are run sequentially.
% To mimick such a procedure in clingo Answer Set Programming, we will use
% a little trick, we will add a dependency for precedent rules in subsequent rules.

% First rule:
% REJECT P if P is a transport, signaling, or synthetic (engineered) pathway. Rationale: Patho-
% Logic predicts only natural metabolic pathways.
ignored_pathway_class("Transport-Pathways"; "Signaling-Pathways"). % Synthetic pathways seems not to be defined as an ontology parent class in pathway-tools model.
pathway_ontology_expanded(Pathway, Type) :- pathway(Pathway) , is_a(Pathway, Type).
pathway_ontology_expanded(Pathway, ParentType) :- is_a(Type, ParentType) , is_a(Pathway, Type).
% achieved: pathway_ontology_expanded contains the ontology flattened.

rule(1, "REJECT", Pathway) :- pathway_ontology_expanded(Pathway, Class) , ignored_pathway_class(Class).

%
% Test data for rule 1:
pathway("PWY-TEST1").
is_a("PWY-TEST1", "My Specific Signal 1").
is_a("My Specific Signal 1", "My generic type 1 Signal").
is_a("My generic type 1 Signal", "Signaling-Pathways").

pathway("PWY-TEST2").
is_a("PWY-TEST1", "My biosynthesis pathway").
is_a("My biosynthesis pathway", "My generic biosynthesis").
is_a("My generic biosynthesis", "Biosynthesis-Pathway").

#show rule/3.
