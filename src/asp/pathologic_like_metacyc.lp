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
ignored_pathway_class("Transport-Pathways" ; "Signaling-Pathways"). % Synthetic pathways seems not to be defined as an ontology parent class in pathway-tools model.
pathway_ontology_expanded(Pathway, Type) :- pathway(Pathway) , is_a(Pathway, Type).
pathway_ontology_expanded(Pathway, ParentType) :- is_a(Type, ParentType) , pathway_ontology_expanded(Pathway, Type).
% achieved: pathway_ontology_expanded contains the ontology flattened.
%!trace_rule {"% is a % which is not handled", Pathway, Class}
rule(1, "REJECT", Pathway) :- pathway_ontology_expanded(Pathway, Class) , ignored_pathway_class(Class) , pathway(Pathway).
rule(1, "PASS", Pathway) :- not rule(1, "REJECT", Pathway) , pathway(Pathway).

% %
% % Test data for rule 1:
% pathway("PWY-TEST1").
% is_a("PWY-TEST1", "My Specific Signal 1").
% is_a("My Specific Signal 1", "My generic type 1 Signal").
% is_a("My generic type 1 Signal", "Signaling-Pathways").

% pathway("PWY-TEST2").
% is_a("PWY-TEST1", "My biosynthesis pathway").
% is_a("My biosynthesis pathway", "My generic biosynthesis").
% is_a("My generic biosynthesis", "Biosynthesis-Pathway").

% Second rule:
% REJECT P if P is an electron transport pathway AND P lacks enzymes for any reaction.
% Rationale: Electron transport pathways are very short and specific, therefore all enzymes
% must be present to include the pathway.
missing_enzyme_for_reaction(Pathway, Reaction) :- is_in_pathway(Reaction, Pathway) , not reactome(Reaction) , not spontaneous(Reaction) , not orphan(Reaction).
missing_enzyme(Pathway) :- missing_enzyme_for_reaction(Pathway, _).
%!trace_rule {"% is an electron transport pathway with a missing enzyme", Pathway}
rule(2, "REJECT", Pathway) :- pathway_ontology_expanded(Pathway, "Electron-Transfer") , missing_enzyme(Pathway).
rule(2, "PASS", Pathway) :- not rule(2, "REJECT", Pathway) , pathway(Pathway).

% % Test data for rule 2:
% % First case is a pathway that is not a electron transfer pathway.
% pathway("PWY-TEST1").
% is_a("Type1", "Catalysis-Pathways").
% is_a("PWY-TEST1", "Type1").

% Second case is a pathway which is an electron transfer pathway,
% % but have all its non orphan, non spontaneous reactions are catalyzed
% % by an enzyme known to be present in the genome, and hence the reaction
% % is expected to be present in the reactome.
% pathway("PWY-TEST2").
% is_a("Type2", "Electron-Transfer").
% is_a("PWY-TEST2", "Type2").
% is_in_pathway("RXN-1", "PWY-TEST2").
% is_in_pathway("RXN-2", "PWY-TEST2").
% is_in_pathway("RXN-3", "PWY-TEST2").
% is_in_pathway("RXN-4", "PWY-TEST2").
% is_in_pathway("RXN-5", "PWY-TEST2").

% reactome("RXN-1").
% reactome("RXN-2").
% reactome("RXN-5").
% orphan("RXN-3").
% spontaneous("RXN-4").

% % Third case is a pathway which is an electron transfer pathway,
% % and has a reaction that has no enzyme catalyzing it in the genome
% % even if it is not orphan.
% pathway("PWY-TEST3").
% is_a("PWY-TEST3", "Type2").
% is_in_pathway("RXN-UNKNOWN", "PWY-TEST3").

% #show rule/3.
