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
% - pathway_reaction_order(Pathway, PredecessorReaction, SuccessorReaction).
%   PredecessorReaction appears before Successor reaction in the pathway reaction graph.
%   Used for a special rule for catabolysis and biosynthesis pathways.
% - pathway_passing_score(Pathway).
% - is_a(Pathway, Type) or is_a(Type, ParentType), define the pathway class ontology


% In the PathoLogic algorithm, the decision rules are run sequentially.
% To mimick such a procedure in clingo Answer Set Programming, we will use
% a little trick, we will add a dependency for precedent rules in subsequent rules.

% If a rule rejected a pathway, subsequent rules cannot  it.
% We ignore the "PASS" decision.
:- rule(RuleNumber, "INCLUDE", Pathway)
  , rule(PreviousRuleNumber, "REJECT", Pathway)
  , PreviousRuleNumber < RuleNumber.

% Similarly, if a rule included a pathway, subsequent rules cannot reject it afterwards.
:- rule(RuleNumber, "REJECT", Pathway)
    , rule(PreviousRuleNumber, "INCLUDE", Pathway)
    , PreviousRuleNumber < RuleNumber.

% Also, we can only have one decision for a pathway, the following is a simple sanity check:
Decision1 = Decision2 :- rule(RuleNumber, Decision1, Pathway), rule(RuleNumber, Decision2, Pathway).

% We can identify rejected pathway and accepted pathway by looking at each of the decisions
rejected_pathway(Pathway) :- rule(_, "REJECT", Pathway).
included_pathway(Pathway) :- rule(_, "INCLUDE", Pathway).
undecided_pathway(Pathway) :- pathway(Pathway) , not rejected_pathway(Pathway) , not included_pathway(Pathway) , reactome(Reaction) , is_in_pathway(Reaction, Pathway).
% undecided pathways are those neither rejected nor included, but having at least one reaction with evidence of presence in the reactome.
% undecided pathways will be subject to a decision rule based on the pathway score (last rule, not implemented here, but in python).


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
pathway_catalysis(Pathway, Reaction) :- is_in_pathway(Reaction, Pathway) , not orphan(Reaction) , not spontaneous(Reaction).
missing_enzyme_for_reaction(Pathway, Reaction) :- pathway_catalysis(Pathway, Reaction) , not reactome(Reaction).
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

% Third rule:
% REJECT P if any key-non-reactions are present.
% Not handled.
% rule(3, "PASS",  Pathway) :- pathway(Pathway).

% Fourth rule:
% INCLUDE P if P has all reactions present (meaning an enzyme is present for each reaction)
% AND if P is outside its taxonomic range, P contains more than 3 reactions. Rationale: We
% keep pathways outside their designated taxonomic range if very strong evidence for the
% pathway exists.

pathway_catalysis_count(Pathway, N) :- N = #count{ Enzyme: pathway_catalysis(Pathway, Enzyme) }, pathway(Pathway).
% achieved: N is the number of reaction that is expected to require an enzyme to be catalyzed (neither spontaneous, nor orphan, but can be absent from the reactome).

rule(4, "INCLUDE", Pathway) :- pathway(Pathway) , not missing_enzyme(Pathway) , in_taxonomic_range(Pathway).
rule(4, "INCLUDE", Pathway) :- pathway(Pathway) , not missing_enzyme(Pathway) , pathway_catalysis_count(Pathway, N) , N > 3.
rule(4, "PASS", Pathway) :- pathway(Pathway) , not rule(4, "INCLUDE", Pathway).


% % Test cases for rule 4:
% reactome("RXN-1").
% reactome("RXN-2").
% reactome("RXN-3").
% reactome("RXN-4").

% % Case 1. A pathway with no missing enzyme, and in taxonomic range
% pathway("PWY-1").
% in_taxonomic_range("PWY-1").
% is_in_pathway("RXN-1", "PWY-1").
% is_in_pathway("RXN-2", "PWY-1").
% is_in_pathway("RXN-3", "PWY-1").

% % Case 2. A pathway with no missing enzyme, outside of taxonomic range, but with more than 3 reactions with an enzyme.
% pathway("PWY-2").
% % not specified atoms are false by default: in_taxonomic_range("PWY-2").
% is_in_pathway("RXN-1", "PWY-2").
% is_in_pathway("RXN-2", "PWY-2").
% is_in_pathway("RXN-3", "PWY-2").
% is_in_pathway("RXN-4", "PWY-2").

% % Case 3. A pathway with no missing enzyme, outside of taxonomic range, but with at most 3 reactions with an enzyme.
% pathway("PWY-3").
% % not specified atoms are false by default: in_taxonomic_range("PWY-2").
% is_in_pathway("RXN-1", "PWY-1").
% is_in_pathway("RXN-2", "PWY-2").
% is_in_pathway("RXN-3", "PWY-3").

% % Case 4. A pathway with missing enzyme, even if in taxonomic range
% pathway("PWY-4").
% in_taxonomic_range("PWY-4").
% is_in_pathway("MISSING-RXN-1", "PWY-4").
% is_in_pathway("RXN-1", "PWY-4").
% is_in_pathway("RXN-2", "PWY-4").

% #show rule/3.

% Fifth rule:
% REJECT P if the score of P is significantly less than the score of a variant pathway of P.
% TODO:

% Sixth rule:
% REJECT P if the score of P is slightly more than a preferred variant pathway of P
% (the “standard” glycolysis and TCA pathways are designated as preferred variants).
% TODO;

% 7th rule:
% REJECT P if P is outside its taxonomic range.
% (applicable, when all above rules are passing, especially the fourth rule. So, we consider this rule to be equivalent to REJECT P if P is outside its taxonomic range and missing an enzyme.)
rule(7, "REJECT", Pathway) :- pathway(Pathway) , missing_pathway(Pathway) , not in_taxonomic_range(Pathway).

% 8th rule:
% REJECT P if P is missing enzymes for any key reaction of P.
has_missing_key_reaction(Pathway) :- not reactome(Reaction) , pathway_key_reaction(Pathway, Reaction).
rule(8, "REJECT", Pathway) :- pathway(Pathway) , has_missing_key_reaction(Pathway).
rule(8, "PASS", Pathway) :- pathway(Pathway) , not rule(8, "REJECT", Pathway).
% % Test cases for rule 8:
% reactome("RXN-1"; "RXN-2"; "RXN-3").
% % Case 1:
% % pathway has a missing key reaction
% pathway("PWY-1").
% is_in_pathway("PWY-1", "MISSING-RXN-1").
% is_in_pathway("PWY-1", "RXN-2").
% pathway_key_reaction("PWY-1", "MISSING-RXN-1").

% %Case 2:
% % pathway has no missing key reaction
% pathway("PWY-2").
% is_in_pathway("RXN-1", "PWY-2").
% is_in_pathway("RXN-2", "PWY-2").
% is_in_pathway("RXN-3", "PWY-2").
% is_in_pathway("RXN-4", "PWY-2").
% pathway_key_reaction("PWY-2", "RXN-1").


% 9th rule:
% REJECT P if P is either:
% – A biosynthetic pathway missing enzymes for its final steps
% – A catabolic pathway missing enzymes for its initial steps
% – An energy pathway missing enzymes for more than half of its steps
is_biosynthesis_pathway(Pathway) :- pathway_ontology_expanded(Pathway, "Biosynthesis").
is_catabolysis_pathway(Pathway) :- pathway_ontology_expanded(Pathway, "Degradation").

first_reaction(Pathway, FirstReaction) :- pathway_reaction_order(Pathway, FirstReaction, _) , is_in_pathway(OtherReaction, Pathway) , not pathway_reaction_order(Pathway, OtherReaction, FirstReaction).
% A first reaction is a reaction with no predecessor
last_reaction(Pathway, FirstReaction) :- pathway_reaction_order(Pathway, _, LastReaction) , is_in_pathway(OtherReaction, Pathway) , not pathway_reaction_order(Pathway, LastReaction, OtherReaction).
% A last reaction is a reaction with no successor
% FIXME: check if a pathway could have multiple extremities, in which case, identifying the reaction of interest for anabolism or catabolism could be not that easy.

%!trace_rule { "% pathway is a biosynthesis pathway and missing the last reaction %.", Pathway, Reaction }
rule(9, "REJECT", Pathway) :-
          pathway(Pathway)
          , last_reaction(Pathway, Reaction)
          , missing_enzyme_for_reaction(Pathway, Reaction)
          , is_biosynthesis_pathway(Pathway).

%!trace_rule { "% pathway is a catabolysis pathway and missing the first reaction %.", Pathway, Reaction }
rule(9, "REJECT", Pathway) :-
          pathway(Pathway)
          , first_reaction(Pathway, Reaction)
          , missing_enzyme_for_reaction(Pathway, Reaction)
          , is_catabolysis_pathway(Pathway).

count_pathway_reactions(Pathway, N) :- N = #count{ Reaction: is_in_pathway(Reaction, Pathway) }, pathway(Pathway).
count_missing_enzyme_reactions(Pathway, N) :- N = #count{ Reaction: missing_enzyme_for_reaction(Pathway, Reaction) }, pathway(Pathway).
rule(9, "REJECT", Pathway) :-
          pathway(Pathway)
            , count_pathway_reactions(Pathway, TotalReactions)
            , count_missing_enzyme_reactions(Pathway, MissingEnzymeReactions)
            , MissingEnzymeReactions > TotalReactions / 2
            , pathway_ontology_expanded(Pathway, "Energy-Metabolism").

rule(9, "PASS", Pathway) :- pathway(Pathway) , not rule(9, "REJECT", Pathway).


% % Test cases for rule 9:
% % Case 1.
% % Test the extraction of the first reaction from the order of the reaction.
% pathway("PWY-1").

% % expect first_reaction("PWY-1", "RXN-1") and last_reaction("PWY-1", "RXN-3").
% reactome("RXN-1"; "RXN-2"; "RXN-3").
% % Case 2.
% % a biosynthesis pathway with a missing last reaction
% pathway("PWY-2").
% pathway_reaction_order("PWY-2", "MISSING-RXN-0", 0).
% pathway_reaction_order("PWY-2", "RXN-1", 1).
% pathway_reaction_order("PWY-2", "RXN-2", 2).
% pathway_reaction_order("PWY-2", "RXN-3", 3).
% pathway_reaction_order("PWY-2", "MISSING-RXN-4", 4).

% is_in_pathway("RXN-1", "PWY-2").
% is_in_pathway("RXN-2", "PWY-2").
% is_in_pathway("RXN-3", "PWY-2").
% is_in_pathway("MISSING-RXN-4", "PWY-2").
% is_a("PWY-2", "Biosynthesis").
% % expect rule 9, alternative 1 to reject the pathway.

% % Case 3.
% % a catabolysis pathway with a missing first reaction
% pathway("PWY-3").

% pathway_reaction_order("PWY-3", "MISSING-RXN-0", 0).
% pathway_reaction_order("PWY-3", "RXN-1", 1).
% pathway_reaction_order("PWY-3", "RXN-2", 2).
% pathway_reaction_order("PWY-3", "RXN-3", 3).

% is_in_pathway("MISSING-RXN-0", "PWY-3").
% is_in_pathway("RXN-1", "PWY-3").
% is_in_pathway("RXN-2", "PWY-3").
% is_in_pathway("RXN-3", "PWY-3").
% is_a("PWY-3", "Degradation").
% % expect rule 9, alternative 2, to reject the pathway.

% % Case 4.
% % a pathway with a missing first reaction, but not a catalysis
% pathway("PWY-4").

% pathway_reaction_order("PWY-4", "MISSING-RXN-0", 0).
% pathway_reaction_order("PWY-4", "RXN-1", 1).
% pathway_reaction_order("PWY-4", "RXN-2", 2).
% pathway_reaction_order("PWY-4", "RXN-3", 3).

% is_in_pathway("MISSING-RXN-0", "PWY-4").
% is_in_pathway("RXN-1", "PWY-4").
% is_in_pathway("RXN-2", "PWY-4").
% is_in_pathway("RXN-3", "PWY-4").
% % expect rule 9 not to reject the pathway.


% 10th rule:
% INCLUDE P if the score of P exceeds the threshold defined by the ptools-init.dat
% parameter PATHWAY-PREDICTION-SCORE-CUTOFF. The value of this parameter can also
% be specified in the graphical user interface for PathoLogic.
% This rule should most probably better be implanted in Python, as it requires computation on floating point numbers.
