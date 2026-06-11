from clyngor import solve


def assert_pass(rule, pathway, answers):
    assert tuple([rule, '"PASS"', f'"{pathway}"']) in answers[0]["rule"]


def assert_reject(rule, pathway, answers):
    assert tuple([rule, '"REJECT"', f'"{pathway}"']) in answers[0]["rule"]


def assert_include(rule, pathway, answers):
    assert tuple([rule, '"INCLUDE"', f'"{pathway}"']) in answers[0]["rule"]


def test_asp_pathologic_rule_1():
    """
    Test rule 1 of pathologic like algorithm implemented in ASP
    """
    asp_code = """
    pathway("PWY-TEST1").
    is_a("PWY-TEST1", "My Specific Signal 1").
    is_a("My Specific Signal 1", "My generic type 1 Signal").
    is_a("My generic type 1 Signal", "Signaling-Pathways").

    pathway("PWY-TEST2").
    is_a("PWY-TEST1", "My biosynthesis pathway").
    is_a("My biosynthesis pathway", "My generic biosynthesis").
    is_a("My generic biosynthesis", "Biosynthesis-Pathway").
    """
    answers = tuple(
        solve(
            ["./src/pan2met/asp/rules/pathologic_like.lp"], inline=asp_code
        ).by_predicate
    )
    assert len(answers) == 1
    assert_reject(1, "PWY-TEST1", answers)
    assert_pass(1, "PWY-TEST2", answers)
    # assert len(answers[0]["rule"]) == 2


def test_asp_pathologic_rule_2():
    """
    Test rule 2 of pathologic like algorithm implemented in ASP, on electron transfer pathway with missing enzyme.
    """
    asp_code = """
    % Test data for rule 2:
    % First case is a pathway that is not a electron transfer pathway.
    pathway("PWY-TEST1").
    is_a("Type1", "Catalysis-Pathways").
    is_a("PWY-TEST1", "Type1").

    % Second case is a pathway which is an electron transfer pathway,
    % but have all its non orphan, non spontaneous reactions are catalyzed
    % by an enzyme known to be present in the genome, and hence the reaction
    % is expected to be present in the reactome.
    pathway("PWY-TEST2").
    is_a("Type2", "Electron-Transfer").
    is_a("PWY-TEST2", "Type2").
    is_in_pathway("RXN-1", "PWY-TEST2").
    is_in_pathway("RXN-2", "PWY-TEST2").
    is_in_pathway("RXN-3", "PWY-TEST2").
    is_in_pathway("RXN-4", "PWY-TEST2").
    is_in_pathway("RXN-5", "PWY-TEST2").

    reactome("RXN-1").
    reactome("RXN-2").
    reactome("RXN-5").
    orphan("RXN-3").
    spontaneous("RXN-4").

    % Third case is a pathway which is an electron transfer pathway,
    % and has a reaction that has no enzyme catalyzing it in the genome
    % even if it is not orphan.
    pathway("PWY-TEST3").
    is_a("PWY-TEST3", "Type2").
    is_in_pathway("RXN-UNKNOWN", "PWY-TEST3").
    """
    answers = tuple(
        solve(
            ["./src/pan2met/asp/rules/pathologic_like.lp"], inline=asp_code
        ).by_predicate
    )
    assert len(answers) == 1
    assert_pass(2, "PWY-TEST1", answers)
    assert_pass(2, "PWY-TEST2", answers)
    assert_reject(2, "PWY-TEST3", answers)
    # assert len(answers[0]["rule"]) == 3


def test_asp_pathologic_rule_4():
    """
    Test rule 4 of pathologic like algorithm implemented in ASP: when all catalyzed reactions' enzymes are represented, include the pathway when in taxonomic range,
    or when not in taxonomic range, when at least 4 reactions of the pathway are catalyzed.
    """
    asp_code = """
    % Test cases for rule 4:
    reactome("RXN-1").
    reactome("RXN-2").
    reactome("RXN-3").
    reactome("RXN-4").

    % Case 1. A pathway with no missing enzyme, and in taxonomic range
    pathway("PWY-1").
    in_taxonomic_range("PWY-1").
    is_in_pathway("RXN-1", "PWY-1").
    is_in_pathway("RXN-2", "PWY-1").
    is_in_pathway("RXN-3", "PWY-1").

    % Case 2. A pathway with no missing enzyme, outside of taxonomic range, but with more than 3 reactions with an enzyme.
    pathway("PWY-2").
    % not specified atoms are false by default: in_taxonomic_range("PWY-2").
    is_in_pathway("RXN-1", "PWY-2").
    is_in_pathway("RXN-2", "PWY-2").
    is_in_pathway("RXN-3", "PWY-2").
    is_in_pathway("RXN-4", "PWY-2").

    % Case 3. A pathway with no missing enzyme, outside of taxonomic range, but with at most 3 reactions with an enzyme.
    pathway("PWY-3").
    % not specified atoms are false by default: in_taxonomic_range("PWY-2").
    is_in_pathway("RXN-1", "PWY-1").
    is_in_pathway("RXN-2", "PWY-2").
    is_in_pathway("RXN-3", "PWY-3").

    % Case 4. A pathway with missing enzyme, even if in taxonomic range
    pathway("PWY-4").
    in_taxonomic_range("PWY-4").
    is_in_pathway("MISSING-RXN-1", "PWY-4").
    is_in_pathway("RXN-1", "PWY-4").
    is_in_pathway("RXN-2", "PWY-4").

    #show rule/3.
    """
    answers = tuple(
        solve(
            ["./src/pan2met/asp/rules/pathologic_like.lp"], inline=asp_code
        ).by_predicate
    )
    assert len(answers) == 1
    assert_include(4, "PWY-1", answers)
    assert_include(4, "PWY-2", answers)
    assert_pass(4, "PWY-3", answers)
    assert_pass(4, "PWY-4", answers)
    # assert len(answers[0]["rule"]) == 3


def test_asp_pathologic_rule_8():
    """
    Test rule 8 of pathologic like algorithm implemented in ASP:

    when a key reaction is missing: reject the pathway.
    """
    asp_code = """
     % Test cases for rule 8:
    reactome("RXN-1"; "RXN-2"; "RXN-3").
    % Case 1:
    % pathway has a missing key reaction
    pathway("PWY-1").
    is_in_pathway("PWY-1", "MISSING-RXN-1").
    is_in_pathway("PWY-1", "RXN-2").
    pathway_key_reaction("PWY-1", "MISSING-RXN-1").

    % Case 2:
    % pathway has no missing key reaction
    pathway("PWY-2").
    is_in_pathway("RXN-1", "PWY-2").
    is_in_pathway("RXN-2", "PWY-2").
    is_in_pathway("RXN-3", "PWY-2").
    is_in_pathway("RXN-4", "PWY-2").
    pathway_key_reaction("PWY-2", "RXN-1").
   """
    answers = tuple(
        solve(
            ["./src/pan2met/asp/rules/pathologic_like.lp"], inline=asp_code
        ).by_predicate
    )
    assert len(answers) == 1
    assert_reject(8, "PWY-1", answers)
    assert_pass(8, "PWY-2", answers)
    # assert len(answers[0]["rule"]) == 3


def test_asp_pathologic_rule_9():
    """
    Test rule 8 of pathologic like algorithm implemented in ASP:

    when a biosynthesis pathway is missing one of its last.
    """
    asp_code = """
        % Test cases for rule 9:
    % Case 1.
    % Test the extraction of the first reaction from the order of the reaction.
    pathway("PWY-1").

    % expect first_reaction("PWY-1", "RXN-1") and last_reaction("PWY-1", "RXN-3").
    reactome("RXN-1"; "RXN-2"; "RXN-3").
    % Case 2.
    % a biosynthesis pathway with a missing last reaction
    pathway("PWY-2").
    pathway_reaction_order("PWY-2", "RXN-1", 1).
    pathway_reaction_order("PWY-2", "RXN-2", 2).
    pathway_reaction_order("PWY-2", "RXN-3", 3).
    pathway_reaction_order("PWY-2", "MISSING-RXN-4", 4).

    is_in_pathway("RXN-1", "PWY-2").
    is_in_pathway("RXN-2", "PWY-2").
    is_in_pathway("RXN-3", "PWY-2").
    is_in_pathway("MISSING-RXN-4", "PWY-2").
    is_a("PWY-2", "Biosynthesis").
    % expect rule 9, alternative 1 to reject the pathway.

    % Case 3.
    % a catabolysis pathway with a missing first reaction
    pathway("PWY-3").

    pathway_reaction_order("PWY-3", "MISSING-RXN-0", 0).
    pathway_reaction_order("PWY-3", "RXN-1", 1).
    pathway_reaction_order("PWY-3", "RXN-2", 2).
    pathway_reaction_order("PWY-3", "RXN-3", 3).

    is_in_pathway("MISSING-RXN-0", "PWY-3").
    is_in_pathway("RXN-1", "PWY-3").
    is_in_pathway("RXN-2", "PWY-3").
    is_in_pathway("RXN-3", "PWY-3").
    is_a("PWY-3", "Degradation").
    % expect rule 9, alternative 2, to reject the pathway.

    % Case 4
    % a pathway with a missing first reaction, but not a catabolysis
    pathway("PWY-4").

    pathway_reaction_order("PWY-4", "MISSING-RXN-0", 0).
    pathway_reaction_order("PWY-4", "RXN-1", 1).
    pathway_reaction_order("PWY-4", "RXN-2", 2).
    pathway_reaction_order("PWY-4", "RXN-3", 3).

    is_in_pathway("MISSING-RXN-0", "PWY-4").
    is_in_pathway("RXN-1", "PWY-4").
    is_in_pathway("RXN-2", "PWY-4").
    is_in_pathway("RXN-3", "PWY-4").
    % expect rule 9 not to reject the pathway.

    #show rule/3.
    """
    answers = tuple(
        solve(
            ["./src/pan2met/asp/rules/pathologic_like.lp"], inline=asp_code
        ).by_predicate
    )
    assert len(answers) == 1
    print(answers)
    assert_reject(9, "PWY-2", answers)
    assert_reject(9, "PWY-3", answers)
    assert_pass(9, "PWY-4", answers)


def test_rule_dependency():
    """
    Test that the rules are applied in the expected order, and that the decision of one rule is available for the next rules.
    For instance, test that rule 9 can use the result of rule 8 to decide whether a reaction is missing or not.
     For instance, in the case of a pathway with a key reaction that is missing, test that rule 9 can use the result of rule 8 to reject the pathway because of the missing key reaction, even if it is missing a first or last reaction.
     In this case, we expect rule 8 to reject the pathway because of the missing key reaction, and then we expect rule 9 to pass the pathway
    """
    asp_code = """
    rule(8, "REJECT", "PWY-TEST").
    rule(7, "INCLUDE", "PWY-TEST").
    """
    models = solve(["./src/pan2met/asp/rules/pathologic_like.lp"], inline=asp_code)
    _model = next(
        models, None
    )  # Important, to trigger the evaluation of parsing of the clingo results.
    assert models.is_unsatisfiable
    assert not models.is_unknown


def test_rule_contradiction():
    """
    Test that no rule return contradictory result for the same pathway.
    """
    asp_code = """
    rule(1, "REJECT", "PWY-TEST").
    rule(1, "INCLUDE", "PWY-TEST").
    """
    models = solve(
        ["./src/pan2met/asp/rules/pathologic_like.lp"], inline=asp_code, stats=False
    )
    _model = next(models, None)
    assert models.is_unsatisfiable
    assert not models.is_unknown
