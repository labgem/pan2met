"""
Reactome inference

Simple inference rule: if a reaction has an enzyme that can catalyze it in an organism,
 simply infer the presence of the reaction in the reactome.

"""

from typing import Iterable
import os

from clyngor import solve
import pythoncyc

from .asp.asp import monomer_asp_rule
from .knowledge_base import KnowledgeBase
from .io import metacyc  # TODO: enable more source of knowledge.
from .utils import logger


def infer_complex_from_monomers(
    monomers: set[str], complex: str, pgdb: pythoncyc.PGDB
) -> bool:
    """
    True if the given complex can be formed by the given set of protein monomers.
    """
    components = metacyc.proteic_complex_subunits(pgdb, complex)
    if len(components) == 0:
        logger.error(f"{complex} complex has no components")
        return False
    for component in components:
        if component not in complex:
            return False
    return True


def infer_complexes_from_monomers(monomers: set[str], pgdb: pythoncyc.PGDB) -> set[str]:
    complexes: set[str] = set()
    for complex in pgdb.all_complexes():
        if infer_complex_from_monomers(monomers, complex, pgdb):
            complexes.add(complex)
    return complexes


def infer_reactome_from_monomers(monomers: set[str], pgdb: pythoncyc.PGDB) -> set[str]:
    """
    Naive inference of a set of reaction.

    Arguments
    ---------

        monomers -- list of monomer identifiers
        pgdb -- PythonCyc PGDB adapter

    Yields
    ------

        reaction identifiers
    """

    # Start by infering all reachable complex
    complexes: set[str] = infer_complexes_from_monomers()
    # Continue, by infering the possible reactions
    reactions: set[str] = set()
    for reaction in pgdb.all_rxns():
        for enzyme in pgdb.enzymes_of_reaction(pgdb, reaction):
            if metacyc.is_proteic_complex(pgdb, enzyme):
                if enzyme in complexes:
                    reactions.add(reaction)
            elif enzyme in monomers:
                reactions.add(reaction)
    return reactions


def infer_reactome_from_monomers_asp(
    monomers: list[str], inference_rules_path: str
) -> Iterable[str]:
    """
    Infer the reactome using Answer Set Programming

    Given a list of 'seed' monomer id,
    infer the list of realized reaction ids.

    Arguments
    ---------

        :monomers: list of monomer id
        :inference_rules_path: Path to a AnsProlog file (i.e., .lp) with reaction inference rules built from the knowledge base

    Yields
    -------

        reaction identifers (e.g., "RXN-1")

    Format of the infered reaction atoms
    ------------------------------------

    This function expects atoms identified with AnsProlog atoms in answer set such as

    .. code:: prolog

      reaction("RXN-1").

    for reaction identifier "RXN-1", when such a reaction is infered to be present in the reactome.

    """
    SHOW_REACTION_DIRECTIVE = "#show reaction/1."
    monomer_asp_rules = "\n".join(map(monomer_asp_rule, monomers))
    monomer_asp_rules += "\n" + SHOW_REACTION_DIRECTIVE
    answers = solve(
        inference_rules_path, inline=monomer_asp_rules, use_clingo_module=False
    )
    answer = next(answers)  # Take the first answer of the clingo output.
    for predicate, value in answer:
        if predicate == "reaction":
            reaction = value[0]
            reaction = reaction.replace('"', "")
            yield reaction  # For all predicate reaction("RXN-1"), yield RXN-1


def infer_reactome_from_ec_numbers(
    ec_numbers: list[str], kb: KnowledgeBase
) -> list[str]:
    reaction_set: set[str] = set()
    for ec_number in ec_numbers:
        for reaction in kb.reactions_by_ec_number(ec_number):
            reaction_set.add(reaction)
    return list(reaction_set)


def minimal_monomer_set(
    reactions: list[str],
    potential_monomer_inference_rule_path: str,
    reaction_inference_rule_path: str,
) -> set[str]:
    """
    Use ASP to identify a minimal set of monomer that is expected to be sufficient to catalyze a set of reactions.

    Arguments
    ---------

    :reactions: a list of reaction identifiers
    :potential_monomer_inference_rule_path: a path to 'potential' involved monomer inference rules
    :reaction_inference_rule_path: a path to inference rules from monomer (to complex) to reaction

    Returns
    -------

    A 'minimal' set of monomer id sufficient to catalyze the given set of reactions
    """

    minimal_set_asp_rule_path = os.path.join(
        os.path.dirname(__file__), "../asp/required_monomer_given_reactions.lp"
    )

    reaction_asp_atoms = "\n".join(
        [f'reaction("{reaction}").' for reaction in reactions]
    )
    target_reaction_asp_atoms = "\n".join(
        [f'target_reaction("{reaction}").' for reaction in reactions]
    )

    print(target_reaction_asp_atoms)
    # First, identify the subset of the whole set of monomer that may be involved in the selected reactions,
    # using the inverse inference rules
    # The output is a set of atom potential_monomer/1.
    answers = solve(potential_monomer_inference_rule_path, inline=reaction_asp_atoms)
    answer = next(answers)
    potential_monomers: set[str] = set()
    for predicate, value in answer:
        if predicate == "potential_monomer":
            identifier = value[0]
            identifier = identifier.replace('"', "")
            potential_monomers.add(identifier)
    logger.debug(potential_monomers)
    potential_monomer_asp_atoms = "\n".join(
        [f'potential_monomer("{monomer}").' for monomer in potential_monomers]
    )

    print(potential_monomer_asp_atoms)
    # Then, find a minimal subset of potential monomer "selected_monomer/1" that satisfies the set of "target_reactions/1"
    answers = solve(
        [minimal_set_asp_rule_path, reaction_inference_rule_path],
        inline=potential_monomer_asp_atoms + "\n" + target_reaction_asp_atoms,
    )

    # Take the first answer set
    selected_monomers: set[str] = set()
    answer = next(answers)
    for predicate, value in answer:
        if predicate == "selected_monomer":
            identifier = value[0]
            identifier = identifier.replace('"', "")
            selected_monomers.add(identifier)
    return selected_monomers
