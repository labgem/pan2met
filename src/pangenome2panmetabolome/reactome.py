"""
Reactome inference

Simple inference rule: if a reaction has an enzyme that can catalyze it in an organism,
 simply infer the presence of the reaction in the reactome.

"""

from typing import Iterable

from clyngor import solve

from .asp import monomer_asp_rule
from .knowledge_base import KnowledgeBase


def infer_reactome_from_monomers(
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
