"""
Reactome inference

Simple inference rule: if a reaction has an enzyme that can catalyze it in an organism,
 simply infer the presence of the reaction in the reactome.

"""

from typing import Iterable

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
