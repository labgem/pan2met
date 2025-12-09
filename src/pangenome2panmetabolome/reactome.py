#!/usr/bin/env python3

from typing import Iterable

from clyngor import solve

from .asp import monomer_asp_rule


def infer_reactome(monomers: list[str], inference_rules_path: str) -> Iterable[str]:
    """
    Infer the reactome using Answer Set Programming

    Given a list of 'seed' monomer id,
    infer the list of realized reaction ids.

    Arguments
    ---------

        monomers: list of monomer id
        inference_rules_path: Path to a AnsProlog file (i.e., .lp) with reaction
        inference rules built from the knowledge base

    Returns
    -------

        The list of reaction identifers (e.g., "RXN-1")

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
    answer = next(answers)
    for predicate, value in answer:
        if predicate == "reaction":
            reaction = value[0]
            reaction = reaction.replace('"', "")
            yield reaction
