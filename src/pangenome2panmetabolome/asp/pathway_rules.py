"""
Generate the AnsProlog rules for the inference of pathway presence based solely on the presence of reactions.

To do so, the chosen heuristic is both simple and naive:
 If all reactions in a pathway is present in a reactome,
 then, the pathway is infered to be present in the metabolome.

Hence, no consideration is done, for the moment for pathways that are incomplete.
If a reaction is "orphan", that is to say, if no known enzyme catalyzes it,
this approach will simply ignore the pathway, arguing that we cannot assure the presence of the reaction in the organism.
On the other hand, if the reaction is known to be "spontaneous", we assume it can occur in any organism, and is assume to be present in the reactome.
Thus, this script will ignore this reaction in the inference rule: it will not appear in the rule.

To launch this script, you will need to launch the pathway-tools python API with:

.. code :: bash

  pathway-tools -lisp -python-local-only-non-strict
"""

import argparse
from typing import Iterable

import pythoncyc

from ..utils import logger
from ..io.metacyc import (
    get_reactions_of_pathway,
    is_spontaneous,
    is_orphan,
    remove_pipes,
)


def pathway_asp_rule(pathway: str, reactions: Iterable[str]) -> str:
    return (
        f'pathway("{pathway}") :- '
        + " , ".join(f'reaction("{reaction}")' for reaction in reactions)
        + "."
    )


def pathway_asp_generator(pgdb, ignore_orphan: bool = False) -> Iterable[str]:
    """
    Generate Answer Set Programming rules for pathway inference rules
    from a BioCyc PGDB (MetaCyc most probably).

    Each pathway will be represented by a AnsProlog rule as follows,
    given that the pathway "PWY-1" contains the reactions RXN-1, RXN-2 and RXN-3:


    .. code :: prolog

      pathway("PWY-1") :- reaction("RXN-1") , reaction("RXN-2") , reaction("RXN-3").

    Arguments
    ---------

        pgdb: the pythoncyc instance (e.g., pythoncyc.select_organism("meta"))

    Yields
    ------

        ASP rules

    """
    for pathway in pgdb.all_pathways():
        required_reactions = []
        for reaction in get_reactions_of_pathway(pgdb, pathway):
            if not is_spontaneous(pgdb, reaction) and (
                ignore_orphan or not is_orphan(pgdb, reaction)
            ):
                required_reactions.append(reaction)
        yield pathway_asp_rule(
            remove_pipes(pathway), map(remove_pipes, required_reactions)
        )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate basic pathway inference ASP rules"
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to output AnsProlog .lp file where ASP rules will be written",
        required=True,
    )
    parser.add_argument(
        "--ignore-orphan",
        action=argparse.BooleanOptionalAction,
        help="If set to --ignore-orphan, orphan reactions will not appear in any pathway ASP rules, thus being considered not required to infer the reaction as present in the metabolism.",
    )
    return parser.parse_args()


def main():
    logger.setLevel("DEBUG")
    args = parse_arguments()
    pgdb = pythoncyc.select_organism("meta")
    with open(args.output, "w") as asp_file:
        for rule in pathway_asp_generator(pgdb):
            asp_file.write(rule + "\n")


if __name__ == "__main__":
    main()
