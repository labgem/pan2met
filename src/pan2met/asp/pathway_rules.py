"""
Generate the AnsProlog rules for the inference of pathway presence based solely on the presence of reactions.

To do so, the chosen heuristic is both simple and naive:
 If all reactions in a pathway is present in a reactome,
 then, the pathway is infered to be present in the metabolism.

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
from collections.abc import Iterable

from pan2met.config import default_config, override_config
from pan2met.io.knowledge_base import KnowledgeBase, select_kb
from pan2met.utils import logger


def pathway_asp_rule(pathway: str, reactions: Iterable[str]) -> str:
    return (
        f'pathway("{pathway}") :- '
        + " , ".join(f'reaction("{reaction}")' for reaction in reactions)
        + "."
    )


def pathway_asp_generator(
    kb: KnowledgeBase, ignore_orphan: bool = False
) -> Iterable[str]:
    """
    Generate Answer Set Programming rules for pathway inference rules

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
    for pathway in kb.pathways():
        if ignore_orphan:
            required_reactions = kb.non_spontaneous_reactions_of_pathway(pathway)
        else:
            required_reactions = kb.non_orphan_non_spontaneous_reactions_of_pathway(
                pathway
            )
        yield pathway_asp_rule(pathway, required_reactions)


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
        default=False,
    )
    return parser.parse_args()


def main():
    logger.setLevel("DEBUG")
    args = parse_arguments()
    if args.config:
        config = override_config(args.config)
    else:
        config = default_config
    kb: KnowledgeBase = select_kb(config)
    with open(args.output, "w") as asp_file:
        asp_file.writelines(
            rule + "\n"
            for rule in pathway_asp_generator(kb, ignore_orphan=args.ignore_orphan)
        )


if __name__ == "__main__":
    main()
