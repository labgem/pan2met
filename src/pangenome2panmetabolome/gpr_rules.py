#!/usr/bin/env python3

"""
Gene-Protein-Reaction rules in Answer Set Programming

Launch Pathway Tools python API with:
    pathway-tools -lisp -python-local-only-non-strict
"""

import logging
import argparse
from typing import Iterable

import pythoncyc

from .io.metacyc import remove_pipes
from .asp import catalysis_asp_rule, protein_complex_asp_rule


logger = logging.getLogger("pangenome2panmetabolome:gpr")


def is_proteic_complex(pgdb, enzyme: str) -> bool:
    """
    Assert True if `enzyme` is a protein complex.
    """
    return pgdb.complex_p(enzyme)


def proteic_complex_subunits(pgdb, enzyme: str) -> Iterable[str]:
    """
    Iter over enzyme protein complex subunits
    """
    enzyme_frame_objects = pgdb.get_frame_objects([enzyme])[0]
    components = enzyme_frame_objects["components"]
    return components


def is_homomeric(pgdb, polymer: str) -> bool:
    components = pgdb.get_frame_objects([polymer])[0]["components"]
    return (
        is_proteic_complex(pgdb, polymer)
        and components is not None
        and len(components) == 1
    )


def gpr_asp_generator(pgdb) -> Iterable[str]:
    """
    Generate Answer Set Programming rules for GPR rules
    from a BioCyc PGDB.

    How are the AnsProlog rules formed?

    Each protein complex composition is represented with 'AND' rules.
    In Clingo the AND symbol is the comma ','.

    So, assuming there exist a complex protein CPLX-1 composed of MONOMER-A and MONOMER-B,
    this function will generate a rule:

    .. code :: prolog

      complex("CPLX-1") :- monomer("MONOMER-A") , monomer("MONOMER-B").

    this will ensure that the atom complex("CPLX-1") will be in the clingo Answer Set whenever
    both monomer("MONOMER-A") and monomer("MONOMER-B") atoms are True.

    Moreover, for each reaction RXN-1 catalyzed by CPLX-1 the following rule will be added:

    .. code :: prolog

      reaction("RXN-1") :- complex("CPLX-1").

    Similarly, if the reaction is catalyzed by MONOMER-C, the following rule will be added:

    .. code :: prolog

       reaction("RXN-1") :- monomer("MONOMER-C").


    Arguments
    ---------

        pgdb: the pythoncyc instance (e.g., pythoncyc.select_organism("meta"))

    Yields
    ------

        ASP rules

    Remark
    ------

    Note that the rules are listing asp rules, coming from the list of known reactions.
    It does not cover exhaustively the set of protein complex: if a protein complex is referenced
    in MetaCyc but is not associated with any referenced reaction, it will simply not be reported
    in any GPR ASP rule.
    """
    # Iter all reactions and generate GPR rules
    for reaction in pgdb.all_rxns():
        logger.debug(f"Exploring {reaction}")
        for enzyme in pgdb.enzymes_of_reaction(pgdb, reaction):
            logger.debug(f"{reaction} catalyzed by {enzyme}")
            if is_proteic_complex(pgdb, enzyme):
                if "MONOMER" in enzyme and not is_homomeric(pgdb, enzyme):
                    logger.info(f"complex MONOMER but not homomeric: {enzyme}")
                yield catalysis_asp_rule(
                    remove_pipes(enzyme), "complex", remove_pipes(reaction)
                )
                logger.debug(f"{enzyme} is a complex")
                components = proteic_complex_subunits(pgdb, enzyme)
                if components is not None:
                    yield protein_complex_asp_rule(
                        remove_pipes(enzyme), map(remove_pipes, components)
                    )
                    for monomer in components:
                        logger.debug(f"{enzyme} contains {monomer}")
                else:
                    logger.info(
                        f"{enzyme} complex does not contain referenced components"
                    )
            else:
                yield catalysis_asp_rule(
                    remove_pipes(enzyme), "monomer", remove_pipes(reaction)
                )
                logger.debug(f"{enzyme} is a monomer")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Export MetaCyc Gene-Protein-Reaction rules in Answer Set Programming AnsProlog format"
    )
    parser.add_argument(
        "-o", "--output", help="Output .lp filename with ASP rules", required=True
    )
    return parser, parser.parse_args()


def main():
    logging.basicConfig()
    logger.setLevel(logging.INFO)
    parser, args = parse_arguments()
    pgdb = pythoncyc.select_organism("meta")
    with open(args.output, "w") as asp_file:
        for rule in gpr_asp_generator(pgdb):
            asp_file.write(rule + "\n")


if __name__ == "__main__":
    main()
