#!/usr/bin/env python3

"""
Command line interface of pangenome2metabolome
"""

import argparse
import logging
import sys


import pythoncyc

from .io import sparql
from .utils import logger, read_list, write_output
from . import reactome
from . import metabolome


def reactome_command(args):
    """
    `pangenome2panmetabolome reactome` subcommand
    """
    monomers: list[str] = read_list(args.input)
    if args.straight:
        if args.monomers is not None:
            try:
                pgdb = pythoncyc.select_organism("meta")
            except Exception as e:
                logger.error(
                    "Make sure you launch PathwayTools python API with "
                    "`pathway-tools -lisp -python-local-only-non-strict` before running this command"
                )

                raise e
            reactions_by_monomers = reactome.infer_reactome_from_monomers(
                monomers, pgdb
            )
            logger.info(
                f"{len(reactions_by_monomers)} reactions identified from monomers."
            )
        if args.ec_numbers is not None:
            kb = sparql.SPARQLBasedKnowledgeBase()
            reactions_by_ec = reactome.infer_reactome_from_ec_numbers(monomers, kb)
            logger.info(
                f"{len(reactions_by_ec)} reactions identified by EC-number matching."
            )
        if args.monomers is not None and args.ec_numbers is not None:
            reactions = reactions_by_monomers | reactions_by_ec

    else:
        reactions = reactome.infer_reactome_from_monomers_asp(
            monomers, inference_rules_path=args.gpr
        )
        write_output(args.output, reactions)


def reverse_reactome_command(args):
    """
    `pangenome2panmetabolome reverse_reactome` subcommand
    """
    reactions: list[str] = read_list(args.reactions)
    monomers: set[str] = reactome.minimal_monomer_set(
        reactions, args.reverse_gpr, args.gpr
    )
    write_output(args.output, monomers)


def metabolome_command(args):
    """
    `pangenome2panmetabolome metabolome` subcommand
    """
    reactome: set[str] = set(read_list(args.input))
    pathways = metabolome.infer_metabolome(reactome)
    write_output(args.output, pathways)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="pangenome2panmetabolome",
        description="Predict metabolic pathway presence from a (pan)genome",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Log level (-v: ERROR, -vv: WARNING, -vvv: INFO, -vvvv: DEBUG)",
    )
    subparsers = parser.add_subparsers(help="commands", dest="command")
    parser_reactome = subparsers.add_parser(
        "reactome",
        help="infer the (pan)reactome from the list of ortholog monomer id",
    )
    parser_reactome.add_argument("--monomers", help="A file listing found monomers")
    parser_reactome.add_argument(
        "-o",
        "--output",
        help="The path of the output file listing the reaction identifiers infered",
        required=True,
    )
    parser_reactome.add_argument(
        "--gpr",
        help="An AnsProlog Gene-Protein-Reaction (GPR) rules reference file",
    )
    parser_reactome.add_argument(
        "--straight",
        action=argparse.BooleanOptionalAction,
        help="Use 'straightforward' inference, using PythonCyc",
    )
    parser_reactome.add_argument("--ec-numbers", help="A file listing EC-numbers")
    parser_metabolome = subparsers.add_parser(
        "metabolome",
        help="infer the (pan)metabolome",
    )
    parser_metabolome.add_argument(
        "-i",
        "--input",
        help="A file listing the reactions found in the (pan)-reactome",
        required=True,
    )
    parser_metabolome.add_argument(
        "-o",
        "--output",
        help="The path of the output file "
        "listing all identifiers of pathway infered to be present",
        required=True,
    )
    parser_reverse_reactome = subparsers.add_parser(
        "reverse_reactome",
        help="Infer a minimal set of monomers required to catalyze a set of reactions",
    )
    parser_reverse_reactome.add_argument(
        "-r",
        "--reactions",
        help="A file containing a list of reaction identifiers",
        required=True,
    )
    parser_reverse_reactome.add_argument(
        "-o",
        "--output",
        help="The path of the output file containing monomer identifiers",
        required=True,
    )
    parser_reverse_reactome.add_argument(
        "--gpr",
        help="An AnsProlog Gene-Protein-Reaction (GPR) rules reference file",
        required=True,
    )
    parser_reverse_reactome.add_argument(
        "--reverse-gpr",
        help="An AnsProlog 'reverse' Gene-Protein-Reaction (GPR) rules reference file",
        required=True,
    )
    parser_reactome.set_defaults(func=reactome_command)
    parser_metabolome.set_defaults(func=metabolome_command)
    parser_reverse_reactome.set_defaults(func=reverse_reactome_command)
    return parser, parser.parse_args()


def set_logging_level(verbose_intensity):
    logging.basicConfig()
    logging_levels = [logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR][::-1]
    logging_level_index = min(verbose_intensity, len(logging_levels) - 1)
    logging_level = logging_levels[logging_level_index]
    logging.getLogger().setLevel(logging_level)


def main():
    parser, args = parse_arguments()
    set_logging_level(args.verbose)
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)
