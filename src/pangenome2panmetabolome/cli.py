#!/usr/bin/env python3

"""
Command line interface of pangenome2metabolome
"""

import argparse
import logging

from .utils import logger, read_list, write_output

from .reactome import infer_reactome_from_monomers_asp, minimal_monomer_set
from .metabolome import infer_metabolome


def reactome(args):
    """
    `pangenome2panmetabolome reactome` subcommand
    """
    monomers: list[str] = read_list(args.input)
    reactions = infer_reactome_from_monomers_asp(
        monomers, inference_rules_path=args.gpr
    )
    write_output(args.output, reactions)


def reverse_reactome(args):
    """
    `pangenome2panmetabolome reverse_reactome` subcommand
    """
    reactions: list[str] = read_list(args.reactions)
    monomers: set[str] = minimal_monomer_set(reactions, args.reverse_gpr, args.gpr)
    write_output(args.output, monomers)


def metabolome(args):
    """
    `pangenome2panmetabolome metabolome` subcommand
    """
    reactome: set[str] = set(read_list(args.input))
    pathways = infer_metabolome(reactome)
    write_output(args.output, pathways)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="pangenome2panmetabolome",
        description="Predict metabolic pathway presence from a (pan)genome",
    )
    parser.add_argument("-v", "--verbose", action="count", default=0)
    subparsers = parser.add_subparsers(help="pangenome2panmetabolome subcommands")
    parser_reactome = subparsers.add_parser(
        "reactome", help="infer the (pan)reactome with Answer Set Programming (clingo)"
    )
    parser_reactome.add_argument(
        "-i", "--input", help="A file listing found monomers", required=True
    )
    parser_reactome.add_argument(
        "-o",
        "--output",
        help="The path of the output file listing the reaction identifiers infered",
        required=True,
    )
    parser_reactome.add_argument(
        "--gpr",
        help="An AnsProlog Gene-Protein-Reaction (GPR) rules reference file",
        required=True,
    )
    parser_reactome.set_defaults(func=reactome)
    parser_metabolome = subparsers.add_parser(
        "metabolome", help="infer the (pan)metabolome"
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
    parser_reverse_reactome.set_defaults(func=reverse_reactome)
    parser_metabolome.set_defaults(func=metabolome)
    return parser, parser.parse_args()


def set_logging_level(verbose_intensity):
    logging.basicConfig()
    logging_levels = [logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR][::-1]
    logging_levels_repr = ["DEBUG", "INFO", "WARNING", "ERROR"][::-1]
    logging_level_index = min(verbose_intensity, len(logging_levels) - 1)
    logging_level = logging_levels[logging_level_index]
    logging.getLogger().setLevel(logging_level)
    logger.debug(f"logging level set to {logging_levels_repr[logging_level_index]}")


def main():
    _parser, args = parse_arguments()
    set_logging_level(args.verbose)
    # Call the subcommand function
    args.func(args)
