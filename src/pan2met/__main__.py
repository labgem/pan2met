#!/usr/bin/env python3

"""
Command line interface of pangenome2metabolism
"""
from typing import List

import argparse
import sys

from .utils import read_list, write_output, set_logging_level
# from .inference import reactome
from .inference.pathologic.pythonic import infer_metabolism
from .config import default_config, override_config
# from .inference.pathologic.aspic import

import pan2met


# def reactome_command(args, config=default_config):
#     """
#     `pan2met reactome` subcommand
#     """
#     monomers: List[str] = read_list(args.input)
#     reactions = reactome.infer_reactome_from_monomers_asp(
#         monomers, inference_rules_path=args.gpr
#     )
#     write_output(args.output, reactions)


# def reverse_reactome_command(args, config=default_config):
#     """
#     `pan2met reverse_reactome` subcommand
#     """
#     reactions: List[str] = read_list(args.reactions)
#     monomers: set[str] = reactome.minimal_monomer_set(
#         reactions, args.reverse_gpr, args.gpr
#     )
#     write_output(args.output, monomers)


def metabolism_command(args, config=default_config):
    """
    `pan2met metabolism` subcommand
    """
    reactome: set[str] = set(read_list(args.reactions))
    taxon_id: int = int(args.taxon_id)
    reason_filename = args.reason
    pathways = infer_metabolism(reactome, taxon_id, reason_filename, config=config)
    write_output(args.output, pathways)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="pan2met",
        description="Reconstruction metabolic network at (pan)-genome scale.",
    )
    # Base options
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {pan2met.__version__}", help="Display the version of pan2met"
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Set log level (-v: ERROR, -vv: WARNING, -vvv: INFO, -vvvv: DEBUG)",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="Override the default configuration with a custom configuration file.",
        default=None,
        required=False
    )
    # Define subcommands
    subparsers = parser.add_subparsers(help="commands", dest="command")
    # parser_reactome = subparsers.add_parser(
    #     "reactome",
    #     help="infer the (pan)reactome from the list of ortholog monomer identifiers",
    # )
    # parser_reactome.add_argument("--monomers", help="A file listing found monomers")
    # parser_reactome.add_argument(
    #     "-o",
    #     "--output",
    #     help="The path of the output file listing the reaction identifiers infered",
    #     required=True,
    # )
    # parser_reactome.add_argument(
    #     "--gpr",
    #     help="An AnsProlog Gene-Protein-Reaction (GPR) rules reference file",
    # )
    parser_metabolism = subparsers.add_parser(
        "metabolism",
        help="infer the (pan)metabolism",
    )
    parser_metabolism.add_argument(
        "-r",
        "--reactions",
        help="A file listing the reactions found in the (pan)-reactome",
        required=True,
    )
    parser_metabolism.add_argument(
        "-o",
        "--output",
        help="The path of the output file "
        "listing all identifiers of pathway infered to be present",
        required=True,
    )
    parser_metabolism.add_argument(
        "--reason",
        help="The path to an output file with a reason log.",
        default=None,
        required=False,
    )
    parser_metabolism.add_argument(
        "-t",
        "--taxon-id",
        help="The NCBI-Taxonomy tax id of the target organism.",
        required=False,
    )
    # # Reverse reactome problem
    # parser_reverse_reactome = subparsers.add_parser(
    #     "reverse_reactome",
    #     help="Infer a minimal set of monomers required to catalyze a set of reactions",
    # )
    # parser_reverse_reactome.add_argument(
    #     "-r",
    #     "--reactions",
    #     help="A file containing a list of reaction identifiers",
    #     required=True,
    # )
    # parser_reverse_reactome.add_argument(
    #     "-o",
    #     "--output",
    #     help="The path of the output file containing monomer identifiers",
    #     required=True,
    # )
    # parser_reverse_reactome.add_argument(
    #     "--gpr",
    #     help="An AnsProlog Gene-Protein-Reaction (GPR) rules reference file",
    #     required=True,
    # )
    # parser_reverse_reactome.add_argument(
    #     "--reverse-gpr",
    #     help="An AnsProlog 'reverse' Gene-Protein-Reaction (GPR) rules reference file",
    #     required=True,
    # )
    # Parser metabolism
    # parser_reactome.set_defaults(func=reactome_command)
    parser_metabolism.set_defaults(func=metabolism_command)
    # parser_reverse_reactome.set_defaults(func=reverse_reactome_command)
    return parser, parser.parse_args()


def main():
    parser, args = parse_arguments()
    set_logging_level(args.verbose)

    if args.config:
        config = override_config(args.config)
    else:
        config = default_config

    if hasattr(args, "func"):
        args.func(args, config=config)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
