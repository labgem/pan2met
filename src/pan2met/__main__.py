#!/usr/bin/env python3

"""
Command line interface of pangenome2metabolism
"""

import argparse
import sys
from typing import List

import pan2met

from .config import default_config, override_config
from .inference import reactome
from .inference.pathologic.pythonic import infer_metabolism
from .inference.proteic_complex import infer_complex
from .io.knowledge_base import KnowledgeBase, select_kb
from .utils import read_list, set_logging_level, write_output


def reactome_command(args, config=default_config):
    """
    `pan2met reactome` subcommand
    """
    monomers: List[str] = read_list(args.input)
    reactions = reactome.infer_reactome_from_monomers_asp(
        monomers, inference_rules_path=args.gpr
    )
    write_output(args.output, reactions)


def reverse_reactome_command(args, config=default_config):
    """
    `pan2met reverse-reactome` subcommand
    """
    reactions: List[str] = read_list(args.reactions)
    monomers: set[str] = reactome.minimal_monomer_set(
        reactions, args.reverse_gpr, args.gpr
    )
    write_output(args.output, monomers)


def metabolism_command(args, config=default_config):
    """
    `pan2met metabolism` subcommand
    """
    reactome: set[str] = set(read_list(args.reactions))
    taxon_id: int = int(args.taxon_id)
    reason_filename = args.reason
    pathways = infer_metabolism(reactome, taxon_id, reason_filename, config=config)
    write_output(args.output, pathways)


def proteic_complex_command(args, config=default_config):
    """
    `pan2met protein-complex` subcommand
    """
    monomers: list[str] = read_list(args.monomers)
    kb: KnowledgeBase = select_kb(config)
    complex = infer_complex(kb, monomers)
    write_output(args.output, complex)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="pan2met",
        description="""Reconstruct metabolic networks at (pan)-genome scale.""",
        epilog=f"""%(prog)s ({pan2met.__version__}) is an open-source bioinformatics tool developed by the LABGeM team, and distributed under the CeCILL Free Sofware License Agreement.""",
    )

    # Base options
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {pan2met.__version__}",
        help="display the version of pan2met",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="set log level (-v: ERROR, -vv: WARNING, -vvv: INFO, -vvvv: DEBUG)",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="override the default configuration with a custom configuration file.",
        default=None,
        required=False,
    )
    # Define subcommands
    subparsers = parser.add_subparsers(help="commands", dest="command")
    parser_metabolism = subparsers.add_parser(
        "metabolism",
        description="infer the (pan)metabolism (i.e., a set of expected metabolic pathways)",
    )
    parser_metabolism.add_argument(
        "-r",
        "--reactions",
        help="a file listing the reactions found in the (pan)-reactome",
        required=True,
    )
    parser_metabolism.add_argument(
        "-o",
        "--output",
        help="the path of the output file "
        "listing all identifiers of pathway infered to be present",
        required=True,
    )
    parser_metabolism.add_argument(
        "--reason",
        help="the path to an output file with a reason log.",
        default=None,
        required=False,
    )
    parser_metabolism.add_argument(
        "-t",
        "--taxon-id",
        help="the NCBI-Taxonomy tax id of the target organism.",
        required=False,
    )
    # Reverse reactome problem
    parser_reverse_reactome = subparsers.add_parser(
        "reverse-reactome",
        help="infer a minimal set of monomers required to catalyze a set of reactions",
    )
    parser_reverse_reactome.add_argument(
        "-r",
        "--reactions",
        help="a file containing a list of reaction identifiers",
        required=True,
    )
    parser_reverse_reactome.add_argument(
        "-o",
        "--output",
        help="the path of the output file containing monomer identifiers",
        required=True,
    )
    parser_reverse_reactome.add_argument(
        "--gpr",
        help="an AnsProlog Gene-Protein-Reaction (GPR) rules reference file",
        required=False,
    )
    parser_reverse_reactome.add_argument(
        "--reverse-gpr",
        help="an AnsProlog 'reverse' Gene-Protein-Reaction (GPR) rules reference file",
        required=False,
    )
    # Parser metabolism

    # Parser proteic complex
    parser_proteic_complex = subparsers.add_parser(
        "proteic-complex",
        description="infer the list of proteic complex constructible from a list of protein monomers",
    )
    parser_proteic_complex.add_argument(
        "-m",
        "--monomers",
        help="input path to a list of protein monomer identifiers",
        required=True,
    )
    parser_proteic_complex.add_argument(
        "-o",
        "--output",
        help="output path to a list of constructible proteic complex",
        required=True,
    )

    # Set default function command function handler
    # parser_reactome.set_defaults(func=reactome_command)
    parser_metabolism.set_defaults(func=metabolism_command)
    parser_reverse_reactome.set_defaults(func=reverse_reactome_command)
    parser_proteic_complex.set_defaults(func=proteic_complex_command)
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
