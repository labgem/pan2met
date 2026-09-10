#!/usr/bin/env python3

"""
Command line interface of pangenome2metabolism
"""

import argparse
import csv
import logging
import sys
from typing import List

import graph_tool as gt

import pan2met

from .config import default_config, override_config
from .inference import reactome
from .inference.pathologic.pythonic import infer_metabolism
from .inference.pathway_operon_filler import pathway_operon_filler
from .inference.proteic_complex import infer_complex
from .io.knowledge_base import KnowledgeBase, select_kb
from .utils import read_list, read_mapping, set_logging_level, write_output

logger = logging.getLogger()


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
    complex = infer_complex(kb, set(monomers))
    write_output(args.output, complex)


def pathway_operon_filler_command(args, config=default_config):
    """
    `pan2met pathway-operon-filler` subcommand
    """
    confident_enzyme_catalyzis = read_mapping(args.confident_enzyme_catalyzis)
    less_confident_enzyme_catalyzis = read_mapping(args.less_confident_enzyme_catalyzis)

    kb: KnowledgeBase = select_kb(config)

    logger.info("Loading the pangenome graph")
    pangenome_graph: gt.Graph = gt.load_graph(args.pangenome_graph)
    logger.info("Pangenome graph loaded.")

    distance = int(config["genomic_context"]["transitive_distance"])
    local_edge_jaccard_threshold = float(
        config["genomic_context"]["local_edge_jaccard_threshold"]
    )
    logger.info(f"Using a gap distance of {distance}")
    logger.info(
        f"Filtering edges with a minimum local Jaccard index of {local_edge_jaccard_threshold}"
    )
    logger.info("Compute pathway operon filler")

    result = pathway_operon_filler(
        pangenome_graph,
        kb,
        confident_enzyme_catalyzis,
        less_confident_enzyme_catalyzis,
        distance,
        local_edge_jaccard_threshold,
    )

    logger.info(f"Writing result to {args.output}")

    with open(args.output, "w") as output_file:
        writer = csv.writer(output_file, delimiter="\t")
        writer.writerow(["gene_family", "reaction", "pathway"])
        for gene_family, reaction_dict in result.items():
            for reaction, assignations in reaction_dict.items():
                for assignation in assignations:
                    writer.writerow([gene_family, reaction, assignation.pathway])


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
        help="infer the (pan)metabolism (i.e., a set of expected metabolic pathways)",
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
        description="infer a minimal set of monomers required to catalyze a set of reactions",
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

    # Parser proteic complex
    parser_proteic_complex = subparsers.add_parser(
        "proteic-complex",
        description="infer the list of proteic complex constructible from a list of protein monomers",
        help="infer the list of proteic complex constructible from a list of protein monomers",
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
    # Add a parser the for pathway-operon-filler subcommand
    parser_pathway_operon_filler = subparsers.add_parser(
        "pathway-operon-filler",
        description="identify potential reaction catalyzis of gene families with uncertain catalyzis based on (pan)-genome genomic context",
        help="identify potential reaction catalyzis of gene families with uncertain catalyzis based on (pan)-genome genomic context",
    )
    parser_pathway_operon_filler.add_argument(
        "--confident-enzyme-catalyzis",
        help="a two-column TSV with confident gene family / catalyzed reaction mapping; serves as anchor for less confident annotations in operonic structures",
        required=True,
    )
    parser_pathway_operon_filler.add_argument(
        "--less-confident-enzyme-catalyzis",
        help="a two-column TSV with less confident gene family / catalyzed reaction mapping: the enzyme we want to assign a reaction with more confidence, using the surrounding pathway reactions to get an hint",
        required=True,
    )
    parser_pathway_operon_filler.add_argument(
        "--pangenome-graph",
        help="a pangenome graph in graph-tool gt format",
        required=True,
    )
    parser_pathway_operon_filler.add_argument(
        "-o",
        "--output",
        help="output file path with a mapping from gene family identifier to pathway, catalyzed reaction and surrounding pathway reaction in the transitive closure",
        required=True,
    )

    # Set default function command function handler
    # parser_reactome.set_defaults(func=reactome_command)
    parser_metabolism.set_defaults(func=metabolism_command)
    parser_reverse_reactome.set_defaults(func=reverse_reactome_command)
    parser_proteic_complex.set_defaults(func=proteic_complex_command)
    parser_pathway_operon_filler.set_defaults(func=pathway_operon_filler_command)
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
