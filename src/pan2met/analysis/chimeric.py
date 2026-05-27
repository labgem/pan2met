"""
Find chimeric pathways among the set of predicted pathways.

A chimeric pathway is defined as a pathway with the set of reactions requiring catalysis
spread across different different genomes of the pangenome.

To identify chimeric pathways, we use the file .Rtab with the matrix of presence-absence of gene families
and a knowledge base with allowing to get the set of reactions in each pathway as well as a set of predicted pathways.
"""

from typing import Dict, Set, List
import argparse
import itertools
import logging
import configparser

from ..io.knowledge_base import KnowledgeBase, select_kb
from ..io.pangenome import read_pangenome_rtab
from ..utils import read_list, read_mapping, reverse_mapping
from ..config import config


def strain_with_reaction(
    reaction: str,
    reaction_to_genes: Dict[str, Set[str]],
    gene_to_strains: Dict[str, Set[str]],
) -> Set[str]:
    """
    Get the list of strains of the pangenomes that is expected to have a enzyme catalyzing the reaction.
    :param reaction: a reaction identifier
    :param reaction_to_genes: a dictionnary with key reaction id and value a set of gene families identifiers
    :param gene_to_strains: a dictionnary with key a gene family identifier and value a set of pangenome strain with a representent of this gene family.
    :return: a set of pangenome strain identifiers.
    """
    if reaction not in reaction_to_genes:
        return []
    genes = reaction_to_genes[reaction]
    return {strain for gene in genes for strain in gene_to_strains[gene]}


def is_chimeric(
    pathway: str,
    reaction_to_genes: Dict[str, Set[str]],
    gene_to_strains: Dict[str, Set[str]],
    kb: KnowledgeBase,
) -> bool:
    """
    Assert whether a pathway is chimeric or not.
    A pathway is chimeric is a part of its reaction is present in a set of strains and another part of the reactions is in another completely distinct set of strains.
    :param pathway: a pathway identifier, available in the knowledge base
    :return: True if the pathway is chimeric, False otherwise.
    """
    pathway_non_orphan_non_spontaneous_reactions: List[str] = (
        kb.non_orphan_non_spontaneous_reactions_of_pathway(pathway)
    )
    reaction_strain_sets: List[Set[str]] = []
    for reaction in pathway_non_orphan_non_spontaneous_reactions:
        reaction_strains: Set[str] = strain_with_reaction(
            reaction, reaction_to_genes, gene_to_strains
        )
        reaction_strain_sets.append(reaction_strains)

    for set_a, set_b in itertools.combinations(reaction_strain_sets):
        if len(set_a) >= 1 and len(set_b) >= 1:
            if len(set_a.intersection(set_b)) == 0:
                return True
    return False


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pathways", help="Input list of pathways", required=True)
    parser.add_argument(
        "--gene-presence-absence",
        help="Input .Rtab file containing a matrix of pangenome gene presence-absence",
        required=True,
    )
    parser.add_argument(
        "--gene-reaction",
        help="TSV file mapping pangenome gene families identifiers to reaction identifiers (column 1: gene family identifier, column 2: reaction identifier)",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file path listing chimeric pathways.",
        required=True,
    )
    parser.add_argument("-c", "--config", help="Config file path", required=False)
    return parser


def main():
    parser = argument_parser()
    args = parser.parse_args()

    strain_to_genes: Dict[str, Set[str]] = read_pangenome_rtab(
        args.gene_presence_absence
    )
    pathways: List[str] = read_list(args.pathways)
    gene_to_reactions: Dict[str, Set[str]] = read_mapping(args.gene_reaction)
    reaction_to_genes: Dict[str, Set[str]] = reverse_mapping(gene_to_reactions)
    gene_to_strains: Dict[str, Set[str]] = reverse_mapping(strain_to_genes)

    # Update config globally overriding default_config with keys from given config filename
    if args.config is not None:
        logging.info(f"Overriding default configuration with {args.config}")
        global config
        config_override = configparser.ConfigParser()
        config_override.read([args.config])
        config.update(config_override)

    kb = select_kb(config["reference"]["source"])

    with open(args.output, "w") as output_file:
        for pathway in pathways:
            if is_chimeric(pathway, reaction_to_genes, gene_to_strains, kb):
                output_file.write(pathway + "\n")


if __name__ == "__main__":
    main()
