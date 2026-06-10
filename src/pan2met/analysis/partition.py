"""
Generate a report of pathway / reaction pangenome gene families partitioning by PPanGGOLiN.

PPanGGOLiN partition gene families of a pangenome into 3 classes: core, shell and cloud genes.

Given a reaction, an enzyme catalyzing this reaction, and the pangenome graph datastructure of PPanGGOLiN,
this script associates a pangenome partition to the reaction.

TODO: summary partition frequency of a pathway.
TODO:
"""

import argparse
import csv
from pathlib import Path
from typing import Dict

from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Gene
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.formats.readBinaries import check_pangenome_info


def ppanggolin_extract_gene_to_partition_mapping(
    pangenome: Pangenome,
) -> Dict[str, str]:
    """
    Construct a mapping between gene identifier and gene family partition

    :param pangenome: a PPanGGOLiN pangenome
    :return: a dictionnary with gene name as key, and corresponding gene partition as value
    """
    mapping: Dict[str, str] = {}
    for contig in pangenome.contigs:
        for gene in contig.genes:
            gene_family = gene.family
            partition = gene_family.named_partition
            mapping[gene.ID] = partition
    return mapping


def ppanggolin_get_gene_by_identifier(
    gene_identifier: str, pangenome: Pangenome
) -> Gene:
    """
    Get a gene object from PPanGGOLiN pangenome by the identifier of the gene.

    :param gene_identifier: a gene identifier in PPanGGOLiN pangenome
    :return: a Gene object from PPanGGOLiN with identifier gene_identifier
    """
    for contig in pangenome.contigs:
        for gene in contig.genes:
            if gene.ID == gene_identifier:
                return gene
    raise ValueError(f"ID {gene_identifier} not found in any contig of the pangenome")


def ppanggolin_partition_of_gene_family(
    gene_identifier: str, pangenome: Pangenome, gene_to_family: Dict[str, str]
) -> str:
    """
    Get the pangenome partition class of a gene as annoted in PPanGGOLiN pangenome.

    :param gene_identifier: a gene identifier in PPanGGOLiN genome
    :param pangenome: a PPanGGOLiN pangenome
    :param gene_to_family: a dictionnary mapping a gene identifier to a gene family partition name
    :return: a PPanGGOLiN partition name (in { 'persistent', 'shell', 'cloud' })
    """
    gene: Gene = ppanggolin_get_gene_by_identifier(gene_identifier, pangenome)
    gene_family: GeneFamily = gene.family
    partition = gene_family.named_partition
    return partition


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pangenome", help="PPanGGOLiN pangenome HDF5 file", required=True
    )
    parser.add_argument("--reactions", help="Reaction list", required=True)
    parser.add_argument(
        "--reaction-enzyme",
        help="Two-column TSV with column 1: reaction identifier and column 2: pangenome monomer identifier catalyzing the reaction.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output, a three-column file with columns 1: reaction identifier, column 2: enzyme gene identifier, column 3: partition (in {persistent, shell, cloud})",
        required=True,
    )
    # TODO: see what to do when a reaction can be catalyzed by several enzymes that may be not on the same pangenome gene family partition
    return parser


def main():
    parser = argument_parser()
    args = parser.parse_args()
    # Open the pangenome
    pangenome = Pangenome()
    pangenome.add_file(Path(args.pangenome))
    check_pangenome_info(
        pangenome, need_families=True, need_annotations=True, disable_bar=True
    )  # do not forget the call to this function, otherwise the generator of contigs and genes will be empty as the pangenome would not be loaded.
    # Prepare a dictionnary mapping a gene identifier to a gene partition
    gene_to_partition: Dict[str, str] = ppanggolin_extract_gene_to_partition_mapping(
        pangenome
    )

    # Read the reference mapping of reaction to enzyme file
    with open(args.reaction_enzyme, "r") as reaction_enzyme_file:
        reaction_to_enzyme_gene: Dict[str, str] = {}
        reader = csv.reader(reaction_enzyme_file, delimiter="\t")
        for row in reader:
            if len(row) == 2:
                reaction = row[0]
                enzyme_gene = row[1]
                reaction_to_enzyme_gene[reaction] = enzyme_gene

    # Write the output file
    with open(args.reactions, "r") as reactions_file:
        with open(args.output, "w") as output_file:
            writer = csv.writer(output_file, delimiter="\t")
            writer.writerow(
                ["reaction", "enzyme_gene", "pangenome_partition"]
            )  # write header
            for row in reactions_file:
                reaction = row.strip()
                if reaction != "" and reaction in reaction_to_enzyme_gene:
                    enzyme_gene = reaction_to_enzyme_gene[reaction]
                    if enzyme_gene in gene_to_partition:
                        pangenome_partition = gene_to_partition[enzyme_gene]
                    else:
                        pangenome_partition = "NA"
                else:
                    enzyme_gene = "NA"
                    pangenome_partition = "NA"
                writer.writerow([reaction, enzyme_gene, pangenome_partition])


if __name__ == "__main__":
    main()
