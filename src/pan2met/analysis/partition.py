"""
Generate a report of pathway / reaction pangenome gene families partitioning by PPanGGOLiN.

PPanGGOLiN partition gene families of a pangenome into 3 classes: core, shell and cloud genes.

Given a reaction, an enzyme catalyzing this reaction, and the pangenome graph datastructure of PPanGGOLiN,
this script associates a pangenome partition to the reaction.

TODO: summary partition frequency of a pathway.
TODO:
"""

from enum import Enum

from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Gene
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.formats.readBinaries import check_pangenome_info


class PangenomePartition(Enum):
    CORE = 1
    SHELL = 2
    CLOUD = 3


def ppanggolin_get_gene_by_identifier(
    gene_identifier: str, pangenome: Pangenome
) -> Gene:
    """
    Get a gene object from PPanGGOLiN pangenome by the identifier of the gene.

    :param gene_identifier: a gene identifier in PPanGGOLiN pangenome
    :return: a Gene object from PPanGGOLiN with identifier gene_identifier
    """
    print(pangenome.organisms)
    for contig in pangenome.contigs:
        for gene in contig.genes:
            if gene.ID == gene_identifier:
                return gene
    raise ValueError(f"ID {gene_identifier} not found in any contig of the pangenome")


def ppanggolin_partition_of_gene_family(
    gene_name: str, pangenome: Pangenome
) -> PangenomePartition:
    """
    Get the pangenome partition class of a gene as annoted in PPanGGOLiN pangenome.

    :param gene_name: a gene_name
    """
    gene: Gene = ppanggolin_get_gene_by_identifier(gene_name, pangenome)
    gene_family: GeneFamily = gene.family
    partition = gene_family.named_partition
    match partition:
        case "persistent":
            return PangenomePartition.CORE
        case "shell":
            return PangenomePartition.SHELL
        case "cloud":
            return PangenomePartition.CLOUD
        case _:
            raise ValueError(
                f"{partition} partition from PPanGGOLiN for family {gene_family.ID} not understood."
            )


def main():
    pangenome = Pangenome()
    pangenome.add_file(
        "/home/sortion/data/home/sortion/data/pangbank/GTDB_refseq_s__Escherichia_coli_id615.h5"
    )
    check_pangenome_info(
        pangenome, need_families=True, need_annotations=True, disable_bar=True
    )
