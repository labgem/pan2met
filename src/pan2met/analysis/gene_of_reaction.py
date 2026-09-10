"""
Find the gene identifiers, of a reaction identifiers, based on a pangenome and gene
"""

from typing import Dict, Set
import csv
import argparse
from pathlib import Path
from collections import defaultdict

from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Gene
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.formats.readBinaries import check_pangenome_info

from ..io.knowledge_base import select_kb
from ..io.pangenome import read_pangenome_rtab
from ..config import default_config as config
from ..utils import read_mapping


def pangenome_strain_family_to_genes(
    pangenome: Pangenome,
) -> Dict[str, Dict[str, Set[str]]]:
    """
    Read a pangenome HDF5 file and return a mapping from each strain to the gene families present in the strain, and for each gene family, the set of gene identifiers of the family in the strain.
    :param pangenome: Path to the pangenome HDF5 file
    :return: A dictionary mapping each strain to a dictionary mapping each gene family identifier to the set of gene identifiers of the family in the strain.
    """
    strain_family_to_genes = defaultdict(lambda: defaultdict(set))
    for contig in pangenome.contigs:
        strain = contig.organism.name
        for gene in contig.genes:
            gene_family = gene.family
            strain_family_to_genes[strain][gene_family.name].add(gene.ID)
    return {
        strain: dict(family_to_genes)
        for strain, family_to_genes in strain_family_to_genes.items()
    }


def main():
    parser = argparse.ArgumentParser(
        description="Find the gene family identifiers of reactions of a pathway"
    )
    parser.add_argument("--pathway", help="Pathway identifier", required=True)
    parser.add_argument(
        "--pangenome-rtab",
        help="Path to the pangenome gene presence/absence table in Rtab format",
        required=True,
    )
    parser.add_argument(
        "--reaction-to-gene-family",
        help="Path to a file containing the mapping from reaction identifiers to gene family identifiers in TSV format, column 1: reaction identifier, column 2: gene family identifier.",
        required=True,
    )
    parser.add_argument(
        "--pangenome",
        help="PPanGGOLiN pangenome HDF5 file, used to find the gene identifiers of the gene families.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to the output file in TSV format, columns: pathway identifier, reaction identifier, gene family identifier",
        required=True,
    )
    args = parser.parse_args()

    pathway = args.pathway


    kb = select_kb(config)
    reactions = kb.non_spontaneous_reactions_of_pathway(pathway)
    reaction_to_gene_families = read_mapping(args.reaction_to_gene_family)



    # Open the pangenome
    pangenome = Pangenome()
    pangenome.add_file(Path(args.pangenome))
    check_pangenome_info(
        pangenome, need_families=True, need_annotations=True, disable_bar=False
    )  # do not forget the call to this function, otherwise the generator of contigs and genes will be empty as the pangenome would not be loaded.
    strain_to_gene_families: Dict[str, Set[str]] = read_pangenome_rtab(args.pangenome_rtab)
    strain_to_gene_families_to_genes: Dict[str, Dict[str, Set[str]]] = pangenome_strain_family_to_genes(
        pangenome
    )

    with open(args.output, "w") as output_file:
        writer = csv.writer(output_file, delimiter="\t")
        # Write header
        strains = sorted(strain_to_gene_families.keys())
        writer.writerow(["pathway", "reaction", "protein_family"] + strains)
        for reaction in reactions:
            for gene_family in reaction_to_gene_families[reaction]:
                writer.writerow(
                    [pathway, reaction, gene_family]
                    + [
                        ",".join(
                            strain_to_gene_families_to_genes[strain][gene_family]
                        )
                            if strain in strain_to_gene_families_to_genes
                            and gene_family in strain_to_gene_families_to_genes[strain]
                            else ""
                            for strain in strains
                    ]
                )


if __name__ == "__main__":
    main()
