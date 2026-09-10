"""
We want to identify reactions that are catalyzed by multiple different enzymes.
"""

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set

from ppanggolin.formats.readBinaries import check_pangenome_info
from ppanggolin.pangenome import Pangenome


def ppanggolin_gene_families_strain_count(
    pangenome: Pangenome,
) -> Dict[str, Dict[str, int]]:
    """
    Count how many genes belongs to each gene family in each strain.
    :param pangenome: a PPanGGOLiN pangenome
    :return: a dictionary mapping each gene family identifier to a dictionnary mapping each strain identifier to the count of genes
    belonging to the gene family in that strain.
    """
    gene_family_strain_count = defaultdict(lambda: defaultdict(int))
    for contig in pangenome.contigs:
        strain = contig.organism.name
        for gene in contig.genes:
            gene_family = gene.family
            gene_family_strain_count[gene_family.name][strain] += 1
    return {
        family_id: dict(strain_counts)
        for family_id, strain_counts in gene_family_strain_count.items()
    }


def count_enzymes_per_reactions(
    enzyme_reaction_mapping: Dict[str, Set[str]],
) -> Dict[str, int]:
    """
    Count the number of enzymes that catalyze each reaction.

    :param enzyme_reaction_mapping: A mapping from enzyme identifiers to sets of reaction identifiers they catalyze.
    :return: A dictionary mapping each reaction identifier to the count of enzymes that catalyze it.
    """
    reaction_enzyme_count = defaultdict(int)
    for _enzyme, catalyzed_reactions in enzyme_reaction_mapping.items():
        for reaction in catalyzed_reactions:
            reaction_enzyme_count[reaction] += 1
    return dict(reaction_enzyme_count)


def count_enzymes_per_reaction_within_strain(
    enzyme_reaction_mapping: Dict[str, Set[str]],
    gene_family_strain_count: Dict[str, Dict[str, int]],
) -> Dict[str, Dict[str, int]]:
    """
    In a reactome predicted at the scale of a pangenome,
    counting the whole number of enzymes that catalyze a reaction is misleading, as some of these enzymes may be present in different strains.
    This function counts the number of enzymes catalyzing a reaction within each strain.
    :param enzyme_reaction_mapping: A mapping from enzyme identifiers to sets of reaction identifiers they catalyze.
    :param gene_family_strain_count: A dictionary mapping each gene family identifier to a dictionary mapping each strain identifier to the count of genes belonging to the gene family in that strain.
    :return: A dictionary mapping each reaction identifier to a dictionary mapping each strain to the count of enzymes that catalyze it within that strain.
    """
    reaction_enzyme_count_by_strain = defaultdict(lambda: defaultdict(int))
    for enzyme, catalyzed_reactions in enzyme_reaction_mapping.items():
        for reaction in catalyzed_reactions:
            if enzyme in gene_family_strain_count:
                for strain, count in gene_family_strain_count[enzyme].items():
                    reaction_enzyme_count_by_strain[reaction][strain] += count
    return {
        reaction: dict(strain_counts)
        for reaction, strain_counts in reaction_enzyme_count_by_strain.items()
    }


def main():
    parser = argparse.ArgumentParser(
        description="Count the number of enzymes that catalyze each reaction."
    )
    parser.add_argument(
        "--enzyme-reaction-mapping",
        type=str,
        required=True,
        help="Path to a file containing the mapping from enzymes to reactions they catalyze in TSV format, column 1: enzyme, column 2: reaction identifier.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to the output file where the reaction enzyme counts will be saved in TSV format, column 1: reaction identifier, column 2: count of enzymes catalyzing the reaction, column 3: strain identifier with the option --by-strain.",
    )
    parser.add_argument(
        "--by-strain",
        action="store_true",
        help="Whether to count the number of enzymes catalyzing each reaction within each strain (if set) or to count the total number of enzymes catalyzing each reaction across the whole pangenome (if not set).",
    )
    parser.add_argument(
        "--pangenome",
        type=str,
        required=False,
        help="Path to the pangenome HDF5 file, required if --by-strain is set.",
    )
    args = parser.parse_args()

    # Load the protein-families to reactions mapping
    input_enzyme_reaction_mapping = args.enzyme_reaction_mapping
    with open(input_enzyme_reaction_mapping, "r") as f:
        enzyme_reaction_mapping = defaultdict(set)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                enzyme, reactions = parts
                for reaction in reactions.split(","):
                    enzyme_reaction_mapping[enzyme].add(reaction)

    if args.by_strain:
        # We need to go back to the pangenome file,
        # to know how many enzymes are present in each strain within a protein family.
        if not args.pangenome:
            raise ValueError(
                "The --pangenome argument is required when using --by-strain option."
            )
        pangenome = Pangenome()
        pangenome.add_file(Path(args.pangenome))
        check_pangenome_info(
            pangenome, need_families=True, need_annotations=True, disable_bar=True
        )  # do not forget the call to this function, otherwise the generator of contigs and genes will be empty as the pangenome would not be loaded.
        gene_family_strain_count = ppanggolin_gene_families_strain_count(pangenome)
        reaction_enzyme_count_by_strain = count_enzymes_per_reaction_within_strain(
            enzyme_reaction_mapping, gene_family_strain_count
        )
        with open(args.output, "w") as output_file:
            for reaction, strain_counts in reaction_enzyme_count_by_strain.items():
                for strain, count in strain_counts.items():
                    output_file.write(f"{reaction}\t{count}\t{strain}\n")
    else:
        reaction_enzyme_count = count_enzymes_per_reactions(enzyme_reaction_mapping)
        with open(args.output, "w") as output_file:
            for reaction, count in reaction_enzyme_count.items():
                output_file.write(f"{reaction}\t{count}\n")


if __name__ == "__main__":
    main()
