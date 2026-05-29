"""
Completion matrix
"""

import argparse
from pathlib import Path
from typing import List, Set, Literal, Dict
import logging
import tqdm

from ..config import config
from ..utils import read_mapping, reverse_mapping, read_list
from ..io.knowledge_base import KnowledgeBase, select_kb
from ..io.pangenome import (
    read_pangenome_rtab,
    read_partition_files,
    read_functional_module_file,
)


logger = logging.getLogger("pan2met:analysis:completion")
logging.basicConfig(level=logging.DEBUG)


def write_reaction_presence_absence_by_strain(
    filename: Path,
    reaction_to_families: Dict[str, Set[str]],
    strain_to_families: Dict[str, Set[str]],
    partition_to_families: Dict[Literal["persistent", "shell", "cloud"], Set[str]],
):
    """
    Write a reaction presence/absence matrix at the strain level.

    The output file is a CSV with the following columns:
    - reaction: an identifier of reaction
    - families: a comma-seperated list of protein families associated to the reaction
    - persistent, shell, cloud - percentage: the percentage of protein families in each of the pangenome partition
    Then, for all strain identifiers, an other column with binary value indicates whether the reaction is predicted to be in the reactome of the strain.

    :param filename: output filename
    :param reaction_to_families: a dictionnary with a reaction identifier key and a set of protein families identifiers value.
    :param strain_to_families: a dictionnary with a strain identifier key and a set of protein families identifiers value.
    :param partition_to_families: a dictionnary with key in {'persistent', 'shell', 'cloud'} and value a set of protein family identifiers.
    """
    logger.info(f"Computing reaction presence absence matrix for {len(reaction_to_families.keys())} reactions in {len(strain_to_families.keys())} strains.")
    with open(filename, "w") as output_file:
        # Write the file header
        header: List[str] = [
            "reaction",
            "families",
            "persistent %",
            "shell %",
            "cloud %",
        ]
        header = header + [strain for strain in strain_to_families.keys()]
        output_file.write("\t".join(header) + "\n")
        for reaction in tqdm.tqdm(reaction_to_families.keys()):
            for strain in strain_to_families.keys():
                families: Set[str] = reaction_to_families[reaction]
                families_str = ",".join(families)
                persistent_percentage: float = (
                    len(families.intersection(partition_to_families["persistent"]))
                    * 100
                    / len(families)
                )
                shell_percentage: float = (
                    len(families.intersection(partition_to_families["shell"]))
                    * 100
                    / len(families)
                )
                cloud_percentage: float = (
                    len(families.intersection(partition_to_families["cloud"]))
                    * 100
                    / len(families)
                )

                record = [
                    reaction,
                    families_str,
                    persistent_percentage,
                    shell_percentage,
                    cloud_percentage,
                ]
                for strain in strain_to_families.keys():
                    if len(families.intersection(strain_to_families[strain])) >= 1:
                        record.append("1")
                    else:
                        record.append("0")
                output_file.write("\t".join(map(str, record)) + "\n")


def write_pathway_completion_by_strain(
    prefix: str,
    pathways: List[str],
    strain_to_families: Dict[str, Set[str]],
    partition_to_families: Dict[str, Set[str]],
    reaction_to_families: Dict[str, Set[str]],
    module_to_families: Dict[str, Set[str]],
    pangenome_reactions: Set[str],
    use_orphan: bool,
    kb: KnowledgeBase,
):
    logger.info(f"Computing completion values for {len(pathways)} pathways in {len(strain_to_families.keys())} strains, {'with orphan reactions' if use_orphan else 'ignoring orphan reactions'}.")
    if use_orphan:
        pathways_to_reactions = {
            pathway: set(kb.non_spontaneous_reactions_of_pathway(pathway))
            for pathway in pathways
        }
        filename = prefix + "_pathway_completion_by_strain.tsv"
    else:
        pathways_to_reactions = {
            pathway: set(kb.non_orphan_non_spontaneous_reactions_of_pathway(pathway))
            for pathway in pathways
        }
        filename = prefix + "_pathway_completion_wo_orphan_by_strain.tsv"

    with open(filename, "w") as output_file:
        logger.info(f"Writing completion values to {filename}.")
        header = "\t".join(
            [
                "pathway",
                "pathway name",
                "nb reactions",
                "families",
                "persistent %",
                "shell %",
                "cloud %",
                "modules (reaction cov.)",
                "global completion",
                "max completion",
            ]
            + [strain for strain in strain_to_families.keys()]
        )
        output_file.write(header + "\n")
        for pathway, pathway_reactions in pathways_to_reactions.items():
            if len(pathway_reactions) == 0:
                continue
            global_completion = len(
                pathway_reactions.intersection(pangenome_reactions)
            ) / len(pathway_reactions)

            families = set()
            for reaction in pathway_reactions:
                if reaction in reaction_to_families:
                    families.update(reaction_to_families[reaction])
            if len(families) == 0:
                continue
            persistent_percentage = len(families.intersection(partition_to_families["persistent"])) * 100 / len(families)
            shell_percentage = len(families.intersection(partition_to_families["shell"])) * 100 / len(families)
            cloud_percentage = len(families.intersection(partition_to_families["cloud"])) * 100 / len(families)

            modules_to_str = ""
            for module in module_to_families:
                nb_common_fam = len(families.intersection(module_to_families[module]))
                if nb_common_fam >= 2:
                    nb_reaction_in_module = 0
                    nb_reaction_with_families = 0
                    for reaction in pathway_reactions:
                        if reaction in reaction_to_families:
                            nb_reaction_with_families += 1
                            if len(reaction_to_families[reaction].intersection(module_to_families[module])) >= 1:
                                nb_reaction_in_module += 1
                    modules_to_str += (
                        module
                        + " ("
                        + str(nb_reaction_in_module / nb_reaction_with_families)
                        + ") ,"
                    )
            modules_to_str = modules_to_str.rstrip(" ,")

            record = [
                    pathway,
                    len(pathway_reactions),
                    ",".join(families),
                    persistent_percentage,
                    shell_percentage,
                    cloud_percentage,
                    modules_to_str,
                    global_completion,
                ]


            strain_completion = dict()
            for strain in strain_to_families:
                nbreactions_in_strain = 0
                for reaction in pathway_reactions:
                    # Check if the reaction is catalyzed by a protein among the gene families of the strains
                    if reaction in reaction_to_families and len(reaction_to_families[reaction].intersection(strain_to_families[strain])) >= 1:
                        nbreactions_in_strain += 1
                strain_completion[strain] = nbreactions_in_strain / len(pathway_reactions)
            max_completion = max(strain_completion.values())
            record.append(max_completion)
            for strain in strain_completion:
                record.append(strain_completion[strain])
            output_file.write("\t".join(list(map(str, record))) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ppanggolin-dir",
        help="PPanGGOLiN directory with files gene_presence_absence.Rtab, functional_modules.tsv, partitions/persistent.txt, partitions/shell.txt and partitions/cloud.txt",
        required=True,
    )
    parser.add_argument(
        "--pathways",
        help="An input list of metabolic pathways to compute completion",
        required=True,
    )
    parser.add_argument("--prefix", help="Output filenames prefix", required=True)
    parser.add_argument(
        "--families-reaction",
        help="A TSV mapping gene families to reaction identifiers",
        required=True,
    )
    parser.add_argument(
        "-c", "--config", help="Path to a config file to override default configuration"
    )

    args = parser.parse_args()
    reaction_presence_absence_by_strain_filename: Path = Path(
        f"{args.prefix}_reaction_presence_absence.Rtab"
    )

    # Input PPanGGOLiN files
    ppanggolin_dir = Path(args.ppanggolin_dir)
    rtab_filename = ppanggolin_dir / "gene_presence_absence.Rtab"
    module_filename = ppanggolin_dir / "functional_modules.tsv"
    persistent_filename = ppanggolin_dir / "partitions" / "persistent.txt"
    shell_filename = ppanggolin_dir / "partitions" / "shell.txt"
    cloud_filename = ppanggolin_dir / "partitions" / "cloud.txt"

    # Strains to families
    strain_to_families: Dict[str, Set[str]] = read_pangenome_rtab(rtab_filename)
    partition_to_families: Dict[Literal["persistent", "shell", "cloud"], Set[str]] = (
        read_partition_files(persistent_filename, shell_filename, cloud_filename)
    )
    # Modules to families
    module_to_families = read_functional_module_file(module_filename)
    # Read a mapping from gene families to reaction identifiers
    families_to_reactions: Dict[str, Set[str]] = read_mapping(args.families_reaction)
    reaction_to_families: Dict[str, Set[str]] = reverse_mapping(families_to_reactions)
    pangenome_reactions: Set[str] = {reaction for reaction in reaction_to_families}

    """write_reaction_presence_absence_by_strain(
        reaction_presence_absence_by_strain_filename,
        reaction_to_families,
        strain_to_families,
        partition_to_families,
    )"""
    pathways: List[str] = read_list(args.pathways)
    kb = select_kb(config["reference"]["source"])
    write_pathway_completion_by_strain(
        args.prefix,
        pathways,
        strain_to_families,
        partition_to_families,
        reaction_to_families,
        module_to_families,
        pangenome_reactions,
        use_orphan=True,
        kb=kb,
    )
    write_pathway_completion_by_strain(
        args.prefix,
        pathways,
        strain_to_families,
        partition_to_families,
        reaction_to_families,
        module_to_families,
        pangenome_reactions,
        use_orphan=False,
        kb=kb,
    )


if __name__ == "__main__":
    main()
