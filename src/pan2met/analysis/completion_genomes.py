"""
Completion matrix for a folder with files with sets of reactions
"""

import argparse
import csv
import logging
from pathlib import Path

from ..config import default_config, override_config
from ..io.knowledge_base import KnowledgeBase, select_kb
from ..utils import read_list

logger = logging.getLogger("pan2met:analysis:completion_genomes")
logging.basicConfig(level=logging.DEBUG)


def write_pathway_completion_by_strain(
    prefix: str,
    pathways: list[str],
    strain_to_reactions: dict[str, set[str]],
    use_orphan: bool,
    kb: KnowledgeBase,
):
    logger.info(
        f"Computing completion values for {len(pathways)} pathways in {len(strain_to_reactions.keys())} strains, {'with orphan reactions' if use_orphan else 'ignoring orphan reactions'}."
    )
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
        logger.info(f"Writing completion values to {filename}")
        writer = csv.writer(output_file, delimiter="\t")
        sorted_strains = sorted(strain_to_reactions.keys())
        header = [
            "pathway",
            "pathway name",
            "nb reactions",
            "max completion",
        ] + sorted_strains
        writer.writerow(header)
        for pathway, pathway_reactions in pathways_to_reactions.items():
            if len(pathway_reactions) == 0:
                logger.warning(f"{pathway} has no reactions.")
                continue

            pathway_name = ""
            try:
                pathway_name = kb.name_of_pathway(pathway)
            except NotImplementedError:
                pass

            record = [
                pathway,
                pathway_name,
                len(pathway_reactions),
            ]
            strain_completion = {}
            for strain, strain_reactions in strain_to_reactions.items():
                strain_completion[strain] = len(
                    strain_reactions.intersection(pathway_reactions)
                ) / len(pathway_reactions)
            max_completion = max(strain_completion.values())
            if max_completion == 0:
                logger.critical(
                    f"pathway {pathway} has a null completion in all strains."
                )
            record.append(max_completion)
            for strain in sorted_strains:
                record.append(strain_completion[strain])
            writer.writerow(record)


def read_folder_sets(folder: Path) -> dict[str, set[str]]:
    folder_sets: dict[str, set[str]] = {}
    for file in folder.iterdir():
        id = file.stem
        file_entries = set(read_list(file))
        folder_sets[id] = file_entries
    return folder_sets


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--genome-reactions",
        help="A folder containing a set of files whose basename is the identifier of the genome and content is a list of reaction identifiers.",
        required=True,
    )
    parser.add_argument(
        "--pathways",
        help="An input list of metabolic pathways to compute completion (e.g., the union of all predicted pathways for the genomes referenced in --genome-reactions folder.)",
        required=True,
    )
    parser.add_argument("--prefix", help="Output filenames prefix", required=True)
    parser.add_argument(
        "-c", "--config", help="Path to a config file to override default configuration"
    )
    args = parser.parse_args()

    if args.config:
        config = override_config(args.config)
    else:
        config = default_config

    pathways: list[str] = read_list(args.pathways)
    strain_to_reactions: dict[str, set[str]] = read_folder_sets(
        Path(args.genome_reactions)
    )

    kb = select_kb(config)
    write_pathway_completion_by_strain(
        args.prefix,
        pathways,
        strain_to_reactions,
        use_orphan=True,
        kb=kb,
    )
    write_pathway_completion_by_strain(
        args.prefix,
        pathways,
        strain_to_reactions,
        use_orphan=False,
        kb=kb,
    )


if __name__ == "__main__":
    main()
