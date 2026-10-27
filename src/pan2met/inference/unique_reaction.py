import argparse
import csv

from ..config import default_config, override_config
from ..io.knowledge_base import KnowledgeBase, select_kb


def reaction_to_pathway_mapping(kb: KnowledgeBase) -> dict[str, set[str]]:
    """

    :param kb:
    :return: a mapping from reaction to pathway
    """
    mapping: dict[str, set[str]] = {}
    for pathway in kb.pathways():
        pathway_reactions = kb.reactions_of_pathway(pathway)
        for reaction in pathway_reactions:
            if reaction in mapping:
                mapping[reaction].add(pathway)
            else:
                mapping[reaction] = {pathway}
    return mapping


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--output",
        help="output path tsv with columns (1) key reaction (2) pathway",
    )
    parser.add_argument("-c", "--config", help="config override")
    return parser


def main():
    parser = argument_parser()
    args = parser.parse_args()

    if args.config:
        config = override_config(args.config, config=default_config)
    else:
        config = default_config

    kb = select_kb(config)

    reaction_to_pathway: dict[str, set[str]] = reaction_to_pathway_mapping(kb)

    with open(args.output, "w") as output_file:
        writer = csv.writer(output_file, delimiter="\t")
        for reaction, pathways in reaction_to_pathway.items():
            if len(pathways) == 1:
                writer.writerow([reaction, next(iter(pathways))])


if __name__ == "__main__":
    main()
