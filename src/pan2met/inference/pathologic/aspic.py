"""
An implementation of a PathoLogic-like pathway inference algorithm,
relying on an Answer Set Programming approach with clingo.

"""

from pathlib import Path
from typing import Dict
import importlib.resources
import configparser
import argparse
import logging

import clyngor

import pan2met
from ...utils import unquote
from ...utils import read_list, write_output, set_logging_level
from ...io.knowledge_base import KnowledgeBase
from ...io.knowledge_base import select_kb
from ...config import config
from ...asp import kb_as_asp as kb_as_asp
from .generic import PathwayInference

logger = logging.getLogger("pan2met:inference:pathologic:aspic")


class AspicPathwayInference(PathwayInference):
    def __init__(
        self,
        kb: KnowledgeBase,
        reactome: set[str],
        taxon_id: int,
        reference_kb_asp: Path,
    ):
        super().__init__(kb, reactome, taxon_id)
        self.reference_kb_asp: Path = reference_kb_asp

    def prepare_clingo_inline_program(self) -> str:
        """
        Prepare reactome/1, in_taxonomic_range/1. inline ASP atoms.
        """
        reactome_atoms = kb_as_asp.target_organism_reactions_as_asp(self.reactome)
        taxonomic_range_atoms = kb_as_asp.pathways_in_taxonomic_range_as_asp(
            self.taxonomic_range_belonging
        )
        return "\n".join([reactome_atoms, taxonomic_range_atoms])

    def solve_asp_problem(
        self, inline_program: str, reference_kb_asp: Path
    ) -> Dict[str, list[str]]:
        """
        Use clingo with clyngor to solve the Answer Set Programming problem to infer pathways.

        :param inline_program: Concatenation of the reactome/1. atoms of the target organism, and the in_taxonomic_range/1. atoms, that depends on the taxon of the target organism.
        :param reference_kb_asp: The path to the ASP representation of the reference knowledge base.
        :return: a list of pathway identifier that are inferred.
        """
        with importlib.resources.path(
            pan2met, "asp/rules/pathologic_like.lp"
        ) as pathologic_like_asp_rule_path:
            answers: clyngor.Answers = clyngor.solve(
                [
                    pathologic_like_asp_rule_path,
                    reference_kb_asp,
                ],
                inline=inline_program,
                nb_model=1,
            ).by_predicate
        answers = tuple(answers)
        results = answers[0]
        # logger.debug("clingo results: %s", results["rule"])
        return {
            "include": [unquote(pathway) for (pathway,) in results["included_pathway"]],
            "reject": [unquote(pathway) for (pathway,) in results["rejected_pathway"]],
            "undecided": [
                unquote(pathway) for (pathway,) in results["undecided_pathway"]
            ],
        }

    def infer_pathways(self) -> list[str]:
        inline_program = self.prepare_clingo_inline_program()
        inferrence_results = self.solve_asp_problem(
            inline_program, self.reference_kb_asp
        )
        undecided_pathways = inferrence_results["undecided"]

        inferred_pathways = inferrence_results["include"]

        if len(undecided_pathways) > 0:
            logger.info(
                f"ASP decision rules left {len(undecided_pathways)} pathway(s) undecided. These pathways will pass the pathway score heuristics decision rule test."
            )

            for pathway in undecided_pathways:
                self.amend_reason(
                    pathway,
                    "ASP decision rules could not decide on this pathway, applying pathway score heuristics decision rule.",
                )
                # TODO: deal with the comparison of variants.
                non_orphan_non_spontaneous_reactions = (
                    self.kb.non_orphan_non_spontaneous_reactions_of_pathway(pathway)
                )
                if len(non_orphan_non_spontaneous_reactions) == 0:
                    logger.warning(
                        f"Pathway {pathway} has no non-orphan non-spontaneous reaction, cannot compute pathway score, skipping pathway score heuristics decision rule, and accepting the pathway by default."
                    )
                    inferred_pathways.append(pathway)
                    if self.record_reason:
                        self.amend_reason(
                            pathway,
                            "ACCEPT: pathway has no non-orphan non-spontaneous reaction, cannot compute pathway score, accepting by default.",
                        )
                    continue
                pathways_keys_reactions = self.kb.key_reactions_of_pathway(pathway)
                pathway_score = self.pathway_score(
                    pathway,
                    non_orphan_non_spontaneous_reactions,
                    pathways_keys_reactions,
                )
                if pathway_score >= PathwayInference.PATHWAY_SCORE_THRESHOLD:
                    inferred_pathways.append(pathway)
                    if self.record_reason:
                        self.amend_reason(
                            pathway,
                            f"ACCEPT: pathway score {pathway_score:.2f} exceeds the minimum value.",
                        )
        return inferred_pathways


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--kb-asp",
        help="filename of the ASP representation of the metabolic knowledge base",
        required=True,
    )
    parser.add_argument(
        "--reactome",
        help="List of reaction identifiers, one per line",
        required=True,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Log level (-v: ERROR, -vv: WARNING, -vvv: INFO, -vvvv: DEBUG)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="The path of the output file listing the pathway identifiers infered",
        required=True,
    )
    parser.add_argument(
        "--version",
        action="version",
        help="Show the version of pan2met and exit",
        version=f"pan2met v{pan2met.__version__}",
    )
    parser.add_argument(
        "--taxon",
        help="The NCBI taxon id of the target organism",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--config",
        help="The path of the configuration file to use to override defaults.",
    )
    args = parser.parse_args()

    set_logging_level(args.verbose)

    if args.config is not None:
        # Update config globally overriding default_config with keys from given config filename
        logging.info(f"Overriding default configuration with {args.config}")
        global config
        config_override = configparser.ConfigParser()
        config_override.read([args.config])
        config.update(config_override)
    logger.info(f"Using config: {(dict(config))}")
    kb = select_kb(config["reference"]["source"])
    reactome = set(read_list(args.reactome))
    pathway_inference = AspicPathwayInference(
        kb, reactome, taxon_id=int(args.taxon), reference_kb_asp=Path(args.kb_asp)
    )
    logger.info("Inferring pathways with ASP")
    pathways = pathway_inference.infer_pathways()
    logger.info(f"Inferred {len(pathways)} pathways")
    write_output(args.output, pathways)
    logger.info("Done writing output to %s", args.output)


if __name__ == "__main__":
    main()
