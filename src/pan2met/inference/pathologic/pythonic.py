"""
Inference rules to infer the presence
of a pathway in an organism
based on its annotated (pan)genome.


Try to reproduce the PathoLogic algorithm
as presented in section 7.3.7 of Pathway-Tools user guide UserGuide.pdf (page 187 for v29.5).


This script relies mainly on a Python-only implementation of the PathoLogic algorithm.
Refer to the `pan2met.inference.pathologic.aspic` module for an alternative implementation,
relying on an Answer Set Programming implementation of a PathoLogic-like inference algorithm.
"""

import logging
import argparse
import importlib.resources
import configparser
from typing import Optional
import io

import graph_tool as gt
import graph_tool.topology


from ...utils import set_logging_level
from ...utils import read_list, write_output
from ...io.knowledge_base import KnowledgeBase, select_kb
from .generic import PathwayInference
import pan2met.conf


logger = logging.getLogger("pan2met:inference")


class PythonicPathwayInference(PathwayInference):
    """
    PathoLogic-like pathway inference algorithm, with decision rules implemented in Python.
    Uses the pathway score heuristics defined in the `pan2met.inference.pathologic.generic.PathwayInference` parent class.
    """

    def __init__(
        self,
        kb: KnowledgeBase,
        reactome: set[str],
        taxon_id: Optional[int]=None,
        record_reason: bool = False,
        config=None
    ):
        super().__init__(kb, reactome, taxon_id=taxon_id, config=config)
        self.record_reason: bool = record_reason
        self.decision_reason: dict[str, str] = {}

    def pathway_presence_decision(
        self,
        pathway_id: str,
        pathway_scores: dict[str, float],
        non_orphan_non_spontaneous_pathway_reactions: list[str],
    ) -> bool:
        """
        Decide whether a pathway is inferred present or not.


        :param pathway_id: the identifier of the pathway
        :param pathway_scores: a dictionnary mapping a pathway identifier to its precomputed pathway score
        :param non_orphan_non_spontaneous_pathway_reactions: the list of reaction identifiers that are both
                - non-orphan (that is to say, they have a known associated enzyme in the knowledge base), and
                - non-spontaneous (that is to say, they do not occur spontaneously in physiological conditions)

        Rules
        -----

        An excerpt of pathway-tools user guide describing the
        decision rules for a pathway in pathologic algorithm:

        REJECT P if P is a transport, signaling, or
        synthetic (engineered) pathway (rule id: pathway_ontology)

        REJECT P if P is an electron transport pathway
        AND P lacks enzymes for any reaction (rule id: pathway_ontology)

        REJECT P if any key-non-reactions are present. (ignored / not implemented)

        INCLUDE P if P has all reactions present
        (meaning an enzyme is present for each reaction)

        AND if P is outside its taxonomic range, P
        contains more than 3 reactions (rule id: taxonomic_range)

        REJECT P if P is missing enzymes for any key reactions of P

        REJECT P if the score of P is significantly less
        than the score of a variant pathway of P (rule id: pathway_variant)

        TODO REJECT P if the score of P is slightly more than a preferred variant pathway of P
        (the "standard" glycolysis and TCA pathways are designated as preferred variants).  (not implemented)

        REJECT P if P is outside its taxonomic range (rule id: taxonomic_range)

        INCLUDE P if the score of P exceeds the
        threshold PATHWAY-PREDICTION-CUTOFF efined in ptools-init.dat or specified in Automated Build dialog

        Default decision: REJECT

        REJECT P if P is either:
        - A biosynthetic pathway missing enzymes for its final steps
        - A catabolic pathway missing enzymes for its initial steps
        - An energy pathway missing enzymes for more than half of its steps

        :return: True if the pathway is inferred present, False otherwise.
        """
        if self.RULES["pathway_ontology"]:
            pathway_ontology_parents: list[str] = (
                self.kb.ontology_parent_class_of_pathway(pathway_id)
            )

            # REJECT P if P is a transport, signaling, or synthetic (engineered) pathway
            # TODO deal with synthetic pathways
            for pathway_class in pathway_ontology_parents:
                if pathway_class == "Transport-Pathways":
                    self.amend_reason(pathway_id, "REJECT: is a transport pathway.")
                    return False
                if (
                    pathway_class == "Signaling-Pathways"
                ):  # INFO: it seems that in MetaCyc, no signaling pathway are reported though
                    self.amend_reason(pathway_id, "REJECT: is a signaling pathway.")
                    return False

        all_reactions_are_present: bool = len(
            non_orphan_non_spontaneous_pathway_reactions
        ) >= 1 and all(
            reaction in self.reactome
            for reaction in non_orphan_non_spontaneous_pathway_reactions
        )
        all_reactions_absent: bool = len(
            non_orphan_non_spontaneous_pathway_reactions
        ) >= 1 and all(
            reaction not in self.reactome
            for reaction in non_orphan_non_spontaneous_pathway_reactions
        )
        if all_reactions_absent:
            self.amend_reason(
                pathway_id, "REJECT: no known catalyzis at all for this pathway."
            )
            return False
        # REJECT P if P is an electron transport pathway AND P lacks enzymes for any reaction
        if self.RULES["pathway_ontology"]:
            if (
                "Electron-Transfer" in pathway_ontology_parents
                and not all_reactions_are_present
            ):
                self.amend_reason(
                    pathway_id,
                    "REJECT: is an electron transport pathway and lacks an enzyme for a reaction.",
                )
                return False

        # INCLUDE P if P has all reactions present (meaning an enzyme is present for each reaction) AND if P is outside its taxonomic range, P contains more than 3 reactions
        if self.RULES["taxonomic_range"]:
            in_taxonomic_range: bool = (
                pathway_id in self.pathway_in_taxonomic_range is not None
                and self.pathway_in_taxonomic_range[pathway_id]
            )

            if all_reactions_are_present:
                if in_taxonomic_range:
                    self.amend_reason(
                        pathway_id,
                        "ACCEPT: all reactions are present and in taxonomic range.",
                    )
                    return True
                elif len(non_orphan_non_spontaneous_pathway_reactions) >= 3:
                    self.amend_reason(
                        pathway_id,
                        details="ACCEPT: all reactions are present, not in taxonomic range, but at least three non-orphan, non-spontaneous reactions.",
                    )
                    return True
                else:
                    # REJECT, otherwise, when, P if P is outside its taxonomic range
                    self.amend_reason(
                        pathway_id,
                        "REJECT: not in taxonomic range, and no more than 3 catalyzed reactions.",
                    )
                    return False
        else:
            if all_reactions_are_present:
                self.amend_reason(
                    pathway_id,
                    details="ACCEPT: all reactions are present and we don't care about the taxonomic range.",
                )
            return True

        # REJECT P if P is missing enzymes for all key reactions of P
        if self.RULES["pathway_key_reaction"]:
            key_reactions = self.kb.key_reactions_of_pathway(pathway_id)
            for key_reaction in key_reactions:
                if key_reaction not in self.reactome:
                    self.amend_reason(
                        pathway_id,
                        f"REJECT: key reaction {key_reaction} not in reactome.",
                    )
                    return False

        # REJECT P if the score of P is significantly less than the score of a variant pathway of P
        if self.RULES["pathway_variant"]:
            variant_pathways = self.kb.variants_of_pathway(pathway_id)
            if variant_pathways is not None:
                for variant_pathway in variant_pathways:
                    if self.consider_pathway_score_to_be_greater(
                        pathway_scores[variant_pathway],
                        pathway_scores[pathway_id],
                        list(pathway_scores.values()),
                    ):
                        self.amend_reason(
                            pathway_id,
                            f"REJECT: variant pathway {pathway_id} has significantly higher pathway score.",
                        )
                        return False

        if self.RULES["pathway_ontology"]:
            one_reaction_graph_ordering = self.reaction_graph_topological_order(
                pathway_id
            )
            if len(one_reaction_graph_ordering) > 2:
                first_reaction_id = one_reaction_graph_ordering[0]
                last_reaction_id = one_reaction_graph_ordering[-1]
                # FIXME: It is possible that a pathway reaction graph has multiple dead ends, so we might miss the correct last reaction id.
                # REJECT P if P is a biosynthetic pathway missing enzymes for its final step
                if "Biosynthesis" in pathway_ontology_parents:
                    if last_reaction_id not in self.reactome:
                        self.amend_reason(
                            pathway_id,
                            f"REJECT: biosynthesis pathway {pathway_id}'s last reaction is not present in known catalyzed reactome.",
                        )
                        return False
                # REJECT P if P is a catabolic pathway missing enzymes for its initial step
                elif "Degradation" in pathway_ontology_parents:
                    if first_reaction_id not in self.reactome:
                        self.amend_reason(
                            pathway_id,
                            f"REJECT: catabolysis pathway {pathway_id}'s first reaction is not present in known catalyzed reactome.",
                        )
                        return False
                # REJECT P if P is an energy metabolism pathway and is missing more than half its catalyzed reactions (understood as both non orphan and non spontaneous reactions)
                elif "Energy-Metabolism" in pathway_ontology_parents:
                    catalyzed_reactions = [
                        reaction
                        for reaction in non_orphan_non_spontaneous_pathway_reactions
                        if reaction in self.reactome
                    ]
                    if (
                        len(catalyzed_reactions)
                        < len(non_orphan_non_spontaneous_pathway_reactions) / 2
                    ):
                        self.amend_reason(
                            pathway_id,
                            f"REJECT: energy metabolism pathway {pathway_id} is missing more than half its catalyzed reactions.",
                        )
                        return False

        # INCLUDE P if the score of P exceeds the threshold PATHWAY-PREDICTION-CUTOFF
        if pathway_scores[pathway_id] > self.PATHWAY_SCORE_THRESHOLD:
            self.amend_reason(
                pathway_id, "ACCEPT: pathway score exceeds the minimum value."
            )
            return True

        # Otherwise, by default, reject the pathway.
        self.amend_reason(
            pathway_id,
            "REJECT: no applicable case to reject or accept the pathway, so reject by default.",
        )
        return False

    def reaction_graph_topological_order(self, pathway):
        """
        Return the topological ordering of a pathway reaction graph.

        Assume that the pathway reaction graph is a directed acyclic graph.
        """
        reaction_order = self.kb.pathway_reaction_order(pathway)
        graph = gt.Graph(directed=True)
        vmap: dict[str, int] = {}
        # Create a graph_tool directed reaction graph
        for predecessor, successor in reaction_order:
            if predecessor in vmap:
                u = vmap[predecessor]
            else:
                u = graph.add_vertex()
                vmap[u] = predecessor
            if successor in vmap:
                v = vmap[successor]
            else:
                v = graph.add_vertex()
                vmap[v] = successor
            graph.add_edge(u, v)
        # Get the topological ordering
        sort = graph_tool.topology.topological_sort(graph)
        sort_named = [vmap[vertex] for vertex in sort]
        return sort_named

    def consider_pathway_score_to_be_greater(
        self,
        pathway_1_score: float,
        pathway_2_score: float,
        pathway_scores: list[float],
    ) -> bool:
        """
        We consider a pathway score $Score(P1)$ to be greater that the pathway score $Score(P2)$ if $Score(P1) > Score(P2) + sigma$

        FIXME: How is this implemented in PathoLogic ?
        It is unclear. It says it test for significance of the pathway score difference, but which test does it use?
        As the pathway score is probably not normally distributed, and we have a single value to compare It is unclear how to deal with this properly.
        It is an heuristic though, so we might simply discard our veleity to do proper statistics, at least for now.
        """
        sigma = 0.1
        return pathway_1_score > pathway_2_score + sigma

    def predict_all_pathway_presence(self) -> dict[str, bool]:
        """
        Predict the presence of all pathway of the template pathway set.
        """
        pathway_reactions: dict[str, list[str]] = {
            pathway: self.kb.non_orphan_non_spontaneous_reactions_of_pathway(pathway)
            for pathway in self.pathways
        }
        pathway_key_reactions: dict[str, list[str]] = {
            pathway: self.kb.key_reactions_of_pathway(pathway)
            for pathway in self.pathways
        }
        pathway_scores: dict[str, float] = {
            pathway: self.pathway_score(
                pathway, pathway_reactions[pathway], set(pathway_key_reactions[pathway])
            )
            for pathway in self.pathways
        }
        pathway_present: dict[str, bool] = {
            pathway: self.pathway_presence_decision(
                pathway, pathway_scores, pathway_reactions[pathway]
            )
            for pathway in self.pathways
        }

        return pathway_present

    def inferred_pathways(self) -> set[str]:
        """
        List all inferred pathways.
        """
        pathway_present: dict[str, bool] = self.predict_all_pathway_presence()
        inferred = {pathway for pathway, present in pathway_present.items() if present}
        return inferred

    def dump_reason(self, filename):
        assert self.record_reason and self.decision_reason is not None, (
            "To dump the reason of the decision, you should have recorded it before."
        )
        with open(filename, "w") as f:
            for pathway_id, decision_reason in self.decision_reason.items():
                f.writelines(["DECISION:" + pathway_id + "\n", decision_reason + "\n"])


def infer_metabolism(
    reactome: set[str], taxon: int, reason_filename: Optional[str] = None, config=None
) -> set[str]:
    kb = select_kb(config["reference"]["source"])
    inference = PythonicPathwayInference(
        kb, reactome, taxon, reason_filename is not None, config=config
    )
    metabolism: set[str] = inference.inferred_pathways()
    if reason_filename is not None:
        inference.dump_reason(reason_filename)
    return metabolism


def main():
    logger.setLevel(logging.DEBUG)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reactome", help="Input list of reaction identifiers.")
    parser.add_argument("-t", "--taxon", help="NCBI Taxonomy tax-id", type=int)
    parser.add_argument(
        "-o", "--output", help="Output list of metabolic pathways.", required=True
    )
    parser.add_argument(
        "--reason",
        help="Give me a reason for your choice and dump it to the file.",
        required=False,
        default=None,
    )
    parser.add_argument("-c", "--config", help="Path to a config file")
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Log level (-v: ERROR, -vv: WARNING, -vvv: INFO, -vvvv: DEBUG)",
    )
    args, remaining_argv = parser.parse_known_args()

    set_logging_level(args.verbose)

    default_config: configparser.ConfigParser = configparser.ConfigParser()
    default_config_str = importlib.resources.read_text(pan2met.conf, "default.ini")
    default_config.read_string(default_config_str)

    config = default_config

    # Update config globally overriding default_config with keys from given config filename
    if args.config:
        logging.info(f"Overriding default configuration with {args.config}")
        config_override = configparser.ConfigParser()
        config_override.read([args.config])
        config.update(config_override)

    # Log the used config
    with io.StringIO() as config_string_stream:
        config.write(config_string_stream)
        logger.info(f"Using config:\n{config_string_stream.getvalue()}")

    # Infer the metabolism
    reactome: set[str] = set(read_list(args.reactome))
    metabolism: set[str] = infer_metabolism(reactome, args.taxon, args.reason, config=config)
    write_output(args.output, metabolism)


if __name__ == "__main__":
    main()
