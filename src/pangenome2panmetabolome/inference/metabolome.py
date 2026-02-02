"""
Inference rules to infer the presence
of a pathway in an organism
based on its annotated (pan)genome.


Try to reproduce the PathoLogic algorithm
as presented in section 7.3.7 of Pathway-Tools user guide UserGuide.pdf (page 187 for v29.5).
"""

import logging
import argparse
from typing import Optional

from ..utils import read_list, write_output
from ..config import config
from ..io.knowledge_base import KnowledgeBase, select_kb
from ..taxonomy import NCBITaxonomyTree

logger = logging.getLogger("pangenome2panmetabolome:inference")
logger.setLevel(logging.INFO)


class PathwayInference:
    PATHWAY_COMPLETION_THRESHOLD: float = float(
        config["inference"]["pathway_completion_threshold"]
    )
    PATHWAY_SCORE_CUTOFF: float = float(config["inference"]["pathway_score_cutoff"])

    def __init__(
        self,
        template: KnowledgeBase,
        reactome: set[str],
        taxon_id: int,
        record_reason: bool = False,
    ):
        self.template: KnowledgeBase = template
        self.reactome: set[str] = reactome
        self.taxonomy: NCBITaxonomyTree = NCBITaxonomyTree(
            config["reference"]["ncbi_taxonomy"]
        )
        self.taxon_id: int = taxon_id
        self.pathways: list[str] = self.template.pathways()
        self.taxonomic_range_belonging_precompute(self.pathways)
        # Trace the decision algorithm if `record_reason` option is set to true.
        self.decision_reason: dict[str, str] = None
        self.record_reason = record_reason
        if self.record_reason:
            self.decision_reason = {}

    def reason(self, pathway_id: str) -> Optional[str]:
        """
        Get a report on why the pathway is retained or not.
        """
        if (
            self.record_reason
            and self.decision_reason is not None
            and pathway_id in self.decision_reason
        ):
            return self.decision_reason[pathway_id]

    def amend_reason(self, pathway_id: str, details: str):
        # Initialize the reason, if no reason exist for the moment
        if pathway_id not in self.decision_reason:
            self.decision_reason[pathway_id] = ""
        # Amend the existing reason
        self.decision_reason[pathway_id] += details

    def pathway_completion(self, template_reactions: list[str]) -> float:
        """
        Compute the completion of a pathway.

        Parameters
        ----------

            template_reactions: the list of reactions in the pathway (excluding orphan reactions and spontaneous reactions)
        """
        count: int = sum(
            1 if reaction in self.reactome else 0 for reaction in template_reactions
        )
        completion: float = count / len(template_reactions)
        return completion

    def pathway_completion_criterion(self, pathway_id: str) -> bool:
        """
        Is the pathway completion criterion satisfied?
        The pathway completion criterion is satisfied whenever the pathway completion
        is greater or equal than the PATHWAY_COMPLETION_THRESHOLD
        """
        non_spontaneous_non_orphan_reactions_of_pathway: list[str] = (
            self.template.non_orphan_non_spontaneous_reactions_of_pathway(pathway_id)
        )
        completion: float = self.pathway_completion(
            non_spontaneous_non_orphan_reactions_of_pathway
        )
        return completion >= PathwayInference.PATHWAY_COMPLETION_THRESHOLD

    def reaction_score(
        self, reaction_id: str, pathway_id: str, pathway_key_reactions: set[str]
    ) -> float:
        r"""
        Compute the reaction score $RS$ following the definition of PathoLogic.

        $$
        RS = P + U + K
        $$

        $P$: presence score

        $$
        P = \begin{cases}
        0.2 & \text{if enzyme catalyzing the reaction is present} \\
        0 & \text{otherwise}
        \end{cases}
        $$

        $U \in [0; 0.6]$: uniqueness score

        U increases if the reaction is present in a single pathway,
          and decreases with the number of pathway it is involved in.

        $K$: key reaction score

        $$
        K = \begin{cases}
        0.5 & \text{if the reaction is a key reaction of the pathway} \\
        0 & \text{otherwise}
        \end{cases}
        $$

        """
        presence_score = self.presence_score(reaction_id)
        uniqueness_score = self.uniqueness_score(reaction_id)
        key_reaction_score = self.key_reaction_score(
            reaction_id, pathway_id, pathway_key_reactions
        )
        score = presence_score + uniqueness_score + key_reaction_score
        if self.record_reason:
            self.amend_reason(
                pathway_id,
                f"- with reaction {reaction_id}: PresenceScore = {presence_score}, UniquenessScore = {uniqueness_score}, KeyReactionScore = {key_reaction_score}\n",
            )
        return score

    def presence_score(self, reaction_id: str) -> float:
        if reaction_id in self.reactome:
            return 0.2
        else:
            return 0

    def uniqueness_score(self, reaction_id: str) -> float:
        r"""
        Unique score $U$ ranges between 0.6 for reaction belonging to a single pathway, and 0 for reaction present in multiple pathways.

        Let $n_r$ be the number of pathway the reaction belongs to.
        Let $\max_r \{n_r\}$ be the maximum number of pathway a reaction can belong to.
        We pose:

        $$
        U = \max \{ 0,  f(n_r) \}
        $$

        where $f(n_r) = (1-x) * a + 0.6$ and where $a = (max_score - min_score) / (1 - max_pathway_uniqueness)$

        It is unclear how PathoLogic compute this score.

        TODO: try to dive deeper in how PathoLogic deals with this computation.
        """
        if reaction_id not in self.reactome:
            return 0
        pathways_having_reaction: int = len(
            self.template.pathways_with_reaction(reaction_id)
        )
        if pathways_having_reaction > 2:
            return 0
        else:
            return 0.6

    def key_reaction_score(
        self, reaction: str, pathway: str, pathway_key_reactions: set[str]
    ) -> float:
        if reaction in pathway_key_reactions:
            return 0.5
        else:
            return 0

    def pathway_score(
        self,
        pathway_id: str,
        non_orphan_non_spontaneous_pathway_reactions: list[str],
        pathway_key_reactions: set[str],
    ) -> float:
        r"""
        Compute the pathway score $PS$ according to the PathoLogic heuristic.

        $RS(r)$ is the reaction score of reaction $r$ of the pathway, where $r$ in both non-spontaneous and non-orphan,

        $T2$ is a boost if the organism belongs the the taxonomic range of the pathway.
        $T2 \geq 1$.
        TODO: deal with T1 and T2.

        $$
        PS = \frac{\sum_{r \in R}RS(r)}{|R|} * T1 * T2
        $$
        """
        n = len(non_orphan_non_spontaneous_pathway_reactions)
        if n == 0:
            return 0
        score: float = (
            sum(
                self.reaction_score(reaction_id, pathway_id, pathway_key_reactions)
                for reaction_id in non_orphan_non_spontaneous_pathway_reactions
            )
            + self.taxonomic_range_boost(pathway_id) / n
        )
        if self.record_reason:
            self.amend_reason(
                pathway_id,
                f"PathwayScore = {score} for n = {n} non-spontaneous, non-orphan reactions.\n",
            )
        return score

    def taxonomic_range_boost(self, pathway: str):
        if self.pathway_in_taxonomic_range[pathway]:
            return 0.1
        else:
            return 0

    def taxonomic_range_belonging_precompute(self, pathways: list[str]):
        """
        Precompute the Boolean saying whether the target organism
        is under the given NCBI-Taxonomy ID for the target organism.
        """
        self.pathway_in_taxonomic_range: dict[str, bool] = {
            pathway: self.taxonomy.is_child_of_parent_tax_id(
                self.template.pathway_taxonomic_range(pathway), self.taxon_id
            )
            for pathway in pathways
        }

    def pathway_presence_decision(
        self,
        pathway_id: str,
        pathway_scores: dict[str, float],
        non_orphan_non_spontaneous_pathway_reactions: list[str],
    ) -> bool:
        """
        REJECT P if P is a transport, signaling, or
        synthetic (engineered) pathway

        REJECT P if P is an electron transport pathway
        AND P lacks enzymes for any reaction

        INCLUDE P if P has all reactions present
        (meaning an enzyme is present for each reaction)
        AND if P is outside its taxonomic range, P
        contains more than 3 reactions

        REJECT P if P is outside its taxonomic range

        REJECT P if P is missing enzymes for all key
        reactions of P

        REJECT P if the score of P is significantly less
        than the score of a variant pathway of P

        INCLUDE P if the score of P exceeds the
        threshold PATHWAY-PREDICTION-CUTOFF efined in ptools-init.dat or specified in Automated Build dialog

        Default decision: REJECT

        Returns True if pathway is predicted, False otherwise.
        """
        # pathway_type: str = self.template.pathway_type(pathway_id) # TODO
        # pathway_type should be in ["transport", "signaling", "pathway", "synthetic"]

        # REJECT P if P is a transport, signaling, or synthetic (engineered) pathway
        # if pathway_type in ["transport", "signaling", "synthetic"]:
        #    return False
        # TODO

        # REJECT P if P is an electron transport pathway AND P lacks enzymes for any reaction
        # TODO

        # INCLUDE P if P has all reactions present (meaning an enzyme is present for each reaction) AND if P is outside its taxonomic range, P contains more than 3 reactions
        all_reactions_are_present: bool = len(
            non_orphan_non_spontaneous_pathway_reactions
        ) >= 1 and all(
            reaction in self.reactome
            for reaction in non_orphan_non_spontaneous_pathway_reactions
        )
        in_taxonomic_range: bool = (
            self.pathway_in_taxonomic_range[pathway_id] is None
            or self.pathway_in_taxonomic_range[pathway_id]
        )
        if all_reactions_are_present:
            if in_taxonomic_range:
                if self.record_reason:
                    self.amend_reason(
                        pathway_id,
                        "ACCEPT: all reactions are present and in taxonomic range.",
                    )
                return True
            elif len(non_orphan_non_spontaneous_pathway_reactions) >= 3:
                if self.record_reason:
                    self.amend_reason(
                        pathway_id,
                        details="ACCEPT: all reactions are present, not in taxonomic range, but at least three non-orphan, non-spontaneous reactions.",
                    )
                return True

        # REJECT P if P is outside its taxonomic range
        if not in_taxonomic_range:
            if self.record_reason:
                self.amend_reason(
                    pathway_id,
                    "REJECT: not all reaction present and not in taxonomic range.",
                )
            return False

        # REJECT P if P is missing enzymes for all key reactions of P
        key_reactions = self.template.key_reactions_of_pathway(pathway_id)
        if key_reactions is not None and all(
            key_reaction not in self.reactome for key_reaction in key_reactions
        ):
            if self.record_reason:
                self.amend_reason(pathway_id, "REJECT: all key reactions are missing.")
            return False

        # REJECT P if the score of P is significantly less than the score of a variant pathway of P
        variant_pathways = self.template.variants_of_pathway(pathway_id)
        if variant_pathways is not None:
            for variant_pathway in variant_pathways:
                if self.consider_pathway_score_to_be_greater(
                    pathway_scores[variant_pathway],
                    pathway_scores[pathway_id],
                    list(pathway_scores.values()),
                ):
                    if self.record_reason:
                        self.amend_reason(
                            pathway_id,
                            f"REJECT: variant pathway {pathway_id} has significantly higher pathway score.",
                        )
                    return False

        # INCLUDE P if the score of P exceeds the threshold PATHWAY-PREDICTION-CUTOFF
        if pathway_scores[pathway_id] > PathwayInference.PATHWAY_SCORE_CUTOFF:
            if self.record_reason:
                self.amend_reason(
                    pathway_id, "ACCEPT: pathway score exceeds the minimum value."
                )
            return True
        # Otherwise, by default, reject the pathway.
        if self.record_reason:
            self.amend_reason(
                pathway_id,
                "REJECT: no applicable case to reject or accept the pathway; reject by default.",
            )
        return False

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
            pathway: self.template.non_orphan_non_spontaneous_reactions_of_pathway(
                pathway
            )
            for pathway in self.pathways
        }
        pathway_key_reactions: dict[str, list[str]] = {
            pathway: self.template.key_reactions_of_pathway(pathway)
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
                f.writelines(["DECISION:" + pathway_id + "\n", decision_reason])


def infer_metabolome(reactome: set[str], taxon: int, reason: str) -> set[str]:
    kb = select_kb(config["reference"]["source"])
    inference = PathwayInference(kb, reactome, taxon, reason is not None)
    metabolome: set[str] = inference.inferred_pathways()
    if reason is not None:
        inference.dump_reason(reason)
    return metabolome


def main():
    logger.setLevel(logging.DEBUG)
    parser = argparse.ArgumentParser()
    parser.add_argument("reactome", help="Inpur list of reaction identifiers.")
    parser.add_argument("-t", "--taxon", help="NCBI Taxonomy tax-id", type=int)
    parser.add_argument(
        "-o", "--output", help="Output list of metabolic pathways.", required=True
    )
    parser.add_argument(
        "--reason",
        help="Give me a reason for your choice.",
        required=False,
        default=None,
    )

    args = parser.parse_args()
    reactome: set[str] = set(read_list(args.reactome))
    metabolome: set[str] = infer_metabolome(reactome, int(args.taxon), args.reason)
    write_output(args.output, metabolome)


if __name__ == "__main__":
    main()
