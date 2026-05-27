import logging
from typing import Optional, Dict

from ...io.knowledge_base import KnowledgeBase
from ...taxonomy import NCBITaxonomyTree
from ...config import config

logger = logging.getLogger("pan2met:inference:pathologic")


class PathwayInference:
    PATHWAY_COMPLETION_THRESHOLD: float = float(
        config["inference"]["pathway_completion_threshold"]
    )
    PATHWAY_SCORE_THRESHOLD: float = float(
        config["inference"]["pathway_score_threshold"]
    )
    AVAILABLE_RULES = [
        "taxonomic_range",
        "pathway_species",
        "pathway_uniqueness",
        "pathway_key_reaction",
        "pathway_ontology",
        "pathway_variant",
    ]
    RULES: dict[str, bool] = {key: True for key in AVAILABLE_RULES}
    RULES.update(
        {
            key: config["inference"]["rules"].getboolean(key)
            for key in AVAILABLE_RULES
            if (
                "inference" in config
                and "rules" in config["inference"]
                and key in config["inference"]["rules"][key]
            )
        }
    )
    # FIXME: rule config does not seem to be correctly taken into account.

    def __init__(
        self,
        kb: KnowledgeBase,
        reactome: set[str],
        taxon_id: int,
        record_reason: bool = False,
    ):
        self.kb: KnowledgeBase = kb
        self.reactome: set[str] = reactome
        self.pathways: list[str] = self.kb.pathways()
        # if PathwayInference.RULES["taxonomic_range"]: # TODO: re-enable the possibility to disable the taxonomic range rule, to support KEGG.
        self.taxonomy: NCBITaxonomyTree = NCBITaxonomyTree(
            config["reference"]["ncbi_taxonomy"]
        )
        self.taxon_id: int = taxon_id
        self.pathway_in_taxonomic_range: dict[str, bool] = (
            self.taxonomic_range_belonging_precompute(self.pathways)
        )
        # Trace the decision algorithm if `record_reason` option is set to true.
        self.decision_reason: dict[str, str] = {}
        self.record_reason = record_reason

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
        """
        Add a string fragment to the reason dictionnary for a pathway.
        """
        # Initialize the reason, if no reason exist for the moment
        if pathway_id not in self.decision_reason:
            self.decision_reason[pathway_id] = ""
        # Amend the existing reason
        self.decision_reason[pathway_id] += details

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

        :param reaction_id: the identifier of the reaction
        :param pathway_id: the identifier of the pathway
        :param pathway_key_reactions: the set of key reaction of the pathway
        :return: the reaction score
        """
        score = 0
        # Add base presence score
        # IDEA: it might be usefull to give more weight to the presence,
        #  when the other score parts are disabled.
        presence_score = self.presence_score(reaction_id)
        score += presence_score

        # Add uniqueness score
        if PathwayInference.RULES["pathway_uniqueness"]:
            uniqueness_score = self.uniqueness_score(reaction_id)
            score += uniqueness_score
        else:
            uniqueness_score = None

        # Add key reaction score
        if PathwayInference.RULES["pathway_key_reaction"]:
            key_reaction_score = self.key_reaction_score(
                reaction_id, pathway_key_reactions
            )
            score += key_reaction_score
        else:
            key_reaction_score = None

        if self.record_reason:
            self.amend_reason(
                pathway_id,
                f"- with reaction {reaction_id}: PresenceScore = {presence_score}, UniquenessScore = {uniqueness_score}, KeyReactionScore = {key_reaction_score}\n",
            )
        return score

    def presence_score(self, reaction_id: str) -> float:
        """
        Pathway presence score of a reaction.
        When reaction is in known to be in the reactome of an organism,
        score is 0.2, otherwise it is 0.

        :param reaction_id: the identifier of the reaction
        :return: the presence score of the reaction
        """
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

        :param reaction_id: the identifier of the reaction
        """
        if reaction_id not in self.reactome:
            return 0
        pathways_having_reaction: int = len(self.kb.pathways_with_reaction(reaction_id))
        if (
            pathways_having_reaction > 5
        ):  # defines "A large number of MetaCyc pathway" (UserGuide p. 188)
            return 0
        else:
            return 0.6

    def key_reaction_score(
        self, reaction: str, pathway_key_reactions: set[str]
    ) -> float:
        """
        Key reaction score $K$ is 0.5 if the reaction is a key reaction of the pathway, and 0 otherwise.

        :param reaction: the identifier of the reaction
        :param pathway_key_reactions: the set of key reaction of the pathway
        :return: the key reaction score of the reaction for the given pathway
        """
        if reaction in self.reactome and reaction in pathway_key_reactions:
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

        If no non-spontaneous and non orphan reaction exists in the pathway, the pathway score is set to 0.

        $$
        PS = \frac{\sum_{r \in R}RS(r)}{|R|} * T1 * T2
        $$

        :param pathway_id: the identifier of the pathway
        :param non_orphan_non_spontaneous_pathway_reactions: the list of reactions of the pathway that are both non-spontaneous and non-orphan.
        :param pathway_key_reactions: the set of key reaction of the pathway
        """
        n = len(non_orphan_non_spontaneous_pathway_reactions)
        if n == 0:
            return 0
        score: float = (
            sum(
                self.reaction_score(reaction_id, pathway_id, pathway_key_reactions)
                for reaction_id in non_orphan_non_spontaneous_pathway_reactions
            )
            / n
        )
        # Species evidence neighborhood
        if PathwayInference.RULES["pathway_species"]:
            pathway_neighbor_species_boost_score = self.taxonomic_neighborhood_boost(
                pathway_id
            )  # T1
            score *= pathway_neighbor_species_boost_score
        else:
            pathway_neighbor_species_boost_score = None
        # Taxonomic range boost
        if PathwayInference.RULES["taxonomic_range"]:
            pathway_taxonomic_range_boost_score = self.taxonomic_range_boost(pathway_id)
            score *= pathway_taxonomic_range_boost_score
        else:
            pathway_taxonomic_range_boost_score = None

        if self.record_reason:
            self.amend_reason(
                pathway_id,
                f"PathwayScore = {score} for n = {n} non-spontaneous, non-orphan reactions with taxonomic range boost T1 = {pathway_taxonomic_range_boost_score} and neighbor species boost T2 = {pathway_neighbor_species_boost_score}.\n",
            )
        return score

    def taxonomic_neighborhood_boost(self, pathway_id: str) -> float:
        """
        "T1 is highest if the two organisms are the same strain, weaker if the same species, weaker still if the same genus." (UserGuide.pdf v29.5 p. 187)

        We do not know the exact weight of this boost in the PathoLogic implementation.
        Hence the choice made below is arbitrary.

        $$
        T1 >= 1.
        $$

        :param pathway_id: the identifier of the pathway
        """
        # FIXME: Warning: in metabiantes, the strain is not taken into account
        for evidence_species_id in self.kb.species_evidence_of_pathway(pathway_id):
            if self.taxonomy.is_under_same_species(self.taxon_id, evidence_species_id):
                return 1.2  # config["inference"]["weights"]["neighbor_species_boost"]
        return 1

    def taxonomic_range_boost(self, pathway_id: str) -> float:
        """
        "PS receives an additional boost T2 if the subject organism is within the expected taxonomic range
        of the pathway as designated within MetaCyc (e.g., if the subject organism is a plant and the path-
        way is designated as a plant pathway). If the subject organism is outside the expected taxonomic
        range of the pathway then T2 < 1." (UserGuide.pdf v29.5 p.187)

        $$
        T2 \geq 1
        $$

        :param pathway_id: the identifier of the pathway
        :return: the taxonomic range boost score of the pathway
        """
        if (
            pathway_id in self.pathway_in_taxonomic_range
            and self.pathway_in_taxonomic_range[pathway_id]
        ):
            return 1.2  # config["inference"]["weights"]["taxonomic_range_boost"]
        else:
            return 1

    def taxonomic_range_belonging_precompute(
        self, pathways: list[str]
    ) -> Dict[str, bool]:
        """
        Precompute the Boolean saying whether the target organism
        is under the given NCBI-Taxonomy ID for the target organism.
        :param pathways: the list of pathway identifier to precompute the taxonomic range belonging for.
        """

        pathway_in_taxonomic_range: dict[str, bool] = {}
        for pathway in pathways:
            taxonomic_range = self.kb.pathway_taxonomic_range(pathway)
            if taxonomic_range is not None:
                pathway_in_taxonomic_range[pathway] = (
                    self.taxonomy.is_child_of_parent_tax_id(
                        taxonomic_range, self.taxon_id
                    )
                )
        return pathway_in_taxonomic_range
