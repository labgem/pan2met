"""
Use the builtin rule kit engine to infer metabolic pathway
"""

from collections.abc import Iterable

from ...io.knowledge_base import KnowledgeBase, select_kb
from ...rulekit.rulekit import Rule, RuleEngine
from ...taxonomy import NCBITaxonomyTree
from .pathway_score import PathwayScore


class MetabolismRuleContext:
    def __init__(
        self,
        config,
        pathway_id: str,
        kb: KnowledgeBase,
        catalyzed_reactions: set[str],
        taxon_id: int | None,
        taxonomy: NCBITaxonomyTree | None,
    ):
        self.config = config
        self.pathway_id = pathway_id
        self.kb = kb
        self.pathway_non_spontaneous_reactions: set[str] = set(
            self.kb.non_spontaneous_reactions_of_pathway(pathway_id)
        )
        self.catalyzed_reactions = catalyzed_reactions
        self.taxonomy = taxonomy
        self.taxon_id = taxon_id

    @property
    def name(self):
        return self.pathway_id

    @property
    def all_reactions_catalyzed(self) -> bool:
        return (
            len(
                self.pathway_non_spontaneous_reactions.difference(
                    self.catalyzed_reactions
                )
            )
            == 0
        )

    @property
    def some_reaction_missing(self) -> bool:
        return (
            len(self.pathway_non_spontaneous_reactions) != 0
            and not self.all_reactions_catalyzed
        )

    @property
    def all_reactions_missing(self) -> bool:
        return (
            len(self.pathway_non_spontaneous_reactions) != 0
            and self.pathway_non_spontaneous_reactions.difference(
                self.catalyzed_reactions
            )
            == self.pathway_non_spontaneous_reactions
        )

    @property
    def in_taxonomic_range(self) -> bool | None:
        if self.taxonomy is not None and self.taxon_id is not None:
            taxonomic_range = self.kb.pathway_taxonomic_range(self.pathway_id)
            if taxonomic_range is not None:
                return self.taxonomy.is_child_of_parent_tax_id(
                    taxonomic_range, self.taxon_id
                )

    @property
    def pathway_ontology_parent_classes(self) -> list[str]:
        return self.kb.ontology_parent_class_of_pathway(self.pathway_id)

    def pathway_is_a(self, pathway_class: str) -> bool:
        for cls in self.pathway_ontology_parent_classes:
            if pathway_class == cls:
                return True
        return False

    def missing_last_reaction(self) -> bool:
        order = self.kb.reaction_graph_topological_order(self.pathway_id)
        if order is None:
            raise ValueError("pathway graph is not a DAG.")
        if len(order) == 0:
            raise ValueError("pathway reaction graph has no reaction.")
        last_reaction = order[-1]
        return (
            last_reaction not in self.catalyzed_reactions
            and last_reaction
            in self.kb.non_spontaneous_reactions_of_pathway(self.pathway_id)
        )

    def missing_first_reaction(self) -> bool:
        order = self.kb.reaction_graph_topological_order(self.pathway_id)
        if order is None:
            raise ValueError("pathway graph is not a DAG.")
        if len(order) == 0:
            raise ValueError("pathway reaction graph has no reaction.")
        first_reaction = order[0]
        return (
            first_reaction not in self.catalyzed_reactions
            and first_reaction
            in self.kb.non_spontaneous_reactions_of_pathway(self.pathway_id)
        )


class RuleNotApplicableException(Exception):
    def __init__(self, rule):
        self.message = f"Rule {rule} not applicable"
        super().__init__(self.message)


class ApplicableRule(Rule):
    ID = "generic"

    def __init__(self, kb):
        pass

    def rule_applicability(self, requirements: list):
        for required_method in requirements:
            if getattr(required_method, "__unsupported", False):
                raise RuleNotApplicableException(self.ID)


class AllReactionsRule(ApplicableRule):
    ID = "all_reactions_catalyzed"

    def __init__(self, kb):
        self.rule_applicability([kb.reactions_of_pathway])

    def evaluate(self, context):
        """
        Applies when all reactions are catalyzed
        """
        return context.all_reactions_catalyzed

    def apply(self, context):
        return True


class NoReactionCatalyzedRule(ApplicableRule):
    ID = "all_reactions_missing"

    def __init__(self, kb):

        self.rule_applicability([kb.reactions_of_pathway])

    def evaluate(self, context):
        """
        Applies if none of the non spontaneous reactions is catalyzed, and there is at least one reaction in this case.
        """

    def apply(self, context):
        """
        Reject the pathway if none of its non spontaneous reactions is catalyzed.
        """
        return False


class TransportPathwayRule(ApplicableRule):
    ID = "transport_pathway"

    def __init__(self, kb):
        # I should probably check the applicacability of this rule before instantiating it
        # this rule requires the ontology of the pathway to be available
        self.rule_applicability([kb.ontology_parent_class_of_pathway])

    def evaluate(self, context):
        """
        Applies if the pathway is a transport pathway
        """
        return context.pathway_is_a("Transport-Pathway")

    def apply(self, context):
        "Reject a transport pathway"
        return False


class ElectronTransferRule(ApplicableRule):
    ID = "electron_transfer_pathway"

    def __init__(self, kb):
        self.rule_applicability([kb.ontology_parent_class_of_pathway])
        # I should probably check the applicacability of this rule before instantiating it
        # this rule requires the ontology of the pathway to be available

    def evaluate(self, context):
        """
        Applies if the pathway is a transport pathway
        """
        return (
            context.pathway_is_a("Electron-Transfer") and context.some_reaction_missing
        )

    def apply(self, context):
        "Reject a transport pathway"
        return False


class PathwayInTaxonomicRangeRule(ApplicableRule):
    ID = "taxonomic_range"

    def __init__(self, kb):
        self.rule_applicability([kb.pathway_taxonomic_range])

    def evaluate(self, context):
        return context.in_taxonomic_range and context.all_reactions_catalyzed

    def apply(self, context):
        return True


class PathwayDegradationMissingFirstRule(ApplicableRule):
    ID = "degradation_pathway_missing_first"

    def __init__(self, kb):
        self.rule_applicability(
            [kb.pathway_reaction_order, kb.ontology_parent_class_of_pathway]
        )

    def evaluate(self, context):
        try:
            return (
                context.pathway_is_a("Degradation") and context.missing_first_reaction()
            )
        except ValueError as _e:
            return False

    def apply(self, context):
        return False


class PathwayBiosynthesisMissingLastRule(ApplicableRule):
    ID = "biosynthesis_pathway_missing_last"

    def __init__(self, kb):
        self.rule_applicability(
            [kb.pathway_reaction_order, kb.ontology_parent_class_of_pathway]
        )

    def evaluate(self, context):
        try:
            return (
                context.pathway_is_a("Biosynthesis") and context.missing_last_reaction()
            )
        except ValueError as _e:
            return False

    def apply(self, context):
        return False


class MissingKeyReactionRule(ApplicableRule):
    ID = "key_reaction_missing"

    def __init__(self, kb):
        self.rule_applicability([kb.key_reactions_of_pathway])

    def evaluate(self, context):
        pathway_key_reactions = context.kb.key_reactions_of_pathway(context.pathway_id)
        return (
            len(set(pathway_key_reactions).difference(context.catalyzed_reactions)) > 0
        )

    def apply(self, context):
        return False


class PathwayScoreRule(ApplicableRule):
    ID = "pathway_score"

    def __init__(self, kb):
        # TODO: Is there some rule applicability to check?
        pass

    def evaluate(self, context):
        scorer = PathwayScore(
            context.config,
            context.pathway_id,
            context.pathway_non_spontaneous_reactions,
            context.kb,
            context.catalyzed_reactions,
            context.taxonomy,
            context.taxon_id,
        )
        pathway_score = scorer.pathway_score()
        return (
            pathway_score >= context.config["inference"]["pathway_score"]["threshold"]
        )


class DefaultRule(ApplicableRule):
    ID = "default"

    def __init__(self, kb):
        pass

    def evaluate(self, context):
        return True

    def apply(self, context):
        return False


AVAILABLE_RULE_CLASSES = [
    AllReactionsRule,
    NoReactionCatalyzedRule,
    TransportPathwayRule,
    ElectronTransferRule,
    PathwayInTaxonomicRangeRule,
    PathwayDegradationMissingFirstRule,
    PathwayBiosynthesisMissingLastRule,
    MissingKeyReactionRule,
    PathwayScoreRule,
]


class MetabolismRuleKit:
    def __init__(
        self,
        config: dict,
        kb: KnowledgeBase,
        catalyzed_reactions: set[str],
        taxon_id: int | None,
        taxonomy: NCBITaxonomyTree | None,
    ):
        self.config = config
        self.kb = kb
        self.catalyzed_reactions = catalyzed_reactions
        self.enabled_rule_identifiers: list[str] = config["inference"]["rules"]
        self.rule_id_to_rule_class = {cl.ID: cl for cl in AVAILABLE_RULE_CLASSES}
        self.rule_class = [
            self.rule_id_to_rule_class[rule_id]
            for rule_id in self.enabled_rule_identifiers
        ]
        self.taxon_id = taxon_id
        self.taxonomy = taxonomy

        self.rules = [rule_class(kb) for rule_class in self.rule_class]
        # Force the rule order
        for i in range(len(self.rules)):
            self.rules[i].priority = i

    def inference(self) -> Iterable[str]:
        for pathway_id in self.kb.pathways():
            engine = RuleEngine(frozenset(self.rules))
            context = MetabolismRuleContext(
                self.config,
                pathway_id,
                self.kb,
                self.catalyzed_reactions,
                self.taxon_id,
                self.taxonomy,
            )
            if engine.evaluate_first(context):
                yield pathway_id


def infer_metabolism(reactome: set[str], taxon_id: int | None, config):
    kb = select_kb(config)
    if taxon_id is not None:
        taxonomy = NCBITaxonomyTree(config["reference"]["ncbi_taxonomy"])
    else:
        taxonomy = None
    rule_kit = MetabolismRuleKit(config, kb, reactome, taxon_id, taxonomy)
    yield from rule_kit.inference()
