"""
Use the builtin rule kit engine to infer metabolic pathway
"""

from ...io.knowledge_base import KnowledgeBase
from ...rulekit.rulekit import Rule
from ...taxonomy import NCBITaxonomyTree


class MetabolismRuleContext:
    def __init__(
        self,
        pathway_id: str,
        kb: KnowledgeBase,
        pathway_non_spontaneous_reactions: set[str],
        catalyzed_reactions: set[str],
        taxon_id: int | None,
        taxonomy: NCBITaxonomyTree | None,
    ):
        self.pathway_id = pathway_id
        self.kb = kb
        self.pathway_non_spontaneous_reactions = pathway_non_spontaneous_reactions
        self.catalyzed_reactions = catalyzed_reactions
        self.taxonomy = taxonomy
        self.taxon_id = taxon_id

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

    def missing_last_reaction(self):
        raise NotImplementedError()

    def missing_first_reaction(self):
        raise NotImplementedError()


class MetabolismRuleKit:
    def __init__(self):
        pass


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
            if required_method.__unsupported == True:
                raise RuleNotApplicableException(self.ID)


class AllReactionsRule(ApplicableRule):
    ID = "all_reactions_catalyzed"

    def __init__(self, kb):
        self.rule_applicability([kb.reactions_of_pathway])

    @property
    def priority(self):
        return 1

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

    @property
    def priority(self):
        return 1

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
        self.rule_applicability([kb.pathway_ontology_parent_classes])

    @property
    def priority(self):
        return 0

    def evaluate(self, context):
        """
        Applies if the pathway is a transport pathway
        """
        return context.pathway_is_a("Transport-Pathway")

    def apply(self, context):
        "Reject a transport pathway"
        return False


class ElectronTransferRule(ApplicableRule):
    ID = "electron_transfer"

    def __init__(self, kb):
        self.rule_applicability([kb.ontology_parent_class_of_pathway])
        # I should probably check the applicacability of this rule before instantiating it
        # this rule requires the ontology of the pathway to be available

    @property
    def priority(self):
        return 1

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

    @property
    def priority(self):
        return 3

    def evaluate(self, context):
        return context.in_taxonomic_range and context.all_reactions_catalyzed

    def apply(self, context):
        return True


class PathwayDegradationMissingFirstRule(ApplicableRule):
    ID = "degradation_pathway_missing_first"

    def __init__(self, kb):
        self.rule_applicability(
            [kb.pathway_reaction_order, kb.pathway_ontology_parent_classes]
        )

    @property
    def priority(self):
        return 4

    def evaluate(self, context):
        return context.pathway_is_a("Degradation") and context.missing_first_reaction()


class PathwayBiosynthesisMissingLastRule(ApplicableRule):
    ID = "biosynthesis_pathway_missing_last"

    def __init__(self, kb):
        self.rule_applicability(
            [kb.pathway_reaction_order, kb.pathway_ontology_parent_classes]
        )

    @property
    def priority(self):
        return 4

    def evaluate(self, context):
        return context.pathway_is_a("Degradation") and context.missing_last_reaction()
