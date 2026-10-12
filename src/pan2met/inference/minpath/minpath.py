import importlib.resources
from collections.abc import Iterable
from typing import Literal

import pan2met
from pan2met.config import default_config
from pan2met.utils import logger

from ...asp import rules
from ...io.knowledge_base import KnowledgeBase, select_kb


def inferred_pathways(
    kb: KnowledgeBase, reactome: set[str], nb_model=1
) -> Iterable[set[str]]:
    """
    List all inferred pathways using Answer Set Programming/ASP inference rules
    with the MinPath approach to select a minimum set of metabolic pathways.


    :param kb: a metabolic pathway knowledgebase
    :param reactome: a set of reactions
    :param nb_model: the max number of model to get (0 for all models)
    :yield: a minimal set of pathway covering the reactions
    :raise: error if the model is unsatisfiable, e.g., if a reaction cannot be coverred.
    """

    # We keep clyngor as an optional dependency
    import clyngor

    with importlib.resources.path(
        pan2met, "asp/rules/minpath_like.lp"
    ) as minpath_rule_path:
        logger.info("Generate ASP atoms for MinPath inference.")
        catalysis_asp = rules.list_to_asp_atoms(
            "catalysis", list(reactome)
        )  # writes reaction/1. atoms : reaction("RXN-1"). if reaction RXN-1 is in the 'known' reactome of the organism
        pathways = kb.pathways()
        pathway_asp = rules.list_to_asp_atoms(
            "pathway", pathways
        )  # writes pathway/1 atoms: pathway("PWY-8089"). for instance.

        is_in_pathway_asp = "\n".join(
            f'is_in_pathway("{reaction}", "{pathway}").'
            for pathway in pathways
            for reaction in kb.reactions_of_pathway(pathway)
        )  # writes is_in_pathway/2 atoms: is_in_pathway(Reaction, Pathway).

        inline_asp = catalysis_asp + pathway_asp + is_in_pathway_asp
        with open("/tmp/debug.lp", "w") as f:
            f.write(inline_asp)

        logger.info("Launch the clingo solver")
        answers = clyngor.solve(
            [minpath_rule_path],
            inline=inline_asp,
            nb_model=nb_model,
        )
        for answer in answers:
            inferred: set[str] = set()
            for predicate, value in answer:
                if predicate == "inferred_pathway":
                    identifier = value[0]
                    identifier = identifier.replace('"', "")
                    inferred.add(identifier)
            yield inferred


def infer_metabolism(
    reactome: set[str], method: Literal["clingo"] = "clingo", config=default_config
) -> set[str]:
    kb = select_kb(config)
    pathways = next(inferred_pathways(kb, reactome, nb_model=1))
    return pathways
