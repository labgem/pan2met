"""
Analysis of the complete metabolic graph
"""

from typing import Tuple, Set

import graph_tool as gt
import padmet


def padmet_to_metabolic_graph(
    padmet_object: padmet.PadmetSpec,
) -> Tuple[Set, Set, Set, gt.Graph]:
    """
    Construct the metabolic multi directed graph from a PADMet format input file.

    :param padmet_object: A PADMet object=

    :return: the set of metabolites, the set of reactions, the set of enzymes
        and a multi-directed graph containing as nodes the both
        the sets of reactions and metabolites and the set of enzymes as edges.
    """


def find_linear_subpathways():
    """
    Find
    """
