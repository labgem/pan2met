#!/usr/bin/env python3

"""
(pan)Metabolome inference
"""


def infer_metabolome(reactome: set[str]) -> list[str]:
    """
    Infer the list of pathways realized by a reactome.

    Arguments
    ---------

        reactome -- a set of reaction identifiers

    Returns
    -------

        a list of pathway identifiers

    Heuristics configuration
    ------------------------

    The heuristics used by the inference algorithm can be configured.
    TODO.
    """
