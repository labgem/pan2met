"""
Compute the set of proteic complex
"""

from typing import Iterable, Set

from pan2met.io.knowledge_base import KnowledgeBase


def is_constructible(complex_id: str, monomers: Set[str], kb: KnowledgeBase) -> bool:
    """
    Check if a complex is constrictible.
    A complex is constructible if all its component is either constructible or a listed monomer.
    :param complex_id: the identifier of the complex in the knowledge base
    :param monomers: the set of monomer identifiers
    :param: a metabolic pathway knowledgebase
    :return: True if the complex is constructible
    """
    for component in kb.components_of_complex(complex_id):
        if component not in monomers and not is_constructible(component, monomers, kb):
            return False
    return True


def infer_complex(kb: KnowledgeBase, monomers: Set[str]) -> Iterable[str]:
    """
    List the set of protein complex that are constructible based on a set of protein monomers.
    """
    potential_complex = kb.complex()
    for complex in potential_complex:
        if is_constructible(complex, monomers, kb):
            yield complex
