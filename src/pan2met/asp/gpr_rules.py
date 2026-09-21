"""
Gene-Protein-Reaction rules in Answer Set Programming

Launch Pathway Tools python API with:
    pathway-tools -lisp -python-local-only-non-strict
"""

from collections.abc import Iterable
from typing import Literal


def catalysis_asp_rule(
    polypeptide: str, polypeptide_type: Literal["monomer", "complex"], reaction: str
) -> str:
    return f'reaction("{reaction}") :- {polypeptide_type}("{polypeptide}").'


def protein_complex_asp_rule(complex: str, components: Iterable[str]) -> str:
    return (
        f'complex("{complex}") :- '
        + " , ".join(f'monomer("{component}")' for component in components)
        + "."
    )
