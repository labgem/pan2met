"""
Format Answer Set Programming AnsProlog rules.
"""

from typing import Literal


def protein_complex_asp_rule(complex: str, components: list[str]) -> str:
    return (
        f'complex("{complex}") :- '
        + " , ".join(f'monomer("{monomer}")' for monomer in components)
        + "."
    )


def catalysis_asp_rule(
    enzyme: str, enzyme_type: Literal["monomer", "complex"], reaction: str
) -> str:
    return f'reaction("{reaction}") :- {enzyme_type}("{enzyme}").'


def monomer_asp_rule(monomer: str) -> str:
    return f'monomer("{monomer}").'


def pathway_asp_rule(pathway: str, reactions: list[str]) -> str:
    return (
        f'pathway("{pathway}") :- '
        + " , ".join(f'reaction("{reaction}")' for reaction in reactions)
        + "."
    )
