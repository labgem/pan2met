from collections.abc import Iterable


def complex_components_asp_rules(complex: str, components: str) -> str:
    return (
        f'complex("{complex}") :-'
        + ".".join(f'monomer("{component}")' for component in components)
        + "."
    )


def monomer_asp_rule(monomer: str):
    return f'monomer("{monomer}").'


def potential_monomer_asp_rule(monomer: str) -> str:
    return f'potential_monomer("{monomer}").'


def complex_asp_rule(complex: str, monomers: Iterable[str]) -> str:
    return (
        f'complex("{complex}") :- '
        + " , ".join('monomer("{monomer}")' for monomer in monomers)
        + "."
    )


def list_to_asp_atoms(predicate: str, literals: Iterable[str]) -> str:
    return "\n".join(f'{predicate}("{literal}").' for literal in literals)


def reactome_to_asp(reactome) -> str:
    return list_to_asp_atoms("reaction", reactome)


def catalyzis_to_asp(reactions) -> str:
    return list_to_asp_atoms("catalyzis", reactions)
