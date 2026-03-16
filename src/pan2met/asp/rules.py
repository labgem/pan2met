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


def complex_asp_rule(complex: str, monomers: list[str]) -> str:
    return (
        'complex("{complex}") :- '
        + " , ".join('monomer("{monomer}")' for monomer in monomers)
        + "."
    )
