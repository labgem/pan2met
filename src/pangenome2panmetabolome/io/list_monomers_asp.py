import pythoncyc

from . import metacyc
from ..utils import write_output


def potential_monomer_asp_rule(monomer: str) -> str:
    return f'potential_monomer("{monomer}").'


pgdb = pythoncyc.select_organism("meta")

write_output(
    "tmp/potential_monomers.lp",
    map(
        potential_monomer_asp_rule,
        map(metacyc.remove_pipes, metacyc.list_all_monomers(pgdb)),
    ),
)
