#!/usr/bin/env python3

"""
Export "reverse" GPR rules
"""

import argparse
from typing import Iterable

from ..utils import read_list, write_output


def reverse_gpr_asp(gpr_rules: Iterable[str]) -> Iterable[str]:
    """
    Convert forward GPR ASP rules to reverse GPR rules.

    .. code:: prolog

      reaction("RXN-1") :- complex("CPLX-1").
      reaction("RXN-2") :- monomer("MONOMER-C").
      complex("CPLX-1") :- monomer("MONOMER-A") , monomer("MONOMER-B").

    into:

    .. code:: prolog

      potential_complex("CPLX-1") :- reaction("RXN-1").
      potential_monomer("MONOMER-C") :- reaction("RXN-2").
      potential_monomer("MONOMER-A") :- potential_complex("CPLX-1").
      potential_monomer("MONOMER-B") :- potential_complex("CPLX-1").

    This is useful to infer the potential set of monomer that may be involved in the catalyzis of a set of reactions.
    """
    for gpr_rule in gpr_rules:
        parts = gpr_rule.split(" :- ")
        for index, part in enumerate(parts):
            parts[index] = part.replace("complex", "potential_complex").replace(
                "monomer", "potential_monomer"
            )
        head, tail = parts
        tail = tail.strip()[
            :-1
        ]  # Remove last char, that sould be the ending dot '.' of an ASP rule.
        for item in tail.split(" , "):
            yield f"{item} :- {head}."


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Convert a GPR rules ASP file into a reverse GPR rule ASP file, where monomer and complex catalyzing reactions can be infered from reactions"
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Input .lp filename with GPR rules to convert",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output .lp filename with 'reverse' ASP rules",
        required=True,
    )
    return parser, parser.parse_args()


def main():
    parser, args = parse_arguments()
    gpr_rules = read_list(args.input)
    reverse_gpr_rules = reverse_gpr_asp(gpr_rules)
    write_output(args.output, reverse_gpr_rules)


if __name__ == "__main__":
    main()
