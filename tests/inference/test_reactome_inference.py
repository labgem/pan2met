#!/usr/bin/env python3

import pytest
import pythoncyc

from pangenome2panmetabolome.utils import read_list
from pangenome2panmetabolome.knowledge_base import KnowledgeBase
from pangenome2panmetabolome.io.sparql import SPARQLBackendKnowledgeBase
from pangenome2panmetabolome.reactome import (
    infer_reactome_from_ec_numbers,
    infer_reactome_from_monomers,
    infer_reactome_from_monomers_asp,
)


REACTIONS = """RXN-5424
ALCOHOL-DEHYDROG-GENERIC-RXN
ALCOHOL-DEHYDROG-RXN
RXN-7693
RXN-7694
RXN-12448
MEVALDATE-REDUCTASE-RXN
RXN3O-4113
RXN-10915
RXN-7657
RXN66-478
RXN-21862
RXN-13198
RXN-7706
RXN-10911
RXN-10781
RXN-7700""".splitlines()
REACTIONS = map(str.strip, REACTIONS)

EC_NUMBERS = ["1.1.1.1"]
"""
PathwayTools function "(MAP-EC-NUMBER-TO-REACTIONS "1.1.1.1")" enables to compute the list of reaction having a the EC-number 1.1.1.1.

To retrieve the list of reactions:

.. code :: lisp

  (loop for rxn in (map-ec-number-to-reactions "1.1.1.1")
        for name (get-frame-name rxn)
        do (print name))
"""


@pytest.fixture
def kb():
    kb: KnowledgeBase = SPARQLBackendKnowledgeBase()
    return kb


def test_infer_reactions_from_ec_number(kb):
    reactions = infer_reactome_from_ec_numbers(EC_NUMBERS, kb)
    assert sorted(reactions) == sorted(REACTIONS)


def test_reactome_inference_from_monomers_sanity_check():
    """
    Ensure the results provided by a ASP Clingo based reactome inference approach (from monomers),
    provides the same results than a direct python approach using the same deduction rule,
    """

    _reactions_test1 = ["RXN-20928", "RXN-20930"]

    _monomer_test1 = []


def test_reactome_inference_from_monomers_compare_asp_and_straightforward():
    """
    Compare the results of reactome given by straightforward python reactome inference
    from the set of monomers, compared with the ASP approach.

    To pass, this test requires a file containing a list of monomers, an active PathwayTools API instance
    and a SPARQL Endpoint instance, as well as the GPR rules.
    """
    monomers = read_list("./tmp/potential_monomers_sample100.list")
    pgdb = pythoncyc.select_organism("meta")
    reactions_from_pythoncyc = infer_reactome_from_monomers(monomers, pgdb)
    reactions_from_asp = infer_reactome_from_monomers_asp(
        monomers, "./tmp/metacyc29.5_gpr.lp"
    )
    assert reactions_from_pythoncyc == reactions_from_asp
