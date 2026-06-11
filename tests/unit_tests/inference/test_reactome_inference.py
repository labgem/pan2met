#!/usr/bin/env python3

import pytest
import pythoncyc

from pan2met.utils import read_list
from pan2met.io.knowledge_base import KnowledgeBase, select_kb
from pan2met.inference.reactome import (
    infer_reactome_from_ec_numbers,
    infer_reactome_from_monomers,
    infer_reactome_from_monomers_asp,
)

from pan2met.config import default_config


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
    kb: KnowledgeBase = select_kb(default_config)
    return kb


def test_infer_reactions_from_ec_number(kb):
    reactions = infer_reactome_from_ec_numbers(EC_NUMBERS, kb)
    assert sorted(reactions) == sorted(REACTIONS)
