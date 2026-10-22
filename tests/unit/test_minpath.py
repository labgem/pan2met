import pytest

from pan2met.config import default_config
from pan2met.inference.minpath import minpath
from pan2met.io.knowledge_base import KnowledgeBase, select_kb


@pytest.fixture
def kb() -> KnowledgeBase:
    return select_kb(default_config)


def test_minpath(kb):
    pathway = "PWY-8089"
    pathway_reactions = [
        "RXN-20930",
        "RXN-20928",
        "RXN-20929",
    ]

    results = list(minpath.inferred_pathways(kb, pathway_reactions))

    assert {pathway} in results


def test_minpath_unsatisfiable(kb):
    pathway_reactions = [
        "RXN-DUMMY",  # this should make the cover unsatisfiable
        "RXN-20930",
        "RXN-20928",
        "RXN-20929",
    ]

    results = list(minpath.inferred_pathways(kb, pathway_reactions))

    assert results == []


def test_minpath_minimal(kb):
    pathway_reactions = [
        "RXN-7567"  # this reaction is present in at least two pathways.
    ]
    # each model should return a set of minimal size of 1 pathway only.

    results = list(minpath.inferred_pathways(kb, pathway_reactions))
    assert len(results) >= 1
    assert all(len(x) == 1 for x in results)
