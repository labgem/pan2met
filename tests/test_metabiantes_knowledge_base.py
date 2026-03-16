import pytest

import pan2met.io.metabiantes.kb


@pytest.fixture
def kb():
    return pan2met.io.metabiantes.kb.MetabiantesKnowledgeBase()


def test_list_pathways(kb):
    assert len(kb.pathways()) > 100
