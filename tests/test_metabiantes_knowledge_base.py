import pytest

import pangenome2panmetabolome.io.metabiantes.kb


@pytest.fixture
def kb():
    return pangenome2panmetabolome.io.metabiantes.kb.MetabiantesKnowledgeBase()


def test_list_pathways(kb):
    assert len(kb.pathways()) > 100
