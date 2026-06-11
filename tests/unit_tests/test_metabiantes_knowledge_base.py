import pytest

from pan2met.io.knowledge_base import select_kb, KnowledgeBase
from pan2met.config import default_config

@pytest.fixture
def kb() -> KnowledgeBase:
    return select_kb(default_config)

def test_list_pathways(kb):
    assert len(kb.pathways()) > 100
