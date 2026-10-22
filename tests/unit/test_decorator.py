from pan2met.config import default_config
from pan2met.io.knowledge_base.biopax_folder import (
    BioPAXFolderKnowledgeBase,
)


def test_unsupported():
    config = default_config
    kb = BioPAXFolderKnowledgeBase(config)

    assert kb.variants_of_pathway.__unsupported
