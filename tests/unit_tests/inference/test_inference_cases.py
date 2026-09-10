import pytest

from pan2met.io.knowledge_base import select_kb
from pan2met.utils import read_list, write_output
from pan2met.config import default_config as config

# @pytest.fixture
# def kb():
#     return select_kb(config)


# def test_d_apiose(kb):
#     """
#     Test a very simple case.
#     A pathway where every reaction is in the verified reactome.
#     """
#     reactome: set[str] = set(read_list("tests/cases/test0/reactome"))
#     expected_pathway_set: set[str] = set(read_list("tests/cases/test0/pathways"))
#     ecoli_tax_id: int = 562
#     inference = pan2met.inference.metabolism.PathwayInference(
#         kb, reactome, ecoli_tax_id
#     )
#     infered_pathway_set: set[str] = inference.inferred_pathways()
#     write_output("/tmp/pathway_list.txt", list(infered_pathway_set))
#     assert expected_pathway_set & infered_pathway_set == expected_pathway_set
