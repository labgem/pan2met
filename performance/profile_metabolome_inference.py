import pan2met.io.metabiantes.kb
import pan2met.inference.metabolome
from pan2met.utils import read_list, write_output


def main():
    kb = pan2met.io.metabiantes.kb.MetabiantesKnowledgeBase()
    reactome: set[str] = set(read_list("tests/cases/test0/reactome"))
    expected_pathway_set: set[str] = set(read_list("tests/cases/test0/pathways"))
    ecoli_tax_id: int = 562
    inference = pan2met.inference.metabolome.PathwayInference(
        kb, reactome, ecoli_tax_id
    )
    infered_pathway_set: set[str] = inference.inferred_pathways()
    write_output("/tmp/pathway_list.txt", list(infered_pathway_set))
    assert expected_pathway_set & infered_pathway_set == expected_pathway_set


if __name__ == "__main__":
    main()
