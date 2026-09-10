import pytest

from pan2met.config import default_config
from pan2met.io.knowledge_base import KnowledgeBase, select_kb


@pytest.fixture
def kb() -> KnowledgeBase:
    # Force the reference source to use metabiantes sql
    default_config["reference"]["knowledge_base"] = "metabiantes"
    return select_kb(default_config)


def test_list_pathways(kb):
    assert len(kb.pathways()) > 100


def test_kb_monomers(kb):
    monomers = kb.monomers()
    assert len(monomers) > 100
    assert "1-PFK-MONOMER" in monomers


def test_kb_pathways(kb):
    pathways = kb.pathways()
    assert len(pathways) > 100
    assert "PWY-8089" in pathways


def test_kb_reactions_of_pathway(kb):
    pathway = "PWY-8089"
    expected_reactions = ["RXN-20930", "RXN-20928", "RXN-20929"]
    reactions = kb.reactions_of_pathway(pathway)
    assert sorted(expected_reactions) == sorted(reactions)


def test_kb_spontaneous_reactions(kb):
    spontaneous_reactions = kb.spontaneous_reactions()
    assert len(spontaneous_reactions) > 100


def test_kb_orphan_reactions(kb):
    orphan_reactions = kb.orphan_reactions()
    assert len(orphan_reactions) > 100


def test_kb_non_spontaneous_reactions_of_pathway(kb):
    pathway = "PWY-8089"
    non_spontanous_reactions = kb.non_spontaneous_reactions_of_pathway(pathway)
    assert len(non_spontanous_reactions) == 3


def test_kb_non_orphan_non_spontaneous_reactions_of_pathway(kb):
    pathway = "PWY-8089"
    non_spontanous_reactions = kb.non_spontaneous_reactions_of_pathway(pathway)
    assert len(non_spontanous_reactions) == 3
    # TODO: test with a pathway with some spontaneous or orphan reactions


def test_kb_enzymes_of_reaction(kb):
    reaction = "3.4.21.22-RXN"
    expected_enzyme = "CPLX66-1735"
    enzymes = kb.enzymes_of_reaction(reaction)
    assert len(enzymes) == 1
    assert enzymes == [expected_enzyme]


def test_kb_reactions_by_ec_number(kb):
    ec_number = "3.4.21.22"
    expected_reactions = ["3.4.21.22-RXN"]
    reactions = kb.reactions_by_ec_number(ec_number)
    assert expected_reactions == reactions


def test_kb_pathway_taxonomic_range(kb):
    pathway = "PWY-8089"
    expected_taxonomic_range = 2
    taxonomic_range = kb.pathway_taxonomic_range(pathway)
    assert taxonomic_range is not None and taxonomic_range == expected_taxonomic_range


def test_kb_reaction_is_key(kb):
    assert kb.reaction_is_key("PWY-8534", "RXN-21500")
    assert not kb.reaction_is_key("PWY-8089", "RXN-21500")


def test_kb_key_reactions_of_pathway(kb):
    pathway = "PWY-8534"
    expected_key_reactions = ["RXN-21500", "RXN-11757", "RXN-24719"]
    key_reactions = kb.key_reactions_of_pathway(pathway)
    assert sorted(expected_key_reactions) == sorted(key_reactions)


def test_kb_pathway_reaction_order(kb):
    pathway = "PWY-8089"
    expected_reaction_order = [("RXN-20928", "RXN-20930"), ("RXN-20929", "RXN-20928")]
    reaction_order = kb.pathway_reaction_order(pathway)
    assert sorted(reaction_order) == sorted(expected_reaction_order)


def test_kb_pathways_with_reaction(kb):
    reaction = "1-ACYLGLYCEROL-3-P-ACYLTRANSFER-RXN"
    pathways_with_reaction = kb.pathways_with_reaction(reaction)
    assert len(pathways_with_reaction) > 5


def test_kb_variants_of_pathway(kb):
    pathway = "PWY-8089"
    expected_variants = ["PWY-8093"]
    variants = kb.variants_of_pathway(pathway)
    assert variants == expected_variants


def test_kb_species_evidence_of_pathway(kb):
    pathway = "PWY-8089"
    expected_species_evidence = [218491, 339671, 435590, 476272, 558270]
    species_evidence = kb.species_evidence_of_pathway(pathway)
    assert sorted(expected_species_evidence) == sorted(species_evidence)


def test_kb_ontology_parent_class_of_pathway(kb):
    pathway = "PWY-8089"
    parent_classes = kb.ontology_parent_class_of_pathway(pathway)
    assert len(parent_classes) > 4


def test_kb_complex(kb):
    complex = kb.complex()
    assert len(complex) > 100


def test_kb_complex_components(kb):
    complex = "ABC-4-CPLX"
    expected_components = [
        "ARTJ-MONOMER",
        "ARTM-MONOMER",
        "ARTP-MONOMER",
        "ARTQ-MONOMER",
    ]
    complex_components = kb.components_of_complex(complex)
    assert sorted(expected_components) == sorted(complex_components)
