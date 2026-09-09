import pytest

from pan2met.config import default_config
from pan2met.io.knowledge_base import KnowledgeBase, select_kb


@pytest.fixture
def kb() -> KnowledgeBase:
    # Force the reference source to use metabiantes sql
    default_config["reference"]["knowledge_base"] = "metabiantes"
    default_config["metabiantes"]["database"] = "tests/test_data/mock_metabiantes.db"
    return select_kb(default_config)


def test_kb_monomers(kb: KnowledgeBase):
    monomers = kb.monomers()
    assert len(monomers) == 2
    assert "MONOMER-1" in monomers


def test_kb_pathways(kb: KnowledgeBase):
    pathways = kb.pathways()
    assert len(pathways) == 3
    assert "PWY-1" in pathways


def test_kb_reactions_of_pathway(kb: KnowledgeBase):
    pathway = "PWY-0"
    expected_reactions = ["RXN-0", "RXN-1", "RXN-2"]
    reactions = kb.reactions_of_pathway(pathway)
    assert sorted(expected_reactions) == sorted(reactions)


def test_kb_spontaneous_reactions(kb: KnowledgeBase):
    spontaneous_reactions = kb.spontaneous_reactions()
    assert len(spontaneous_reactions) == 3


def test_kb_orphan_reactions(kb: KnowledgeBase):
    orphan_reactions = kb.orphan_reactions()
    assert len(orphan_reactions) == 2
    assert sorted(orphan_reactions) == sorted(["RXN-4", "RXN-3"])


def test_kb_non_spontaneous_reactions_of_pathway(kb: KnowledgeBase):
    pathway = "PWY-0"
    non_spontanous_reactions = kb.non_spontaneous_reactions_of_pathway(pathway)
    assert len(non_spontanous_reactions) == 1
    assert non_spontanous_reactions == ["RXN-1"]


def test_kb_non_orphan_non_spontaneous_reactions_of_pathway(kb: KnowledgeBase):
    pathway = "PWY-0"
    non_spontanous_reactions = kb.non_spontaneous_reactions_of_pathway(pathway)
    assert len(non_spontanous_reactions) == 1
    assert non_spontanous_reactions == ["RXN-1"]
    # TODO: test with a pathway with some spontaneous or orphan reactions


def test_kb_enzymes_of_reaction(kb: KnowledgeBase):
    reaction = "RXN-2"
    expected_enzyme = "CPLX-1"
    enzymes = kb.enzymes_of_reaction(reaction)
    assert len(enzymes) == 1
    assert enzymes == [expected_enzyme]


def test_kb_reactions_by_ec_number(kb: KnowledgeBase):
    ec_number = "1.1.1.1"
    expected_reactions = ["RXN-0"]
    reactions = kb.reactions_by_ec_number(ec_number)
    assert expected_reactions == reactions


def test_kb_pathway_taxonomic_range(kb: KnowledgeBase):
    pathway = "PWY-0"
    expected_taxonomic_range = 1234
    taxonomic_range = kb.pathway_taxonomic_range(pathway)
    assert taxonomic_range is not None and taxonomic_range == expected_taxonomic_range


def test_kb_reaction_is_key(kb: KnowledgeBase):
    assert kb.reaction_is_key("PWY-0", "RXN-0")
    assert not kb.reaction_is_key("PWY-1", "RXN-0")


def test_kb_key_reactions_of_pathway(kb: KnowledgeBase):
    pathway = "PWY-0"
    expected_key_reactions = ["RXN-0"]
    key_reactions = kb.key_reactions_of_pathway(pathway)
    assert sorted(expected_key_reactions) == sorted(key_reactions)


def test_kb_pathway_reaction_order(kb: KnowledgeBase):
    pathway = "PWY-0"
    expected_reaction_order = [("RXN-0", "RXN-1"), ("RXN-1", "RXN-2")]
    reaction_order = kb.pathway_reaction_order(pathway)
    assert sorted(reaction_order) == sorted(expected_reaction_order)


def test_kb_pathways_with_reaction(kb: KnowledgeBase):
    reaction = "RXN-0"
    pathways_with_reaction = kb.pathways_with_reaction(reaction)
    assert len(pathways_with_reaction) == 2


def test_kb_variants_of_pathway(kb: KnowledgeBase):
    pathway = "PWY-0"
    expected_variants = ["PWY-2"]
    variants = kb.variants_of_pathway(pathway)
    assert variants == expected_variants


def test_kb_species_evidence_of_pathway(kb: KnowledgeBase):
    pathway = "PWY-0"
    expected_species_evidence = [1234]
    species_evidence = kb.species_evidence_of_pathway(pathway)
    assert sorted(expected_species_evidence) == sorted(species_evidence)


def test_kb_ontology_parent_class_of_pathway(kb: KnowledgeBase):
    pathway = "PWY-0"
    parent_classes = kb.ontology_parent_class_of_pathway(pathway)
    assert sorted(parent_classes) == sorted(
        ["SyntheticMetabolism", "Synthetic", "Root"]
    )


def test_kb_complex(kb: KnowledgeBase):
    complex = kb.complex()
    assert complex == ["CPLX-1"]


def test_kb_complex_components(kb: KnowledgeBase):
    complex = "CPLX-1"
    expected_components = ["MONOMER-1", "MONOMER-3"]
    complex_components = kb.components_of_complex(complex)
    assert sorted(expected_components) == sorted(complex_components)
