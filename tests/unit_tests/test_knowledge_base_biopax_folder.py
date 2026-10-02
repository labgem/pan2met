import pytest

from pan2met.config import default_config
from pan2met.io.knowledge_base import KnowledgeBase, select_kb
from pan2met.io.knowledge_base.biopax_folder import BioPAXFolderKnowledgeBase

BIOPAX_XML_BASE = "http://www.biopax.org/examples/myExample#"


def prefixed(id: str) -> str:
    return f"{BIOPAX_XML_BASE}{id}"


@pytest.fixture
def kb() -> KnowledgeBase:
    # Force the reference source to use metabiantes sql
    default_config["reference"]["knowledge_base"] = "biopax_folder"
    return select_kb(default_config)


def test_list_pathways(kb: BioPAXFolderKnowledgeBase):
    assert len(kb.pathways()) > 0


def test_kb_monomers(kb: BioPAXFolderKnowledgeBase):
    monomers = kb.monomers()
    assert len(monomers) > 0
    assert prefixed("#Protein_54") in monomers


def test_kb_pathways(kb: BioPAXFolderKnowledgeBase):
    pathways = kb.pathways()
    assert len(pathways) > 0
    assert prefixed("#Pathway50") in pathways


def test_kb_reactions_of_pathway(kb: BioPAXFolderKnowledgeBase):
    pathway = prefixed("#Pathway50")
    expected_reactions = ["#glucokinase", "#phosphoglucoisomerase"]
    expected_prefixed_reactions = list(map(prefixed, expected_reactions))
    reactions = kb.reactions_of_pathway(pathway)
    assert sorted(expected_prefixed_reactions) == sorted(reactions)


def test_kb_spontaneous_reactions(kb: BioPAXFolderKnowledgeBase):
    spontaneous_reactions = kb.spontaneous_reactions()
    assert len(spontaneous_reactions) >= 1


def test_kb_orphan_reactions(kb: BioPAXFolderKnowledgeBase):
    orphan_reactions = kb.orphan_reactions()
    assert len(orphan_reactions) >= 1


def test_kb_non_spontaneous_reactions_of_pathway(kb: BioPAXFolderKnowledgeBase):
    pathway = "#Pathway50"
    prefixed_pathway = prefixed(pathway)
    non_spontanous_reactions = kb.non_spontaneous_reactions_of_pathway(prefixed_pathway)
    assert len(non_spontanous_reactions) == 2


def test_kb_non_orphan_non_spontaneous_reactions_of_pathway(
    kb: BioPAXFolderKnowledgeBase,
):
    pathway = "#Pathway50"
    prefixed_pathway = prefixed(pathway)
    non_spontanous_reactions = kb.non_orphan_non_spontaneous_reactions_of_pathway(
        prefixed_pathway
    )
    assert len(non_spontanous_reactions) == 0


def test_kb_enzymes_of_reaction(kb: BioPAXFolderKnowledgeBase):
    reaction = prefixed("#glucokinase")
    expected_enzyme = prefixed("#Protein_54")
    # corresponds to BioPAX bp:Catalyzis' bp:controlled/bp:controller
    enzymes = kb.enzymes_of_reaction(reaction)
    assert len(enzymes) == 1
    assert enzymes == [expected_enzyme]


def test_kb_reactions_by_ec_number(kb: BioPAXFolderKnowledgeBase):
    ec_number = "5.3.1.9"
    expected_reactions = ["#phosphoglucoisomerase"]
    expected_reactions_prefixed = list(map(prefixed, expected_reactions))
    reactions = kb.reactions_by_ec_number(ec_number)
    assert expected_reactions_prefixed == reactions


def test_kb_pathway_taxonomic_range(kb: BioPAXFolderKnowledgeBase):
    pathway = prefixed("#Pathway50")
    with pytest.raises(NotImplementedError):
        kb.pathway_taxonomic_range(pathway=pathway)


# def test_kb_reaction_is_key(kb: BioPAXFolderKnowledgeBase):
#     dummy_reaction = "nope"
#     dummy_pathway = "nope"
#     with pytest.raises(NotImplementedError):
#         kb.reaction_is_key(dummy_pathway, dummy_reaction)

# def test_kb_key_reactions_of_pathway(kb: BioPAXFolderKnowledgeBase):
#     pathway = "PWY-8534"
#     expected_key_reactions = ["RXN-21500", "RXN-11757", "RXN-24719"]
#     key_reactions = kb.key_reactions_of_pathway(pathway)
#     assert sorted(expected_key_reactions) == sorted(key_reactions)


# def test_kb_pathway_reaction_order(kb: BioPAXFolderKnowledgeBase):
#     pathway = "PWY-8089"
#     expected_reaction_order = [("RXN-20928", "RXN-20930"), ("RXN-20929", "RXN-20928")]
#     reaction_order = kb.pathway_reaction_order(pathway)
#     assert sorted(reaction_order) == sorted(expected_reaction_order)


def test_kb_pathways_with_reaction(kb: BioPAXFolderKnowledgeBase):
    reaction = prefixed("#glucokinase")
    pathways_with_reaction = kb.pathways_with_reaction(reaction)
    expected_pathways = [prefixed("#Pathway50")]
    assert len(pathways_with_reaction) == 1
    assert expected_pathways == pathways_with_reaction


# def test_kb_variants_of_pathway(kb: BioPAXFolderKnowledgeBase):
#     pathway = "PWY-8089"
#     expected_variants = ["PWY-8093"]
#     variants = kb.variants_of_pathway(pathway)
#     assert variants == expected_variants


# def test_kb_species_evidence_of_pathway(kb: BioPAXFolderKnowledgeBase):
#     pathway = "PWY-8089"
#     expected_species_evidence = [218491, 339671, 435590, 476272, 558270]
#     species_evidence = kb.species_evidence_of_pathway(pathway)
#     assert sorted(expected_species_evidence) == sorted(species_evidence)


# def test_kb_ontology_parent_class_of_pathway(kb: BioPAXFolderKnowledgeBase):
#     pathway = "PWY-8089"
#     parent_classes = kb.ontology_parent_class_of_pathway(pathway)
#     assert len(parent_classes) > 4


def test_kb_complex(kb: BioPAXFolderKnowledgeBase):
    complex = kb.complex()
    assert len(complex) == 1
    assert complex == [prefixed("#Complex_37")]


def test_kb_complex_components(kb: BioPAXFolderKnowledgeBase):
    complex = prefixed("#Complex_37")
    expected_components = [
        "#Protein_2",
        "#Protein_38",
    ]
    prefixed_expected_components = list(map(prefixed, expected_components))
    complex_components = kb.components_of_complex(complex)
    print(prefixed_expected_components)
    print(complex_components)
    assert sorted(prefixed_expected_components) == sorted(complex_components)
