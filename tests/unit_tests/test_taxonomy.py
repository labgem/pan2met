#!/usr/bin/env python3

from pan2met import taxonomy
from pan2met.config import config


def test_taxonomy_parent_tax_id():
    tree = taxonomy.NCBITaxonomyTree(config["reference"]["ncbi_taxonomy"])

    # Escherichia coli is a bacteria
    ecoli_tax_id = 562
    bacteria_tax_id = 2
    assert tree.is_child_of_parent_tax_id(bacteria_tax_id, ecoli_tax_id)

    # Homo sapiens is not a flowering plant
    homosapiens_tax_id = 9606
    plant_tax_id = 3398
    assert not tree.is_child_of_parent_tax_id(plant_tax_id, homosapiens_tax_id)

    # Angiospermae is not an orchidaceae
    orchidaceae_tax_id = 4747
    assert not tree.is_child_of_parent_tax_id(orchidaceae_tax_id, plant_tax_id)
