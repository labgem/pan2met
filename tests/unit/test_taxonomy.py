from pathlib import Path

import pytest

from pan2met import taxonomy
from pan2met.config import default_config as config


def _require_ncbi_taxonomy():
    has_taxonomy = False
    try:
        assert Path(config["reference"]["ncbi_taxonomy"]).exists()
        has_taxonomy = True
    except AssertionError:
        pass

    return pytest.mark.skipif(
        not has_taxonomy, reason="access to the NCBI Taxonomy dump required"
    )


require_ncbi_taxonomy = _require_ncbi_taxonomy()


@require_ncbi_taxonomy
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
