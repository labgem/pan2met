from pathlib import Path
import pytest

from ppanggolin.pangenome import Pangenome
from ppanggolin.formats.readBinaries import check_pangenome_info


import pan2met
import pan2met.analysis
from pan2met.analysis.partition import PangenomePartition


@pytest.fixture
def pangenome():
    nils_ecoli_pangenome = Pangenome()
    nils_ecoli_pangenome.add_file(
        Path(
            "/home/sortion/Documents/projects/analysis/NILSmetabolism/results/ppanggolin/pangenome.h5"
        )
    )
    check_pangenome_info(
        nils_ecoli_pangenome,
        need_families=True,
        need_annotations=True,
        need_graph=True,
        disable_bar=True,
    )
    return nils_ecoli_pangenome


def test_find_partition_of_a_gene_family_by_gene_name(pangenome):
    gene_name: str = "NILS01_CDS_0212"
    gene_family_partition: PangenomePartition = (
        pan2met.analysis.partition.ppanggolin_partition_of_gene_family(
            gene_name, pangenome
        )
    )
    assert gene_family_partition == PangenomePartition.CORE
