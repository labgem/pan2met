"""
pan2met pathway-operon-filler --confident-enzyme-catalyzis confident_catalyzis.tsv --less-confident-enzyme-catalyzis less_confident_catalyzis.tsv --pangenome-graph pangenome.gt --output report.tsv
"""

import subprocess
import tempfile
from pathlib import Path


def run(cmd):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def test_metabolism_inference():
    with tempfile.TemporaryDirectory() as tmpdirname:
        outdir = Path(tmpdirname)

        confident_filename = "./tests/test_data/operon_filler/confident_reactions.tsv"
        less_confident_filename = (
            "./tests/test_data/operon_filler/confident_reactions.tsv"
        )
        pangenome_filename = "./tests/test_data/mock_pangenome_graph.gt"
        output_filename = outdir / "operon_filler_report.tsv"

        cmd = (
            f"pan2met pathway-operon-filler"
            f' --confident-enzyme-catalyzis "{confident_filename}"'
            f' --less-confident-enzyme-catalyzis "{less_confident_filename}"'
            f' --pangenome-graph "{pangenome_filename}" --output "{output_filename}"'
        )

        run(cmd)

        expected_outputs = [
            output_filename,
        ]

        for file in expected_outputs:
            assert file.exists(), f"Expected file {file} not found after `{cmd}`."
            assert file.stat().st_size > 0, f"File {file} is empty after `{cmd}`"
