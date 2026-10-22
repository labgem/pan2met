"""
pan2met proteic-complex --monomers monomer.list --output complex.list
"""

import subprocess
import tempfile
from pathlib import Path

from pan2met.utils import write_output

REACTIONS = [
    "RXN-20930",
    "RXN-20928",
    "RXN-20929",
]


def run(cmd):
    print(f"Running {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def test_metabolism_inference():
    with tempfile.TemporaryDirectory() as tmpdirname:
        outdir = Path(tmpdirname)

        reactions_filename = outdir / "reactions.list"
        output_filename = outdir / "pathways.list"

        # Write a sample list of protein monomers
        write_output(reactions_filename, REACTIONS)

        cmd = f'pan2met minpath --reactions "{reactions_filename}" --output "{output_filename}"'

        run(cmd)

        expected_outputs = [
            output_filename,
        ]

        for file in expected_outputs:
            assert file.exists(), f"Expected file {file} not found after `{cmd}`."
            assert file.stat().st_size > 0, f"File {file} is empty after `{cmd}`"
