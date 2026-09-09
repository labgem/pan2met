"""
pan2met proteic-complex --monomers monomer.list --output complex.list
"""

import subprocess
import tempfile
from pathlib import Path

from pan2met.utils import write_output

MONOMERS = [
    "ARTJ-MONOMER",
    "ARTM-MONOMER",
    "ARTP-MONOMER",
    "ARTQ-MONOMER",
    "LEUC-MONOMER",
]


def run(cmd):
    subprocess.run(cmd, shell=True, check=True)


def test_metabolism_inference():
    with tempfile.TemporaryDirectory() as tmpdirname:
        outdir = Path(tmpdirname)

        monomers_filename = outdir / "monomers.list"
        output_filename = outdir / "complex.list"

        # Write a sample list of protein monomers
        write_output(monomers_filename, MONOMERS)

        cmd = f'pan2met proteic-complex --monomers "{monomers_filename}" --output "{output_filename}"'

        run(cmd)

        expected_outputs = [
            output_filename,
        ]

        for file in expected_outputs:
            assert file.exists(), f"Expected file {file} not found after `{cmd}`."
            assert file.stat().st_size > 0, f"File {file} is empty after `{cmd}`"
