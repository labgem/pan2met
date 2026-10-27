"""
pan2met metabolism --reactions reactions.list --output pathways.list \
    --taxon-id 562 --reason pathways_pan2met_reason.log
"""

import subprocess
import tempfile
from pathlib import Path

import pytest

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


def run(cmd):
    subprocess.run(cmd, shell=True, check=True)


@require_ncbi_taxonomy
def test_metabolism_inference():
    with tempfile.TemporaryDirectory() as tmpdirname:
        outdir = Path(tmpdirname)
        reactions_filename = "./tests/data/cases/test0/reactions.list"
        output_filename = outdir / "pathways.list"
        # output_log_filename = outdir / "reason.log"
        taxon_id = 562  # E. coli

        cmd = (
            f'pan2met metabolism --reactions "{reactions_filename}" --output "{output_filename}"'
            f" --taxon-id {taxon_id} "
            # --reason "{output_log_filename}"'
        )

        run(cmd)

        expected_outputs = [
            output_filename,
            # output_log_filename,
        ]

        for file in expected_outputs:
            assert file.exists(), f"Expected file {file} not found after `{cmd}`."
            assert file.stat().st_size > 0, f"File {file} is empty after `{cmd}`"
