


"""
pan2met metabolism --reactions reactions.list --output pathways.list \
    --taxon-id 562 --reason pathways_pan2met_reason.log
"""
import os
import logging
import tempfile
import importlib
from pathlib import Path

import pytest


def test_metabolism_inference():

    with tempfile.TemporaryDirectory() as tmpdirname:
        outdir = Path(tmpdirname)
        reactions_filename = "./tests/test_data/reactions.list"
        output_filename = outdir / "pathways.list"
        output_log_filename = outdir / "reason.log"
        taxon_id = 562 # E. coli

        cmd = (
            f"pan2met metabolism --reactions {reactions_filename} --output pathways.list"
            f" --taxon-id {taxon_id} --reason {output_log_filename}"
        )


        expected_outputs = [
            output_filename,
            output_log_filename,
        ]

        for file in expected_outputs:
            assert filename.exists(), f"Expected file {file} not found after `{cmd}`."
            assert file.stat().st_size() > 0, f"File {file} is empty after `{cmd}`"
