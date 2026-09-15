"""
Default configuration
"""

import importlib.resources
import io
import logging
from pathlib import Path

import yaml

import pan2met.conf
import pan2met.conf.schema

logger = logging.getLogger("pan2met:config")

# Load default config
default_config_str = importlib.resources.read_text(pan2met.conf, "default.yaml")
default_config = yaml.safe_load(default_config_str)


def override_config(filename: Path) -> dict:
    """
    Read default config and override it with values from config of a file

    :param filename: path to a config file.
    :return: the config
    """

    config = default_config

    # Update config globally overriding default_config with keys from given config filename
    with open(filename, "r") as file:
        config_override = yaml.safe_load(file)
        if pan2met.conf.schema.validate(config_override):
            config.update(config_override)
        else:
            raise ValueError(
                "The provided configuration does not match the expected schema."
            )

    # Log the used config
    with io.StringIO() as config_string_stream:
        config.write(config_string_stream)
        logger.info(f"Using config:\n{config_string_stream.getvalue()}")

    return config
