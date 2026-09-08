"""
Default configuration
"""

import configparser
import importlib.resources
import io
import logging
from pathlib import Path

import pan2met.conf

logger = logging.getLogger("pan2met:config")

# Load default config
default_config: configparser.ConfigParser = configparser.ConfigParser()
default_config_str = importlib.resources.read_text(pan2met.conf, "default.ini")
default_config.read_string(default_config_str)


def override_config(filename: Path) -> configparser.ConfigParser:
    """
    Read default config and override it with values from config of a file

    :param filename: path to a config file.
    :return: the config
    """

    config = default_config

    # Update config globally overriding default_config with keys from given config filename
    config_override = configparser.ConfigParser()
    config_override.read([filename])
    config.update(config_override)

    # Log the used config
    with io.StringIO() as config_string_stream:
        config.write(config_string_stream)
        logger.info(f"Using config:\n{config_string_stream.getvalue()}")

    return config
