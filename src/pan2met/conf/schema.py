"""
Configuration schema validation
"""

import importlib.resources

import yaml
from cerberus import Validator

import pan2met
import pan2met.conf

from ..utils import logger

schema_definition = yaml.safe_load(
    importlib.resources.read_text(pan2met.conf, "schema.yaml")
)


def validate(config: dict) -> bool:
    """
    Validate a configuration dictionnary.

    :param config: the configuration dictionnary.
    :return: True if the configuration respects the schema.
    """
    validator = Validator(schema_definition)
    valid = validator.validate(config, schema_definition)
    if not valid:
        logger.error(validator.errors)
    return valid
