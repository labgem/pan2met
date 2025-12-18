"""
Default configuration
"""

import configparser
import importlib.resources

from dotenv import load_dotenv

import pangenome2panmetabolome
from .utils import logger

load_dotenv()

config = configparser.ConfigParser()
default_config = importlib.resources.read_text(
    pangenome2panmetabolome, "conf/default.ini"
)
config.read_string(default_config)

logger.debug(config)
