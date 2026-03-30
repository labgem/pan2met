"""
Default configuration
"""

import configparser
import importlib.resources

from dotenv import load_dotenv

import pan2met.conf

load_dotenv()

config: configparser.ConfigParser = configparser.ConfigParser()
default_config = importlib.resources.read_text(pan2met.conf, "default.ini")
config.read_string(default_config)
