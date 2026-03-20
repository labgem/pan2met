"""
(pan)genome to (pan)metabolome
"""

from . import config
import importlib.metadata

__all__ = ["config"]
__version__ = importlib.metadata.version(__name__)
