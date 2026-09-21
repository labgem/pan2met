import logging
from collections import defaultdict
from collections.abc import Iterable

logger = logging.getLogger("pan2met")


def read_list(filename: str) -> list[str]:
    """
    Read a list of strings from a file, one per line
    """
    with open(filename, "r") as f:
        return f.read().splitlines()


def read_mapping(filename: str, sep="\t") -> dict[str, set[str]]:
    """
    Read a dictionary from a file.
    """
    mapping: dict[str, set[str]] = defaultdict(set)
    with open(filename, "r") as f:
        for row in f:
            parts = row.strip().split("\t")
            if len(parts) == 2:
                key, value = parts
                for item in value.split(","):
                    mapping[key].add(item)
    return dict(mapping)


def reverse_mapping(mapping: dict[str, set[str]]) -> dict[str, set[str]]:
    """
    Reverse a dictionnary. Set values as keys and keys as values.
    """
    reversed: dict[str, set[str]] = defaultdict(set)
    for key, values in mapping.items():
        for value in values:
            reversed[value].add(key)
    return dict(reversed)


def write_output(filename: str, content: Iterable[str]):
    """
    Write a list of strings to a file, one per line
    """
    with open(filename, "w") as output_file:
        output_file.writelines(line + "\n" for line in content)


def static_vars(**kwargs):
    """Decorate a function with local attributes"""

    def decorate(func):
        for key, value in kwargs.items():
            setattr(func, key, value)
        return func

    return decorate


def set_logging_level(verbose_intensity: int, logger=None):
    """
    Set the logging level based on the number of -v flags provided to the CLI.
    The more -v flags, the more verbose the logging level (up to DEBUG).

    :param verbose_intensity: The number of -v flags provided to the CLI

    verbose_intensity is mapped to logging levels as follows:
    0: ERROR
    1: WARNING
    2: INFO
    3 or more: DEBUG
    """
    logging.basicConfig()
    logging_levels = [logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR][::-1]
    logging_level_index = min(verbose_intensity, len(logging_levels) - 1)
    logging_level = logging_levels[logging_level_index]
    logging.getLogger().setLevel(logging_level)
    if logger is not None:
        logger.setLevel(logging_level)


def unquote(text: str) -> str:
    """
    Remove leading and trailing quotes from a string, if present
    """
    if text.startswith('"') and text.endswith('"'):
        return text[1:-1]
    return text
