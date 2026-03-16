from typing import Iterable

import logging

logger = logging.getLogger("pan2met")


def read_list(filename: str) -> list[str]:
    with open(filename, "r") as f:
        return f.read().splitlines()


def write_output(filename: str, content: Iterable[str]):
    with open(filename, "w") as output_file:
        output_file.writelines(line + "\n" for line in content)


def static_vars(**kwargs):
    """Decorate a function with local attributes"""

    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func

    return decorate
