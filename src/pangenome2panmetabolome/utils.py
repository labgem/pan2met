#!/usr/bin/env python3

import logging

logger = logging.getLogger("pangenome2panmetabolome")
logger.setLevel(logging.DEBUG)
logging.basicConfig(level=logging.DEBUG)


def static_vars(**kwargs):
    """Decorate a function with local attributes"""

    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func

    return decorate
