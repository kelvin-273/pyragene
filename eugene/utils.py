from random import randint, randrange, sample, choice
from typing import List, Optional, NewType, Callable, Tuple
from collections import namedtuple


def count_ones(x: int):
    """
    Returns the number of one bits in the integer
    """
    if x < 0:
        raise ValueError("Gave negative int and don't know how to deal with this yet")
    out = 0
    while x > 0:
        x, r = x >> 1, x & 1
        out += r
    return out
