from dataclasses import dataclass
from typing import Tuple, Optional
import typing


"""Custom classes"""

TreeA = None
TreeB = None


@dataclass
class TreeA:
    a: typing.Any
    history: Optional[Tuple[TreeB, TreeB]]

    def prob(self):
        if self.history is None:
            return [1]
        else:
            upper, lower = self.history
            return [p1 * p2 for p1 in upper.prob() for p2 in lower.prob()]


@dataclass
class TreeB:
    b: typing.Any
    history: Optional[Tuple[TreeA, TreeA]]

    def prob(self):
        pass
