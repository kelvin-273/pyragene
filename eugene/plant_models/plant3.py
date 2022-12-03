"""
A set of data structures for plant breeding.

Functional patterns:
    pure :: Gamete -> TreeGamete
    pure :: Genotype -> TreeGenotype

    map :: (gamete -> gamete) -> TreeGamete -> TreeGamete
    map :: (genotype -> genotype) -> TreeGenotype -> TreeGenotype

    map :: (genotype -> genotype) -> TreeGenotype -> TreeGenotype
"""


class Gamete:
    def __init__(self, gamete):
        self._gamete = gamete


class Genotype:
    def __init__(self, chromosomes):
        self._chromosomes = chromosomes


class TreeGamete:
    def __init__(self, gamete, history):
        self._gamete = gamete
        self.history = history


class TreeGenotype:
    def __init__(self, genotype, history):
        self._genotype = genotype
        self.history = history


class Segment:
    def __init__(self, start, end, gamete):
        self.s = start
        self.e = end
        self.g = gamete
