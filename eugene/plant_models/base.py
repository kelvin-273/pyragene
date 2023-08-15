from bitarray import frozenbitarray
from dataclasses import dataclass
from .abstract_plants import Genotype, Dominance, Diploid


@dataclass
class BaseGenotype(Genotype, Diploid):
    n_loci: int
    chrom1: frozenbitarray
    chrom2: frozenbitarray


@dataclass
class BaseGamete:
    n_loci: int
    gamete: frozenbitarray


class DomBaseGenotype(Dominance):
    @staticmethod
    def dom(x: BaseGenotype, y: BaseGenotype):
        return (
            DomBaseGamete.dom(x.upper(), y.upper())
            and DomBaseGamete.dom(x.lower(), y.lower())
            or DomBaseGamete.dom(x.upper(), y.lower())
            and DomBaseGamete.dom(x.lower(), y.upper())
        )


class DomBaseGamete(Dominance):
    @staticmethod
    def dom(gx: BaseGamete, gy: BaseGamete):
        return (gx.gamete | ~gy.gamete).all()
