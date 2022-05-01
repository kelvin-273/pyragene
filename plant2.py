from collections import Counter
from dataclasses import dataclass
from random import randint, randrange, sample
from typing import NewType
from abc import ABC, abstractmethod, abstractclassmethod


Chromosome = NewType("Chromosome", int)


@dataclass
class Plant(ABC):

    """Abstract class for plants"""

    @abstractmethod
    def crosspoints(self):
        pass

    @abstractmethod
    def gamete_specified(self, crosspoint):
        pass

    @abstractmethod
    def cross_specified(self, other, crosspoints):
        pass

    @abstractmethod
    def random_crosspoint(self):
        pass

    def gamete_random(self):
        crosspoint = self.random_crosspoint()
        return self.gamete_specified(crosspoint)

    @abstractmethod
    def cross_random(self, other):
        pass


@dataclass(order=True)
class PlantSPC(Plant):

    """Plant with single-point crossover"""

    n_loci: int
    chrom1: int
    chrom2: int

    def cross_specified(self, other, crosspoints):
        k1, k2 = crosspoints
        g1 = self.gamete_specified(k1)
        g2 = other.gamete_specified(k2)
        return PlantSPC(self.n_loci, g1, g2)

    def random_crosspoint(self):
        """
        Generates a random crosspoint as a tuple (start, split) where
        start in {0, 1} is the chromosome whose prefix forms the prefix of the
        gamete and split in {0..n-1} (inclusive) is the locus up to which the prefix of
        the start chromosome is created from.
        Mainly for use with random gamete creation.
        """
        return (randint(0, 1), randrange(0, self.n_loci))

    def gamete_specified(self, crosspoint):
        """
        Creates a gamete with crosspoint = (start, split)
        """
        assert type(crosspoint) is tuple and len(crosspoint) == 2
        assert type(crosspoint[0]) is int and type(crosspoint[1]) is int
        start, split = crosspoint
        # bit-string representing the source of the ith allele ∀i∈[n]
        mask = (1 << self.n_loci + 1) - 1
        crossing = (mask >> split) << split
        if not start:
            crossing ^= mask

        return (crossing & self.chrom1) | (~crossing & self.chrom2)

    def cross_random(self, other):
        return PlantSPC(
            self.n_loci,
            self.gamete_random(),
            other.gamete_random()
        )

    def dom_strong(self, other):
        """
        Returns True iff every potential gamete of self dominates all
        potential gametes of other.
        """
        assert self.n_loci == other.n_loci
        return all(
            self.chrom1 | other.chrom1 == self.chrom1,
            self.chrom1 | other.chrom2 == self.chrom1,
            self.chrom2 | other.chrom1 == self.chrom2,
            self.chrom2 | other.chrom2 == self.chrom2,
        )

    def dom_weak(self, other: "Plant"):
        """
        Returns True iff the upper and lower chromosomes of other do not have
        favourable alleles that are not already in the respective chromosome of
        self.
        """
        assert self.n_loci == other.n_loci
        return (
            self.chrom1 | other.chrom1 == self.chrom1
            and self.chrom2 | other.chrom2 == self.chrom2
        )

    def reachable_gametes(self):
        """
        Returns the set of gametes that can be created from this plant.
        """
        d = Counter(
            self.create_gamete_with_crosspoint(i, j)
            for i in range(2)
            for j in range(self.n_loci)
        )
        return set(d.keys())

    def reachable_gametes_with_counts(self):
        d = Counter(
            self.create_gamete_with_crosspoint(i, j)
            for i in range(2)
            for j in range(self.n_loci)
        )
        return list(d.items())

    def crosspoints(self):
        for i in range(2):
            for j in range(self.n_loci):
                yield (i, j)

    def __hash__(self):
        return (self.chrom1 << self.n_loci) | self.chrom2

    @staticmethod
    def initial_pop_trait_introgression(n_loci: int, n_holes: int, n_elite=1, n_donor=1):
        """
        Creates an initial population for trait introgression.
        Output list is divided into [elite_pop : donor_pop].
        Elite and donor populations each have a complimentary set of homozygous
        loci and loci not in that set are hetrozygous or homozygous
        unfavourable.
        """
        assert n_holes <= n_loci // 2
        holes = sample(range(n_loci), n_holes)

        bitbed = (1 << n_loci) - 1
        elite_mask = bitbed
        for i in holes:
            elite_mask ^= 1 << i

        def gen_elite(mask_elite: int):
            xu = xl = (1 << n_loci) - 1
            # get bits that need fixing
            d = (~mask_elite) & xu & xl
            while d:
                u = randrange(1 << n_loci)
                l = randrange(1 << n_loci)
                xu = ((~d) & xu) | (d & u)
                xl = ((~d) & xl) | (d & l)
                d = (~mask_elite) & xu & xl
            return PlantSPC(n_loci, xu, xl)

        elite_pop = [gen_elite(elite_mask) for _ in range(n_elite)]
        donor_pop = [gen_elite(bitbed ^ elite_mask) for _ in range(n_elite)]

        return elite_pop + donor_pop

    @staticmethod
    def initial_pop_random(n_loci: int, n_individuals: int):
        """
        Generates a population of uniform randomly generated plants such that
        the ideotype is in its span.
        """
        out = [PlantSPC(
            n_loci,
            randrange(1 << n_loci),
            randrange(1 << n_loci)
        ) for _ in range(n_individuals)]
        while PlantSPC.union(out) != (1 << n_loci) - 1:
            out = [PlantSPC(
                n_loci,
                randrange(1 << n_loci),
                randrange(1 << n_loci)
            ) for _ in range(n_individuals)]
        return out

    @staticmethod
    def initial_pop_singles_homo(n_loci: int):
        return [PlantSPC(n_loci, 1 << i, 1 << i) for i in range(n_loci)]

    @staticmethod
    def initial_pop_singles_hetero(n_loci: int):
        return [PlantSPC(n_loci, 1 << i, 0) for i in range(n_loci)] + \
            [PlantSPC(n_loci, 0, 1 << i) for i in range(n_loci)]

    @staticmethod
    def union(pop) -> int:
        """
        Takes the union of all chromosomes in the population.
        Useful for testing feasibility.

        :pop: List of individual PlantSPC's
        :returns: integer containing the bits after union of all chromosomes in the population
        """
        out = 0
        for x in pop:
            out |= x.chrom1 | x.chrom2
        return out


if __name__ == "__main__":
    x = PlantSPC(4, 5, 10)
