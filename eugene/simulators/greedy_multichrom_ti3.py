from dataclasses import dataclass
from itertools import chain, product, starmap, zip_longest
from random import choice, randrange
from typing import List
from functools import reduce

from eugene.simulators import greedy_time as gt
from eugene.simulators.abc import BreedingProgram
from eugene.plant_models.plant2 import PlantSPC, Crossable, DomStrong


@dataclass(order=True)
class PlantMCh(Crossable, DomStrong):

    """Diploid plant with multiple homologous chromosome pairs and single-point
    crossover.

    n_loci[i] is the number of loci in the ith chromosome.
    """

    n_loci: List[int]
    chrom1: List[int]
    chrom2: List[int]

    def crosspoints(self):
        return product(
            *(product(range(2), range(n_loci_i)) for n_loci_i in self.n_loci)
        )

    def random_crosspoint(self):
        return [(randrange(2), randrange(n_loci_i)) for n_loci_i in self.n_loci]

    def gamete_specified(self, crosspoint):
        pass

    def cross_specified(self, other, crosspoints):
        pass

    def from_gametes(*args):
        pass

    def dom_weak(other):
        pass

    def dom_strong(other):
        pass


class BreedingProgramMultiChromTI(BreedingProgram):

    """Greedy algorithm to solve multi-chromosome trait introgression with
    single-point crossover

    This is done by extracting segments for each chromosome of the elite and
    donor, sorting them and ...
    TODO: Continue this train of thought
    """

    def __init__(self, n_loci: List[int], plant_type):
        """Constructor for BreedingProgramMultiChromTI

        Here n_loci is a list of ints representing the number of loci in either
        chromosome in each homologous pair."""
        BreedingProgram.__init__(self, n_loci, plant_type)

    def run(self):
        # TODO: test to see if initial population is a TI population
        # This should maybe go in to the init_pop
        # but for now go on the assumption that it is
        n_plants = len(self._pop_0)
        assert n_plants == 2
        assert is_ti_pop(self._pop_0)

        # Extract a table `segments` where segments[i][x] is a list of segments
        # in the ith chromosome of plant x in sorted order.
        segments = [
            [
                gt.segments_from_genotype(
                    n_loci_i, PlantSPC(n_loci_i, plant.chrom1[i], plant.chrom2[i])
                )
                for plant in self._pop_0
            ]
            for i, n_loci_i in enumerate(self._n_loci)
        ]

        # create zipped segments
        def f(n_loci_i, seg_a, seg_b):
            if seg_a is None and seg_b is None:
                raise ValueError("segments are both None")

            if seg_a is None:
                return seg_b
            if seg_b is None:
                return seg_a

            # under the assumption that the segments are disjoint, order them
            # cross them
            return segment_join(n_loci_i, seg_a, seg_b)

        segments_next_gen = [
            [f(n_loci_i, seg_a, seg_b) for seg_a, seg_b in zip_longest(*segments_i)]
            for segments_i, n_loci_i in zip(segments, self._n_loci)
        ]


def segment_join(n_loci, cx, cy):
    (sx, ex, gx) = cx
    (sy, ey, gy) = cy

    if sx > sy:
        return segment_join(cy, cx)
    elif sx == sy and ex < ey:
        return segment_join(cy, cx)
    elif sx == sy and ex >= ey:
        return cx

    assert sx < sy <= ex + 1 < ey
    return (sx, ey, PlantSPC(n_loci, gx, gy).gamete_specified((0, ex)))


def segment_is_joinable(cx, cy):
    (sx, ex, gx) = cx
    (sy, ey, gy) = cy
    return sx < sy <= ex + 1 < ey


def stack_trees(n_loci, trees_individual):
    return reduce(mappend, map(pure, trees_individual), mempty)


def single_chromsome_breeding_program_TI(n_loci, pop_0):
    from eugene.simulators.greedy_time import BreedingProgramGreedyTime

    bp = BreedingProgramGreedyTime(n_loci, PlantSPC)
    bp.set_ideotype(PlantSPC(n_loci, (1 << n_loci) - 1, (1 << n_loci) - 1,))
    bp.set_init_pop(pop_0)
    bp.run()


def is_ti_pop(pop_0):
    return True
