import argparse
import sys
import numpy as np
import pandas as pd
from bitarray import frozenbitarray
from bitarray.util import zeros
from tqdm import tqdm
from random import seed, choices, randrange, shuffle
from time import time
from eugene.plant_models.plant2 import PlantSPC, PlantSPCBitarray
from eugene.simulators import greedy_time
from eugene.simulators import enumerators

"""
We want to benchmark the greedy algorithm against the enumerative approach.
There will come the question of why we don't benchmark against CANZAR.
- CANZAR assumes that the number of crossovers is probabilistic
- for us the number of crossovers is constant

How do we show that greedy is better?
- average relative difference
- average solving times
- worst case examples
    - singles
"""


def benchmark_random_greedy_bigint(n_loci, n_pop, n_iterations=100):
    print(n_loci, file=sys.stderr)
    results = [0] * n_iterations
    seed(0)
    # for i in range(n_iterations):
    for i in tqdm(range(n_iterations)):
        bp = greedy_time.BreedingProgramGreedyTime(n_loci, PlantSPC)
        bp.set_init_pop(PlantSPC.initial_pop_random(n_loci, n_pop))
        bp.set_ideotype(PlantSPC(n_loci, (1 << n_loci) - 1, (1 << n_loci) - 1))
        start = time()
        bp.run()
        end = time()
        delta = end - start
        results[i] = delta
    print(f"Res\tGreedy\t{n_loci}\t{sum(results) / n_iterations}\t{np.std(results)}")


def benchmark_random_greedy(n_loci, n_pop, n_iterations=100):
    print(n_loci, file=sys.stderr)
    results = [0] * n_iterations
    seed(0)
    # for i in range(n_iterations):
    for i in tqdm(range(n_iterations)):
        bp = greedy_time.BreedingProgramGreedyTime(n_loci, PlantSPCBitarray)
        bp.set_init_pop(PlantSPCBitarray.initial_pop_random(n_loci, n_pop))
        bp.set_ideotype(
            PlantSPCBitarray(
                n_loci, frozenbitarray(~zeros(n_loci)), frozenbitarray(~zeros(n_loci))
            )
        )
        start = time()
        bp.run()
        end = time()
        delta = end - start
        results[i] = delta
    print(f"Res\tGreedy\t{n_loci}\t{sum(results) / n_iterations}\t{np.std(results)}")


def benchmark_random_enumdom(n_loci, n_pop, n_iterations=100):
    results = [0] * n_iterations
    seed(0)
    # for i in tqdm(range(n_iterations)):
    for i in range(n_iterations):
        bp = enumerators.BreedingProgramDom(n_loci, PlantSPC)
        bp.set_init_pop(PlantSPC.initial_pop_random(n_loci, n_pop))
        bp.set_ideotype(PlantSPC(n_loci, (1 << n_loci) - 1, (1 << n_loci) - 1))
        start = time()
        bp.run()
        end = time()
        delta = end - start
        results[i] = delta
    print(f"Res\tEnumDom\t{n_loci}\t{sum(results) / n_iterations}\t{np.std(results)}")


def benchmark_random_enum(n_loci, n_pop, n_iterations=100):
    results = [0] * n_iterations
    seed(0)
    # for i in tqdm(range(n_iterations)):
    for i in range(n_iterations):
        bp = enumerators.BreedingProgram(n_loci, PlantSPC)
        bp.set_init_pop(PlantSPC.initial_pop_random(n_loci, n_pop))
        bp.set_ideotype(PlantSPC(n_loci, (1 << n_loci) - 1, (1 << n_loci) - 1))
        start = time()
        bp.run()
        end = time()
        delta = end - start
        results[i] = delta
    print(f"Res\tEnum\t{n_loci}\t{sum(results) / n_iterations}\t{np.std(results)}")


# for n_loci in range(1, 11):
#     benchmark_random_enum(n_loci, 6, 1000)
# for n_loci in range(1, 16):
#     benchmark_random_enumdom(n_loci, 6, 1000)
# for n_loci in range(101, 10001):
#     benchmark_random_greedy(n_loci, 6, 1000)
# benchmark_random_greedy(10_000, 6, 1000)


class RangeSet:

    """Data structure for maintaining a compressed set of elements"""

    def __init__(self, lbound, ubound, remove_list=None):
        """Initialises the set to be a complete set from lbound to ubound inclusive.

        :lbound: Upper bound of the set
        :ubound: Lower bound of the set
        """
        self._lbound = lbound
        self._ubound = ubound
        self._card = max(0, ubound - lbound)
        if remove_list is None:
            self.ranges = {lbound: (lbound, ubound)}
        else:
            self.ranges = {}
            remove_list = sorted(remove_list)
            lb = lbound
            for rx in remove_list:
                if lb < rx:
                    self.ranges[lb] = (lb, rx)
                lb = max(lb, rx + 1)
            if lb < ubound:
                self.ranges[lb] = (lb, ubound)

    def sample(self):
        if len(self) == 0:
            return None
        lb, ub = choices(
            list(self.ranges.values()),
            weights=[ub - lb for (lb, ub) in self.ranges.values()],
        )[0]
        return randrange(lb, ub)

    def remove(self, x):
        r = self.find(x)
        if r is None:
            return None
        lb, ub = r
        if len(r) == 1:
            self.ranges.pop(lb)
        elif x == lb:
            self.ranges.pop(lb)
            self.ranges[lb + 1] = (lb + 1, ub)
        elif x == ub - 1:
            self.ranges[lb] = (lb, ub - 1)
        else:
            self.ranges[lb] = (lb, x)
            self.ranges[x + 1] = (x + 1, ub)
        self._card -= 1

    def find(self, x):
        for lb, ub in self.ranges.values():
            if lb <= x < ub:
                return (lb, ub)

    def __len__(self):
        return self._card

    def sample_remove(self):
        x = self.sample()
        if x is not None:
            self.remove(x)
        return x


class RangeSet2:
    """Keeps a shuffled list of elements"""

    def __init__(self, lbound, ubound, remove_list=None):
        self.lbound = lbound
        self.ubound = ubound
        if remove_list is not None:
            self.range_set = set(range(lbound, ubound))
            self.range_set.difference_update(set(remove_list))
            self.range_set = list(self.range_set)
        else:
            self.range_set = list(range(lbound, ubound))
        shuffle(self.range_set)
        self.idx = 0

    def __len__(self):
        return len(self.range_set) - self.idx

    def sample_remove(self):
        if self.idx == len(self.range_set):
            return None
        res = self.range_set[self.idx]
        self.idx += 1
        return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some command line arguments")

    parser.add_argument("-b", "--bounds", nargs=2, type=int, help="range of numbers")
    parser.add_argument("-f", "--files", nargs="+", type=str, help="list of file names")

    args = parser.parse_args()

    if not args.bounds and not args.files:
        benchmark_random_greedy(20_000, 6, 1000)
        exit(0)

    remove_list = set()

    if args.files:
        # print(f"Filenames: {args.files}")
        for filename in args.files:
            df = pd.read_csv(filename, header=None, sep="\t")
            try:
                column2 = df.iloc[:, 2]
                remove_list.update(column2)
            except IndexError:
                pass

    if args.bounds:
        # print(f"Range of numbers: {range(args.bounds[0], args.bounds[1])}")
        lbound = args.bounds[0]
        ubound = args.bounds[1] + 1
    else:
        lbound = min(remove_list)
        ubound = max(remove_list)

    range_set = RangeSet2(lbound, ubound, remove_list=remove_list)

    # define iterator
    def gen_samples():
        while len(range_set) > 0:
            yield range_set.sample_remove()

    # define pickleable function
    def f(n_loci):
        return benchmark_random_greedy(n_loci, 6, 1_000)

    # create a pool of 5 processes
    # import multiprocessing
    # with multiprocessing.Pool(processes=5) as pool:
    #     # apply the function to each element of the iterable using imap
    #     result_iterator = pool.imap(f, gen_samples())
    #     # process the results as they become available
    #     for _ in result_iterator:
    #         pass
    for n_loci in gen_samples():
        f(n_loci)
