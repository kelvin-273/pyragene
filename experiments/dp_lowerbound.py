from functools import lru_cache
from random import random, randrange
from collections import defaultdict
from itertools import takewhile, dropwhile


def count_ones(n):
    assert n >= 0
    count = 0
    while n > 0:
        n, r = n >> 1, n & 1
        count += r
    return count


def time_lowerbound_naive(n_loci, ranges) -> int:
    """Computes the lower bound time required for the breeding program given
    the runs of favourable alleles.

    @param ranges: A list tuples [(start, end)]
    @return: the minimum number of ranges to cover {0..n-1}

    >>> time_lowerbound_naive(16, [])
    inf
    >>> time_lowerbound_naive(16, [(0, 15)])
    inf
    >>> time_lowerbound_naive(16, [(0, 16)])
    1
    >>> time_lowerbound_naive(16, [\
        (6,8), (3, 9), (10, 13), (1, 5), (7, 11), (5, 12), (7, 11), (9, 14),\
        (9, 16), (15, 16), (0, 1)])
    4
    """
    if not all(0 <= s < e <= n_loci for s, e in ranges):
        raise ValueError("One of the ranges is invalid")

    smallest_subset_size = float("inf")
    for selection_int in range(1 << len(ranges)):
        # does some stuff
        if all(
            any(
                (selection_int >> j) & 1 and i in range(s, e)
                for j, (s, e) in enumerate(ranges)
            )
            for i in range(n_loci)
        ):
            smallest_subset_size = min(smallest_subset_size, count_ones(selection_int))
    return smallest_subset_size


def time_lowerbound_dp(n_loci, ranges):
    """Computes the lower bound time required for the breeding program given
    the runs of favourable alleles.

    @param ranges: A list tuples [(start, end)]
    @return: the minimum number of ranges to cover {0..n-1}

    >>> time_lowerbound_dp(16, [])
    inf
    >>> time_lowerbound_dp(16, [(0, 15)])
    inf
    >>> time_lowerbound_dp(16, [(0, 16)])
    1
    >>> time_lowerbound_dp(16, [\
        (6,8), (3, 9), (10, 13), (1, 5), (7, 11), (5, 12), (7, 11), (9, 14),\
        (9, 16), (15, 16), (0, 1)])
    4
    """
    if not all(0 <= s < e <= n_loci for s, e in ranges):
        raise ValueError("One of the ranges is invalid")

    ranges.sort(key=lambda x: (x[0], -x[1]))

    @lru_cache
    def aux(covered_loci, i):
        if covered_loci == n_loci:
            return 0

        if i == len(ranges) and covered_loci < n_loci:
            return float("inf")

        s, e = ranges[i]
        if s > covered_loci:
            return float("inf")

        if e <= covered_loci:
            return aux(covered_loci, i+1)

        return min(
            1 + aux(e, i+1),
            aux(covered_loci, i+1)
        )

    return aux(0, 0)


def time_lowerbound_shortest_path_pierre(n_loci, ranges):
    g = defaultdict(lambda: set())


def time_lowerbound_sorting(n_loci, ranges):
    if not all(0 <= s < e <= n_loci for s, e in ranges):
        raise ValueError("One of the ranges is invalid")

    if n_loci == 0:
        return 0

    n = len(ranges)
    if n == 0:
        return float('inf')

    ranges.sort(key=lambda x: (x[0], -x[1]))

    i = j = 0
    count = 0
    covered_loci = 0

    while i < n and covered_loci < n_loci:
        next_endpoint = covered_loci
        while j < n and ranges[j][0] <= covered_loci:
            _, e = ranges[j]
            if e > next_endpoint:
                next_endpoint = e
            j += 1

        if j == i:
            return float('inf')
        count += 1
        covered_loci = next_endpoint
        i = j

    if covered_loci == n_loci:
        return count
    else:
        return float('inf')


def random_testcase(n_loci, f_range_generator):
    out = [f_range_generator(n_loci) for _ in range(n_loci)]
    return out


def f_range_generator_geometric(p):
    def aux(n_loci):
        run_len = 1
        while random() < p:
            run_len += 1
        start = randrange(n_loci - run_len + 1)
        end = start + run_len
        return (start, end)
    return aux


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=False)

    n_loci = 16
    for i in range(100):
        test_case = random_testcase(n_loci, f_range_generator_geometric(0.5))
        print(i)
        # assert time_lowerbound_naive(n_loci, test_case) == \
            # time_lowerbound_dp(n_loci, test_case)
        # assert time_lowerbound_dp(n_loci, test_case) == \
            # time_lowerbound_sorting(n_loci, test_case)
        assert time_lowerbound_dp(n_loci, test_case) == \
            time_lowerbound_sorting(n_loci, test_case), \
            f"{time_lowerbound_dp(n_loci, test_case)} != {time_lowerbound_sorting(n_loci, test_case)}\n" + \
            f"n_loci: {n_loci}, test_case:\n{test_case}"
