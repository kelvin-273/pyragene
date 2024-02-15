"""
High-Low solving

Enumerates gamete subsets (high level) and solves the distribute instance
that each gamete subset corresponds to (low level).
"""


from math import ceil
from dataclasses import dataclass
from collections import namedtuple
from eugene.simulators.greedy_time import segments_from_genotype, min_segment_cover


def gen_covering_subsets(n_loci: int, segments: list):
    """
    Give an integer `n_loci` and a list of `segments`, returns a generator that
    yields subsets of segments that cover the loci 0..n_loci-1.
    """
    n = len(segments)
    segments.sort()
    selection = [False] * n

    def aux(i, current_end):
        if i == n or current_end == n_loci:
            yield [seg for seg, b in zip(segments, selection) if b]
        elif segments[i][0] > current_end:
            return
        else:
            yield from aux(i + 1, current_end)
            if segments[i][1] > current_end:
                selection[i] = True
                yield from aux(i + 1, max(current_end, segments[i][1]))
                selection[i] = False

    return aux(0, 0)


def breeding_program(n_loci: int, pop0: list) -> int:
    # generate segments
    segments = []
    for x in pop0:
        segments.extend(segments_from_genotype(n_loci, x))
    n = len(segments)
    segments.sort()
    selection_seg = [False] * n

    # array for counting the number of times each parent gets used
    map_s = {g for (_, _, g) in segments}
    map_d = {g: i for (i, g) in enumerate(map_s)}
    map_l = list(map_s)
    m = len(map_s)
    segments = [(s, e, map_d[g]) for (s, e, g) in segments]
    selection_gam = [0] * m
    # struct to hold the number of segments and gametes used
    used = [0, 0]
    # struct to hold the lower and upper bound
    bounds = [-float("inf"), float("inf")]

    # enumerate covering subsets
    def aux(i, current_end):
        if i == n or current_end == n_loci:
            yield [seg for seg, b in zip(segments, selection_seg) if b]
            bounds[1] = min(bounds[1], used[0]*2)
        elif segments[i][0] > current_end:
            return
        else:
            yield from aux(i + 1, current_end)
            (s, e, gi) = segments[i]
            if e > current_end \
                    and ceil(sum(used) / 2) <= bounds[1]:
                selection_seg[i] = True
                selection_gam[gi] += 1
                used[0] += 1
                used[1] += selection_gam[gi] == 1

                yield from aux(i + 1, max(current_end, segments[i][1]))

                selection_seg[i] = False
                selection_gam[gi] -= 1
                used[0] -= 1
                used[1] -= selection_gam[gi] == 0

    return aux(0, 0)
