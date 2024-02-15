from dataclasses import dataclass
from collections import defaultdict

from .abc import BreedingProgram
from ..plant_models.plant2 import PlantSPC, PlantSPCBitarray


class BreedingProgramGreedyTime(BreedingProgram):
    """
    Breeding program solver that uses the favourable segments of alleles in
    gametes to determine the set of crossovers that achieve the target in the
    fewest number of generations
    """

    def __init__(self, n_loci, plant_type):
        BreedingProgram.__init__(self, n_loci, plant_type)

    def run(self):
        @dataclass
        class Genotype:
            plant: self._plant_type
            histories: tuple

            def crosspoints(self):
                return self.plant.crosspoints()

            def gamete_specified(self, crosspoint):
                return self.plant.gamete_specified(crosspoint)

            def format(self):
                return str(self.plant)

        @dataclass
        class Gamete:
            gamete: int
            histories: list

            def __str__(g):
                return format(g.gamete, f"0{self._n_loci}")

        self._Genotype = Genotype
        self._Gamete = Gamete

        # Special case: target is in pop_0
        if self._ideotype in self._pop_0:
            return Genotype(self._ideotype, None)

        # Wrap genotypes
        pop = [Genotype(x, None) for x in self._pop_0]

        # TODO: for histories
        # # construct the set of segments the initial population
        # wrapped_gamete_store = {}
        # segments = []
        # for x in pop:
        #     for s, e, g in segments_from_genotype(self._n_loci, x.plant):
        #         if g not in wrapped_gamete_store:
        #             wrapped_gamete_store[g] = Gamete(g, [])
        #         wg = wrapped_gamete_store[g]
        #         wg.histories.append(x)
        #         segments.append((s, e, wg))

        segments_best = min_segment_cover(
            self._n_loci,
            (c for x in pop for c in segments_from_genotype(self._n_loci, x)),
        )
        assert segments_best == sorted(segments_best)

        # construct dag
        def aux(segments_best):
            if len(segments_best) == 1:
                return segments_best[0]

            segments_best_l = segments_best[: len(segments_best) // 2]
            segments_best_r = segments_best[len(segments_best) // 2 :]

            sl, el, wgl = aux(segments_best_l)
            sr, er, wgr = aux(segments_best_r)

            new_gen = Genotype(
                self._plant_type(self._n_loci, wgl.gamete, wgr.gamete), (wgl, wgr)
            )
            new_gam = Gamete(new_gen.gamete_specified((0, sr)), [new_gen])
            new_seg = (sl, er, new_gam)
            return new_seg

        # # Shortcut tree evaluation
        # return Gamete(
        #     PlantSPC(
        #         self._n_loci,
        #         (1 << self._n_loci) - 1,
        #         (1 << self._n_loci) - 1,
        #     ),
        #     None
        # )

        _, _, wgam = aux(segments_best)
        ideotype = Genotype(self._ideotype, (wgam, wgam))
        return ideotype


def segments_from_gamete(n_loci: int, gamete: int) -> list:
    """
    Extracts a list of segments from a raw gamete.
    Each segment is represented as a tuple (s, e, g) where:
        s - start point
        e - end point
        g - gamete of origin
    Note: start and end points are 0-indexed

    >>> extract_segments_from_gamete(10, int('0111000110', base=2))
    [(1, 3, 454), (7, 8, 454)]
    """
    # assert 0 <= gamete, f"gamete is a negative number {gamete}"
    # assert gamete < 1 << n_loci, f"gamete has more than {n_loci} loci {gamete}"
    out = [None] * ((n_loci + 1) // 2)
    out_i = 0
    in_segment = False
    s = e = None
    # traverse alleles in reverse order
    for i in range(n_loci):
        allele = (gamete >> i) & 1  # for bigints
        # allele = gamete[i]  # for bitarrays
        if allele:
            if not in_segment:
                e = n_loci - 1 - i
            s = n_loci - 1 - i
        elif not allele and in_segment:
            out[out_i] = (s, e, gamete)
            out_i += 1
        in_segment = bool(allele)
    if in_segment:
        out[out_i] = (s, e, gamete)
        out_i += 1
    return [out[i] for i in reversed(range(out_i))]


def segments_from_genotype(n_loci: int, genotype: PlantSPC) -> list:
    """
    Returns a list of all segments extractable from a given genotype.
    The segments are sorted by start point.
    All segments are maximal.
    All segments are unique in their ranges i.e. no two segments share both
    start and end point.
    """
    q1 = segments_from_gamete(n_loci, genotype.chrom1)
    q2 = segments_from_gamete(n_loci, genotype.chrom2)

    used1 = [False] * len(q1)
    used2 = [False] * len(q2)

    i = j = 0
    out = [None] * max(len(q1) + len(q1), n_loci - 1)
    out_i = 0

    # This works.
    while i < len(q1) and j < len(q2):
        c1 = (s1, e1, g1) = q1[i]
        c2 = (s2, e2, g2) = q2[j]

        if s1 > s2 or s1 == s2 and e1 < e2:
            q1, q2 = q2, q1
            used1, used2 = used2, used1
            i, j = j, i
        elif s1 <= s2 and e1 >= e2:
            j += 1
        elif e1 + 1 < s2:
            # assert s1 < e1 + 1 < s2 <= e2
            if not used1[i]:
                out[out_i] = c1
                out_i += 1
                used1[i] = True
            i += 1
        elif e2 == n_loci - 1:
            out[out_i] = (
                s1,
                e2,
                PlantSPC(n_loci, g1, g2).gamete_specified((0, e1 + 1)),
            )
            # out[out_i] = (s1, e2, 0)
            out_i += 1
            used1[i] = used2[j] = True
            # break and don't use any more segments
            i = len(q1)
            j = len(q2)
        else:
            # assert e2 < n_loci - 1
            out[out_i] = (
                s1,
                e2,
                PlantSPC(n_loci, g1, g2).gamete_specified((0, e1 + 1)),
            )
            # out[out_i] = (s1, e2, 0)
            out_i += 1
            used1[i] = used2[j] = True
            i += 1

    out = out[:out_i]
    return out + q1[i:] + q2[j:]  # works because i == len(q1) or j == len(q2)


def min_segment_cover(n_loci, segments):
    """
    Given a list of segments, represented by triples [(s_x, e_x, g^x), ...],
    returns a minimum cardinality subset of segments that covers [n_loci].
    """
    out = [None] * n_loci

    for (s, e, g) in segments:
        if out[s] is None or e > out[s][1]:
            out[s] = (s, e, g)

    e_covered = -1
    j = 0
    i = 0
    while i < n_loci and e_covered < n_loci - 1:
        # INV: j == e_covered + 1
        c_next = s_next, e_next, g_next = out[i]
        while i < n_loci and out[i][0] <= e_covered + 1:
            c = s, e, g = out[i]
            if s <= e_covered + 1 <= e and e > e_next:
                c_next = s_next, e_next, g_next = c
            i += 1
        out[j] = c_next
        j += 1
        e_covered = e_next
    return out[:j]


def min_segment_cover_key(n_loci, segments, key):
    """
    Given a list of segments, represented by triples [(s_x, e_x, g^x), ...],
    returns a minimum cardinality subset of segments that covers [n_loci].
    """
    out = [None] * n_loci

    for c in segments:
        (s, e, g) = key(c)
        if out[s] is None or e > key(out[s])[1]:
            out[s] = c

    e_covered = -1
    j = 0
    i = 0
    while i < n_loci and e_covered < n_loci - 1:
        # INV: j == e_covered + 1
        c_next = out[i]
        s_next, e_next, g_next = key(out[i])
        while i < n_loci and (out[i] is None or key(out[i])[0] <= e_covered + 1):
            if out[i] is None:
                i += 1
                continue
            c = out[i]
            s, e, g = key(out[i])
            if s <= e_covered + 1 <= e and e > e_next:
                c_next = c
                s_next, e_next, g_next = key(c)
            i += 1
        out[j] = c_next
        j += 1
        e_covered = e_next
    return out[:j]


def min_segment_cover_old(n_loci, segments):
    segments.sort(key=lambda seg: (seg[0], -seg[1]))
    n_segments = len(segments)
    e_covered = -1
    if n_segments == 0:
        raise ValueError(f"Segments list does not cover [n_loci], namely {e_covered+1}")

    i = 0
    out = []
    while i < n_segments and e_covered < n_loci - 1:
        # INV loci 0..e are covered by `out`
        c_next = s_next, e_next, g_next = segments[i]
        if e_covered + 1 < s_next:
            raise ValueError(
                f"Segments list does not cover [n_loci], namely {e_covered+1}"
            )

        # under the assumption that s_next <= e_covered + 1, choose c where
        # c = segment with largest end point that contains `e_covered + 1`
        while i < n_segments and segments[i][0] <= e_covered + 1:
            c = s, e, g = segments[i]
            if s <= e_covered + 1 <= e and e > e_next:
                c_next = s_next, e_next, g_next = c
            i += 1
        if e_next <= e_covered:
            raise ValueError(
                f"Segments list does not cover [n_loci], namely {e_covered+1}"
            )
        out.append(c_next)
        e_covered = e_next
    return out
