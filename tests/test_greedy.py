import unittest
from random import randrange, choice, random
from itertools import combinations
from eugene.simulators import greedy_time
from eugene.plant_models.plant2 import PlantSPC


class TestGreedy(unittest.TestCase):

    """Test case docstring."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_segments_from_gamete(self):
        g = int("1010110111", base=2)
        self.assertEqual(
            greedy_time.segments_from_gamete(10, g),
            [(0, 0, g), (2, 2, g), (4, 5, g), (7, 9, g)],
        )

    def test_segments_from_genotype(self):
        self.assertEqual(
            list(
                map(
                    lambda x: (x[0], x[1]),
                    greedy_time.segments_from_genotype(
                        10,
                        PlantSPC(
                            10, int("1001100111", base=2), int("0111001010", base=2),
                        ),
                    ),
                )
            ),
            [(0, 3), (1, 4), (6, 9)],
        )

        self.assertEqual(
            list(
                map(
                    lambda x: (x[0], x[1]),
                    greedy_time.segments_from_genotype(
                        10,
                        PlantSPC(
                            10, int("1011010111", base=2), int("0111111110", base=2),
                        ),
                    ),
                )
            ),
            [(0, 8), (1, 9)],
        )

    def test_min_segment_cover(self):
        def generate_gamete_for_segment_with_random_background(n_loci, s, e):
            assert 0 <= s <= e < n_loci
            x = randrange(1 << n_loci)
            # construct the favourable alleles in the segment
            fav_alleles = (1 << (e - s)) - 1
            # need to shift the first bit of fav_alleles to the correct position
            # i.e. e - s + Δ = n_loci - 1 - s
            #   -> Δ = n_loci - 1 - e
            x |= fav_alleles << (n_loci - 1 - e)
            # make unfavourable alleles surrounding the segment
            mask = (1 << n_loci) - 1
            if s > 0:
                mask ^= 1 << (n_loci - 1 - (s - 1))
            if e < n_loci - 1:
                mask ^= 1 << (n_loci - 1 - (e + 1))
            return (s, e, x & mask)

        def generate_gamete_for_segment_no_background(n_loci, s, e):
            assert 0 <= s <= e < n_loci
            # construct the favourable alleles in the segment
            fav_alleles = (1 << (e - s + 1)) - 1
            # need to shift the first bit of fav_alleles to the correct position
            # i.e. e - s + Δ = n_loci - 1 - s
            #   -> Δ = n_loci - 1 - e
            x = fav_alleles << (n_loci - 1 - e)
            return (s, e, x)

        def gen_random_segments(n_loci):
            in_seg = False
            for i in range(n_loci):
                b = random() < 0.5
                if not in_seg and b:
                    s = i
                if b:
                    e = i
                if in_seg and not b:
                    yield (s, e)
                in_seg = b
            if in_seg:
                yield (s, e)

        non_empty = lambda f: (lambda res: res if res else non_empty(f))(f())

        n_loci = 14

        # generate instance
        def generate_instance(n_loci):
            segments = []
            collected_loci = 0
            while collected_loci != (1 << n_loci) - 1:
                s, e = choice(non_empty(lambda: list(gen_random_segments(n_loci))))
                c = s, e, g = generate_gamete_for_segment_no_background(n_loci, s, e)
                collected_loci |= g
                segments.append(c)
            return segments

        def solve_naive(segments):
            n_segs = len(segments)
            # feasibility check
            z = 0
            for s, e, g in segments:
                z |= g
            for k in range(1, n_segs):
                for selection in combinations(segments, k):
                    # test return k if selection covers
                    z = 0
                    # assuming g is created without background
                    for s, e, g in selection:
                        z |= g
                    if z == (1 << n_loci) - 1:
                        return list(selection)

        for _ in range(20):
            inst = generate_instance(n_loci)
            check = solve_naive(inst)
            guess = greedy_time.min_segment_cover(n_loci, inst)
            self.assertEqual(len(guess), len(check))
