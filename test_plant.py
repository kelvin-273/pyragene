import unittest
import random

import plant
from plant import Plant


class TestPlant(unittest.TestCase):

    """Test case docstring."""

    def test_random_plant_generates_gamete_with_alleles_from_self(self):
        """
        Checking that all alleles are taken from their respective loci.
        """
        for _ in range(1000):
            n_loci = random.randrange(1, 10)

            x = plant.generate_random_plant(n_loci)
            su = format(x.chrom1, f"0{n_loci}b")
            sl = format(x.chrom2, f"0{n_loci}b")
            cs = format(x.create_gamete(), f"0{n_loci}b")
            for i in range(n_loci):
                assert cs[i] == su[i] or cs[i] == sl[i]

    def test_generates_all_gametes_simple(self):
        a = Plant(4, 5, 10)
        # self.assertEqual(a.create_gamete_with_crosspoint())
        self.assertEqual(a.create_gamete_with_crosspoint(0, 0), 10)
        self.assertEqual(a.create_gamete_with_crosspoint(0, 1), 11)
        self.assertEqual(a.create_gamete_with_crosspoint(0, 2), 9)
        self.assertEqual(a.create_gamete_with_crosspoint(0, 3), 13)
        self.assertEqual(a.create_gamete_with_crosspoint(1, 0), 5)
        self.assertEqual(a.create_gamete_with_crosspoint(1, 1), 4)
        self.assertEqual(a.create_gamete_with_crosspoint(1, 2), 6)
        self.assertEqual(a.create_gamete_with_crosspoint(1, 3), 2)
        

    def test_prob_z_given_xy_fast(self):
        """
        The two algorithms for calculating the probability of z given x and y
        are the same.
        """
        for _ in range(1000):
            n_loci = random.randrange(1, 10)
            x = plant.generate_random_plant(n_loci)
            y = plant.generate_random_plant(n_loci)
            z = plant.generate_random_plant(n_loci)
            p1 = plant.prob_z_given_xy(z, x, y)
            p2 = plant.prob_z_given_xy_fast(z, x, y)
            assert p1 == p2

    def test_dom_weak_one_allele(self):
        a = Plant(4, 13, 13)
        b = Plant(4, 13, 15)
        c = Plant(4, 15, 13)
        d = Plant(4, 15, 15)
        assert b.dom_weak(b)
        assert c.dom_weak(c)
        assert not c.dom_weak(b)
        assert not b.dom_weak(c)
        assert b.dom_weak(a)
        assert d.dom_weak(a)
        assert d.dom_weak(b)
        assert d.dom_weak(c)

    def test_all_gametes(self):
        # create plant
        # create set of all gametes
        a = Plant(4, 5, 10)
        assert sorted(a.reachable_gametes()) == sorted(
            [
                5,
                4,
                6,
                2,
                10,
                11,
                9,
                13
            ]
        ), f"reachable_gametes = {a.reachable_gametes()}"

    def test_all_intermediates(self):
        pass

