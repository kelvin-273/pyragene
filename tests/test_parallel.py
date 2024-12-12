import unittest
import math

from eugene.plant_models.plant2 import PlantSPC
from eugene.solvers.enumerators import BreedingProgram


class TestClass(unittest.TestCase):

    """Test case docstring."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_trivial(self):
        """Tests that populations with the goal terminate in 0 generations"""
        def aux(n_loci):
            runner = BreedingProgram(n_loci, PlantSPC)
            runner.set_init_pop(
                PlantSPC.initial_pop_trait_introgression(n_loci, n_holes=0)
            )
            runner.set_ideotype(PlantSPC(
                n_loci,
                (1 << n_loci) - 1,
                (1 << n_loci) - 1),
            )
            runner.run()
            res = runner.get_results()
            self.assertEqual(
                res.n_generations, 0,
                msg="Didn't find the goal in the initial population"
            )
        for n_loci in range(1, 10):
            aux(n_loci)

    def test_worstcase_homo(self):
        """Tests that populations with the goal terminate in 0 generations"""
        def aux(n_loci):
            runner = BreedingProgram(n_loci, PlantSPC)
            runner.set_init_pop(
                PlantSPC.initial_pop_singles_homo(n_loci)
            )
            runner.set_ideotype(PlantSPC(
                n_loci,
                (1 << n_loci) - 1,
                (1 << n_loci) - 1),
            )
            runner.run()
            res = runner.get_results()
            self.assertEqual(
                res.n_generations, [0, 2, 3, 3, 4, 4, 4, 4, 5][n_loci-1],
                msg="Incorrect worst-case for homozygous singles"
            )
        for n_loci in range(1, 8):
            aux(n_loci)
