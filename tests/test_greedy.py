import unittest
from eugene.simulators import greedy_time
from eugene.plant_models.plant2 import PlantSPC 


class TestGreedy(unittest.TestCase):

    """Test case docstring."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_chunks_from_gamete(self):
        g = int("1010110111", base=2)
        self.assertEqual(
            greedy_time.extract_chunks_from_gamete(10, g),
            [(0, 0, g), (2, 2, g), (4, 5, g), (7, 9, g)],
        )

    def test_chunks_from_genotype(self):
        self.assertEqual(
            list(map(lambda x: (x[0], x[1]), greedy_time.extract_chunks_from_genotype(
                10, PlantSPC(
                    10,
                    int('1001100111', base=2),
                    int('0111001010', base=2),
                )
            ))),
            [(0, 3), (1, 4), (6, 9)]
        )
