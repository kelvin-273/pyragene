import unittest
import eugene.solvers.base_min_crossings_wedges as ew


class TestClass(unittest.TestCase):

    """Test class for max_repeated_wedges."""

    def test_max_repeated_wedges(self):
        f = lambda s: [int(c) for c in s]
        self.assertEqual(ew.max_repeated_wedges(f("0")), (0, []))
        self.assertEqual(ew.max_repeated_wedges(f("01")), (0, [True]))
        self.assertEqual(ew.max_repeated_wedges(f("010")), (0, [True, False]))
        self.assertEqual(ew.max_repeated_wedges(f("0101")), (1, [True, False, True]))
        self.assertEqual(ew.max_repeated_wedges(f("01201")), (1, [True, False, False, True]))
        self.assertEqual(ew.max_repeated_wedges(f("01210")), (1, [True, False, False, True]))
        self.assertEqual(ew.max_repeated_wedges(f("01212")), (1, [False, True, False, True]))
        res_n, res_wedges = ew.max_repeated_wedges(f("01234563210"))
        self.assertEqual(res_n, 2)
        self.assertEqual(len(res_wedges), 10)
        self.assertEqual(res_wedges[-4:], [False, True, False, True])
        self.assertEqual(res_wedges[:4], [True, False, True, False])
        self.assertEqual(
            ew.max_repeated_wedges(f("010101010120")),
            (4, [True, False, True, False, True, False, True, False, True, False, True])
        )

        [0, 1, 2, 1, 3, 4, 3, 5, 4, 6, 3, 4]
        [0, 1, 0, 2, 0, 3, 0, 4, 0, 2]
