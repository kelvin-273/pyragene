import unittest
import eugene.solvers.base_min_crossings_distribute_astar as eda


class TestDistAstar(unittest.TestCase):

    """Test case docstring."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_branching_unsimplified(self):
        def f(case):
            return sorted(
                [
                    sol.xs
                    for sol in eda.branching(
                        eda.Node.from_dist_array(case),
                        eda.BRANCHING_CTX,
                        simplify_results=False,
                    )
                ]
            )

        self.assertEqual(f([0, 1]), [[0, 0]])
        self.assertEqual(f([0, 1, 0]), [[0, 0, 1], [0, 1, 1]])
        self.assertEqual(f([0, 1, 2]), [[0, 0, 1], [0, 1, 0], [0, 1, 1]])
        self.assertEqual(f([0, 1, 0, 1]), [[0, 1], [0, 1, 0]])
        self.assertEqual(f([0, 1, 0, 2]), [[0, 1, 0], [0, 1, 0, 1], [0, 1, 2]])
        self.assertEqual(f([0, 1, 2, 0]), [[0, 1, 0], [0, 1, 2]])
        self.assertEqual(f([0, 1, 2, 1]), [[0, 1, 0], [0, 1, 2]])
        self.assertEqual(f([0, 1, 2, 3]), [[1, 2, 3]])
        self.assertEqual(f([0, 1, 0, 1, 0]), [])
        self.assertEqual(f([0, 1, 0, 1, 2]), [])
        self.assertEqual(f([0, 1, 0, 2, 0]), [])
        self.assertEqual(f([0, 1, 0, 2, 1]), [])
        self.assertEqual(f([0, 1, 0, 2, 3]), [])
        self.assertEqual(f([0, 1, 2, 0, 1]), [])
        self.assertEqual(f([0, 1, 2, 0, 2]), [])
        self.assertEqual(f([0, 1, 2, 0, 3]), [])
        self.assertEqual(f([0, 1, 2, 1, 0]), [])
        self.assertEqual(f([0, 1, 2, 1, 2]), [])
        self.assertEqual(f([0, 1, 2, 1, 3]), [])
        self.assertEqual(f([0, 1, 2, 3, 0]), [])
        self.assertEqual(f([0, 1, 2, 3, 1]), [])
        self.assertEqual(f([0, 1, 2, 3, 2]), [])
        self.assertEqual(f([0, 1, 2, 3, 4]), [])

    def test_branching_simplified(self):
        def f(case):
            return sorted(
                [
                    sol.xs
                    for sol in eda.branching(
                        eda.Node.from_dist_array(case), eda.BRANCHING_CTX
                    )
                ]
            )

        self.assertEqual(f([0, 1]), [[0]])
        self.assertEqual(f([0, 1, 0]), [[0, 1]])
        self.assertEqual(f([0, 1, 2]), [[0, 1]])
        self.assertEqual(f([0, 1, 0, 1]), [[0, 1], [0, 1, 0]])
        self.assertEqual(f([0, 1, 0, 2]), [[0, 1, 0], [0, 1, 0, 1], [0, 1, 2]])
        self.assertEqual(f([0, 1, 2, 0]), [[0, 1, 0], [0, 1, 2]])
        self.assertEqual(f([0, 1, 2, 1]), [[0, 1, 0], [0, 1, 2]])
        self.assertEqual(f([0, 1, 2, 3]), [[1, 2, 3]])
        self.assertEqual(f([0, 1, 0, 1, 0]), [])
        self.assertEqual(f([0, 1, 0, 1, 2]), [])
        self.assertEqual(f([0, 1, 0, 2, 0]), [])
        self.assertEqual(f([0, 1, 0, 2, 1]), [])
        self.assertEqual(f([0, 1, 0, 2, 3]), [])
        self.assertEqual(f([0, 1, 2, 0, 1]), [])
        self.assertEqual(f([0, 1, 2, 0, 2]), [])
        self.assertEqual(f([0, 1, 2, 0, 3]), [])
        self.assertEqual(f([0, 1, 2, 1, 0]), [])
        self.assertEqual(f([0, 1, 2, 1, 2]), [])
        self.assertEqual(f([0, 1, 2, 1, 3]), [])
        self.assertEqual(f([0, 1, 2, 3, 0]), [])
        self.assertEqual(f([0, 1, 2, 3, 1]), [])
        self.assertEqual(f([0, 1, 2, 3, 2]), [])
        self.assertEqual(f([0, 1, 2, 3, 4]), [])
