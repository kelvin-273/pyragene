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

        """
        In the current mininzinc (2024-12-13) the length of the output array is
        always the same of the length of the input array. As a consequence,
        many of the solutions that the minizinc model returns are  identical.
        For example, given the distribute array `[1, 2, 3]`, the minizinc model
        will return `[[0, 0, 2], [0, 1, 0], [0, 1, 1]]`. However, `[0, 0, 2]`
        is structurally identical to `[0, 1, 1]` in that the any optimal
        solution to `[0, 0, 2]` exhibits the same structure to some opitmal
        solution to `[0, 1, 1]`.

        This necessitates a second layer of symmetry breaking,
        left as an exercise to a future computer scientist.

        Wait (2024-12-13), this symmetry breaking is meant to be performed when
        the `simplify_results` option is set to `True`.
        See `test_branching_simplified` below.
        """
        # TODO: Break symmetry on equivalent output arrays <13-12-24> #
        self.assertEqual(f([0, 1]), [[0, 0]])
        self.assertEqual(f([0, 1, 0]), [[0, 0, 1], [0, 1, 1]])
        self.assertEqual(f([0, 1, 2]), [[0, 0, 2], [0, 1, 0], [0, 1, 1]])
        self.assertEqual(f([0, 1, 0, 1]), [[0, 0, 1, 0], [0, 0,  1, 1], [0, 1, 0, 0], [0, 1, 1, 0]])
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
