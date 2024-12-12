import unittest

import eugene.utils as eut
import eugene.db as edb
import eugene.solvers.base_min_crossings_minizinc as emz
import eugene.solvers.base_min_crossings_distribute as eds


class TestDistributeDominance(unittest.TestCase):

    """Tests to ensure the dominance decomposition works"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_distribute_dominance(self):
        with open("./distribute_data_2.json") as f:
            db = edb.DistributeDB.from_json_file(f)
        solver_mzn = edb.DBSolver(
            lambda case: emz.breeding_program_distribute(
                len(case), case, emz.DEFAULT_CTX
            ),
            db,
        )
        solver_dom = lambda case: eds.breeding_program_distribute(
            len(case), case, lambda _, case: solver_mzn(case)
        )
        for n_loci in range(2, 9):
            for case in eut.gen_distribute_instances(n_loci):
                self.assertEqual(solver_dom(case).objective, solver_mzn(case).objective)
