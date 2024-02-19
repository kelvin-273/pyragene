import unittest
import eugene.db as ed


class TestDB(unittest.TestCase):

    """testing the functions of the DistributeDB"""

    def test_DistributeDB(self):
        db = ed.DistributeDB()
        db[[0, 1, 0]] = 3
        self.assertTrue(db[[0, 1, 0]] == 3)
        self.assertNotIn([0], db)

    def test_caching(self):
        import os
        test_file = "tests/test_distribute_db.json"

        if os.path.exists(test_file):
            os.remove(test_file)

        with ed.CachedDB(
            test_file,
            ed.DistributeDB,
        ) as db:
            db[[0, 1, 2]] = 3
            self.assertEqual(db[[0, 1, 2]], 3)
            self.assertNotIn([0], db)

        with ed.CachedDB(
            test_file,
            ed.DistributeDB,
        ) as db:
            db[[0, 1, 0, 1]] = 3
            db[[0, 1, 0, 2]] = 4

            self.assertIn([0, 1, 2], db)
            self.assertEqual(db[[0, 1, 2]], 3)
            self.assertEqual(db[[0, 1, 0, 1]], 3)
            self.assertEqual(db[[0, 1, 0, 2]], 4)

        os.remove(test_file)
