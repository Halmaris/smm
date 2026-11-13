"""Check TPE and pruning using artificial scores, without model fitting."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'python'))
import tpe
import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)


class SearchTests(unittest.TestCase):
    def test_anchors(self):
        for method, count in zip(('M00', 'M10', 'M01', 'M11'), (4, 12, 12, 36)):
            rows = tpe.anchors(method)
            self.assertEqual(len(rows), count)
            search = tpe.Search(method, 10, count)
            for row in rows:
                trial = search.ask()
                self.assertEqual(trial['parameters'], row)
                self.assertTrue(trial['anchor'])
                self.assertFalse(search.report(.1, .1, 1))
                search.tell(.7)
            self.assertIsNone(search.ask())
            self.assertEqual(search.best()['parameters'], rows[0])

    def test_bound_and_final_fold(self):
        search = tpe.Search('M00', 2, 6)
        for _ in range(4):
            search.ask()
            search.tell(.8)
        search.ask()
        self.assertFalse(search.report(.7, .8, 1))
        self.assertTrue(search.report(.5, .7, 2))
        search.tell(pruned=True)
        search.ask()
        self.assertFalse(search.report(.5, .5, 5))
        search.tell(.5)

    def test_reproducibility(self):
        searches = [tpe.Search('M11', 91001, 40) for _ in range(2)]
        for j in range(40):
            left, right = [search.ask() for search in searches]
            self.assertEqual(left, right)
            for search in searches:
                search.tell(.5 + j / 1000)


if __name__ == '__main__':
    unittest.main()
