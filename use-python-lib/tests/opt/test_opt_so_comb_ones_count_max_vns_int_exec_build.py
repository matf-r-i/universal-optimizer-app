import unittest
import sys
import os
import sys
from pathlib import Path
directory = Path(__file__).resolve()
root_dir = directory.parent.parent.parent.parent
sys.path.append(str(root_dir))

from opt_so_comb_ones_count_max_vns_int_exec_build import (
    OnesCountMaxProblem2,
    OnesCountMaxProblemBinaryIntSolution,
    OnesCountMaxProblemBinaryIntSolutionVnsSupport
)

class TestOnesCountMaxProblem2(unittest.TestCase):
    def test_dimension_valid(self):
        p = OnesCountMaxProblem2(5)
        self.assertEqual(p.dimension, 5)

    def test_dimension_invalid_type(self):
        with self.assertRaises(TypeError):
            OnesCountMaxProblem2('a')

    def test_dimension_invalid_value(self):
        with self.assertRaises(ValueError):
            OnesCountMaxProblem2(0)
        with self.assertRaises(ValueError):
            OnesCountMaxProblem2(32)

class TestOnesCountMaxProblemBinaryIntSolution(unittest.TestCase):
    def setUp(self):
        self.problem = OnesCountMaxProblem2(5)
        self.solution = OnesCountMaxProblemBinaryIntSolution()

    def test_init_random(self):
        self.solution.init_random(self.problem)
        self.assertTrue(0 <= self.solution.representation < 2**self.problem.dimension)

    def test_init_from(self):
        self.solution.init_from(7, self.problem)
        self.assertEqual(self.solution.representation, 7)

    def test_argument(self):
        self.solution.init_from(7, self.problem)
        self.assertEqual(self.solution.argument(self.problem), bin(7))

    def test_calculate_quality_directly(self):
        q = self.solution.calculate_quality_directly(7, self.problem)
        self.assertEqual(q.fitness_value, bin(7).count('1'))

    def test_native_representation(self):
        self.assertEqual(self.solution.native_representation('101'), 5)

    def test_representation_distance_directly(self):
        d = self.solution.representation_distance_directly('101', '111')
        self.assertEqual(d, 1)

class TestOnesCountMaxProblemBinaryIntSolutionVnsSupport(unittest.TestCase):
    def setUp(self):
        self.problem = OnesCountMaxProblem2(5)
        self.solution = OnesCountMaxProblemBinaryIntSolution()
        self.solution.init_random(self.problem)
        self.optimizer = type('DummyOptimizer', (), {
            'finish_control': type('DummyFinish', (), {'is_finished': lambda *a, **kw: False})(),
            'evaluation': 0,
            'iteration': 0,
            'elapsed_seconds': lambda self=None: 0,
            'is_first_better': lambda self, s1, s2, p: True
        })()
        self.vns_support = OnesCountMaxProblemBinaryIntSolutionVnsSupport()

    def test_shaking(self):
        result = self.vns_support.shaking(1, self.problem, self.solution, self.optimizer)
        self.assertTrue(result)

    def test_local_search_best_improvement(self):
        result = self.vns_support.local_search_best_improvement(1, self.problem, self.solution, self.optimizer)
        self.assertTrue(result)

    def test_local_search_first_improvement(self):
        result = self.vns_support.local_search_first_improvement(1, self.problem, self.solution, self.optimizer)
        self.assertTrue(result)

if __name__ == '__main__':
    unittest.main()
