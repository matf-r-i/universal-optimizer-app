import unittest
from opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem_binary_bit_array_solution import OnesCountMaxProblemBinaryBitArraySolution
from opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem import OnesCountMaxProblem
from uo.algorithm.metaheuristic.variable_neighborhood_search.vns_ls_support_standard_bi_bit_array import VnsLocalSearchSupportStandardBestImprovementBitArray
from uo.algorithm.metaheuristic.variable_neighborhood_search.vns_shaking_support_standard_bit_array import VnsShakingSupportStandardBitArray 
from uo.algorithm.metaheuristic.variable_neighborhood_search.vns_optimizer import VnsOptimizer
from uo.algorithm.metaheuristic.variable_neighborhood_search.vns_optimizer import VnsOptimizerConstructionParameters

from uo.algorithm.metaheuristic.additional_statistics_control import AdditionalStatisticsControl
from uo.algorithm.metaheuristic.finish_control import FinishControl
from uo.algorithm.output_control import OutputControl
from bitstring import Bits, BitArray, BitStream, pack
from random import choice
from random import randint
from copy import deepcopy
import unittest.mock as mocker
import sys
import os
from pathlib import Path
directory = Path(__file__).resolve()
root_dir = directory.parent.parent.parent.parent.parent.parent.parent
sys.path.append(str(root_dir))

directory = Path(__file__).resolve()
sys.path.append(directory.parent)


class TestOnesCountMaxProblemVnsBinaryBitArraySolutionLsbi(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print(
            "setUpClass TestIntegrationOnesCountMaxProblemVnsBinaryBitArraySolutionLsbi\n")

    def setUp(self):
        self.output_control = OutputControl()
        dimension: int = 24
        self.problem_to_solve: OnesCountMaxProblem = OnesCountMaxProblem.from_dimension(
            dimension=dimension)
        self.solution: OnesCountMaxProblemBinaryBitArraySolution = OnesCountMaxProblemBinaryBitArraySolution(
            random_seed=43434343)
        self.finish_control: FinishControl = FinishControl(
            criteria='evaluations', evaluations_max=1000)
        self.vns_ls_support: VnsLocalSearchSupportStandardBestImprovementBitArray[int] = VnsLocalSearchSupportStandardBestImprovementBitArray[int](
            dimension=dimension)
        self.vns_shaking_support: VnsShakingSupportStandardBitArray[int] = VnsShakingSupportStandardBitArray[int](dimension=dimension)
        self.additional_stat = AdditionalStatisticsControl(keep='')
        vns_construction_params: VnsOptimizerConstructionParameters = VnsOptimizerConstructionParameters()
        vns_construction_params.output_control = self.output_control
        vns_construction_params.problem = self.problem_to_solve
        vns_construction_params.solution_template = self.solution
        vns_construction_params.vns_shaking_support = self.vns_shaking_support
        vns_construction_params.vns_ls_support = self.vns_ls_support
        vns_construction_params.finish_control = self.finish_control
        vns_construction_params.random_seed = 43434343
        vns_construction_params.additional_statistics_control = self.additional_stat
        vns_construction_params.k_min = 1
        vns_construction_params.k_max = 3
        vns_construction_params.local_search_type = 'localSearchBestImprovement'
        self.optimizer: VnsOptimizer = VnsOptimizer.from_construction_tuple(
            vns_construction_params)
        self.bs = self.optimizer.optimize()

    def test_best_solution_after_optimization_should_be_all_optimal(self):
        result: str = '111111111111111111111111'
        self.assertEqual(
            self.optimizer.best_solution.string_representation(), result)

    def test_best_solution_after_optimization_should_be_optimal_2(self):
        self.assertEqual(len(self.optimizer.best_solution.string_representation(
        )), self.problem_to_solve.dimension)

    def test_best_solution_after_optimization_should_have_optimal_fitness(self):
        self.assertEqual(self.optimizer.best_solution.fitness_value,
                         self.problem_to_solve.dimension)

    def test_best_solution_after_optimization_should_have_optimal_objective_value(self):
        self.assertEqual(
            self.optimizer.best_solution.objective_value, self.problem_to_solve.dimension)

    def tearDown(self):
        return

    @classmethod
    def tearDownClass(cls):
        print("\ntearDownClass TestIntegrationOnesCountMaxProblemVnsBinaryBitArraySolutionLsbi")


if __name__ == '__main__':
    unittest.main()
