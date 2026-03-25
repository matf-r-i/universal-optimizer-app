import unittest

from bitstring import BitArray

from opt.single_objective.comb.knapsack_problem.knapsack_problem import KnapsackProblem
from opt.single_objective.comb.knapsack_problem.knapsack_problem_bit_array_solution import (
    KnapsackProblemBitArraySolution,
)


class TestKnapsackProblemBitArraySolution(unittest.TestCase):

    def test_initialize_instance_with_default_parameters(self):
        solution = KnapsackProblemBitArraySolution()

        self.assertIsNone(solution.fitness_value)
        self.assertIsNone(solution.fitness_values)
        self.assertIsNone(solution.objective_value)
        self.assertIsNone(solution.objective_values)
        self.assertFalse(solution.is_feasible)

    def test_init_random_method_with_problem(self):
        problem = KnapsackProblem(10, [2, 3, 4, 5], [5, 6, 7, 8])
        solution = KnapsackProblemBitArraySolution()

        solution.init_random(problem)

        self.assertIsInstance(solution.representation, BitArray)
        self.assertEqual(len(solution.representation), problem.dimension)

    def test_init_random_raises_when_problem_has_no_weights(self):
        solution = KnapsackProblemBitArraySolution()

        class DummyProblem:
            weights = None
            values = [1, 2, 3]

        with self.assertRaises(ValueError):
            solution.init_random(DummyProblem())

    def test_init_random_raises_when_problem_has_no_values(self):
        solution = KnapsackProblemBitArraySolution()

        class DummyProblem:
            weights = [1, 2, 3]
            values = None

        with self.assertRaises(ValueError):
            solution.init_random(DummyProblem())

    def test_init_from_method_with_bitarray_and_problem(self):
        problem = KnapsackProblem(10, [2, 3, 4], [5, 6, 7])
        representation = BitArray(bin="101")
        solution = KnapsackProblemBitArraySolution()

        solution.init_from(representation, problem)

        self.assertEqual(solution.representation.bin, "101")

    def test_init_from_raises_for_invalid_type(self):
        problem = KnapsackProblem(10, [2, 3, 4], [5, 6, 7])
        solution = KnapsackProblemBitArraySolution()

        with self.assertRaises(TypeError):
            solution.init_from("101", problem)

    def test_init_from_raises_for_empty_representation(self):
        problem = KnapsackProblem(10, [2, 3, 4], [5, 6, 7])
        solution = KnapsackProblemBitArraySolution()

        with self.assertRaises(ValueError):
            solution.init_from(BitArray(bin=""), problem)

    def test_is_feasible_sol_true(self):
        solution = KnapsackProblemBitArraySolution()
        representation = BitArray(bin="1100")

        self.assertTrue(solution.is_feasible_sol(representation, 10, [2, 3, 4, 5]))

    def test_is_feasible_sol_false(self):
        solution = KnapsackProblemBitArraySolution()
        representation = BitArray(bin="111")

        self.assertFalse(solution.is_feasible_sol(representation, 4, [2, 3, 4]))

    def test_native_representation_method_with_string_representation(self):
        representation_str = "101010"
        solution = KnapsackProblemBitArraySolution()

        native_representation = solution.native_representation(representation_str)

        self.assertEqual(native_representation.bin, representation_str)

    def test_native_representation_raises_for_invalid_type(self):
        solution = KnapsackProblemBitArraySolution()

        with self.assertRaises(TypeError):
            solution.native_representation(101010)

    def test_native_representation_raises_for_invalid_characters(self):
        solution = KnapsackProblemBitArraySolution()

        with self.assertRaises(ValueError):
            solution.native_representation("101201")

    def test_representation_distance_directly_method_with_string_representations(self):
        solution = KnapsackProblemBitArraySolution()

        distance = solution.representation_distance_directly("101010", "111000")

        self.assertEqual(distance, 2)

    def test_representation_distance_directly_raises_for_invalid_first_type(self):
        solution = KnapsackProblemBitArraySolution()

        with self.assertRaises(TypeError):
            solution.representation_distance_directly(1010, "1010")

    def test_representation_distance_directly_raises_for_invalid_second_type(self):
        solution = KnapsackProblemBitArraySolution()

        with self.assertRaises(TypeError):
            solution.representation_distance_directly("1010", 1010)

    def test_representation_distance_directly_raises_value_error(self):
        solution = KnapsackProblemBitArraySolution()

        with self.assertRaises(ValueError):
            solution.representation_distance_directly("1010", "101")

    def test_calculate_quality_directly_method_with_feasible_solution(self):
        problem = KnapsackProblem(10, [2, 3, 4, 5], [5, 6, 7, 8])
        solution = KnapsackProblemBitArraySolution()
        representation = BitArray(bin="1100")

        quality = solution.calculate_quality_directly(representation, problem)

        self.assertEqual(quality.objective_value, 11.0)
        self.assertEqual(quality.fitness_value, 10.995)
        self.assertTrue(quality.is_feasible)

    def test_calculate_quality_directly_method_with_infeasible_solution(self):
        problem = KnapsackProblem(4, [2, 3, 4], [5, 6, 7])
        solution = KnapsackProblemBitArraySolution()
        representation = BitArray(bin="111")

        quality = solution.calculate_quality_directly(representation, problem)

        self.assertEqual(quality.objective_value, float("-inf"))
        self.assertEqual(quality.fitness_value, -5000009.0)
        self.assertFalse(quality.is_feasible)

    def test_calculate_quality_directly_raises_for_invalid_representation_type(self):
        problem = KnapsackProblem(10, [2, 3, 4], [5, 6, 7])
        solution = KnapsackProblemBitArraySolution()

        with self.assertRaises(TypeError):
            solution.calculate_quality_directly("101", problem)

    def test_calculate_quality_directly_raises_for_invalid_problem_type(self):
        solution = KnapsackProblemBitArraySolution()
        representation = BitArray(bin="101")

        with self.assertRaises(TypeError):
            solution.calculate_quality_directly(representation, "not_a_problem")

    def test_calculate_quality_directly_raises_for_wrong_length(self):
        problem = KnapsackProblem(10, [2, 3, 4], [5, 6, 7])
        solution = KnapsackProblemBitArraySolution()
        representation = BitArray(bin="1010")

        with self.assertRaises(ValueError):
            solution.calculate_quality_directly(representation, problem)

    def test_copy_method_returns_deep_copy(self):
        problem = KnapsackProblem(10, [2, 3, 4], [5, 6, 7])
        solution = KnapsackProblemBitArraySolution()
        solution.init_from(BitArray(bin="101"), problem)

        copy_solution = solution.copy()

        self.assertIsNot(solution, copy_solution)
        self.assertEqual(solution.representation.bin, copy_solution.representation.bin)

    def test_argument_method_returns_correct_string_representation(self):
        solution = KnapsackProblemBitArraySolution()
        representation = BitArray(bin="101010")

        argument = solution.argument(representation)

        self.assertEqual(argument, "101010")


if __name__ == "__main__":
    unittest.main()