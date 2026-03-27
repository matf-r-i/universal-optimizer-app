import tempfile
import unittest
from pathlib import Path

from opt.single_objective.comb.knapsack_problem.knapsack_problem import KnapsackProblem


class TestKnapsackProblem(unittest.TestCase):

    def test_initialize_instance_with_valid_parameters(self):
        problem = KnapsackProblem(10, [2, 3, 4], [5, 6, 7])

        self.assertEqual(problem.capacity, 10)
        self.assertEqual(problem.weights, [2, 3, 4])
        self.assertEqual(problem.values, [5, 6, 7])
        self.assertEqual(problem.dimension, 3)
        self.assertFalse(problem.is_minimization)
        self.assertFalse(problem.is_multi_objective)

    def test_copy_returns_independent_copy(self):
        problem = KnapsackProblem(10, [2, 3, 4], [5, 6, 7])

        copied = problem.copy()

        self.assertIsNot(problem, copied)
        self.assertEqual(problem.capacity, copied.capacity)
        self.assertEqual(problem.weights, copied.weights)
        self.assertEqual(problem.values, copied.values)

        copied.weights[0] = 99
        self.assertEqual(problem.weights[0], 2)

    def test_from_capacity_weights_and_values(self):
        problem = KnapsackProblem.from_capacity_weights_and_values(
            capacity=15,
            weights=[1, 2, 3],
            values=[4, 5, 6],
        )

        self.assertEqual(problem.capacity, 15)
        self.assertEqual(problem.weights, [1, 2, 3])
        self.assertEqual(problem.values, [4, 5, 6])
        self.assertEqual(problem.dimension, 3)

    def test_capacity_type_error(self):
        with self.assertRaises(TypeError):
            KnapsackProblem("10", [2, 3, 4], [5, 6, 7])

    def test_weights_type_error(self):
        with self.assertRaises(TypeError):
            KnapsackProblem(10, "234", [5, 6, 7])

    def test_values_type_error(self):
        with self.assertRaises(TypeError):
            KnapsackProblem(10, [2, 3, 4], "567")

    def test_capacity_must_be_positive(self):
        with self.assertRaises(ValueError):
            KnapsackProblem(0, [2, 3, 4], [5, 6, 7])

    def test_weights_must_be_non_empty(self):
        with self.assertRaises(ValueError):
            KnapsackProblem(10, [], [])

    def test_weights_and_values_must_have_same_length(self):
        with self.assertRaises(ValueError):
            KnapsackProblem(10, [2, 3, 4], [5, 6])

    def test_weights_must_be_positive_integers(self):
        with self.assertRaises(ValueError):
            KnapsackProblem(10, [2, 0, 4], [5, 6, 7])

    def test_values_must_be_non_negative_integers(self):
        with self.assertRaises(ValueError):
            KnapsackProblem(10, [2, 3, 4], [5, -1, 7])

    def test_from_input_file_with_valid_file(self):
        content = "\n".join([
            "15",
            "2 10",
            "3 5",
            "5 15",
            "7 7",
        ])

        with tempfile.TemporaryDirectory() as tmp_dir:
            input_path = Path(tmp_dir) / "knapsack.txt"
            input_path.write_text(content, encoding="utf-8")

            problem = KnapsackProblem.from_input_file(str(input_path))

        self.assertEqual(problem.capacity, 15)
        self.assertEqual(problem.weights, [2, 3, 5, 7])
        self.assertEqual(problem.values, [10, 5, 15, 7])
        self.assertEqual(problem.dimension, 4)

    def test_from_input_file_raises_for_missing_file(self):
        with self.assertRaises(FileNotFoundError):
            KnapsackProblem.from_input_file("definitely_missing_file.txt")

    def test_from_input_file_raises_for_too_short_file(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_path = Path(tmp_dir) / "bad.txt"
            input_path.write_text("15\n", encoding="utf-8")

            with self.assertRaises(ValueError):
                KnapsackProblem.from_input_file(str(input_path))

    def test_from_input_file_raises_for_bad_item_line(self):
        content = "\n".join([
            "15",
            "2 10",
            "3",
            "5 15",
        ])

        with tempfile.TemporaryDirectory() as tmp_dir:
            input_path = Path(tmp_dir) / "bad.txt"
            input_path.write_text(content, encoding="utf-8")

            with self.assertRaises(ValueError):
                KnapsackProblem.from_input_file(str(input_path))


if __name__ == "__main__":
    unittest.main()