import tempfile
import unittest
from pathlib import Path

from opt.single_objective.comb.knapsack_problem.knapsack_problem import KnapsackProblem
from opt.single_objective.comb.knapsack_problem.solver import (
    build_parser,
    solve_ga,
    solve_vns,
)


class TestKnapsackProblemSolver(unittest.TestCase):

    def test_build_parser_accepts_ga_method(self):
        parser = build_parser()

        args = parser.parse_args([
            "--input-file", "dummy.txt",
            "--method", "ga",
        ])

        self.assertEqual(args.input_file, "dummy.txt")
        self.assertEqual(args.method, "ga")

    def test_build_parser_accepts_vns_method(self):
        parser = build_parser()

        args = parser.parse_args([
            "--input-file", "dummy.txt",
            "--method", "vns",
        ])

        self.assertEqual(args.input_file, "dummy.txt")
        self.assertEqual(args.method, "vns")

    def test_build_parser_rejects_invalid_method(self):
        parser = build_parser()

        with self.assertRaises(SystemExit):
            parser.parse_args([
                "--input-file", "dummy.txt",
                "--method", "invalid",
            ])

    def test_solve_ga_returns_optimizer_and_solution(self):
        problem = KnapsackProblem(
            capacity=15,
            weights=[2, 3, 5, 7, 1, 4],
            values=[10, 5, 15, 7, 6, 18],
        )

        optimizer, best_solution = solve_ga(
            problem=problem,
            random_seed=123,
            evaluations_max=300,
            population_size=20,
            elite_count=2,
            crossover_probability=0.9,
            mutation_probability=0.05,
        )

        self.assertIsNotNone(optimizer)
        self.assertIsNotNone(best_solution)
        self.assertIsNotNone(best_solution.representation)
        self.assertEqual(len(best_solution.representation), problem.dimension)

    def test_solve_vns_returns_optimizer_and_solution(self):
        problem = KnapsackProblem(
            capacity=15,
            weights=[2, 3, 5, 7, 1, 4],
            values=[10, 5, 15, 7, 6, 18],
        )

        optimizer, best_solution = solve_vns(
            problem=problem,
            random_seed=123,
            evaluations_max=300,
            k_min=1,
            k_max=3,
        )

        self.assertIsNotNone(optimizer)
        self.assertIsNotNone(best_solution)
        self.assertIsNotNone(best_solution.representation)
        self.assertEqual(len(best_solution.representation), problem.dimension)

    def test_problem_can_be_loaded_from_temp_file(self):
        content = "\n".join([
            "15",
            "2 10",
            "3 5",
            "5 15",
            "7 7",
            "1 6",
            "4 18",
        ])

        with tempfile.TemporaryDirectory() as tmp_dir:
            input_path = Path(tmp_dir) / "knapsack.txt"
            input_path.write_text(content, encoding="utf-8")

            problem = KnapsackProblem.from_input_file(str(input_path))

        self.assertEqual(problem.capacity, 15)
        self.assertEqual(problem.dimension, 6)


if __name__ == "__main__":
    unittest.main()