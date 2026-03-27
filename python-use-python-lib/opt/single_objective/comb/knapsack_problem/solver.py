"""
.. _py_knapsack_problem_solver:
"""

import argparse
from random import seed

from uo.algorithm.metaheuristic.finish_control import FinishControl

from uo.algorithm.metaheuristic.genetic_algorithm.ga_selection_roulette import (
    GaSelectionRoulette,
)
from uo.algorithm.metaheuristic.genetic_algorithm.ga_crossover_support_one_point_bit_array import (
    GaCrossoverSupportOnePointBitArray,
)
from uo.algorithm.metaheuristic.genetic_algorithm.ga_mutation_support_one_point_bit_array import (
    GaMutationSupportOnePointBitArray,
)
from uo.algorithm.metaheuristic.genetic_algorithm.ga_optimizer_gen import (
    GaOptimizerGenerationalConstructionParameters,
    GaOptimizerGenerational,
)

from uo.algorithm.metaheuristic.variable_neighborhood_search.vns_shaking_support_standard_bit_array import (
    VnsShakingSupportStandardBitArray,
)
from uo.algorithm.metaheuristic.variable_neighborhood_search.vns_ls_support_standard_bi_bit_array import (
    VnsLocalSearchSupportStandardBestImprovementBitArray,
)
from uo.algorithm.metaheuristic.variable_neighborhood_search.vns_optimizer import (
    VnsOptimizerConstructionParameters,
    VnsOptimizer,
)

from opt.single_objective.comb.knapsack_problem.knapsack_problem import KnapsackProblem
from opt.single_objective.comb.knapsack_problem.knapsack_problem_bit_array_solution import (
    KnapsackProblemBitArraySolution,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Solve 0/1 Knapsack problem.")

    parser.add_argument(
        "--input-file",
        type=str,
        required=True,
        help="Path to input file describing the problem instance.",
    )
    parser.add_argument(
        "--method",
        type=str,
        required=True,
        choices=["ga", "vns"],
        help="Optimization method to use.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=43434343,
        help="Random seed.",
    )
    parser.add_argument(
        "--evaluations-max",
        type=int,
        default=5000,
        help="Maximum number of evaluations.",
    )

    # GA params
    parser.add_argument(
        "--population-size",
        type=int,
        default=100,
        help="Population size for GA.",
    )
    parser.add_argument(
        "--elite-count",
        type=int,
        default=10,
        help="Elite count for GA.",
    )
    parser.add_argument(
        "--crossover-probability",
        type=float,
        default=0.95,
        help="Crossover probability for GA.",
    )
    parser.add_argument(
        "--mutation-probability",
        type=float,
        default=0.01,
        help="Mutation probability for GA.",
    )

    # VNS params
    parser.add_argument(
        "--k-min",
        type=int,
        default=1,
        help="Minimal neighborhood size for VNS.",
    )
    parser.add_argument(
        "--k-max",
        type=int,
        default=3,
        help="Maximal neighborhood size for VNS.",
    )

    return parser


def solve_ga(
    problem: KnapsackProblem,
    random_seed: int,
    evaluations_max: int,
    population_size: int,
    elite_count: int,
    crossover_probability: float,
    mutation_probability: float,
) -> tuple[GaOptimizerGenerational, KnapsackProblemBitArraySolution]:
    solution = KnapsackProblemBitArraySolution()

    finish = FinishControl(
        criteria="evaluations",
        evaluations_max=evaluations_max,
    )

    ga_selection = GaSelectionRoulette()
    ga_crossover_support = GaCrossoverSupportOnePointBitArray[str](
        crossover_probability=crossover_probability
    )
    ga_mutation_support = GaMutationSupportOnePointBitArray[str](
        mutation_probability=mutation_probability
    )

    params = GaOptimizerGenerationalConstructionParameters()
    params.problem = problem
    params.solution_template = solution
    params.finish_control = finish
    params.ga_selection = ga_selection
    params.ga_crossover_support = ga_crossover_support
    params.ga_mutation_support = ga_mutation_support
    params.random_seed = random_seed
    params.population_size = population_size
    params.elite_count = elite_count

    seed(random_seed)

    optimizer = GaOptimizerGenerational.from_construction_tuple(params)
    best_solution = optimizer.optimize()
    return optimizer, best_solution


def solve_vns(
    problem: KnapsackProblem,
    random_seed: int,
    evaluations_max: int,
    k_min: int,
    k_max: int,
) -> tuple[VnsOptimizer, KnapsackProblemBitArraySolution]:
    solution = KnapsackProblemBitArraySolution()

    finish = FinishControl(
        criteria="evaluations",
        evaluations_max=evaluations_max,
    )

    dimension = problem.dimension
    vns_shaking_support = VnsShakingSupportStandardBitArray[str](
        dimension=dimension
    )
    vns_ls_support = VnsLocalSearchSupportStandardBestImprovementBitArray[str](
        dimension=dimension
    )

    params = VnsOptimizerConstructionParameters()
    params.problem = problem
    params.solution_template = solution
    params.finish_control = finish
    params.random_seed = random_seed
    params.vns_shaking_support = vns_shaking_support
    params.vns_ls_support = vns_ls_support
    params.k_min = k_min
    params.k_max = k_max

    seed(random_seed)

    optimizer = VnsOptimizer.from_construction_tuple(params)
    best_solution = optimizer.optimize()
    return optimizer, best_solution


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    problem = KnapsackProblem.from_input_file(args.input_file)

    if args.method == "ga":
        optimizer, best_solution = solve_ga(
            problem=problem,
            random_seed=args.seed,
            evaluations_max=args.evaluations_max,
            population_size=args.population_size,
            elite_count=args.elite_count,
            crossover_probability=args.crossover_probability,
            mutation_probability=args.mutation_probability,
        )
    else:
        optimizer, best_solution = solve_vns(
            problem=problem,
            random_seed=args.seed,
            evaluations_max=args.evaluations_max,
            k_min=args.k_min,
            k_max=args.k_max,
        )

    print("Best solution representation: {}".format(best_solution.representation.bin))
    print("Best solution code: {}".format(best_solution.string_representation()))
    print("Best solution objective: {}".format(best_solution.objective_value))
    print("Best solution fitness: {}".format(best_solution.fitness_value))
    print("Best solution feasible: {}".format(best_solution.is_feasible))
    print("Number of iterations: {}".format(optimizer.iteration))
    print("Number of evaluations: {}".format(optimizer.evaluation))


if __name__ == "__main__":
    main()