from random import seed

from uo.algorithm.metaheuristic.finish_control import FinishControl
from uo.algorithm.metaheuristic.genetic_algorithm.ga_selection_roulette import GaSelectionRoulette
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

from opt.single_objective.comb.knapsack_problem.knapsack_problem import KnapsackProblem
from opt.single_objective.comb.knapsack_problem.knapsack_problem_bit_array_solution import (
    KnapsackProblemBitArraySolution,
)


def main():
    problem_to_solve = KnapsackProblem.from_capacity_weights_and_values(
        capacity=15,
        weights=[2, 3, 5, 7, 1, 4],
        values=[10, 5, 15, 7, 6, 18],
    )

    solution = KnapsackProblemBitArraySolution()

    finish = FinishControl(criteria="evaluations", evaluations_max=5000)

    ga_selection = GaSelectionRoulette()
    ga_crossover_support = GaCrossoverSupportOnePointBitArray[str](
        crossover_probability=0.95
    )
    ga_mutation_support = GaMutationSupportOnePointBitArray[str](
        mutation_probability=0.01
    )

    ga_construction_params = GaOptimizerGenerationalConstructionParameters()
    ga_construction_params.problem = problem_to_solve
    ga_construction_params.solution_template = solution
    ga_construction_params.finish_control = finish
    ga_construction_params.ga_selection = ga_selection
    ga_construction_params.ga_crossover_support = ga_crossover_support
    ga_construction_params.ga_mutation_support = ga_mutation_support
    ga_construction_params.random_seed = 43434343
    ga_construction_params.population_size = 100
    ga_construction_params.elite_count = 10

    seed(ga_construction_params.random_seed)

    optimizer = GaOptimizerGenerational.from_construction_tuple(
        ga_construction_params
    )
    best_solution = optimizer.optimize()

    print("Best solution representation: {}".format(best_solution.representation.bin))
    print("Best solution code: {}".format(best_solution.string_representation()))
    print("Best solution objective: {}".format(best_solution.objective_value))
    print("Best solution fitness: {}".format(best_solution.fitness_value))
    print("Best solution feasible: {}".format(best_solution.is_feasible))
    print("Number of iterations: {}".format(optimizer.iteration))
    print("Number of evaluations: {}".format(optimizer.evaluation))


if __name__ == "__main__":
    main()