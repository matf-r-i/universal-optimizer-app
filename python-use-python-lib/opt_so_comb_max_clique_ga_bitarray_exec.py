import sys
import os
from pathlib import Path

from random import randint
from random import choice
from random import randint
import networkx as nx

from uo.algorithm.output_control import OutputControl
from uo.algorithm.metaheuristic.finish_control import FinishControl

from uo.algorithm.metaheuristic.genetic_algorithm.ga_selection_roulette import GaSelectionRoulette
from uo.algorithm.metaheuristic.genetic_algorithm.ga_crossover_support_one_point_bit_array import GaCrossoverSupportOnePointBitArray
from uo.algorithm.metaheuristic.genetic_algorithm.ga_mutation_support_one_point_bit_array import GaMutationSupportOnePointBitArray
from uo.algorithm.metaheuristic.genetic_algorithm.ga_optimizer_gen import GaOptimizerGenerationalConstructionParameters
from uo.algorithm.metaheuristic.genetic_algorithm.ga_optimizer_gen import GaOptimizerGenerational

from opt.single_objective.comb.max_clique_problem.max_clique_problem import MaxCliqueProblem
from opt.single_objective.comb.max_clique_problem.max_clique_problem_binary_bit_array_solution import MaxCliqueProblemBitArraySolution

def main():
        input_path:str = 'opt/single_objective/comb/max_clique_problem/inputs/graph_01.txt'
        problem_to_solve:MaxCliqueProblem = MaxCliqueProblem.from_input_file(input_path)
        solution:MaxCliqueProblemBitArraySolution = MaxCliqueProblemBitArraySolution()
        finish:FinishControl = FinishControl(criteria='iterations', iterations_max=100)
        select:GaSelectionRoulette = GaSelectionRoulette()
        ga_cross_support:GaCrossoverSupportOnePointBitArray[str] = GaCrossoverSupportOnePointBitArray[str](crossover_probability=0.999)
        ga_mut_support:GaMutationSupportOnePointBitArray[str] = GaMutationSupportOnePointBitArray(mutation_probability=0.05)
        ga_construction_params:GaOptimizerGenerationalConstructionParameters = GaOptimizerGenerationalConstructionParameters()
        ga_construction_params.problem = problem_to_solve
        ga_construction_params.solution_template = solution
        ga_construction_params.finish_control = finish
        ga_construction_params.ga_selection = select
        ga_construction_params.ga_crossover_support = ga_cross_support
        ga_construction_params.ga_mutation_support = ga_mut_support
        ga_construction_params.random_seed = 43434343
        ga_construction_params.population_size = 30
        ga_construction_params.elite_count = 3
        optimizer:GaOptimizerGenerational = GaOptimizerGenerational.from_construction_tuple(ga_construction_params)
        optimizer.optimize()

        #print('Best solution representation: {}'.format(optimizer.best_solution.representation.bin))
        print('Best solution code: {}'.format(optimizer.best_solution.string_representation()))
        print('Best solution objective: {}'.format(optimizer.best_solution.objective_value))
        print('Best solution fitness: {}'.format(optimizer.best_solution.fitness_value))
        print('Number of iterations: {}'.format(optimizer.iteration))
        print('Number of evaluations: {}'.format(optimizer.evaluation))

if __name__ == '__main__':
        main()