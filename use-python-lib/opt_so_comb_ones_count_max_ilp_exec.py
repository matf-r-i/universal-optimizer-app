from random import randint
from random import choice

from uo.algorithm.output_control import OutputControl

from opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem import OnesCountMaxProblem
from opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem_ilp_solution import OnesCountMaxProblemIntegerLinearProgrammingSolution
from opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem_ilp_linopy import OnesCountMaxProblemIntegerLinearProgrammingSolver


def main():
    problem_to_solve: OnesCountMaxProblem = OnesCountMaxProblem.from_dimension(
        dimension=10)
    solver: OnesCountMaxProblemIntegerLinearProgrammingSolver = OnesCountMaxProblemIntegerLinearProgrammingSolver(
        problem=problem_to_solve)
    bs:OnesCountMaxProblemIntegerLinearProgrammingSolution = solver.optimize()
    print('Best solution:' + bs.string_representation())
    print(f'Best solution: {bs}')



if __name__ == '__main__':
    main()
