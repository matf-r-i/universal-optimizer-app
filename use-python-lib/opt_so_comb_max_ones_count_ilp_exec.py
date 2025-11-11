from random import randint
from random import choice

import sys
import os
from pathlib import Path
directory = Path(__file__).resolve()
root_dir = directory.parent
sys.path.append(str(root_dir))
if 'LIB_SOURCE' in os.environ and os.environ['LIB_SOURCE']=='CODE':
        sys.path.append(str(root_dir/ "lib"))

from uo.algorithm.output_control import OutputControl

from opt.single_objective.comb.max_ones_count_problem.max_ones_count_problem import MaxOnesCountProblem
from opt.single_objective.comb.max_ones_count_problem.max_ones_count_problem_ilp_linopy import \
                MaxOnesCountProblemIntegerLinearProgrammingSolver

def main():
        problem_to_solve:MaxOnesCountProblem = MaxOnesCountProblem.from_dimension(dimension=10)
        solver:MaxOnesCountProblemIntegerLinearProgrammingSolver = MaxOnesCountProblemIntegerLinearProgrammingSolver(
                        problem=problem_to_solve)
        bs = solver.optimize()
        print('Best solution: {}'.format(bs.string_representation()))
        #print('Best solution code: {}'.format(solver.model.solution.x))            

if __name__ == '__main__':
        main()
