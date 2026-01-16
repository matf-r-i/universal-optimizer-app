

import sys
import os
from pathlib import Path
directory = Path(__file__).resolve()
sys.path.append(directory.parent)
sys.path.append(directory.parent.parent)
sys.path.append(directory.parent.parent.parent)
sys.path.append(directory.parent.parent.parent.parent)
root_dir = directory.parent.parent.parent.parent.parent
sys.path.append(str(root_dir))

from typing import Optional
from copy import deepcopy
from random import choice
from random import randint

from uo.utils.logger import logger
from uo.utils.complex_counter_uniform_ascending import ComplexCounterUniformAscending

from uo.solution.quality_of_solution import QualityOfSolution
from uo.algorithm.algorithm import Algorithm
from uo.algorithm.metaheuristic.variable_neighborhood_search.vns_shaking_support import \
        VnsShakingSupport

from opt.single_objective.glob.function_one_variable_max_problem.function_one_variable_max_problem import \
        FunctionOneVariableMaxProblemMax
from opt.single_objective.glob.function_one_variable_max_problem.function_one_variable_max_problem_binary_int_solution \
        import FunctionOneVariableMaxProblemBinaryIntSolution

class FunctionOneVariableMaxProblemBinaryIntSolutionVnsShakingSupport(VnsShakingSupport[int,float]):
    
    def __init__(self)->None:
        return

    def __copy__(self):
        sup = deepcopy(self)
        return sup

    def copy(self):
        return self.__copy__()
        
    def shaking(self, k:int, problem:FunctionOneVariableMaxProblemMax, solution:FunctionOneVariableMaxProblemBinaryIntSolution, 
            optimizer:Algorithm)->bool:
        if k <= 0:
            return False
        if optimizer.finish_control.is_finished(optimizer.evaluation, optimizer.iteration, optimizer.elapsed_seconds()):
            return False
        tries:int = 0
        limit:int = 10000
        representation_length:int = 32
        while tries < limit:
            positions:list[int] = []
            for _ in range(0,k):
                positions.append(choice(range(representation_length)))
            mask:int = 0
            for p in positions:
                mask |= 1 << p
            solution.representation ^= mask
            all_ok:bool = True
            if solution.representation.bit_count() > representation_length:
                all_ok = False
            if all_ok:
                break
        if tries < limit:
            solution.representation = solution.obtain_feasible_representation(problem)
            if optimizer.finish_control.is_finished(optimizer.evaluation, optimizer.iteration, optimizer.elapsed_seconds()):
                return False
            optimizer.write_output_values_if_needed("before_evaluation", "b_e")
            optimizer.evaluation += 1
            solution.evaluate(problem)
            optimizer.write_output_values_if_needed("after_evaluation", "a_e")
            return True
        else:
            return False 

 
    def string_rep(self, delimiter:str, indentation:int=0, indentation_symbol:str='', group_start:str ='{', 
        group_end:str ='}')->str:
        return 'FunctionOneVariableMaxProblemBinaryIntSolutionVnsShakingSupport'

    def __str__(self)->str:
        return self.string_rep('|')

    def __repr__(self)->str:
        return self.string_rep('\n')


    def __format__(self, spec:str)->str:
        return self.string_rep('|')



