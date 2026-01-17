""" 
..  _py_ones_count_max_problem_bit_array_solution_vns_support:

The :mod:`~opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem_binary_bit_array_solution_vns_support` 
contains class :class:`~opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem_binary_bit_array_solution_vns_support.OnesCountMaxProblemBinaryBitArraySolutionVnsSupport`, 
that represents supporting parts of the `VNS` algorithm, where solution of the :ref:`Problem_Ones_Count_Max` have `BitArray` 
representation.
"""

import sys
import os
from pathlib import Path
directory = Path(__file__).resolve()
sys.path.append(directory)
sys.path.append(directory.parent)
sys.path.append(directory.parent.parent.parent)
sys.path.append(directory.parent.parent.parent.parent)
root_dir = directory.parent.parent.parent.parent.parent
sys.path.append(str(root_dir))

from copy import deepcopy
from random import choice
from random import random

from bitstring import Bits, BitArray, BitStream, pack


from uo.utils.logger import logger
from uo.utils.complex_counter_uniform_ascending import ComplexCounterUniformAscending

from uo.solution.quality_of_solution import QualityOfSolution
from uo.algorithm.algorithm import Algorithm
from uo.algorithm.metaheuristic.variable_neighborhood_search.vns_shaking_support import VnsShakingSupport

from opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem import OnesCountMaxProblem
from opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem_binary_bit_array_solution import OnesCountMaxProblemBinaryBitArraySolution

class OnesCountMaxProblemBinaryBitArraySolutionVnsShakingSupport(VnsShakingSupport[BitArray,str]):
    
    def __init__(self)->None:
        """
        Create new `OnesCountMaxProblemBinaryBitArraySolutionVnsShakingSupport` instance
        """

    def __copy__(self):
        """
        Internal copy of the `OnesCountMaxProblemBinaryBitArraySolutionVnsShakingSupport`

        :return: new `OnesCountMaxProblemBinaryBitArraySolutionVnsShakingSupport` instance with the same properties
        :rtype: `OnesCountMaxProblemBinaryBitArraySolutionVnsShakingSupport`
        """
        sol = deepcopy(self)
        return sol

    def copy(self):
        """
        Copy the `OnesCountMaxProblemBinaryBitArraySolutionVnsShakingSupport` instance

        :return: new `OnesCountMaxProblemBinaryBitArraySolutionVnsShakingSupport` instance with the same properties
        :rtype: `OnesCountMaxProblemBinaryBitArraySolutionVnsShakingSupport`
        """
        return self.__copy__()

    def shaking(self, k:int, problem:OnesCountMaxProblem, solution:OnesCountMaxProblemBinaryBitArraySolution, 
            optimizer:Algorithm)->bool:
        """
        Random shaking of k parts such that new solution code does not differ more than k from all solution codes 
        inside shakingPoints 

        :param int k: int parameter for VNS
        :param `OnesCountMaxProblem` problem: problem that is solved
        :param `OnesCountMaxProblemBinaryBitArraySolution` solution: solution used for the problem that is solved
        :param `Algorithm` optimizer: optimizer that is executed
        :return: if randomization is successful
        :rtype: bool
        """    
        if optimizer.finish_control.is_finished(optimizer.evaluation, optimizer.iteration, optimizer.elapsed_seconds()):
            return False
        tries:int = 0
        limit:int = 10000
        while tries < limit:
            repres:BitArray = BitArray(solution.representation)
            positions:list[int] = []
            for _ in range(0,k):
                positions.append(choice(range(len(repres))))
            for pos in positions:
                repres.invert(pos)
            solution.representation = repres
            all_ok:bool = True
            if solution.representation.count(value=1) > problem.dimension:
                all_ok = False
            if all_ok:
                break
        if tries < limit:
            if optimizer.finish_control.is_finished(optimizer.evaluation, optimizer.iteration, optimizer.elapsed_seconds()):
                return False
            optimizer.write_output_values_if_needed("before_evaluation", "b_e")
            optimizer.evaluation += 1
            solution.evaluate(problem)
            optimizer.write_output_values_if_needed("after_evaluation", "a_e")
            optimizer.write_output_values_if_needed("after_step_in_iteration", "shaking")
            return True
        else:
            return False 

    def string_rep(self, delimiter:str, indentation:int=0, indentation_symbol:str='', group_start:str ='{', 
        group_end:str ='}')->str:
        """
        String representation of the vns support structure

        :param delimiter: delimiter between fields
        :type delimiter: str
        :param indentation: level of indentation
        :type indentation: int, optional, default value 0
        :param indentation_symbol: indentation symbol
        :type indentation_symbol: str, optional, default value ''
        :param group_start: group start string 
        :type group_start: str, optional, default value '{'
        :param group_end: group end string 
        :type group_end: str, optional, default value '}'
        :return: string representation of vns support instance
        :rtype: str
        """        
        return 'OnesCountMaxProblemBinaryBitArraySolutionVnsShakingSupport'

    def __str__(self)->str:
        """
        String representation of the vns support instance

        :return: string representation of the vns support instance
        :rtype: str
        """
        return self.string_rep('|')

    def __repr__(self)->str:
        """
        Representation of the vns support instance

        :return: string representation of the vns support instance
        :rtype: str
        """
        return self.string_rep('\n')


    def __format__(self, spec:str)->str:
        """
        Formatted the vns support instance

        :param str spec: format specification
        :return: formatted vns support instance
        :rtype: str
        """
        return self.string_rep('|')


