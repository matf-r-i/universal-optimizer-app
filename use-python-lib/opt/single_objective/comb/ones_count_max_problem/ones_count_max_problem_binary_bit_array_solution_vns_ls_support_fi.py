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
from uo.algorithm.metaheuristic.variable_neighborhood_search.vns_ls_support import VnsLocalSearchSupport

from opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem import OnesCountMaxProblem
from opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem_binary_bit_array_solution import OnesCountMaxProblemBinaryBitArraySolution

class OnesCountMaxProblemBinaryBitArraySolutionVnsLocalSearchSupportFI(VnsLocalSearchSupport[BitArray,str]):
    
    def __init__(self)->None:
        """
        Create new `OnesCountMaxProblemBinaryBitArraySolutionVnsLocalSearchSupportFI` instance
        """

    def __copy__(self):
        """
        Internal copy of the `OnesCountMaxProblemBinaryBitArraySolutionVnsLocalSearchSupportFI`

        :return: new `OnesCountMaxProblemBinaryBitArraySolutionVnsLocalSearchSupportFI` instance with the same properties
        :rtype: `OnesCountMaxProblemBinaryBitArraySolutionVnsLocalSearchSupportFI`
        """
        sol = deepcopy(self)
        return sol

    def copy(self):
        """
        Copy the `OnesCountMaxProblemBinaryBitArraySolutionVnsLocalSearchSupportFI` instance

        :return: new `OnesCountMaxProblemBinaryBitArraySolutionVnsLocalSearchSupportFI` instance with the same properties
        :rtype: `OnesCountMaxProblemBinaryBitArraySolutionVnsLocalSearchSupportFI`
        """
        return self.__copy__()
    
    def local_search(self, k:int, problem:OnesCountMaxProblem, solution:OnesCountMaxProblemBinaryBitArraySolution, 
            optimizer: Algorithm)->bool:
        """
        Executes "first improvement" variant of the local search procedure 
        
        :param int k: int parameter for VNS
        :param `OnesCountMaxProblem` problem: problem that is solved
        :param `OnesCountMaxProblemBinaryBitArraySolution` solution: solution used for the problem that is solved
        :param `Algorithm` optimizer: optimizer that is executed
        :return: result of the local search procedure 
        :rtype: if local search is successful
        """
        if optimizer.finish_control.is_finished(optimizer.evaluation, optimizer.iteration, optimizer.elapsed_seconds()):
            return False
        if k < 1 or k > problem.dimension:
            return False
        start_sol:OnesCountMaxProblemBinaryBitArraySolution = solution.copy()
        # initialize indexes
        indexes:ComplexCounterUniformAscending = ComplexCounterUniformAscending(k, problem.dimension)
        in_loop:bool = indexes.reset()
        while in_loop:
            # collect positions for inversion from indexes
            positions:list[int] = indexes.current_state()
            # invert and compare, switch and exit if new is better
            solution.representation.invert(positions) 
            if optimizer.finish_control.is_finished(optimizer.evaluation, optimizer.iteration, optimizer.elapsed_seconds()):
                solution.copy_from(start_sol)
                return False
            optimizer.write_output_values_if_needed("before_evaluation", "b_e")
            optimizer.evaluation += 1
            solution.evaluate(problem)
            optimizer.write_output_values_if_needed("after_evaluation", "a_e")
            if optimizer.is_first_better(solution, start_sol, problem):
                return True
            solution.representation.invert(positions)
            # increment indexes and set in_loop accordingly
            in_loop = indexes.progress()
        solution.copy_from(start_sol)
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
        return 'OnesCountMaxProblemBinaryBitArraySolutionVnsLocalSearchSupportFI'

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


