"""
..  _py_knapsack_problem_bit_array_solution:
"""

import sys
from pathlib import Path

directory = Path(__file__).resolve()
sys.path.append(str(directory))
sys.path.append(str(directory.parent))
sys.path.append(str(directory.parent.parent.parent))
sys.path.append(str(directory.parent.parent.parent.parent))
root_dir = directory.parent.parent.parent.parent.parent
sys.path.append(str(root_dir))

from typing import Optional
from random import random

from bitstring import BitArray

from uo.problem.problem import Problem
from uo.solution.quality_of_solution import QualityOfSolution
from uo.solution.solution import Solution

from opt.single_objective.comb.knapsack_problem.knapsack_problem import KnapsackProblem


class KnapsackProblemBitArraySolution(Solution[BitArray, str]):
    """
    BitArray-based solution for the 0/1 Knapsack Problem.
    """

    def __init__(
        self,
        random_seed: Optional[int] = None,
        evaluation_cache_is_used: bool = False,
        evaluation_cache_max_size: int = 0,
        distance_calculation_cache_is_used: bool = False,
        distance_calculation_cache_max_size: int = 0
    ) -> None:
        """
        Create new `KnapsackProblemBitArraySolution` instance.
        """
        if not isinstance(random_seed, int) and random_seed is not None:
            raise TypeError("Parameter 'random_seed' must be 'int' or 'None'.")

        super().__init__(
            random_seed=random_seed,
            fitness_value=None,
            fitness_values=None,
            objective_value=None,
            objective_values=None,
            is_feasible=False,
            evaluation_cache_is_used=evaluation_cache_is_used,
            evaluation_cache_max_size=evaluation_cache_max_size,
            distance_calculation_cache_is_used=distance_calculation_cache_is_used,
            distance_calculation_cache_max_size=distance_calculation_cache_max_size
        )
        self.is_minimization = False

    def copy(self) -> "KnapsackProblemBitArraySolution":
        """
        Internal copy of the solution.
        """
        sol = KnapsackProblemBitArraySolution(self.random_seed)
        sol.copy_from(self)
        return sol

    def copy_from(self, original) -> None:
        """
        Copy all data from the original target solution.
        """
        super().copy_from(original)

    def argument(self, representation: BitArray) -> str:
        """
        Convert internal representation to solution code.
        """
        return representation.bin

    def init_random(self, problem: Problem) -> None:
        """
        Random initialization of the solution.
        """
        if not hasattr(problem, "weights") or problem.weights is None:
            raise ValueError("Can not randomly initialize solution without item weights.")
        if not hasattr(problem, "values") or problem.values is None:
            raise ValueError("Can not randomly initialize solution without item values.")

        self.representation = BitArray(len(problem.weights))
        for i in range(len(self.representation)):
            if random() > 0.5:
                self.representation[i] = True

    def init_from(self, representation: BitArray, problem: Problem) -> None:
        """
        Initialization of the solution by setting its native representation.
        """
        if not isinstance(representation, BitArray):
            raise TypeError("Parameter 'representation' must have type 'BitArray'.")
        if len(representation) == 0:
            raise ValueError("Representation must have positive length.")
        self.representation = BitArray(bin=representation.bin)

    def is_feasible_sol(
        self,
        representation: BitArray,
        capacity: int,
        weights: list[int]
    ) -> bool:
        total_weight = 0
        for i in range(representation.len):
            if representation[i]:
                total_weight += weights[i]
        return total_weight <= capacity

    def calc_fitness(
        self,
        representation: BitArray,
        capacity: int,
        weights: list[int],
        values: list[int]
    ) -> tuple[bool, float, float]:
        total_weight = 0
        total_value = 0

        for i in range(representation.len):
            if representation[i]:
                total_weight += weights[i]
                total_value += values[i]

        if total_weight > capacity:
            overflow = total_weight - capacity
            objective_value = float("-inf")
            fitness_value = float(-(1000000 * overflow + total_weight))
            return (False, objective_value, fitness_value)

        objective_value = float(total_value)
        fitness_value = float(total_value) - 0.001 * float(total_weight)
        return (True, objective_value, fitness_value)

    def calculate_quality_directly(
        self,
        representation: BitArray,
        problem: KnapsackProblem
    ) -> QualityOfSolution:
        """
        Fitness calculation of the knapsack BitArray solution.
        """
        if not isinstance(representation, BitArray):
            raise TypeError("Parameter 'representation' must have type 'BitArray'.")
        if not isinstance(problem, KnapsackProblem):
            raise TypeError("Parameter 'problem' must have type 'KnapsackProblem'.")

        if len(representation) != problem.dimension:
            raise ValueError("Representation length must match problem dimension.")

        feasible, objective_value, fitness_value = self.calc_fitness(
            representation, problem.capacity, problem.weights, problem.values
        )

        return QualityOfSolution(
            fitness_value=fitness_value,
            fitness_values=None,
            objective_value=objective_value,
            objective_values=None,
            is_feasible=feasible
        )

    def native_representation(self, representation_str: str) -> BitArray:
        """
        Native solution representation from its string representation.
        """
        if not isinstance(representation_str, str):
            raise TypeError("Parameter 'representation_str' must be 'str'.")
        if any(ch not in ("0", "1") for ch in representation_str):
            raise ValueError("Representation string should contain only '0' and '1'.")
        return BitArray(bin=representation_str)

    def representation_distance_directly(
        self,
        solution_code_1: str,
        solution_code_2: str
    ) -> int:
        """
        Hamming distance between two solution codes.
        """
        if not isinstance(solution_code_1, str):
            raise TypeError("Parameter 'solution_code_1' should be 'str'.")
        if not isinstance(solution_code_2, str):
            raise TypeError("Parameter 'solution_code_2' should be 'str'.")
        if len(solution_code_1) != len(solution_code_2):
            raise ValueError("Representations should have the same length.")

        return sum(1 for a, b in zip(solution_code_1, solution_code_2) if a != b)

    def string_rep(
        self,
        delimiter: str = "\n",
        indentation: int = 0,
        indentation_symbol: str = "   ",
        group_start: str = "{",
        group_end: str = "}"
    ) -> str:
        """
        String representation of the solution instance.
        """
        s = group_start
        s += super().string_rep(
            delimiter=delimiter,
            indentation=indentation,
            indentation_symbol=indentation_symbol,
            group_start="",
            group_end=""
        )
        s += delimiter + "string_representation()="
        s += "" if self.representation is None else self.representation.bin
        s += group_end
        return s

    def __str__(self) -> str:
        return self.string_rep("|", 0, "", "{", "}")

    def __repr__(self) -> str:
        return self.string_rep("\n", 0, "   ", "{", "}")

    def __format__(self, spec: str = "") -> str:
        return self.string_rep("|")