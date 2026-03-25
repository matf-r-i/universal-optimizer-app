"""
..  _py_knapsack_problem:
"""

import sys
from pathlib import Path

directory = Path(__file__).resolve()
sys.path.append(str(directory.parent))
sys.path.append(str(directory.parent.parent))
sys.path.append(str(directory.parent.parent.parent))
sys.path.append(str(directory.parent.parent.parent.parent))
root_dir = directory.parent.parent.parent.parent.parent
sys.path.append(str(root_dir))

from uo.problem.problem import Problem
from uo.utils.logger import logger


class KnapsackProblem(Problem):
    """
    Class representing the 0/1 Knapsack Problem.

    The problem is defined by:
    - knapsack capacity
    - item weights
    - item values

    A feasible solution is a subset of items whose total weight does not exceed
    the knapsack capacity. The goal is to maximize total value.
    """

    def __init__(self, capacity: int, weights: list[int], values: list[int]) -> None:
        """
        Create new `KnapsackProblem` instance.

        :param int capacity: knapsack capacity
        :param list[int] weights: list of item weights
        :param list[int] values: list of item values
        """
        if not isinstance(capacity, int):
            raise TypeError("Parameter 'capacity' for KnapsackProblem should be 'int'.")
        if not isinstance(weights, list):
            raise TypeError("Parameter 'weights' for KnapsackProblem should be 'list'.")
        if not isinstance(values, list):
            raise TypeError("Parameter 'values' for KnapsackProblem should be 'list'.")
        if capacity <= 0:
            raise ValueError("Parameter 'capacity' must be positive.")
        if len(weights) == 0:
            raise ValueError("Parameter 'weights' must be non-empty.")
        if len(weights) != len(values):
            raise ValueError("Parameters 'weights' and 'values' must have the same length.")
        if any((not isinstance(w, int) or w <= 0) for w in weights):
            raise ValueError("All item weights must be positive integers.")
        if any((not isinstance(v, int) or v < 0) for v in values):
            raise ValueError("All item values must be non-negative integers.")

        super().__init__(
            name="KnapsackProblem",
            is_minimization=False,
            is_multi_objective=False
        )

        self.__capacity = capacity
        self.__weights = weights
        self.__values = values
        self.__dimension = len(weights)

    def copy(self) -> "KnapsackProblem":
        """
        Copy the target problem.
        """
        return KnapsackProblem(
            capacity=self.capacity,
            weights=self.weights.copy(),
            values=self.values.copy()
        )

    @classmethod
    def from_capacity_weights_and_values(
        cls,
        capacity: int,
        weights: list[int],
        values: list[int]
    ) -> "KnapsackProblem":
        """
        Additional constructor when capacity, weights and values are specified.
        """
        return cls(capacity, weights, values)

    @classmethod
    def __load_from_file__(cls, input_file_path: str) -> tuple[int, list[int], list[int]]:
        """
        Static function that reads problem data from specified file.

        Expected file format:
            first line: capacity
            each next line: weight value

        Example:
            15
            2 10
            3 5
            5 15

        :param str input_file_path: path to the input file
        :return: capacity, weights, values
        :rtype: tuple[int, list[int], list[int]]
        """
        logger.debug("Load parameters: input file path=" + str(input_file_path))

        with open(input_file_path, "r", encoding="utf-8") as file:
            lines = [line.strip() for line in file if line.strip()]

        if len(lines) < 2:
            raise ValueError("Input file must contain capacity and at least one item.")

        capacity = int(lines[0])
        weights: list[int] = []
        values: list[int] = []

        for line in lines[1:]:
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(
                    "Each item line must contain exactly two integers: weight and value."
                )
            weight = int(parts[0])
            value = int(parts[1])
            weights.append(weight)
            values.append(value)

        return capacity, weights, values

    @classmethod
    def from_input_file(cls, input_file_path: str) -> "KnapsackProblem":
        """
        Additional constructor. Create new `KnapsackProblem` instance
        when input file and input format are specified.

        :param str input_file_path: path to the input file
        :return: class instance
        :rtype: KnapsackProblem
        """
        capacity, weights, values = cls.__load_from_file__(input_file_path)
        return cls(capacity=capacity, weights=weights, values=values)

    @property
    def capacity(self) -> int:
        """
        Property getter for knapsack capacity.
        """
        return self.__capacity

    @property
    def weights(self) -> list[int]:
        """
        Property getter for item weights.
        """
        return self.__weights

    @property
    def values(self) -> list[int]:
        """
        Property getter for item values.
        """
        return self.__values

    @property
    def dimension(self) -> int:
        """
        Property getter for problem dimension.
        """
        return self.__dimension

    def string_rep(
        self,
        delimiter: str,
        indentation: int = 0,
        indentation_symbol: str = "",
        group_start: str = "{",
        group_end: str = "}"
    ) -> str:
        """
        String representation of the `KnapsackProblem` instance.
        """
        s = delimiter
        for _ in range(0, indentation):
            s += indentation_symbol
        s += group_start
        s += super().string_rep(delimiter, indentation, indentation_symbol, "", "")
        s += delimiter
        s += "capacity=" + str(self.__capacity)
        s += delimiter
        s += "weights=" + str(self.__weights)
        s += delimiter
        s += "values=" + str(self.__values)
        s += delimiter
        s += "dimension=" + str(self.__dimension)
        s += group_end
        return s

    def __str__(self) -> str:
        return self.string_rep("|", 0, "", "{", "}")

    def __repr__(self) -> str:
        return self.string_rep("\n", 0, "   ", "{", "}")

    def __format__(self, spec: str = "") -> str:
        return self.string_rep("|")