from uo.solution.quality_of_solution import QualityOfSolution
from uo.solution.solution import Solution
from uo.problem.problem import Problem
from uo.utils.logger import logger
from linopy import Model
import xarray as xr
from datetime import datetime
from copy import deepcopy
from dataclasses import dataclass
import typing
from typing import Optional

import sys
import os
from pathlib import Path
directory = Path(__file__).resolve()
sys.path.append(directory)
sys.path.append(directory.parent)
sys.path.append(directory.parent.parent)
sys.path.append(directory.parent.parent.parent)
sys.path.append(directory.parent.parent.parent.parent)
root_dir = directory.parent.parent.parent.parent.parent
sys.path.append(str(root_dir))


class OnesCountMaxProblemIntegerLinearProgrammingSolution(Solution[xr.DataArray, str]):

    def __init__(self, random_seed=None, fitness_value=None, fitness_values=None,                        objective_value=None, objective_values=None, is_feasible=True) -> None:
        """
        Create new `OnesCountMaxProblemIntegerLinearProgrammingSolution` instance
        """
        super().__init__(
            random_seed, fitness_value, fitness_values, objective_value, objective_values, is_feasible)

    def __init__(self, sol: xr.DataArray, random_seed=None, fitness_value=None, fitness_values=None, objective_value=None, objective_values=None, is_feasible=True) -> None:
        super().__init__(
            random_seed, fitness_value, fitness_values, objective_value, objective_values, is_feasible)
        try:
            self.representation = deepcopy(sol)
        except TypeError:
            self.representation = sol

    def copy(self) -> 'OnesCountMaxProblemIntegerLinearProgrammingSolution':
        obj: OnesCountMaxProblemIntegerLinearProgrammingSolution = OnesCountMaxProblemIntegerLinearProgrammingSolution(
            self.representation, self.random_seed, self.fitness_value, self.fitness_values, self.objective_value, self.objective_values, self.is_feasible)
        return obj

    def copy_from(self, original: Solution) -> None:
        self.__init__(
            original.representation, original.random_seed, original.fitness_value, original.fitness_values, original.objective_value, original.objective_values, original.is_feasible)

    def argument(self, representation: object) -> str:
        return str(representation)

    def init_random(self, problem: Problem) -> None:
        pass

    def init_from(self, representation: object, problem: Problem) -> None:
        if not isinstance(problem, Problem):
            raise TypeError(
                'Parameter \'problem\' must be of type \'Problem\'.')
        self.representation = representation

    def native_representation(self, representation_str: str) -> object:
        return representation_str

    def calculate_quality_directly(self, representation: object, problem: Problem) -> QualityOfSolution:
        return QualityOfSolution(0, None, 0, None, True)

    def representation_distance_directly(self, solution_code_1: str, solution_code_2: str) -> float:
        return 0

    def string_representation(self):
        x: str = 'OnesCountMaxProblemIntegerLinearProgrammingSolution - ['
        for v in self.representation.values:
            x += str(int(v)) + ', '
        x = x[:-2] + '] \n'
        x += super().string_representation()
        return x

    def string_rep(self, separator: str) -> str:
        return separator.join(str(int(v)) for v in self.representation.values)

    def __str__(self) -> str:
        return self.string_rep("|")

    def __repr__(self) -> str:
        return self.string_rep("|")

    def __format__(self, spec: str) -> str:
        return self.string_rep(spec)
