
from datetime import datetime
import unittest   
import unittest.mock as mocker

from uo.problem.problem_void_min_so import ProblemVoidMinSO

from uo.algorithm.output_control import OutputControl
from uo.utils import logger

from opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem import OnesCountMaxProblem
from opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem_ilp_linopy import OnesCountMaxProblemIntegerLinearProgrammingSolver
from opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem_ilp_linopy import OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters
from uo.solution.solution import Solution
from uo.solution.solution_void_representation_int import SolutionVoidInt
from uo.solution.solution_void_representation_object import SolutionVoidObject

from opt.single_objective.comb.ones_count_max_problem.ones_count_max_problem_ilp_solution import OnesCountMaxProblemIntegerLinearProgrammingSolution

class TestOnesCountMaxProblemIlpLinopy(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        print("setUpClass TestIntegrationOnesCountMaxProblemIlpLinopy\n")

    def setUp(self):
        self.problem:OnesCountMaxProblem = OnesCountMaxProblem.from_dimension(dimension=5) 
        self.solver: OnesCountMaxProblemIntegerLinearProgrammingSolver = OnesCountMaxProblemIntegerLinearProgrammingSolver(problem=self.problem)
        self.bs: OnesCountMaxProblemIntegerLinearProgrammingSolution = self.solver.optimize()

            
    def test_best_solution_after_optimization_should_be_optimal(self):
        #Assert
        result = str(self.bs)
        expected = ''
        for i in range(self.bs.representation.size):
            expected += '1'
        self.assertEqual(expected, result)

    # creating an instance of OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters with valid OutputControl and Problem parameters should return an instance of the class with the same parameters
    def test_valid_parameters(self):
        # Arrange
        output_control = OutputControl(moments="after_algorithm")
        problem = ProblemVoidMinSO('problem_name', False)
        # Act
        construction_params = OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters(
            output_control=output_control, problem=problem)
        # Assert
        self.assertIsInstance(construction_params, OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters)
        self.assertEqual(construction_params.output_control, output_control)
        self.assertEqual(construction_params.problem, problem)

    # creating an instance of OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters with default parameters should return an instance of the class with None parameters
    def test_default_parameters(self):
        # Arrange
        problem = ProblemVoidMinSO('problem_name', False)
        # Act & Assert
        with self.assertRaises(TypeError):
            construction_params = OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters()

    # creating an instance of OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters with invalid OutputControl parameter should raise a TypeError
    def test_invalid_output_control(self):
        # Arrange
        output_control = "invalid_output_control"
        problem = ProblemVoidMinSO('problem_name', False)
        # Act & Assert
        with self.assertRaises(TypeError):
            OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters(output_control, problem)

    # creating an instance of OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters with invalid Problem parameter should raise a TypeError
    def test_invalid_problem(self):
        # Arrange
        problem = "invalid_problem"
        # Act & Assert
        with self.assertRaises(TypeError):
            OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters(problem)

    # creating an instance of OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters with OutputControl and Problem parameters of different types should raise a TypeError
    def test_different_types(self):
        # Arrange
        problem = SolutionVoidInt(42, None, None, False)
        # Act & Assert
        with self.assertRaises(TypeError):
            OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters(problem)

    # creating an instance of OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters with OutputControl and Problem parameters of the same type but different from OutputControl and Problem should raise a TypeError
    def test_same_types_different_classes(self):
        # Arrange
        problem = SolutionVoidObject()
        # Act & Assert
        with self.assertRaises(TypeError):
            OnesCountMaxProblemIntegerLinearProgrammingSolverConstructionParameters(problem)

    def tearDown(self):
        return

    @classmethod
    def tearDownClass(cls):
        print("\ntearDownClass TestIntegrationOnesCountMaxProblemIlpLinopy")
    
if __name__ == '__main__':
    unittest.main()


class TestIntegration(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("setUpClass TestIntegration\n")

    def setUp(self):
        self.problem:OnesCountMaxProblem = OnesCountMaxProblem.from_dimension(dimension=5) 
        self.solver: OnesCountMaxProblemIntegerLinearProgrammingSolver = OnesCountMaxProblemIntegerLinearProgrammingSolver(problem=self.problem)
        self.bs: OnesCountMaxProblemIntegerLinearProgrammingSolution = self.solver.optimize()

    def tearDown(self):
        return super().tearDown()
    
    # The method creates an instance of Model and adds variables to it.
    def test_model_instance_and_variables_added(self):
        # Assert
        self.assertGreater(len(self.solver.model.variables), 0)

    # The method sets the objective function of the model to minimize or maximize the sum of the variables.
    def test_invalid_instance_type_error(self):
        # Arrange
        output_control = "invalid_output_control"
        problem = OnesCountMaxProblem(dim=5)
        # Act & Assert
        with self.assertRaises(TypeError):
            solver = OnesCountMaxProblemIntegerLinearProgrammingSolver(output_control, problem)
            bs = solver.optimize()
        output_control = OutputControl()
        problem = "invalid_problem"
        # Act & Assert
        with self.assertRaises(TypeError):
            solver = OnesCountMaxProblemIntegerLinearProgrammingSolver(output_control, problem)
            bs = solver.optimize()

    # The method solves the model and sets the best solution to the solution of the model.
    def test_model_solved_and_best_solution_set(self):
        # Arrange
        problem = OnesCountMaxProblem(dim=5)
        solver = OnesCountMaxProblemIntegerLinearProgrammingSolver(problem=problem)
        # Act
        bs = solver.optimize()
        # Assert
        self.assertIsNotNone(bs)
        self.assertIsNotNone(solver.best_solution)
        

class TestStringRep(unittest.TestCase):

    # Returns a string representation of the 'OnesCountMaxProblemIntegerLinearProgrammingSolver' instance
    def test_returns_string_representation(self):
        # Arrange
        problem = OnesCountMaxProblem(dim=5)
        solver = OnesCountMaxProblemIntegerLinearProgrammingSolver(problem=problem)
        solver.execution_started = datetime.now()
        self.best_solution_mock = mocker.MagicMock(spec=Solution)
        self.best_solution_mock.copy = mocker.Mock(return_value=self.best_solution_mock)
        self.best_solution_mock.string_rep = mocker.Mock(return_value="solution mock")
        solver.best_solution = self.best_solution_mock
        # Act
        result = solver.string_rep("|")
        # Assert
        self.assertIsInstance(result, str)

    # The string representation contains the name of the class and its properties
    def test_contains_class_name_and_properties(self):
        # Arrange
        output_control = OutputControl()
        problem = OnesCountMaxProblem(dim=5)
        solver = OnesCountMaxProblemIntegerLinearProgrammingSolver(output_control, problem)
        solver.execution_started = datetime.now()
        self.best_solution_mock = mocker.MagicMock(spec=Solution)
        self.best_solution_mock.copy = mocker.Mock(return_value=self.best_solution_mock)
        self.best_solution_mock.string_rep = mocker.Mock(return_value="solution mock")
        solver.best_solution = self.best_solution_mock    
        # Act
        result = solver.string_rep("|")
        # Assert
        self.assertIn("OnesCountMaxProblemIntegerLinearProgrammingSolver", result)
        self.assertIn("output_control", result)
        self.assertIn("problem", result)

    # The string representation is properly formatted with indentation and grouping symbols
    def test_properly_formatted_with_indentation_and_grouping_symbols(self):
        # Arrange
        output_control = OutputControl()
        problem = OnesCountMaxProblem(dim=5)
        solver = OnesCountMaxProblemIntegerLinearProgrammingSolver(output_control, problem)
        solver.execution_started = datetime.now()
        self.best_solution_mock = mocker.MagicMock(spec=Solution)
        self.best_solution_mock.copy = mocker.Mock(return_value=self.best_solution_mock)
        self.best_solution_mock.string_rep = mocker.Mock(return_value="solution mock")
        solver.best_solution = self.best_solution_mock    
        # Act
        result = solver.string_rep("|", indentation=2, indentation_symbol="-", group_start="[", group_end="]")    
        # Assert
        self.assertIn( "output_control", result)
        self.assertIn( "problem", result)

 