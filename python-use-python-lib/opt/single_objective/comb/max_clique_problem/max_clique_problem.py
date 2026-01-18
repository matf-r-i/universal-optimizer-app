""" 
..  _py_minimum_multi_cut_problem:

"""

import sys
import os
from pathlib import Path
from typing import Optional
directory = Path(__file__).resolve()
sys.path.append(directory.parent)
sys.path.append(directory.parent.parent)
sys.path.append(directory.parent.parent.parent)
sys.path.append(directory.parent.parent.parent.parent)
root_dir = directory.parent.parent.parent.parent.parent
sys.path.append(str(root_dir))

import networkx as nx
import json

from uo.problem.problem import Problem
from uo.utils.logger import logger

class MaxCliqueProblem(Problem):
    """
    Class representing the Max Clique Problem.

    This class inherits from the `Problem` class and is used to define and solve the Max Clique Problem. The problem is defined by a graph.

    Attributes:
        graph (nx.Graph): Graph.

    Methods:
        __init__(graph: nx.Graph, source_terminal_pairs: list): Initializes a new instance of the MaxCliqueProblem class.
        
        __load_from_files__(graph_file_path: str) -> tuple: Static function that reads problem data from files.
        from_input_files(graph_file_path: str): Creates a new MaxCliqueProblem instance when the input file and input format are specified.
        __copy__() -> MaxCliqueProblem: Internal copy of the MaxCliqueProblem problem.
        copy() -> MaxCliqueProblem: Copy the MaxCliqueProblem problem.
        graph() -> nx.Graph: Property getter for the graph of the target problem.
        string_rep(delimiter: str, indentation: int = 0, indentation_symbol: str = '', group_start: str = '{', group_end: str = '}') -> str: String representation of the MaxCliqueProblem instance.
        __str__() -> str: String representation of the MaxCliqueProblem structure.
        __repr__() -> str: Representation of the MaxCliqueProblem instance.
        __format__() -> str: Formatted MaxCliqueProblem instance.
    """
    
    def __init__(self, graph:nx.Graph)->None:
        """
        Create new `MaxCliqueProblem` instance

        :param nx.Graph graph: graph of the problem
        """
        if not isinstance(graph, nx.Graph):
            raise TypeError('Parameter \'graph\' for  MaxCliqueProblem should be \'nx.Graph\'.')
        super().__init__(name="MaxCliqueProblem", is_minimization=False, is_multi_objective=False)
        self.__graph = graph

    def copy(self):
        obj:MaxCliqueProblem = MaxCliqueProblem( self.graph.copy())
        return obj
            
    @classmethod
    def from_graph(cls, graph:nx.Graph):
        """
        Additional constructor. Create new `MaxCliqueProblem` instance when graph and source_terminal_pairs are specified

        :param nx.Graph graph: graph of the problem
        """
        return cls(graph)

    @classmethod
    def __load_from_file__(cls,file_path:str, verbosity: bool=False)->nx.Graph:

        edges = []
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('c'):  # graph description
                    if verbosity:
                        print(*line.split()[1:])
                # first line: p name num_of_vertices num_of_edges
                elif line.startswith('p'):
                    _, name, vertices_num, edges_num = line.split()
                    if verbosity:
                        print('{0} {1} {2}'.format(name, vertices_num, edges_num))
                elif line.startswith('e'):
                    _, v1, v2 = line.split()
                    edges.append((v1, v2))
                else:
                    continue
        return nx.Graph(edges)


    @classmethod
    def from_input_file(cls, graph_file_path: str, verbosity: bool=False)->'MaxCliqueProblem':
        """
        Additional constructor. Create new `MaxCliqueProblem` instance when input file and input format are specified

        :param str input_file_path: path of the input file with problem graph

        :return: class instance
        :rtype: MaxCliqueProblem
        """
        result:nx.Graph = MaxCliqueProblem.__load_from_file__(graph_file_path, verbosity=verbosity)
        graph:nx.Graph = result
        return cls(graph=graph)

    @property
    def graph(self)->nx.Graph:
        """
        Property getter for graph of the target problem

        :return: graph of the target problem instance 
        :rtype: nx.Graph
        """
        return self.__graph
    
    def string_rep(self, delimiter:str, indentation:int=0, indentation_symbol:str='', group_start:str ='{', 
        group_end:str ='}')->str:
        """
        String representation of the `MaxCliqueProblem` instance

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
        :return: string representation of instance that controls output
        :rtype: str
        """          
        s = delimiter
        for _ in range(0, indentation):
            s += indentation_symbol  
        s += group_start
        s+= super().string_rep(delimiter, indentation, indentation_symbol, '', '')
        s+= delimiter
        for _ in range(0, indentation):
            s += indentation_symbol  
        s += 'graph=' + str(self.__graph) + delimiter
        s += group_end 
        return s

    def __str__(self)->str:
        """
        String representation of the minimum multi cut problem structure

        :return: string representation of the minimum multi cut problem structure
        :rtype: str
        """
        return self.string_rep('|', 0, '', '{', '}')


    def __repr__(self)->str:
        """
        Representation of the minimum multi cut problem instance
        :return: str -- string representation of the minimum multi cut problem instance
        """
        return self.string_rep('\n', 0, '   ', '{', '}')

    def __format__(self)->str:
        """
        Formatted the minimum multi cut problem instance
        :return: str -- formatted minimum multi cut problem instance
        """
        return self.string_rep('|')


