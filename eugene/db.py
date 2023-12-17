from abc import ABC, abstractmethod
from typing import Callable


class DB(ABC):

    @abstractmethod
    def __getitem__(self, instance):
        pass

    @abstractmethod
    def __setitem__(self, instance, solution):
        pass


class DistributeDB(DB):
    """
    Database for storing distribute instances
    """

    def __init__(self):
        self.db = {}

    def __getitem__(self, instance):
        try:
            return self.db[str(instance)]
        except ValueError:
            return None

    def __setitem__(self, instance, solution):
        self.db[str(instance)] = solution


class DBSolver:

    def __init__(self, solver, db: DB):
        self.solver = solver
        self.db = db

    def solve(self, instance):
        solution = self.db[instance]
        if solution is None:
            solution = self.solver(instance)
            self.db[instance] = solution
        return solution
