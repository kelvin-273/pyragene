import json
from abc import ABC, abstractmethod
from typing import Callable, List
from eugene.utils import cmp_le_distribute_arrays


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

    def __init__(self, db={}):
        self.db = db

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


def distribute_db_from_json(filename: str):
    with open(filename) as f:
        return DistributeDB(json.load(f))


def parse_dist_array(s: str) -> List[int]:
    return eval(s)


def dic_to_list(d: dict):
    class A:
        def __init__(self, arr):
            self.arr = arr
        def __lt__(self, other):
            return cmp_le_distribute_arrays(self.arr, other.arr) and self != other
    dist_arrays = [A(parse_dist_array(s)) for s in d.keys()]
    dist_arrays.sort()
    return [{"instance": a.arr, "solution": d[str(a.arr)]} for a in dist_arrays]


if __name__ == "__main__":
    db = distribute_db_from_json("./distribute_data.json")
    # __import__('pprint').pprint(dic_to_list(db.db))
    print(json.dumps(dic_to_list(db.db)))
