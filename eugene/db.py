import json
from abc import ABC, abstractmethod
from typing import List, Optional
from eugene.utils import cmp_le_distribute_arrays
from eugene.solution import BaseSolution


class DB(ABC):
    @property
    @abstractmethod
    def db(self):
        pass

    @abstractmethod
    def __getitem__(self, instance):
        pass

    @abstractmethod
    def __setitem__(self, instance, solution):
        pass

    @abstractmethod
    def __contains__(self, instance, solution):
        pass


class DeserializableABC(ABC):

    @abstractmethod
    def to_json_file(self, fp):
        pass

    @staticmethod
    @abstractmethod
    def from_json_file(fp):
        pass

    @staticmethod
    @abstractmethod
    def default():
        pass


class DistributeDB(DB, DeserializableABC):
    """
    Database for storing distribute instances
    """

    def __init__(self, db={}):
        self._db = db

    @property
    def db(self):
        return self._db

    def __getitem__(self, instance):
        return self.db[str(instance)]

    def __setitem__(self, instance, solution):
        self.db[str(instance)] = solution

    def __contains__(self, instance):
        return self[instance] is not None

    def default():
        return DistributeDB()

    def get_base_solution(self, instance):
        try:
            return BaseSolution.from_dict(self.db[str(instance)])
        except KeyError:
            return None

    def to_json_file(self, fp):
        json.dump(self._db, fp)

    def from_json_file(fp):
        return DistributeDB(db=json.load(fp))


class DBSolver:
    def __init__(self, solver, db: DB, filename: Optional[str] = None):
        self.solver = solver
        self.db = db

    def solve(self, instance) -> BaseSolution:
        try:
            solution = BaseSolution.from_dict(self.db[instance])
        except KeyError:
            solution = self.solver(instance)
            self.db[instance] = solution.to_dict()
        return solution

    def __call__(self, instance):
        return self.solve(instance)


def distribute_db_from_json(filename: str):
    with open(filename) as f:
        return DistributeDB(json.load(f))


def parse_dist_array(s: str) -> List[int]:
    return eval(s)


def distribute_db_dict_to_list(d: dict):
    class A:
        def __init__(self, arr):
            self.arr = arr

        def __lt__(self, other):
            return (
                cmp_le_distribute_arrays(self.arr, other.arr) and self.arr != other.arr
            )

    dist_arrays = [A(parse_dist_array(s)) for s in d.keys()]
    dist_arrays.sort()
    return [{"instance": a.arr, "solution": d[str(a.arr)]} for a in dist_arrays]


class CachedDB(DB):
    def __init__(
        self,
        filename: str,
        db_cls: DB,
    ):
        self._filename = filename
        self._db_cls = db_cls

    @property
    def db(self):
        return self._db

    def __enter__(self):
        try:
            with open(self._filename) as f:
                self._db = self._db_cls.from_json_file(f)
        except FileNotFoundError:
            self._db = self._db_cls.default()
        finally:
            return self

    def __exit__(self, *args):
        with open(self._filename, "w") as f:
            self.db.to_json_file(f)

    def __contains__(self, item):
        return self.db.__contains__(item)

    def __getitem__(self, item):
        return self.db.__getitem__(item)

    def __setitem__(self, key, value):
        self.db[key] = value


if __name__ == "__main__":
    db = distribute_db_from_json("./distribute_data.json")
    print(json.dumps(distribute_db_dict_to_list(db.db)))
