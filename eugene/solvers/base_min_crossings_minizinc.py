import contextlib
import minizinc
from math import ceil
from typing import List

import eugene.utils as eu
import eugene.plant_models.plant2 as ep2
from eugene.solution import BaseSolution


def breeding_program_distribute(n_loci, dist_array, ctx=None) -> BaseSolution:
    """
    Solves a distribute instance using the minizinc model.
    A MinizincContext can be passed in as an optional parameter,
    otherwise a MinizincContext is constructed using cp-sat from OR-Tools.
    """
    items = instance_array_genotype_homo(
        instance_array=dist_array, max_crossovers=1
    ).items()

    if ctx is None:
        ctx = MinizincContext.from_solver_and_model_file(
            "sat", "./eugene/solvers/minizinc/mincross.mzn",
        )
        instance = ctx.instance

        for k, v in items:
            instance[k] = v

        result = instance.solve(free_search=True)
        return BaseSolution(
            tree_data=result["xs"],
            tree_type=result["treeType"],
            tree_left=result["treeLeft"],
            tree_right=result["treeRight"],
            objective=result.objective,
        )
    else:
        with ctx.instance.branch() as instance:

            for k, v in items:
                instance[k] = v

            result = instance.solve(free_search=True)
            return BaseSolution(
                tree_data=result["xs"],
                tree_type=result["treeType"],
                tree_left=result["treeLeft"],
                tree_right=result["treeRight"],
                objective=result.objective,
            )


def breeding_program(n_loci, pop_0, ctx=None) -> BaseSolution:
    """
    Solves a distribute instance using the minizinc model.
    A MinizincContext can be passed in as an optional parameter,
    otherwise a MinizincContext is constructed using cp-sat from OR-Tools.
    """
    items = instance_array_genotype(pop_0, max_crossovers=1).items()

    if ctx is None:
        ctx = MinizincContext.from_solver_and_model_file(
            "sat", "./eugene/solvers/minizinc/mincross.mzn",
        )
    with ctx.instance.branch() as instance:

        for k, v in items:
            instance[k] = v

        result = instance.solve(free_search=True)
        return BaseSolution(
            tree_data=result["xs"],
            tree_type=result["treeType"],
            tree_left=result["treeLeft"],
            tree_right=result["treeRight"],
            objective=result.objective,
        )


def breeding_program_distribute_optimised(n_loci, dist_array, ctx=None) -> BaseSolution:
    raise NotADirectoryError("still to decide what optimisations will go here")
    # bounds
    n_pop = max(dist_array) + 1
    bound_lower = ceil((n_loci + n_pop) / 2)
    bound_upper = ceil(n_loci / 2) + n_pop * n_pop - n_pop

    leading_mins, mid, trailing_maxes = eu.distribute_sington_decomposition(dist_array)

    # minizinc model
    return breeding_program_distribute(n_loci, dist_array, ctx=ctx)


class MinizincContext:
    """
    A container for the model, solver, and instance configurations for
    minizinc models.
    """

    def __init__(
        self,
        model: minizinc.Model,
        solver: minizinc.Solver,
        instance: minizinc.Instance,
    ):
        self.model = model
        self.solver = solver
        self.instance = instance

    @staticmethod
    def from_solver_and_model_file(solver: str, model_file: str):
        model = minizinc.Model(model_file)
        solver = minizinc.Solver.lookup(solver)
        instance = minizinc.Instance(solver, model)
        return MinizincContext(model, solver, instance)

    @contextlib.contextmanager
    def branch(self):
        with self.instance.branch() as child:
            yield MinizincContext(self.model, self.solver, child)


def instance_array_genotype_homo(instance_array, max_crossovers):
    """
    Creates an instance as a dictionary of an instance of homozygous genotypes
    from an instance array and an integer representing the maximum allowed
    crossovers.
    """
    n_loci = len(instance_array)
    n_gametes = max(instance_array) + 1
    return {
        "maxCrossovers": max_crossovers,
        "nLoci": n_loci,
        "nGenotypes": n_gametes,
        "nTreeCells": n_loci + n_gametes,
        "genotypes": [
            [[1 if j == i else 0 for j in instance_array]] * 2 for i in range(n_gametes)
        ],
    }


def instance_array_genotype(pop_0: List[ep2.PlantSPC], max_crossovers):
    """
    Creates an instance as a dictionary of an instance of homozygous genotypes
    from an instance array and an integer representing the maximum allowed
    crossovers.
    """
    if len(pop_0) == 0:
        raise NotImplementedError("Need a way to account for empty instances")

    n_loci = pop_0[0].n_loci
    n_pop = len(pop_0)
    return {
        "maxCrossovers": max_crossovers,
        "nLoci": n_loci,
        "nGenotypes": n_pop,
        "nTreeCells": n_loci + n_pop,
        "genotypes": [x.to_bitlist() for x in pop_0],
    }


DEFAULT_CTX = MinizincContext.from_solver_and_model_file(
    "sat", "../minizinc/genotype/modelGenotypes.mzn"
)


if __name__ == "__main__":
    from eugene.utils import gen_distribute_instances
    from eugene.db import DistributeDB, distribute_db_from_json, DBSolver
    import signal
    import json
    import eugene.utils as eu

    # db = DistributeDB()
    db = distribute_db_from_json("distribute_data.json")

    def handle_interupt(signum, frame, ask=True):
        with open("distribute_data_2.json", "w") as f:
            json.dump(db.db, f)
        exit(0)

    signal.signal(signal.SIGINT, handle_interupt)

    for n_loci in range(1, 10):
        for dist_array in gen_distribute_instances(n_loci):
            print(dist_array)
            if dist_array not in db:
                solution = breeding_program_distribute(n_loci, dist_array)
                db[dist_array] = solution

    handle_interupt(None, None)
