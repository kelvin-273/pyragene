import minizinc
from math import ceil
import eugene.utils as eu


def breeding_program_distribute(n_loci, dist_array, ctx=None):
    if ctx is None:
        ctx = MinizincContext.from_solver_and_model_file(
            "sat", "../minizinc/genotype/modelGenotypes.mzn",
        )
    instance = ctx.instance

    items = instance_array_genotype_homo(
        instance_array=dist_array, max_crossovers=1
    ).items()

    for k, v in items:
        instance[k] = v

    result = instance.solve(free_search=True)
    return {
        "treeData": result["xs"],
        "treeType": result["treeType"],
        "treeLeft": result["treeLeft"],
        "treeRight": result["treeRight"],
        "objective": result.objective,
    }


def breeding_program_distribute_optimised(n_loci, dist_array, ctx=None):
    # bounds
    n_pop = max(dist_array) + 1
    bound_lower = ceil((n_loci + n_pop) / 2)
    bound_upper = ceil(n_loci / 2) + n_pop * n_pop - n_pop

    leading_mins, mid, trailing_maxes = eu.distribute_sington_decomposition(dist_array)

    # minizinc model
    return breeding_program_distribute(n_loci, dist_array, ctx=ctx)


class MinizincContext:
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

    def branch(self):
        return MinizincContext(self.model, self.solver, self.instance.branch())


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
        "nTreeCells": n_loci * 2,
        "genotypes": [
            [[1 if j == i else 0 for j in instance_array]] * 2 for i in range(n_gametes)
        ],
    }


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
