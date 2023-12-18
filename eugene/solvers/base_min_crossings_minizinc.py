import minizinc


def breeding_program_distribute(n_loci, dist_array, ctx=None):
    if ctx is None:
        ctx = MinizincContext.from_solver_and_model_file(
            # "or-tools", "../minizinc/genotype/modelGenotypes.mzn"
            "sat",
            "../minizinc/genotype/modelGenotypes.mzn",
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
    from eugene.db import DistributeDB
    import signal
    import json

    dist_array = [0, 1, 0, 2, 0, 1, 3, 2, 3, 4, 1]
    __import__("pprint").pprint(
        (breeding_program_distribute(len(dist_array), dist_array))
    )

    db = DistributeDB()

    def handle_interupt(signum, frame, ask=True):
        with open("distribute_data.json", "w") as f:
            json.dump(db.db, f)
        exit(0)

    signal.signal(signal.SIGINT, handle_interupt)

    for n_loci in range(1, 10):
        for dist_array in gen_distribute_instances(n_loci):
            print(dist_array)
            solution = breeding_program_distribute(n_loci, dist_array)
            db[dist_array] = solution

    handle_interupt(None, None)
