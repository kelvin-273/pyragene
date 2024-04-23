"""
In this experiment, we ask how different solvers perform on randomly generated
distribute instances.

We'll measure the expected the run-time for randomly generated distribute
instances. As the number of distribute instances for n_loci = 1..6 is less
than 100, we shall simply take the expected value over all loci. The number
instances where n_loci = 7..8 is less than 1000, so these can be chosen from
the full set. Instances where n_loci ≥ 9can sampled by generating and checking
against a list of previously observed instances.
"""

import eugene.utils as eu
import eugene.solvers.base_min_crossings_astar as east
import eugene.solvers.base_min_crossings_mip as emip
import eugene.solvers.base_min_crossings_minizinc as emzn

import time
import matplotlib.pyplot as plt


class TimeOut:
    pass


def measure_performance(foo, n_loci):
    instances = [eu.random_distribute_instance(n_loci)]
    total_time = 0
    while total_time < 10:
        start = time.time()
        for instance in instances:
            foo(instance)
        end = time.time()
        total_time = end - start
        result_time = total_time / len(instances)
        instances.extend(
            [
                eu.random_distribute_instance(n_loci)
                for _ in range(len(instances))
            ]
        )
    return result_time


def solver_astar(dist_array):
    return east.breeding_program_distribute(dist_array)


CTX_SAT = emzn.MinizincContext.from_solver_and_model_file(
    "sat", "./eugene/solvers/minizinc/modelGenotypes.mzn"
)


def solver_cp_sat(dist_array):
    return emzn.breeding_program_distribute(
        len(dist_array), dist_array, ctx=CTX_SAT
    )


CTX_MIP = emzn.MinizincContext.from_solver_and_model_file(
    "gurobi", "./eugene/solvers/minizinc/modelGenotypes.mzn"
)


def solver_cp_mip(dist_array):
    return emzn.breeding_program_distribute(
        len(dist_array), dist_array, ctx=CTX_MIP
    )


def solver_mip(dist_array):
    return emip.breeding_program_distribute(dist_array)


SOLVERS = {
    # "CANZAR": None,
    "CP-SAT": solver_cp_sat,
    "CP-MIP": solver_cp_mip,
    "MIP": solver_mip,
    "A*": solver_astar,
}

RESULTS = {
    "CP-SAT": [None] * 8,
    "CP-MIP": [None] * 8,
    "MIP": [None] * 8,
    "A*": [None] * 8,
}

# for n_loci in range(1, 5):
#     print(n_loci)
#     print(RESULTS)
#     RESULTS["A*"][n_loci - 1] = measure_performance(SOLVERS["A*"], n_loci)

# for n_loci in range(1, 5):
#     print(n_loci)
#     print(RESULTS)
#     RESULTS["MIP"][n_loci - 1] = measure_performance(SOLVERS["MIP"], n_loci)

RESULTS = {
    "CP-SAT": [None, None, None, None, None, None, None, None],
    "CP-MIP": [None, None, None, None, None, None, None, None],
    "MIP": [
        0.0005094455336802639,
        0.003369288519024849,
        0.023645638953894377,
        0.19293437525629997,
        None,
        None,
        None,
        None,
    ],
    "A*": [
        2.3788530597812496e-05,
        0.0001384845036227489,
        0.45167065411806107,
        109.64520788192749,
        None,
        None,
        None,
        None,
    ],
}

for n_loci in range(1, 9):
    print(n_loci)
    print(RESULTS)
    RESULTS["CP-SAT"][n_loci - 1] = measure_performance(
        SOLVERS["CP-SAT"], n_loci
    )
    RESULTS["CP-MIP"][n_loci - 1] = measure_performance(
        SOLVERS["CP-MIP"], n_loci
    )


for name, results in RESULTS.items():
    plt.scatter(list(range(1, 9)), RESULTS[name], label=name)
plt.title("Runtime vs #genes")
plt.xlabel("#genes")
plt.ylabel("Runtime (s)")
plt.yscale("log")
plt.show()
