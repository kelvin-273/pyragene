"""
In this experiment, we as how different solvers perform on randomly generated
distribute instances.
"""

import os
import json
from random import seed

import eugene.solvers.base_min_crossings_astar as east
import eugene.solvers.base_min_crossings_mip as emip
import eugene.solvers.base_min_crossings_minizinc as emzn
from eugene.plant_models.plant2 import PlantSPC

seed(0)

N_LOCI = list(range(2, 11))
N_POP = [2, 4, 6, 8]
N_INST = 100
TIMEOUT = 300

INSTANCES = {
    n_loci: {
        n_pop: [
            PlantSPC.initial_pop_random(
                n_loci=n_loci, n_individuals=n_pop, p=0.5
            )
            for _ in range(N_INST)
        ]
        for n_pop in N_POP
    }
    for n_loci in N_LOCI
}


def solver_astar(n_loci, pop_0):
    return east.breeding_program(n_loci, pop_0, timeout=TIMEOUT)


CTX_SAT = emzn.MinizincContext.from_solver_and_model_file(
    "sat", "./eugene/solvers/minizinc/mincross.mzn"
)


def solver_cp_sat(n_loci, pop_0):
    return emzn.breeding_program(
        n_loci, pop_0, ctx=CTX_SAT, timeout=TIMEOUT
    )


CTX_MIP = emzn.MinizincContext.from_solver_and_model_file(
    "gurobi", "./eugene/solvers/minizinc/mincross.mzn"
)


def solver_cp_mip(n_loci, pop_0):
    return emzn.breeding_program(
        n_loci, pop_0, ctx=CTX_MIP, timeout=TIMEOUT
    )


def solver_mip(n_loci, pop_0):
    return emip.breeding_program(n_loci, pop_0, timeout=TIMEOUT)


SOLVERS = {
    # "CANZAR": None,
    "CP-SAT": solver_cp_sat,
    "CP-MIP": solver_cp_mip,
    "MIP": solver_mip,
    "A*": solver_astar,
}

if __name__ == "__main__":
    import argparse
    import time

    parser = argparse.ArgumentParser(
        "Experiment runner for MinCross algorithms on random instances"
    )
    parser.add_argument("-s", "--solver")
    parser.add_argument("-o", "--output_file")
    parser.add_argument("-l", "--loci")
    parser.add_argument("-p", "--pop")
    args = parser.parse_args()
    print(args)

    solver = args.solver
    output_file = args.output_file
    loci_string = args.loci
    pop_string = args.pop

    if output_file is None:
        output_file = ""

    if os.path.exists(output_file):
        print(f"output_file exists: {output_file}")
        exit()

    output_dir = os.path.dirname(output_file)
    if output_dir != "" and not os.path.exists(output_dir):
        # make directory
        os.makedirs(output_dir)
        # then we can assume that output_file does not exist

    if solver is None:
        solvers = SOLVERS.items()
    elif solver.upper() == "ASTAR":
        solvers = [("A*", SOLVERS["A*"])]
    elif solver.upper() not in SOLVERS:
        raise ValueError("solver not available")
    else:
        solvers = [(solver.upper(), SOLVERS[solver.upper()])]

    if loci_string is None:
        loci_string = "2-8"

    if loci_string.isnumeric():
        N_LOCI = [eval(loci_string)]
    elif loci_string.count("-") == 1 and all(
        s.isnumeric() for s in loci_string.split("-")
    ):
        s, e = (int(c) for c in loci_string.split("-"))
        N_LOCI = list(range(s, e + 1))
    elif loci_string.count("-") == 0 and all(
        s.isnumeric() for s in loci_string.split(",")
    ):
        N_LOCI = [eval(s) for s in loci_string.split(",")]
    else:
        raise ValueError("incorrect format for loci set")

    if pop_string is None:
        pop_string = "2,4,6,8"

    if pop_string.isnumeric():
        N_POP = [eval(pop_string)]
    elif all(
        s.isnumeric() and s in "2468" for s in pop_string.split(",")
    ):
        N_POP = [eval(s) for s in pop_string.split(",")]
    else:
        raise ValueError("incorrect format for pop set")

    for name, solver in solvers:
        for n_loci in N_LOCI:
            for n_pop in N_POP:
                instances = INSTANCES[n_loci][n_pop]
                objectives = [None] * N_INST
                times_l1 = [0] * N_INST
                start = time.time()
                for i, inst in enumerate(instances):
                    start = time.time()
                    result = solver(n_loci, inst)
                    time_res = time.time() - start
                    objectives[i] = None if result is None else result.objective
                    times_l1[i] = time_res
                if output_file == "":
                    print(
                        json.dumps(
                            {
                                "solver": name,
                                "n_loci": n_loci,
                                "n_pop": n_pop,
                                "times_l1": times_l1,
                                "objectives": objectives,
                            }
                        )
                    )
                else:
                    with open(output_file, "a") as f:
                        print(
                            json.dumps(
                                {
                                    "solver": name,
                                    "n_loci": n_loci,
                                    "n_pop": n_pop,
                                    "times_l1": times_l1,
                                    "objectives": objectives,
                                }
                            ),
                            file=f
                        )
