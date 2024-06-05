"""
In this experiment, we as how different solvers perform on randomly generated
distribute instances.
"""

import os
import json
import sys
from random import seed

import eugene.solvers.base_min_generations_enumerator as enum
import eugene.solvers.base_min_generations_enumerator_dominance as edom
import eugene.solvers.base_min_generations_segment as eseg
from eugene.plant_models.plant2 import PlantSPC

TIMEOOUT = 300


def solver_enum(n_loci, pop_0):
    return enum.breeding_program(n_loci, pop_0, timeout=TIMEOOUT)


def solver_edom(n_loci, pop_0):
    return edom.breeding_program(n_loci, pop_0, timeout=TIMEOOUT)


def solver_eseg(n_loci, pop_0):
    return eseg.breeding_program(n_loci, pop_0, timeout=TIMEOOUT)


SOLVERS = {
    # "CANZAR": None,
    "ENUM": solver_enum,
    "EDOM": solver_edom,
    "ESEG": solver_eseg,
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
    elif solver.upper() not in SOLVERS:
        raise ValueError("solver not available")
    else:
        solvers = [(solver.upper(), SOLVERS[solver.upper()])]

    # Parse loci_string
    if loci_string is None:
        loci_string = "2-1000000000"

    if loci_string.isnumeric():
        N_LOCI = [eval(loci_string)]
    elif loci_string.count("-") == 1 and all(
        s.isnumeric() for s in loci_string.split("-")
    ):
        s, e = (int(c) for c in loci_string.split("-"))
        N_LOCI = range(s, e + 1)
    elif loci_string.count("-") == 0 and all(
        s.isnumeric() for s in loci_string.split(",")
    ):
        N_LOCI = [eval(s) for s in loci_string.split(",")]
    else:
        raise ValueError("incorrect format for loci set")

    # Parse pop_string
    if pop_string is None:
        pop_string = "2,4,6,8"

    if pop_string.isnumeric():
        N_POP = [eval(pop_string)]
    elif all(
        s.isnumeric() for s in pop_string.split(",")
    ):
        N_POP = [eval(s) for s in pop_string.split(",")]
    else:
        raise ValueError("incorrect format for pop set")

    for line in sys.stdin:
        try:
            obj = json.loads(line)
            n_loci = obj["n_loci"]
            n_pop = obj["n_pop"]

            if n_loci not in N_LOCI or n_pop not in N_POP:
                continue

            inst = eval(obj["pop_0"])
        except json.JSONDecodeError as e:
            print(f"Error decoding JSON: {e}", file=sys.stderr)
        except Exception as e:
            print(f"Unexpected error: {e}", file=sys.stderr)

        for name, solver in solvers:
            start = time.time()
            result = solver(n_loci, inst)
            time_res = time.time() - start
            if output_file == "":
                print(
                    json.dumps(
                        {
                            "solver": name,
                            "n_loci": n_loci,
                            "n_pop": n_pop,
                            "time_l1": time_res,
                            "objective": result.objective
                            if result is not None
                            else float("inf"),
                            "success": True if result is not None else False,
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
                                "time_l1": time_res,
                                "objective": result.objective
                                if result is not None
                                else float("inf"),
                                "success": True if result is not None else False,
                            }
                        ),
                        file=f,
                    )
