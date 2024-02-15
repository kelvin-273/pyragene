"""
Here we run distribute instances until we find one where there is a mismatch
between generations and crossings
"""

import json
from random import sample, choices
import eugene.utils as eu
import eugene.solvers.base_min_crossings_wedges as ew
import eugene.solvers.base_min_crossings_minizinc as em
import eugene.db as ed

ctx_c = em.MinizincContext.from_solver_and_model_file(
    "sat", "eugene/solvers/minizinc/mincross_with_gen.mzn"
)

ctx_g = em.MinizincContext.from_solver_and_model_file(
    "sat", "eugene/solvers/minizinc/mingen_with_cross.mzn"
)

db_cross = ed.distribute_db_from_json("./distribute_data_2.json")

# solver_c = ed.DBSolver(
#     lambda instance: em.breeding_program_distribute(len(instance), instance), db_cross
# )
solver_c = lambda instance: em.breeding_program_distribute(len(instance), instance, ctx_c)

# solver_g = lambda instance: ew.breeding_program(
#     len(instance),
#     instance,
#     wedges_selector=lambda instance: ew.max_repeated_wedges(
#         instance, strict_mingen=True
#     ),
# )
solver_g = lambda instance: em.breeding_program_distribute(len(instance), instance, ctx_g)


def experiment1():
    solutions = []
    for n_loci in range(2, 10):
        print(f"n_loci: {n_loci}")
        n_trials = 100
        cases = [eu.random_distribute_instance(n_loci) for _ in range(n_trials)]
        sol_opt = [None] * n_trials
        sol_heu = [None] * n_trials
        for i, case in enumerate(cases):
            opt = solver_c(case)
            heu = solver_g(case)
            og = opt.generations()
            hg = heu.generations()
            oc = opt.crossings()
            hc = heu.crossings()
            if og > hg and oc < hc:
                print(case, (oc, hc, hc - oc), (og, hg, hg - og), sep="\t")
            sol_opt[i] = opt
            sol_heu[i] = heu
        solutions.append((sol_opt, sol_heu))

    import matplotlib.pyplot as plt

    xs = list(range(2, 10))
    ys_diff = [
        sum(he.objective - op.objective for op, he in zip(sols_opt, sols_heu))
        / len(sols_opt)
        for sols_opt, sols_heu in solutions
    ]
    ys_diff_rel = [
        sum(he.objective / op.objective for op, he in zip(sols_opt, sols_heu))
        / len(sols_opt)
        for sols_opt, sols_heu in solutions
    ]


def experiment2():
    n_trials = 100
    for n_loci in range(4, 6):
        print(n_loci)
        cases = list(eu.gen_distribute_instances(n_loci))
        for dist_arr in cases:
            res_c = solver_c(dist_arr)
            res_g = solver_g(dist_arr)
            gen_c = res_c.generations()
            gen_g = res_g.generations()
            cro_c = res_c.crossings()
            cro_g = res_g.crossings()
            if cro_c < cro_g and gen_c > gen_g:
                print(dist_arr, (cro_c, cro_g, cro_g - cro_c),
                      (gen_c, gen_g, gen_g - gen_c), sep='\t')

    for n_loci in range(6, 9):
        print(n_loci)
        cases = sorted(sample(list(eu.gen_distribute_instances(n_loci)), k=n_trials))
        for dist_arr in cases:
            res_c = solver_c(dist_arr)
            res_g = solver_g(dist_arr)
            gen_c = res_c.generations()
            gen_g = res_g.generations()
            cro_c = res_c.crossings()
            cro_g = res_g.crossings()
            if cro_c < cro_g and gen_c > gen_g:
                print(dist_arr, (cro_c, cro_g, cro_g - cro_c),
                      (gen_c, gen_g, gen_g - gen_c), sep='\t')

    for n_loci in range(9, 14):
        print(n_loci)
        cases = sorted([eu.random_distribute_instance(n_loci) for _ in range(n_trials)])
        for dist_arr in cases:
            res_c = solver_c(dist_arr)
            res_g = solver_g(dist_arr)
            gen_c = res_c.generations()
            gen_g = res_g.generations()
            cro_c = res_c.crossings()
            cro_g = res_g.crossings()
            if cro_c < cro_g and gen_c > gen_g:
                print(dist_arr, (cro_c, cro_g, cro_g - cro_c),
                      (gen_c, gen_g, gen_g - gen_c), sep='\t')


def experiment3():
    n_trials = 100
    d = {}

    for n_loci in range(1, 7):
        cases = list(eu.gen_distribute_instances(n_loci))
        d[n_loci] = cases

    for n_loci in range(7, 11):
        cases = sorted(sample(list(eu.gen_distribute_instances(n_loci)), k=n_trials))
        d[n_loci] = cases

    for n_loci in range(11, 16):
        cases = [None] * n_trials
        for i in range(n_trials):
            while cases[i] is None:
                xs = eu.random_distribute_instance(n_loci)
                if xs not in cases[:i]:
                    cases[i] = xs
        d[n_loci] = cases

    print(json.dumps(d))


if __name__ == "__main__":
    experiment3()
