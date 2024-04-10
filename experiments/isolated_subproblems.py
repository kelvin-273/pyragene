import itertools
from random import seed
from math import ceil

import eugene.utils as eu
import eugene.db as ed
import eugene.solvers.base_min_crossings_minizinc as em
import eugene.solvers.base_min_crossings_distribute as edist


seed(0)
minizinc_ctx = em.DEFAULT_CTX


def solver(case):
    return em.breeding_program_distribute(
        len(case), case, ctx=minizinc_ctx, processes=4
    )


def lb(case):
    return ceil((len(case) + max(case) + 1) / 2)


with ed.CachedDB("distribute_data_2.json", ed.DistributeDB) as db:
    solver1 = ed.DBSolver(solver, db)

    def solver2(case): return edist.breeding_program_distribute(
        len(case), case, sub_solver=lambda _, case: solver1(case)
    )
    solver = ed.DBSolver(solver2, db)
    for n_loci in range(2, 16):
        for case in eu.gen_distribute_instances(n_loci):
            case_subs = eu.distribute_to_isolated_subproblems(case)
            if len(case_subs) > 1:
                case_lb = lb(case)
                obj_subs = [
                    solver(eu.sanitise_distribute_array(case[s:e+1])).objective if e - s > 1 else 1
                    for s, e in case_subs
                ]
                obj_heu = sum(obj_subs)
                if obj_heu != case_lb:
                    sol_opt = solver(case)
                    obj_opt = sol_opt.objective
                    print(case)
                    print(sol_opt)
                    print(case_subs)
                    print(obj_subs)
                    print(case_lb, obj_opt, obj_heu)
