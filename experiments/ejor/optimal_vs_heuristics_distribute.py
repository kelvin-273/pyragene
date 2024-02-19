"""
In this experiment, we will explore the gap between optimal optimal methods
and the wedges heuristic on randomly generated distribute instances.
"""

import json

import eugene.utils as eu
import eugene.solvers.base_min_crossings_wedges as ew
import eugene.solvers.base_min_crossings_minizinc as em
import eugene.db as ed

db = ed.distribute_db_from_json("./distribute_data_2.json")
solver_opt = ed.DBSolver(lambda instance: em.breeding_program_distribute(len(instance), instance), db)
solver_heu = lambda instance: ew.breeding_program(len(instance), instance)

# for n_loci in range(2, 10):
#     print()
#     cases = [eu.random_distribute_instance(n_loci) for _ in range(100)]
#     obj_opt = [solver_opt(arr).objective for arr in cases]
#     obj_heu = [solver_heu(arr).objective for arr in cases]
#     for case, opt, heu in zip(cases, obj_opt, obj_heu):
#         print(case, opt, heu, heu - opt, sep='\t')

with open("./data/optgap_instances.json") as f:
    case_dict = json.load(f)

for k, cases in case_dict.items():
    n_loci = str(k)
    obj_opt = [solver_opt(arr).objective for arr in cases]
    obj_heu = [solver_heu(arr).objective for arr in cases]
    for case, opt, heu in zip(cases, obj_opt, obj_heu):
        print(case, opt, heu, heu - opt, sep='\t')
