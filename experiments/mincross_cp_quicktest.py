from time import time
import eugene.solvers.base_min_crossings_minizinc as emzn
import eugene.utils as eutl

MIP_CTX = emzn.MinizincContext.from_solver_and_model_file(
    model_file="./eugene/solvers/minizinc/mincross_distribute.mzn",
    solver="sat"
)

for n_loci in range(3, 10):
    for i in range(10):
        case = eutl.random_distribute_instance(n_loci)
        with MIP_CTX.branch() as ctx:
            time_start = time()
            res = emzn.breeding_program_distribute(n_loci, case, ctx)
            time_delta = time() - time_start
            print(f"n_loci: {n_loci}\tinstance: {case}\tsolve-time: {time_delta}\tn_instances: {i+1}")
