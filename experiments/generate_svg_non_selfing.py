import eugene.solvers.base_min_crossings_minizinc as em
import eugene.plant_models.plant2 as ep2

ctx = em.MinizincContext.from_solver_and_model_file(
    "sat", "eugene/solvers/minizinc/modelGenotypesNonSelfed.mzn"
)

n_loci = 8
n_pop = 4


def trial():
    pop_0 = ep2.PlantSPC.initial_pop_random(n_loci, n_pop)
    solution = em.breeding_program(n_loci, pop_0, ctx)
    return solution


if __name__ == "__main__":
    print(repr(trial()))
