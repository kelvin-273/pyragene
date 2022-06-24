from eugene.plant_models.plant2 import PlantSPC
from eugene.simulators.enumerators import BreedingProgram


def all_good_run(n_loci, pop_0):
    runner = BreedingProgram(n_loci, PlantSPC)
    runner.set_ideotype(PlantSPC(
        n_loci,
        (1 << n_loci) - 1,
        (1 << n_loci) - 1,
    ))
    runner.set_init_pop(pop_0)
    runner.run()
    runner.print_path_to_goal()


if __name__ == "__main__":
    n_loci = 5
    all_good_run(n_loci, [PlantSPC(n_loci, 27, 27), PlantSPC(n_loci, 4, 4)])
    print()
    # all_good_run(n_loci, [PlantSPC(n_loci, 10, 21)])
    # print()
    # all_good_run(n_loci, [PlantSPC(n_loci, 10, 10), PlantSPC(n_loci, 21, 21)])
