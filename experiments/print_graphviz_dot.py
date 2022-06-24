from eugene.plant_models.plant2 import PlantSPC
from eugene.simulators.enumerators import BreedingProgram


def print_dot_viz_experiment(n_loci, pop_0):
    runner = BreedingProgram(n_loci, PlantSPC)
    runner.set_ideotype(
        PlantSPC(
            n_loci,
            (1 << n_loci) - 1,
            (1 << n_loci) - 1
        )
    )
    runner.set_init_pop(pop_0)
    runner.run()
    print(runner.to_dot(runner._pop_current[0]))


if __name__ == "__main__":
    # print_dot_viz_experiment(3, PlantSPC.initial_pop_singles_hetero(3))
    # print_dot_viz_experiment(8, PlantSPC.initial_pop_singles_hetero(8))
    print_dot_viz_experiment(8, PlantSPC.initial_pop_trait_introgression(8, 3, 1, 1))
