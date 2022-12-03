import matplotlib.pyplot as plt

from eugene.plant_models.plant2 import PlantSPC
from eugene.simulators.greedy_time import BreedingProgramGreedyTime


def homo(n_loci: int):
    pop_0 = PlantSPC.initial_pop_singles_homo(n_loci)
    runner = BreedingProgramGreedyTime(n_loci, PlantSPC)
    runner.set_ideotype(PlantSPC(n_loci, (1 << n_loci) - 1, (1 << n_loci) - 1,))
    runner.set_init_pop(pop_0)
    res = runner.run()
    return height(res)


def hetero(n_loci: int):
    pop_0 = PlantSPC.initial_pop_singles_hetero(n_loci)
    runner = BreedingProgramGreedyTime(n_loci, PlantSPC)
    runner.set_ideotype(PlantSPC(n_loci, (1 << n_loci) - 1, (1 << n_loci) - 1,))
    runner.set_init_pop(pop_0)
    res = runner.run()
    return height(res)


def height(x, depth=0, plant_type=None):
    if plant_type is None:
        plant_type = type(x)

    if x.histories is None:
        return depth

    assert isinstance(x.histories, (list, tuple))

    return max(
        depth,
        *(
            height(y, depth=depth + (type(y) is plant_type), plant_type=plant_type)
            for y in x.histories
        )
    )


if __name__ == "__main__":
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True)

    ax2.scatter(
        [n for n in range(1, 16)],
        [homo(n) for n in range(1, 16)]
    )
    ax2.set_title("$t_{\max}$ for homozygous_singles")
    ax2.set_ylabel("number of generations")
    ax2.set_xlabel("number of loci")
    ax2.set_xticks(range(1, 16))
    ax2.set_yticks(range(0, 6))

    ax1.scatter(
        [n for n in range(1, 16)],
        [hetero(n) for n in range(1, 16)]
    )
    ax1.set_title("$t_{\max}$ for heterozygous_singles")
    ax1.set_ylabel("number of generations")
    ax1.set_xlabel("number of loci")
    ax1.set_xticks(range(1, 16))
    ax1.set_yticks(range(0, 6))

    plt.show()
