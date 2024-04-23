import eugene_rs as eurs
from eugene.solution import BaseSolution
from eugene.plant_models.plant2 import PlantSPC
from eugene.utils import distribute_to_plants
from typing import List


def breeding_program_distribute(dist_array: List[int]) -> BaseSolution:
    return breeding_program(len(dist_array), distribute_to_plants(dist_array))


def breeding_program(n_loci: int, pop_0: List[PlantSPC]) -> BaseSolution:
    pop_1 = [[[bool(b) for b in c] for c in x.to_bitlist()] for x in pop_0]
    res = eurs.breeding_program_python(n_loci, pop_1)
    return BaseSolution(*res)


if __name__ == "__main__":
    n_loci = 2
    print(breeding_program(n_loci, [
        PlantSPC(
            n_loci,
            0b10,
            0b10,
        ),
        PlantSPC(
            n_loci,
            0b01,
            0b01,
        ),
    ]))

    print("yo")

    n_loci = 3
    print(breeding_program(n_loci, [
        PlantSPC(
            n_loci,
            0b101,
            0b101,
        ),
        PlantSPC(
            n_loci,
            0b010,
            0b010,
        ),
    ]))

    n_loci = 3
    print(breeding_program(n_loci, [
        PlantSPC(
            n_loci,
            0b100,
            0b100,
        ),
        PlantSPC(
            n_loci,
            0b010,
            0b010,
        ),
        PlantSPC(
            n_loci,
            0b001,
            0b001,
        ),
    ]))
