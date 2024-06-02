import eugene_rs
import eugene.utils as eu
from eugene.plant_models.plant2 import PlantSPC
from eugene.solution import BaseSolution
from typing import List


def breeding_program_distribute(dist_array: List[int]) -> BaseSolution:
    """
    Solves a distribute instance using the minizinc model.
    A MinizincContext can be passed in as an optional parameter,
    otherwise a MinizincContext is constructed using cp-sat from OR-Tools.
    """
    n_loci = len(dist_array)
    pop_0 = eu.distribute_to_plants(dist_array)

    return breeding_program(n_loci, pop_0)


def breeding_program(n_loci: int, pop_0: List[PlantSPC]) -> BaseSolution:
    """
    Solves a distribute instance using the minizinc model.
    A MinizincContext can be passed in as an optional parameter,
    otherwise a MinizincContext is constructed using cp-sat from OR-Tools.
    """

    result = eugene_rs.min_gen.segment.breeding_program_python(
        n_loci,
        [
            [[bool(allele) for allele in row] for row in x.to_bitlist()]
            for x in pop_0
        ],
    )
    return BaseSolution(*result)
