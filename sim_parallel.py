from collections import namedtuple
from plant import Population, generate_goal
from plant import prob_z_given_xy_fast, number_of_trials_to_create
from plant2 import PlantSPC
from sim import filter_non_dominating


def breeding_par(
    n_loci=10, pop_0=[], pruning=True, constraint_time=None, constraint_resource=None
):
    Results = namedtuple(
        "Results", ["n_generations", "n_plants_max", "n_plants_tot", "success"]
    )
    # check feasibility
    if not union(pop_0) == union([generate_goal(n_loci)]):
        return Results(0, 0, 0, False)

    pop = pop_0.copy()
    ideotype = generate_goal(n_loci)
    t = 0
    n_tot = 0
    n_max = 0

    if pruning:
        pop = filter_non_dominating(pop)

    while ideotype not in pop and (not constraint_time or t < constraint_time):
        # generate reachable genotypes
        pop = generate_all_progeny(n_loci, pop)
        # prune reachable set
        # pop = filter_non_dominating(pop)
        t += 1

    return Results(
        t, 0, 0, True
    )


def generate_goal(n_loci) -> PlantSPC:
    return PlantSPC(n_loci, (1 << n_loci) - 1, (1 << n_loci) - 1,)


def filter_non_dominating_gametes(gametes):
    n = len(gametes)
    to_keep = [True] * n
    for i, x in enumerate(gametes):
        for j in range(i):
            if not to_keep[j]:
                continue
            y = gametes[j]
            if dom_gamete(x, y):
                to_keep[i] = False
                break
            if dom_gamete(y, x):
                to_keep[j] = False
    return [gametes[i] for i in range(n) if to_keep[i]]


def dom_gamete(x, y):
    """
    Returns tru if y dominates x i.e. x ≤ y.
    """
    return not (x & ~y)


def generate_all_progeny(n_loci: int, pop: Population):
    gametes = list({
        x.gamete_specified(k1)
        for x in pop
        for k1 in x.crosspoints()
    })
    gametes = filter_non_dominating_gametes(gametes)
    n_gametes = len(gametes)
    # TODO: can we break the symmetry on gametes?
    # return [PlantSPC(n_loci, gx, gy) for gx in gametes for gy in gametes]
    # gonna assume yes, we can
    return [PlantSPC(n_loci, gametes[i], gametes[j])
            for i in range(n_gametes)
            for j in range(i, n_gametes)
            ]


def generate_all_progeny_1(pop: Population):
    return list({
        x.cross_specified(y, (k1, k2))
        for x in pop
        for y in pop
        for k1 in x.crosspoints()
        for k2 in y.crosspoints()
    })


def union(pop: Population):
    out = 0
    for x in pop:
        out |= x.chrom1
        out |= x.chrom2
    return out


if __name__ == "__main__":
    from pprint import pprint
    from random import seed
    seed(1)

    # pprint(generate_all_progeny([PlantSPC(4, 5, 10)]))
    # print(len(generate_all_progeny([PlantSPC(4, 5, 10)])))
    # print()
    print(breeding_par(4, [PlantSPC(4, 5, 10)]))
    print(breeding_par(10, PlantSPC.initial_pop_trait_introgression(10, 1)))
    print(PlantSPC.initial_pop_singles_homo(1))
    print(PlantSPC.initial_pop_singles_homo(2))
    print(PlantSPC.initial_pop_singles_homo(3))
    print()
    n_loci = 32
    # # for i in range(1, (n_loci >> 1) + 1):
    #     # print(f"i: {i}")
    #     # for _ in range(10):
    #         # print(breeding_par(
    #             # n_loci,
    #             # PlantSPC.initial_pop_singles_homo(n_loci)
    #         # ))
    #     # print()
    for n_loci in range(1, 17):
        print(f"n_loci: {n_loci}")
        res_hom = breeding_par(
            n_loci,
            PlantSPC.initial_pop_singles_homo(n_loci)
        )
        print("hom:", res_hom)
        res_het = breeding_par(
            n_loci,
            PlantSPC.initial_pop_singles_hetero(n_loci)
        )
        print("het:", res_het)
        res_tin = breeding_par(
            n_loci,
            PlantSPC.initial_pop_trait_introgression(n_loci, n_loci)
        )
        print("tin:", res_tin)
