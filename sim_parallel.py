from collections import namedtuple, Counter
from sortedcollections import SortedDict
from plant import Population, generate_goal
from plant import prob_z_given_xy_fast, number_of_trials_to_create
from plant2 import PlantSPC, WDataG, WDataP
from sim import filter_non_dominating


def breeding_par(
    n_loci=10, pop_0=[], pruning=True, constraint_time=None, constraint_resource=None
):
    Results = namedtuple(
        "Results", ["success", "n_generations", "n_plants_max", "n_plants_tot",],
    )
    # check feasibility
    if not union(pop_0) == union([generate_goal(n_loci)]):
        return Results(False, 0, 0, 0)

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
        # remove random subset
        t += 1

    return Results(True, t, 0, 0)


def generate_goal(n_loci) -> PlantSPC:
    return PlantSPC(n_loci, (1 << n_loci) - 1, (1 << n_loci) - 1,)


def filter_non_dominating_gametes(gametes):
    n = len(gametes)
    to_keep = [True] * n
    if n > 0 and isinstance(gametes[0], int):
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
    elif n > 0 and isinstance(gametes[0], WDataG):
        for i, x in enumerate(gametes):
            for j in range(i):
                if not to_keep[j]:
                    continue
                y = gametes[j]
                if y.dom_gamete(x):
                    to_keep[i] = False
                    break
                if x.dom_gamete(y):
                    to_keep[j] = False
    return [gametes[i] for i in range(n) if to_keep[i]]


def dom_gamete(x, y):
    """
    Returns tru if y dominates x i.e. x ≤ y.
    """
    return not (x & ~y)


def generate_all_progeny(n_loci: int, pop: Population):
    """
    Returns a list of all unique progeny that can be produced from the given
    population.
    """
    plant_type = type(pop[0])
    gametes = list({x.gamete_specified(k1) for x in pop for k1 in x.crosspoints()})
    # n_gametes = len(gametes)

    gametes_fil = filter_non_dominating_gametes(gametes)
    n_gametes_fil = len(gametes_fil)

    # TODO: can we break the symmetry on gametes?
    # return [PlantSPC(n_loci, gx, gy) for gx in gametes for gy in gametes]
    # gonna assume yes, we can
    # progeny = [
    #     PlantSPC(n_loci, gametes[i], gametes[j])
    #     for i in range(n_gametes)
    #     # for j in range(i, n_gametes)
    #     for j in range(n_gametes)
    # ]
    progeny_fil = [
        plant_type.from_gametes(n_loci, gametes_fil[i], gametes_fil[j])
        for i in range(n_gametes_fil)
        # for j in range(i, n_gametes_fil)
        for j in range(n_gametes_fil)
    ]

    # assert sorted(filter_non_dominating(progeny)) == sorted(progeny_fil)
    return progeny_fil


def generate_all_progeny_wdata(n_loci: int, pop: Population):
    """
    Returns a list of all unique progeny that can be produced from the given
    population.
    """
    gametes = list({x.gamete_specified(k1) for x in pop for k1 in x.crosspoints()})
    # n_gametes = len(gametes)

    gametes_fil = filter_non_dominating_gametes(gametes)
    n_gametes_fil = len(gametes_fil)

    # TODO: can we break the symmetry on gametes?
    # return [PlantSPC(n_loci, gx, gy) for gx in gametes for gy in gametes]
    # gonna assume yes, we can
    # progeny = [
    #     PlantSPC(n_loci, gametes[i], gametes[j])
    #     for i in range(n_gametes)
    #     # for j in range(i, n_gametes)
    #     for j in range(n_gametes)
    # ]
    progeny_fil = [
        PlantSPC(n_loci, gametes_fil[i], gametes_fil[j])
        for i in range(n_gametes_fil)
        # for j in range(i, n_gametes_fil)
        for j in range(n_gametes_fil)
    ]

    # assert sorted(filter_non_dominating(progeny)) == sorted(progeny_fil)
    return progeny_fil


def generate_all_progeny_1(pop: Population):
    return list(
        {
            x.cross_specified(y, (k1, k2))
            for x in pop
            for y in pop
            for k1 in x.crosspoints()
            for k2 in y.crosspoints()
        }
    )


def union(pop: Population):
    out = 0
    for x in pop:
        out |= x.chrom1
        out |= x.chrom2
    return out


class Experiments:
    @staticmethod
    def basic_tests():
        """
        :returns: TODO

        """
        pass

    @staticmethod
    def enumerations():
        """
        TODO: Docstring for enumerations.
        """
        for n_loci in range(1, 17):
            print(f"n_loci: {n_loci}")
            res_hom = breeding_par(n_loci, PlantSPC.initial_pop_singles_homo(n_loci))
            print("hom:", res_hom)
            res_het = breeding_par(n_loci, PlantSPC.initial_pop_singles_hetero(n_loci))
            print("het:", res_het)
            print(
                [
                    breeding_par(
                        n_loci, PlantSPC.initial_pop_random(n_loci, n_loci)
                    ).n_generations
                    for _ in range(10)
                ]
            )
            # res_tin = breeding_par(
            # n_loci,
            # PlantSPC.initial_pop_trait_introgression(n_loci, n_loci >> 1)
            # )
            # print("tin:", res_tin)

    @staticmethod
    def trait_introgression_dispersion(n_loci):
        """
        How do the number of donor loci affect the run times?

        :n_loci: number of loci per individual plant
        :returns: TODO

        """
        for donor_loci in range(n_loci // 2 + 1):
            results = []
            for _ in range(10_000):
                pop_0 = PlantSPC.initial_pop_trait_introgression(n_loci, donor_loci)
                res = breeding_par(n_loci, pop_0)
                results.append(res)
            print(f"{donor_loci} / {n_loci}")
            print(SortedDict(Counter([x.n_generations for x in results])))

    @staticmethod
    def tree_extraction(n_loci):
        """
        Extract a tree to the target from an enumeration run.

        :n_loci: TODO
        :returns: TODO

        """
        # pop_0 = PlantSPC.initial_pop_singles_hetero(n_loci)
        pop_0_unfiltered = PlantSPC.initial_pop_trait_introgression(n_loci, n_loci >> 1)
        print(pop_0_unfiltered)

        pop_0 = [WDataP(x, history=i) for i, x in enumerate(pop_0_unfiltered)]

        Results = namedtuple(
            "Results", ["success", "n_generations", "n_plants_max", "n_plants_tot",],
        )
        # check feasibility
        if not union(pop_0) == union([generate_goal(n_loci)]):
            return Results(False, 0, 0, 0)

        pop = pop_0.copy()
        ideotype = generate_goal(n_loci)
        t = 0

        pop = filter_non_dominating(pop)

        while not any(ideotype == plant.x for plant in pop):
            # generate reachable genotypes
            pop = generate_all_progeny(n_loci, pop)
            # remove random subset
            t += 1

        # find target
        target = None
        for plant in pop:
            if plant.x == ideotype:
                target = plant
                break

        print(target.history)
        # print(solve_crosspoints(pop_0_unfiltered))

        return Results(True, t, 0, 0)


if __name__ == "__main__":
    from pprint import pprint
    from random import seed

    # seed(0)

    # Experiments.tree_extraction(8)
    Experiments.trait_introgression_dispersion(8)
