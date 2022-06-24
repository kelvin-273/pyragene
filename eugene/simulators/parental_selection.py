from collections import Counter, defaultdict, namedtuple


class BreedingProgram:
    """
    Class to represent an enumerative breeding program
    """

    def __init__(self, n_loci, plant_type):
        self._n_loci = n_loci
        self._plant_type = plant_type
        self._pruning = True
        self._pop_0 = None
        self._ideotype = None
        self._constraint_time = None

        self._Genotype = namedtuple("Genotype", [
            "plant",
            "histories",
        ])

        self._Gamete = namedtuple("Gametes", [
            "gamete",
            "parents",
            "counts",
        ])

    def set_init_pop(self, pop_0):
        """Set the initial population unwrapped"""

        if not type(pop_0) is list:
            raise TypeError("pop_0 must be a list")

        if not all(isinstance(x, self._plant_type) for x in pop_0):
            raise TypeError("Element in pop_0 mismatches the plant type")

        self._pop_0 = pop_0
        self._pop_current = self._pop_0

    def set_ideotype(self, ideotype):
        self._ideotype = ideotype

    def run(self):
        # Pre-flight checks
        if self._pop_0 is None:
            raise ValueError("Initial population not set")

        if self._ideotype is None:
            raise ValueError("ideotype must be set")

        if not isinstance(self._ideotype, self._plant_type):
            raise ValueError("ideotype must be a plant of the same type")

        # create tuple types for genotypes, gametes, and Results
        Genotype = self._Genotype
        Genotype.crosspoints = lambda self: self.plant.crosspoints()
        Genotype.gamete_specified = lambda self, crosspoint: self.plant.gamete_specified(crosspoint)
        Genotype.format = lambda self: str(self.plant)

        Gamete = self._Gamete
        Gamete.__str__ = lambda g: format(g.gamete, f"0{self._n_loci}b")

        Results = namedtuple(
            "Results", ["success", "n_generations", "n_plants_max", "n_plants_tot"],
        )

        # check feasibility
        if not self._plant_type.union(self._pop_0) == self._plant_type.union(
            [self._ideotype]
        ):
            return Results(False, 0, 0, 0)

        # initialise starting variables
        pop = [Genotype(x, None) for i, x in enumerate(self._pop_0)]
        t = 0
        n_tot = 0
        n_max = 0

        if self._pruning:
            # pop = filter_non_dominating(pop)
            pass

        while not any(x.plant == self._ideotype for x in pop) and (
            self._constraint_time is None or t < self._constraint_time
        ):
            # generate reachable genotypes
            gametes = {}
            for plant in pop:
                gametes_per_individual = Counter(
                    plant.gamete_specified(k) for k in plant.crosspoints()
                )
                # update gametes with paths to g from plant
                for g, c in gametes_per_individual.items():
                    if g in gametes:
                        _g, ps, cs = gametes[g]
                        ps.append(plant)
                        gametes[g] = Gamete(_g, ps, cs + c)
                    else:
                        gametes[g] = Gamete(g, [plant], c)
            gametes = list(gametes.values())

            # filter non-dominating gametes
            gametes_fil = filter_non_dominating_gametes(gametes)
            n_gametes_fil = len(gametes_fil)

            # TODO: can we break the symmetry on gametes?
            # construct set of producible progeny
            progeny = {}
            for i in range(n_gametes_fil):
                for j in range(i, n_gametes_fil):
                    g1 = gametes_fil[i]
                    g2 = gametes_fil[j]
                    p = self._plant_type.from_gametes(
                        self._n_loci, g1.gamete, g2.gamete
                    )
                    if p in progeny:
                        progeny[p].histories.append((g1, g2))
                    else:
                        progeny[p] = Genotype(p, [(g1, g2)])
            progeny = list(progeny.values())

            progeny_fil_old = [
                Genotype(
                    plant=self._plant_type.from_gametes(
                        self._n_loci,
                        gametes_fil[i].gamete,
                        gametes_fil[j].gamete
                    ),
                    histories=(
                        gametes_fil[i],
                        gametes_fil[j]
                    ),
                )
                for i in range(n_gametes_fil)
                # for j in range(i, n_gametes_fil)
                for j in range(n_gametes_fil)
            ]

            # assert sorted(filter_non_dominating(progeny)) == sorted(progeny_fil)
            pop = progeny
            self._pop_current = pop

            # remove random subset
            t += 1

        self._results = Results(True, t, None, None)

    def get_results(self):
        return self._results

    def print_path_to_goal(self):
        # Extract the ideotype from the population
        for x in self._pop_current:
            if x.plant == self._ideotype:
                goal = x
                break
        else:
            raise ValueError("goal is not found yet")

        # print via post-order traversal of tree choosing only the first
        # element of each list of histories / parents
        def post_order_traversal(node, depth=0):
            if type(node) is self._Genotype:
                if node.histories is not None:
                    g1, g2 = node.histories[0]
                    post_order_traversal(g1, depth=depth+1)
                    post_order_traversal(g2, depth=depth+1)
                print(f"Genotype: {node.format()}\t generation: {self._results.n_generations - depth}")
            elif type(node) is self._Gamete:
                post_order_traversal(node.parents[0], depth=depth)
                print(node)
            else:
                raise TypeError(f"Got something funky fresh {node}")

        post_order_traversal(goal)


def filter_non_dominating_gametes(gametes):
    n = len(gametes)
    gametes_naked = [x.gamete for x in gametes]
    to_keep = [True] * n
    for i, x in enumerate(gametes_naked):
        for j in range(i):
            if not to_keep[j]:
                continue
            y = gametes_naked[j]
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

