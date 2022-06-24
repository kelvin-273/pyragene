from itertools import product
from pprint import pprint

from eugene.plant_models.plant2 import PlantSPC
from eugene.simulators.enumerators import BreedingProgram

# self._Genotype = namedtuple("Genotype", ["plant", "histories",])

# self._Gamete = namedtuple("Gametes", ["gamete", "parents", "counts",])


def achievable_by_parental_selection(runner: BreedingProgram, node):
    """
    Assertions:
    - node is the ideotype
      - node_u == node_l == 111..111
        - there can only be one history of node but many parents of node_u
      - we can iterate through each parent of node_u

    """
    d_aux = {}

    recurrent, _ = runner._pop_0

    def f(x, y):
        px, py = x.plant, y.plant
        return px is recurrent and py is not recurrent and aux(y) or \
            py is recurrent and px is not recurrent and aux(x)

    def aux(z):
        print(type(z))
        if id(z) in d_aux:
            return d_aux[id(z)]
        if z.histories is None:
            d_aux[id(z)] = True
            return True
        else:
            for gx, gy in z.histories:
                for x, y in product(gx.parents, gy.parents):
                    if f(x, y):
                        d_aux[id(z)] = True
                        return True

        d_aux[id(z)] = False
        return False

    assert len(node.histories) == 1, "Genotype has more than one history"
    assert node.histories[0][0] == node.histories[0][0],\
        "different gametes are being generated"
    return any(aux(p) for p in node.histories[0][0].parents)


if __name__ == "__main__":
    n_loci = 5
    runner = BreedingProgram(n_loci, PlantSPC)
    runner.set_ideotype(PlantSPC(n_loci, (1 << n_loci) - 1, (1 << n_loci) - 1))
    for n_holes in range(1, (n_loci >> 1) + 1):
        runner.set_init_pop(
            PlantSPC.initial_pop_trait_introgression_homo(n_loci, n_holes, 1, 1)
        )
        runner.run()
        # print(runner._pop_current)
        assert achievable_by_parental_selection(runner, runner._pop_current[0]),\
            f"Parental selection cannot be used to achieve the target optimally for\n {runner._pop_0}"
