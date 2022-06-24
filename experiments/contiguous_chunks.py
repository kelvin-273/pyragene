from itertools import groupby, tee
from collections import Counter, namedtuple
from eugene.simulators.enumerators import BreedingProgram
from eugene.plant_models.plant2 import PlantSPC

def count_contiguous_chunks(plant: PlantSPC) -> int:
    """
    Counts the number of contiguous chunks in the given genotype.
    Assumes that the input plant is homozygous.

    Example:
    11110001101110011 -> 7 contiguous chunks
    11110001101110011

    >>> f = lambda n: count_contiguous_chunks(PlantSPC(8, n, n))
    >>> f(0)
    1
    >>> f(1)
    2
    >>> f(2)
    3
    >>> f(3)
    2
    >>> f(5)
    4
    >>> f(251)
    3
    >>> f(255)
    1
    """
    if not plant.chrom1 == plant.chrom2:
        raise ValueError("Plant must be homozygous")

    x = plant.chrom1
    out = 1
    for i in range(plant.n_loci-1):
        out += (x & 1) ^ ((x >> 1) & 1)
        x >>= 1
    return out


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=False)
    # 
    n_loci = 10
    data = []
    for n_holes in range((n_loci >> 1) + 1):
        for _ in range(1000):
            runner = BreedingProgram(n_loci, PlantSPC)
            pop_0 = PlantSPC.initial_pop_trait_introgression(n_loci, n_holes)
            # make the population homozygous
            pop_0 = [PlantSPC(
                n_loci,
                x.chrom1 & x.chrom2,
                x.chrom1 & x.chrom2,
            ) for x in pop_0]
            n_chunks = count_contiguous_chunks(pop_0[0])
            runner.set_init_pop(pop_0)
            runner.set_ideotype(PlantSPC(
                n_loci,
                (1 << n_loci) - 1,
                (1 << n_loci) - 1,
            ))
            runner.run()
            res = runner.get_results()
            print('\t'.join([
                f"n_loci: {n_loci}",
                f"n_holes: {n_holes}",
                f"n_chunks: {n_chunks}",
                f"n_generations: {res.n_generations}",
                f"pop_0: {[str(x) for x in pop_0]}",
            ]))
            data.append({
                "n_loci": n_loci,
                "n_holes": n_holes,
                "n_chunks": n_chunks,
                "n_generations": res.n_generations,
                "pop_0": [str(x) for x in pop_0],
            })

    import matplotlib.pyplot as plt
    # plot the number of generations against the number of chunks

    data.sort(key=lambda x: (x["n_holes"], x["n_generations"]))

    def f(accessor):
        for nh, column in groupby(data, key=lambda x: x[accessor]):
            l, r = tee(column)
            total = len(list(l))
            for t, xs in groupby(r, key=lambda x: x["n_generations"]):
                yield (nh, t, len(list(xs)) / total)

    data_holes = list(f("n_holes"))

    plt.scatter(
        [nh for nh, _, _ in data_holes],
        [t for _, t, _ in data_holes]
    )
    for nh, t, c in data_holes:
        plt.annotate(c, (nh, t))
    plt.title("distribution of lower-bound times to introgress n holes")
    plt.ylabel("number of generations")
    plt.xlabel("number of holes")
    plt.show()

    # plot the number generations against the number of chunks
    data_chunks = list(f("n_chunks"))
    plt.scatter(
        [nc for nc, _, _ in data_chunks],
        [t for _, t, _ in data_chunks]
    )
    for nc, t, c in data_chunks:
        plt.annotate(c, (nc, t))
    plt.title("distribution of lower-bound times to combine n non-overlapping chunks")
    plt.ylabel("number of generations")
    plt.xlabel("number of chunks")
    plt.show()
