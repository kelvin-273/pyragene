from eugene.simulators.enumerators import BreedingProgram
from eugene.plant_models.plant2 import PlantSPC

def count_contiguous_chunks(plant: PlantSPC) -> int:
    """
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
        for _ in range(100):
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
    plt.scatter(
        [x["n_holes"] for x in data],
        [x["n_generations"] for x in data]
    )
    plt.title("number of generations to introgress n holes")
    plt.ylabel("number of generations")
    plt.xlabel("number of holes")
    plt.show()

    plt.scatter(
        [x["n_chunks"] for x in data],
        [x["n_generations"] for x in data]
    )
    plt.title("number of generations to introgress n chunks")
    plt.ylabel("number of generations")
    plt.xlabel("number of chunks")
    plt.show()
