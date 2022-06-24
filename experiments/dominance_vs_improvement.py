from itertools import count, product
from eugene.plant_models.plant2 import PlantSPC


def imp_weak(x, y):
    """Returns True if and only if x is a weak improvement over y.

    Defined by De Beukelaer et al. 2015, x is a weak improvement over y if
    (a) one of the chromosomes of x contains a desired stretch where there is
    not in either chromosome of y or (b) x contains a desired allele in a
    homozygous state which does not occur in y in a homozygous state.
    """
    n_loci = x.n_loci
    mask = (1 << n_loci) - 1
    clause1 = ((x.chrom1 | x.chrom2) & (y.chrom1 ^ mask) & (y.chrom2 ^ mask))
    clause2 = x.chrom1 & x.chrom2 & ((y.chrom1 ^ mask) | (y.chrom2 ^ mask))
    return clause1 | clause2 > 0


def bi_imp_weak(x, y):
    imp_weak(x, y) and imp_weak(y, x)


def non_dom(x, y):
    not (x.dom_weak(y) or y.dom_weak(y))


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=False)

    for n_loci in count(1):
        print(f"n_loci: {n_loci}")
        for (p1c1, p1c2), (p2c1, p2c2) in product(
            product(range(1 << n_loci), range(1 << n_loci)),
            product(range(1 << n_loci), range(1 << n_loci))
        ):
            p1 = PlantSPC(n_loci, p1c1, p1c2)
            p2 = PlantSPC(n_loci, p2c1, p2c2)
            # assert (p1.dom_weak(p2) and p1 != p2) == imp_weak(p1, p2), \
                # f"{p1}, {p2}, {p1.dom_weak(p2)}, {imp_weak(p1, p2)}"
            assert non_dom(p1, p2) == bi_imp_weak(p1, p2), \
                f"{p1}, {p2}, {(p1.dom_weak(p2), p2.dom_weak(p1))}, {(imp_weak(p1, p2), imp_weak(p2, p1))}"
