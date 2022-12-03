import random
from eugene.simulators.greedy_multichrom_ti import BreedingProgramMultiChromTI
import eugene.simulators.greedy_multichrom_ti as mcti
from eugene.utils import count_ones


def random_instance(n_loci):
    """
    Create a diploid multi-chrom TI instance.

    This returns a list of two complimentary homozygous genotypes where the
    first genotype has at least as many favourable alleles as the second loci.
    """
    chrom1 = [random.randrange(1 << n_loci_i) for n_loci_i in n_loci]
    chrom2 = [c ^ ((1 << n_loci_i) - 1) for c, n_loci_i in zip(chrom1, n_loci)]

    if sum(count_ones(c) for c in chrom1) < sum(count_ones(c) for c in chrom2):
        chrom1, chrom2 = chrom2, chrom1

    return [
        mcti.PlantMCh(n_loci, chrom1, chrom1),
        mcti.PlantMCh(n_loci, chrom2, chrom2),
    ]


print(random_instance([2, 3]))
