from .abc import BreedingProgram
from ..plant_models.plant2 import PlantSPC


class BreedingProgramGreedyTime(BreedingProgram):

    """Breeding program solver that uses the favourable chunks of alleles in gametes to determine the set of crossovers that achieve the target in the fewest number of generations"""

    def __init__(self, n_loci, plant_type):
        """TODO: to be defined. """
        BreedingProgram.__init__(self, n_loci, plant_type)

    def run(self):
        # construct the set of chunks the initial population in (x, -y) order
        # choose the chunks to use by largest of extension
        # wrap chunks
        # construct dag
        pass


Gamete = int


def extract_chunks_from_gamete(n_loci, gamete) -> list:
    assert gamete < 1 << n_loci, f"gamete has more than {n_loci} loci {gamete}"
    out = []
    in_chunk = False
    s = e = None
    for i in range(n_loci):
        allele = (gamete >> i) & 1
        if allele:
            if not in_chunk:
                e = n_loci - i - 1
            s = n_loci - i - 1
        elif not allele and in_chunk:
            out.append((s, e, gamete))
        in_chunk = bool(allele)
    if in_chunk:
        out.append((s, e, gamete))
    return list(reversed(out))


def extract_chunks_from_genotype(n_loci, genotype: PlantSPC):
    q1 = extract_chunks_from_gamete(n_loci, genotype.chrom1)
    q2 = extract_chunks_from_gamete(n_loci, genotype.chrom2)

    used1 = [False] * len(q1)
    used2 = [False] * len(q1)

    i = j = 0
    out = []

    while i < len(q1) and j < len(q2):
        c1 = (s1, e1, g1) = q1[i]
        c2 = (s2, e2, g2) = q2[j]

        if s1 > s2 or s1 == s2 and e1 < e2:
            q1, q2 = q2, q1
            used1, used2 = used2, used1
            i, j = j, i
        elif s1 <= s2 and e1 >= e2:
            j += 1
        elif e1 + 1 < s2:
            assert s1 < e1 + 1 < s2 <= e2
            if not used1[i]:
                out.append(c1)
                used1[i] = True
            i += 1
        elif e2 == n_loci - 1:
            out.append((s1, e2, PlantSPC(n_loci, g1, g2).gamete_specified((0, e1 + 1))))
            used1[i] = used2[j] = True
            # break and don't use any more chunks
            i = len(q1)
            j = len(q2)
        else:
            assert e2 < n_loci - 1
            out.append((s1, e2, PlantSPC(n_loci, g1, g2).gamete_specified((0, e1 + 1))))
            used1[i] = used2[j] = True
            i += 1

    return out + q1[i:] + q2[j:]  # works because i == len(q1) or j == len(q2)
