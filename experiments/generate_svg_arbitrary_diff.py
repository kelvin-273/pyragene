from random import seed
import svg
import eugene.plant_models.plant2 as ep2
import eugene.solvers.greedy_time as eg
import eugene.solvers.base_min_crossings_wedges as ew
import eugene.solvers.base_min_crossings_minizinc as em
import eugene.solvers.enumerators as ed

# seed(7)
seed(0)
n_loci = 16
n_pop = 5
# n_loci = 8
# n_pop = 4
SIDE_LENGTH = 40
STROKE_WIDTH = 2
STROKE_COLOR = "black"
ALLELE_0 = "#4444FF"
ALLELE_1 = "#FFFF89"
ALLELE_0_FADED = "#000044"
ALLELE_1_FADED = "#444424"
GAP_WITH = STROKE_WIDTH
PLOIDY = 2

# initialise random population
pop_0 = ep2.PlantSPC.initial_pop_random(n_loci, n_pop, p=0.3)

# extract gametes
gop_0 = []
for (i, x) in enumerate(pop_0):
    gop_0.extend([(i, x.gamete_specified(cp)) for cp in x.crosspoints()])
print(len(gop_0))

# extract segments
sop_0 = []
for (j, (i, gx)) in enumerate(gop_0):
    sop_0.extend([(j, c) for c in eg.segments_from_gamete(n_loci, gx)])

sop_0.sort(key=lambda jc: jc[1][0])
print(len(sop_0))

# sort and remove dominated gametes
gop_1 = sorted(gop_0, key=lambda igx: igx[1], reverse=True)


class G:
    def __init__(self, i, gamete):
        self.i = i
        self.gamete = gamete


gop_2 = [G(i, gx) for i, gx in gop_1]
gop_2 = ed.filter_non_dominating_gametes(gop_2)
gop_2 = [(g.i, g.gamete) for g in gop_2]
print(len(gop_2))

# extract segments
sop_1 = []
for (j, (i, gx)) in enumerate(gop_2):
    sop_1.extend([(j, c) for c in eg.segments_from_gamete(n_loci, gx)])

print(len(sop_1))
print(sop_1)

sop_final = eg.min_segment_cover_key(n_loci, sop_1, lambda jc: jc[1])
print(len(sop_final))

Allele = lambda xy, a: svg.Rect(
    x=xy[0],
    y=xy[1],
    stroke=STROKE_COLOR,
    stroke_width=STROKE_WIDTH,
    fill=ALLELE_1 if a else ALLELE_0,
    width=SIDE_LENGTH,
    height=SIDE_LENGTH,
)

AlleleFaded = lambda xy, a: svg.Rect(
    x=xy[0],
    y=xy[1],
    stroke=STROKE_COLOR,
    stroke_width=STROKE_WIDTH,
    fill=ALLELE_1_FADED if a else ALLELE_0_FADED,
    width=SIDE_LENGTH,
    height=SIDE_LENGTH,
)

Gamete = lambda xy, bits: svg.G(
    elements=[Allele((xy[0] + SIDE_LENGTH * i, xy[1]), a) for i, a in enumerate(bits)]
)

Genotype = lambda xy, chromosomes: svg.G(
    elements=[
        Gamete((xy[0], xy[1] + SIDE_LENGTH * i), chrom)
        for i, chrom in enumerate(chromosomes)
    ]
)

Segment = lambda xy, seg: svg.G(
    elements=[
        Allele((xy[0] + SIDE_LENGTH * i, xy[1]), a)
        if i in range(seg[0], seg[1]+1)
        else AlleleFaded((xy[0] + SIDE_LENGTH * i, xy[1]), a)
        for i, a in enumerate(seg[2])
    ]
)


def gamete_to_bits(n_loci, gx):
    out = [0] * n_loci
    i = n_loci - 1
    while i >= 0 and gx > 0:
        gx, out[i] = gx >> 1, gx & 1
        i -= 1
    return out


# draw svgs
with open("/tmp/pop_0.svg", "w") as f:
    print(
        svg.G(
            elements=[
                Genotype((0, 3 * SIDE_LENGTH * i), x.to_bitlist())
                for i, x in enumerate(pop_0)
            ],
        ),
        file=f,
    )

with open("/tmp/gop_0.svg", "w") as f:
    print(
        svg.G(
            elements=[
                Gamete((0, 2 * SIDE_LENGTH * i), gamete_to_bits(n_loci, gx))
                for i, (_, gx) in enumerate(gop_0)
            ],
        ),
        file=f,
    )

with open("/tmp/gop_1.svg", "w") as f:
    print(
        svg.G(
            elements=[
                Gamete((0, 2 * SIDE_LENGTH * i), gamete_to_bits(n_loci, gx))
                for i, (_, gx) in enumerate(gop_1)
            ],
        ),
        file=f,
    )


with open("/tmp/gop_2.svg", "w") as f:
    print(
        svg.G(
            elements=[
                Gamete((0, 2 * SIDE_LENGTH * i), gamete_to_bits(n_loci, gx))
                for i, (_, gx) in enumerate(gop_2)
            ],
        ),
        file=f,
    )


with open("/tmp/sop_1.svg", "w") as f:
    print(
        svg.G(
            elements=[
                Segment((0, 2 * SIDE_LENGTH * i), (s, e, gamete_to_bits(n_loci, gx)))
                for i, (_, (s, e, gx)) in enumerate(sop_1)
            ],
        ),
        file=f,
    )

with open("/tmp/sop_final.svg", "w") as f:
    print(
        svg.G(
            elements=[
                Segment((0, 2 * SIDE_LENGTH * i), (s, e, gamete_to_bits(n_loci, gx)))
                for i, (_, (s, e, gx)) in enumerate(sop_final)
            ],
        ),
        file=f,
    )
