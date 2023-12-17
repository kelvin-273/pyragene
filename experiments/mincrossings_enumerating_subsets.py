"""
This was an attempt to enumerate subsets of segments.
We found that most covering subsets of segments even with bounding yield distribute instances
of the form [1, ..., n_loci].

A potential future idea could be to keep track of the distribute instances we've seen and avoid repeating such instances.
instances can be stored in a suffix tree
"""
from eugene.plant_models.plant2 import PlantSPC
from eugene.simulators.base_min_crossings import gen_covering_subsets
from eugene.simulators.greedy_time import segments_from_genotype
from eugene.simulators.base_min_crossings import breeding_program

n_loci = 50
pop0 = PlantSPC.initial_pop_random(n_loci, 6)
print(pop0)

segments = []
for x in pop0:
    segments.extend(segments_from_genotype(n_loci, x))

for i, x in enumerate(breeding_program(n_loci, pop0)):
    print(i)
