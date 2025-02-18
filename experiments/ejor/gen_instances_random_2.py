"""
This script reads individual loci from stdin and outputs a set of instances
(passed by the -i option)
"""

import argparse
import random
import json
import sys
from eugene.plant_models.plant2 import PlantSPCBitarray

# Create the parser
parser = argparse.ArgumentParser(description="Process some integers.")

# Add arguments
parser.add_argument('-i', '--n_instances', type=int, required=True, help='Number of instances')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
n_instances = args.n_instances

random.seed(1)

for line in sys.stdin:
    n_loci = int(line)
    for n_pop in [2, 4, 6, 8]:
        for i in range(n_instances):
            pop_0 = PlantSPCBitarray.initial_pop_random(n_loci, n_pop)
            print(json.dumps({
                "n_loci": n_loci,
                "n_pop": n_pop,
                "pop_0": repr(pop_0),
            }))
