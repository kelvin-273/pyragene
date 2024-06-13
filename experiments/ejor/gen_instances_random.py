import argparse
import random
import json
from eugene.plant_models.plant2 import PlantSPCBitarray

# Create the parser
parser = argparse.ArgumentParser(description="Process some integers.")

# Add arguments
parser.add_argument('-l', '--n_loci_max', type=int, required=True, help='Maximum number of loci')
parser.add_argument('-i', '--n_instances', type=int, required=True, help='Number of instances')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
n_loci_max = args.n_loci_max
n_instances = args.n_instances

random.seed(0)

for n_loci in range(2, n_loci_max+1):
    for n_pop in [2, 4, 6, 8]:
        for _ in range(n_instances):
            pop_0 = PlantSPCBitarray.initial_pop_random(n_loci, n_pop)
            print(json.dumps({
                "n_loci": n_loci,
                "n_pop": n_pop,
                "pop_0": repr(pop_0),
            }))
