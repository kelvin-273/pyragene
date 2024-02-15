"""
In this experiment, we ask how different solvers perform on randomly generated
distribute instances.

We'll measure the expected the run-time for randomly generated distribute
instances. As the number of distribute instances for n_loci = 1..6 is less
than 100, we shall simply take the expected value over all loci. The number
instances where n_loci = 7..8 is less than 1000, so these can be chosen from
the full set. Instances where n_loci ≥ 9can sampled by generating and checking
against a list of previously observed instances.
"""

import eugene.utils as eu

SOLVERS = [
    ("CANZAR", None),
    ("CP-SAT", None),
    ("CP-MIP", None),
    ("A*", None),
]
