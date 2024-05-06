from typing import List
from eugene.solution import BaseSolution
from eugene.plant_models.plant2 import PlantSPC
from eugene.utils import distribute_to_plants
import gurobipy as gp


def breeding_program_distribute(dist_array: List[int]) -> BaseSolution:
    return breeding_program(len(dist_array), distribute_to_plants(dist_array))


def breeding_program(n_loci: int, pop_0: List[PlantSPC]) -> BaseSolution:
    m = gp.Model()
    T = upper_bound(n_loci, pop_0)

    # Binary variables denoting whether genotype (gx, gy) is AVAILABLE in
    # generation t
    # u = m.addMVar((1 << n_loci, 1 << n_loci, T + 1), vtype="B")
    # u = [
    #     [
    #         [m.addVar(vtype="B") for t in range(T + 1)]
    #         for gy in range(1 << n_loci)
    #     ]
    #     for gx in range(1 << n_loci)
    # ]
    u = [
        [
            [m.addVar(vtype="B") if gy <= gx else None for t in range(T + 1)]
            for gy in range(1 << n_loci)
        ]
        for gx in range(1 << n_loci)
    ]

    # Binary variables denoting whether gamete gx is AVAILABLE in generation t
    # v = m.addMVar((1 << n_loci, T), vtype="B")
    v = [[m.addVar(vtype="B") for t in range(T)] for gx in range(1 << n_loci)]

    # # Binary variables denoting whether gametes gx and gy are used to CREATE
    # # a new genotype in generation t
    # # w = m.addMVar((1 << n_loci, 1 << n_loci, T), vtype="B")
    # w = [
    #     [[m.addVar(vtype="B") for t in range(T)] for gy in range(1 << n_loci)]
    #     for gx in range(1 << n_loci)
    # ]

    # Minimise the number of crossings
    # m.setObjective(sum(xt for r in w for c in r for xt in c), gp.GRB.MINIMIZE)
    # m.setObjective(
    #     sum(
    #         u[gx][gy][T] - u[gx][gy][0]
    #         for gx in range(1 << n_loci)
    #         for gy in range(1 << n_loci)
    #     ),
    #     gp.GRB.MINIMIZE,
    # )
    m.setObjective(
        sum(
            u[gx][gy][T] - u[gx][gy][0]
            for gx in range(1 << n_loci)
            for gy in range(gx + 1)
        ),
        gp.GRB.MINIMIZE,
    )

    # A genotype is available in generation 0 if and only if it is
    # in the initial population
    gen_0 = [[0] * (1 << n_loci) for _ in range(1 << n_loci)]
    for x in pop_0:
        gen_0[x.chrom1][x.chrom2] = 1
        gen_0[x.chrom2][x.chrom1] = 1
    for gx in range(1 << n_loci):
        # for gy in range(1 << n_loci):
        for gy in range(gx + 1):
            m.addConstr(u[gx][gy][0] == gen_0[gx][gy])

    # Target x* available in generation T+1
    m.addConstr(u[(1 << n_loci) - 1][(1 << n_loci) - 1][T] == 1)

    # If a genotype x is available in generation t-1 then it is available
    # in generation t
    for gx in range(1 << n_loci):
        for gy in range(gx + 1):
            for t in range(T):
                m.addConstr(u[gx][gy][t] <= u[gx][gy][t + 1])

    # # A genotype can only be available if it is created
    # for gx in range(1 << n_loci):
    #     for gy in range(1 << n_loci):
    #         for t in range(T):
    #             m.addConstr(w[gx][gy][t] == u[gx][gy][t + 1] - u[gx][gy][t])

    # # If a gamete gx is available in generation t-1 then it is available
    # # in generation t
    # for gx in range(1 << n_loci):
    #     for t in range(T - 1):
    #         m.addConstr(v[gx][t] <= v[gx][t + 1])
    #         # NOTE: This should be redundant given the next constraints

    # gx is available if and only if there is some x that can produce gx
    # v[gx][t] <-> exists(u[x] : gx in x)
    # ≡ -v[gx][t] <-> forall(-x : gx in x)

    # -v[gx][t] -> forall(-u[x] : gx in x)
    # ≡ v[gx][t] \/ forall(-u[x] : gx in x)
    # ≡ forall(v[gx][t] \/ -u[x] : gx in x)
    for gx in range(1 << n_loci):
        # for gy in range(1 << n_loci):
        for gy in range(gx + 1):
            for gz in gen_recombined_gametes(n_loci, gx, gy):
                for t in range(T):
                    m.addConstr(v[gz][t] >= u[gx][gy][t])

    # v[gx][t] -> exists(u[x] : gx in x)
    # ≡ -v[gx][t] \/ exists(u[x] : gx in x)
    parent_gamete_masks = list(gen_parent_gametes_fast_masks(n_loci))
    for gz in range(1 << n_loci):
        parent_gametes = parent_gametes_from_masks(
            n_loci, gz, parent_gamete_masks
        )
        for t in range(T):
            m.addConstr(
                (1 - v[gz][t])
                + sum(u[gx][gy][t] for gx, gy in parent_gametes if gy <= gx)
                >= 1
            )

    # The gametes required to create x must be available
    # in order for x to be created
    for gx in range(1 << n_loci):
        for gy in range(gx + 1):
            for t in range(T):
                m.addConstr(
                    2 * u[gx][gy][t + 1] - 2 * u[gx][gy][t]
                    <= v[gx][t] + v[gy][t]
                )
                # m.addConstr(w[gx][gy][t] <= v[gx][t])
                # m.addConstr(w[gx][gy][t] <= v[gy][t])

    # # Symmetry breaking
    # for gx in range(1 << n_loci):
    #     for gy in range(1 << n_loci):
    #         for t in range(T-1):
    #             m.addConstr(v[gx][t] + v[gy][t] + u[gx][gy][t + 2] - u[gx][gy][t + 1] <= 2)

    m.optimize()

    # for t in range(T):
    #     # __import__('pprint').pprint(w[:, :, t].X)
    #     # __import__('pprint').pprint(v[:, t].X)
    #     __import__("pprint").pprint(
    #         [
    #             [w[gx][gy][t].X for gy in range(1 << n_loci)]
    #             for gx in range(1 << n_loci)
    #         ]
    #     )
    #     print()

    return BaseSolution(
        tree_data=[],
        tree_type=[],
        tree_left=[],
        tree_right=[],
        objective=m.objVal,
    )


def upper_bound(n_loci: int, pop_0: List[PlantSPC]) -> int:
    return n_loci


def is_child_gamete(n_loci: int, gx: int, gy: int, gz: int) -> bool:
    mask_0 = (1 << n_loci) - 1
    for k in range(n_loci):
        mask_l = (mask_0 >> k) << k
        mask_r = mask_0 ^ mask_l
        if ((gx & mask_l) | (gy & mask_r)) == gz:
            return True
        if ((gy & mask_l) | (gx & mask_r)) == gz:
            return True
    return False


def gen_parent_gametes(n_loci: int, gz: int):
    for gx in range(1 << n_loci):
        for gy in range(1 << n_loci):
            if is_child_gamete(n_loci, gx, gy, gz):
                yield (gx, gy)


def gen_parent_gametes_fast_masks(n_loci: int):
    mask_0 = (1 << n_loci) - 1

    # No targetG
    for len_tail in range(1, n_loci):
        len_head = n_loci - len_tail
        base_tail = (1 << len_tail) - 1
        base_head = mask_0 ^ base_tail
        for comp_head in range((1 << len_head) - 1):
            for comp_tail in range(1 << (len_tail - 1)):
                c1 = base_head | comp_tail
                c2 = (comp_head << len_tail) | base_tail
                yield c1, c2
                yield c2, c1

    # One targetG
    for comp in range((1 << n_loci) - 1):
        yield mask_0, comp
        yield comp, mask_0

    # Both targetG
    yield (mask_0, mask_0)


def parent_gametes_from_masks(
    n_loci: int, gz: int, parent_gamete_masks: List[int]
):
    mask_0 = (1 << n_loci) - 1
    return [
        (
            mask_0 ^ mask_l ^ gz,
            mask_0 ^ mask_r ^ gz
        )
        for mask_l, mask_r in parent_gamete_masks
    ]


def gen_recombined_gametes(n_loci, gx, gy):
    mask_0 = (1 << n_loci) - 1
    for k in range(n_loci):
        mask_l = (mask_0 >> k) << k
        mask_r = mask_0 ^ mask_l
        yield (gx & mask_l) | (gy & mask_r)
        yield (gy & mask_l) | (gx & mask_r)


if __name__ == "__main__":
    n_loci = 2
    print(
        breeding_program(
            n_loci,
            [PlantSPC(n_loci, 0b10, 0b10,), PlantSPC(n_loci, 0b01, 0b01,),],
        )
    )
    n_loci = 3
    print(
        breeding_program(
            n_loci,
            [PlantSPC(n_loci, 0b101, 0b101,), PlantSPC(n_loci, 0b010, 0b010,),],
        )
    )
    n_loci = 3
    print(
        breeding_program(
            n_loci,
            [
                PlantSPC(n_loci, 0b100, 0b100,),
                PlantSPC(n_loci, 0b010, 0b010,),
                PlantSPC(n_loci, 0b001, 0b001,),
            ],
        )
    )
    n_loci = 4
    print(
        breeding_program(
            n_loci,
            [
                PlantSPC(n_loci, 0b1010, 0b1010,),
                PlantSPC(n_loci, 0b0101, 0b0101,),
            ],
        )
    )
    n_loci = 4
    print(
        breeding_program(
            n_loci,
            [
                PlantSPC(n_loci, 0b1000, 0b1000,),
                PlantSPC(n_loci, 0b0101, 0b0101,),
                PlantSPC(n_loci, 0b0010, 0b0010,),
            ],
        )
    )
    n_loci = 5
    print(
        breeding_program(
            n_loci,
            [
                PlantSPC(n_loci, 0b10100, 0b10100,),
                PlantSPC(n_loci, 0b01001, 0b01001,),
                PlantSPC(n_loci, 0b00010, 0b00010,),
            ],
        )
    )

    from eugene.solvers.base_min_crossings_minizinc import (
        breeding_program as g_mzn,
    )
    from random import seed
    seed(0)

    g_mip = breeding_program
    # for n_loci in range(6, 7):
    for n_loci in range(2, 6):
        for _ in range(100):
            case = PlantSPC.initial_pop_random(n_loci, 3, p=0.2)
            res_mzn = g_mzn(n_loci, case)
            res_mip = g_mip(n_loci, case)
            assert res_mzn.objective == res_mip.objective
