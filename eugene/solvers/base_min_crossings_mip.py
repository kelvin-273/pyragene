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
    u = [
        [
            [m.addVar(vtype="B") for t in range(T + 1)]
            for gy in range(1 << n_loci)
        ]
        for gx in range(1 << n_loci)
    ]

    # Binary variables denoting whether gamete gx is AVAILABLE in generation t
    # v = m.addMVar((1 << n_loci, T), vtype="B")
    v = [[m.addVar(vtype="B") for t in range(T)] for gx in range(1 << n_loci)]

    # Binary variables denoting whether gametes gx and gy are used to CREATE
    # a new genotype in generation t
    # w = m.addMVar((1 << n_loci, 1 << n_loci, T), vtype="B")
    w = [
        [[m.addVar(vtype="B") for t in range(T)] for gy in range(1 << n_loci)]
        for gx in range(1 << n_loci)
    ]

    # Minimise the number of crossings
    m.setObjective(sum(xt for r in w for c in r for xt in c), gp.GRB.MINIMIZE)

    # A genotype is available in generation 0 if and only if it is
    # in the initial population
    gen_0 = [[0] * (1 << n_loci) for _ in range(1 << n_loci)]
    for x in pop_0:
        gen_0[x.chrom1][x.chrom2] = 1
    for gx in range(1 << n_loci):
        for gy in range(1 << n_loci):
            m.addConstr(u[gx][gy][0] == gen_0[gx][gy])

    # Target x* available in generation T+1
    m.addConstr(u[(1 << n_loci) - 1][(1 << n_loci) - 1][T] == 1)

    # # If a genotype x is available in generation t-1 then it is available
    # # in generation t
    # for gx in range(1 << n_loci):
    #     for gy in range(1 << n_loci):
    #         for t in range(T):
    #             m.addConstr(u[gx][gy][t] <= u[gx][gy][t + 1])

    # A genotype can only be available if it is created
    for gx in range(1 << n_loci):
        for gy in range(1 << n_loci):
            for t in range(T):
                m.addConstr(w[gx][gy][t] == u[gx][gy][t + 1] - u[gx][gy][t])

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
        for gy in range(1 << n_loci):
            for t in range(T):
                for gz in gen_recombined_gametes(n_loci, gx, gy):
                    m.addConstr(v[gz][t] >= u[gx][gy][t])

    # v[gx][t] -> exists(u[x] : gx in x)
    # ≡ -v[gx][t] \/ exists(u[x] : gx in x)
    for gz in range(1 << n_loci):
        for t in range(T):
            pass
            m.addConstr(
                (1 - v[gz][t])
                + sum(u[gx][gy][t] for gx, gy in gen_parent_gametes(n_loci, gz))
                >= 1
            )

    # The gametes required to create x must be available
    # in order for x to be created
    for gx in range(1 << n_loci):
        for gy in range(1 << n_loci):
            for t in range(T):
                m.addConstr(w[gx][gy][t] <= v[gx][t])
                m.addConstr(w[gx][gy][t] <= v[gy][t])

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


def gen_parent_gametes_fast(gz: int):
    # Far out, this is gonna be a backtracking algorithm
    gx_arr = [0] * n_loci
    gy_arr = [0] * n_loci
    gz_arr = [0] * n_loci
    gz_tmp = gz
    for i in reversed(n_loci):
        b = gz_tmp & 1
        gz_tmp >> 1
        gz_arr[i] = b

    # At the start of an iteration at locus p, k_state indicates which parent
    # is responsible for the prefix up to locus p-1.
    # k_state = 0 -> both parents could be responsible for the prefix
    # k_state = 1 -> only gx could be responsible for the prefix
    # k_state = 2 -> only gy could be responsible for the prefix
    #
    # Furthermore, if k != 0 then the suffix is set by the other parent
    k_state = 0
    k_point = 0

    def intify(arr: List[int]) -> int:
        """
        Returns the integer representation of the bit.
        Assumes every element in the bit array is either 0 or 1.
        """
        out = 0
        for b in arr:
            out <<= 1
            out &= b
        return out

    def aux(p: int):
        if p == n_loci:
            yield (intify(gx_arr), intify(gy_arr))
        else:
            if k_state == 0:
                if gz_arr[p] == 0:
                    pass

                pass

                pass

                if gz_arr[p] == 1:
                    pass

            elif k_state == 1:
                gy_arr[p] = 0
                yield from aux(p)
                gy_arr[p] = 1
                yield from aux(p)

            elif k_state == 2:
                gx_arr[p] = 0
                yield from aux(p)
                gx_arr[p] = 1
                yield from aux(p)

            else:
                raise ValueError("k_state takes an undefined value")


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
            [PlantSPC(n_loci, 0b1010, 0b1010,), PlantSPC(n_loci, 0b1010, 0b1010,),],
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
