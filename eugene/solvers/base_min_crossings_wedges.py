from typing import List


def max_repeated_wedges(dist_arr: List[int]) -> (int, List[bool]):
    n_loci = len(dist_arr)
    if n_loci <= 3:
        return 0

    # construct graph
    n_diff = n_loci - 1
    classlist = [
        tuple(sorted((dist_arr[i], dist_arr[i + 1]))) for i in range(n_diff)
    ]
    adjacent = [
        [classlist[i] == classlist[j] and i != j for j in range(n_diff)]
        for i in range(n_diff)
    ]

    xs = [False] * n_diff

    def aux(i, prev_skip=False, obj_curr=0):
        if i >= n_diff:
            # list is fully allocated
            return obj_curr, xs.copy()
        else:
            out = obj_curr
            out_xs = None

            # choose current and continue
            xs[i] = True
            repeated = any(adjacent[i][j] and xs[j] for j in range(i))
            # print(i, repeated, sep='\t')
            obj_new = obj_curr + repeated
            res, res_xs = aux(i + 2, prev_skip=False, obj_curr=obj_new)
            if res > out:
                out = res
                out_xs = res_xs

            # remove current from choice
            xs[i] = False
            # skip
            if not prev_skip:
                res, res_xs = aux(i + 1, prev_skip=True, obj_curr=obj_curr)
                if res > out:
                    out = res
                    out_xs = res_xs
            # print(f"return {out}")
            return out, out_xs

    return aux(0)


def min_crossings_wedges_heuristic(dist_arr: List[int]) -> int:
    res, res_xs = max_repeated_wedges(dist_arr)
    return len(dist_arr) - res


def min_crossings_wedges_heuristic_w_solution(dist_arr: List[int]):
    pass


def min_crossings_wedges_heuristic_w_bounded_objective(dist_arr: List[int], bound: int):
    return min(min_crossings_wedges_heuristic(dist_arr), bound)
