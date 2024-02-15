from math import ceil, log2
from typing import List
from eugene.solution import BaseSolution
from eugene.utils import sanitise_distribute_array


ONE_INDEX = 1


def max_repeated_wedges(dist_arr: List[int], strict_mingen=False, **kwargs) -> (int, List[bool]):
    """
    Returns a tuple containing the largest number `out_obj` of repeated wedges
    in any selection of wedges, and the selection `out_selection` that yields
    `out_obj` repeated wedges.
    """
    n_loci = len(dist_arr)
    if n_loci == 0:
        raise ValueError("Empty dist_arr")
    if n_loci == 1:
        return (0, [])
    if n_loci == 2:
        return (0, [True])
    if n_loci == 3:
        return (0, [True, False])

    # determine bound on the number of singles
    if strict_mingen:
        k = ceil(log2(n_loci))
        bound_singles = (1 << k) - n_loci
    else:
        bound_singles = float("inf")

    # construct graph
    n_diff = n_loci - 1
    classlist = [tuple(sorted((dist_arr[i], dist_arr[i + 1]))) for i in range(n_diff)]
    adjacent = [
        [classlist[i] == classlist[j] and i != j for j in range(n_diff)]
        for i in range(n_diff)
    ]

    selected = [False] * n_diff

    def aux(i, prev_skip=False, obj_curr=0, singles=0):
        if i >= n_diff:
            # list is fully allocated
            return obj_curr, selected.copy()
        else:
            out_obj = obj_curr
            out_selection = selected

            # choose current wedge and continue
            selected[i] = True
            repeated = any(selected[j] and adjacent[i][j] for j in range(i))
            obj_new = obj_curr + repeated
            res_obj, res_selection = aux(
                i + 2, prev_skip=False, obj_curr=obj_new, singles=singles
            )
            if res_obj > out_obj:
                out_obj = res_obj
                out_selection = res_selection

            # remove current wedge from choice
            selected[i] = False
            # skip
            if not prev_skip and (singles < bound_singles or not strict_mingen):
                res_obj, res_selection = aux(
                    i + 1, prev_skip=True, obj_curr=obj_curr, singles=singles + 1
                )
                if res_obj > out_obj:
                    out_obj = res_obj
                    out_selection = res_selection
            return 0 if out_selection is None else out_obj, out_selection

    return aux(0)


def default_wedges(dist_arr: List[int], *args, **kwargs) -> (int, List[bool]):
    n_loci = len(dist_arr)
    out_selection = [i % 2 == 0 for i in range(n_loci - 1)]
    n_repeats = sum(
        any(
            (dist_arr[i] == dist_arr[j] and dist_arr[i + 1] == dist_arr[j + 1])
            or (dist_arr[i] == dist_arr[j + 1] and dist_arr[i + 1] == dist_arr[j])
            for j in range(0, i, 2)
        )
        for i in range(2, n_loci - 1, 2)
    )
    return n_repeats, out_selection


def breeding_program(
    n_loci: int, dist_arr: List[int],
    wedges_selector=None
) -> BaseSolution:
    #####################################################
    #         Heap             #  Wedges #      P_0     #
    #####################################################
    n_loci = len(dist_arr)  # this makes n_loci redundant
    dist_arr = sanitise_distribute_array(dist_arr)
    if wedges_selector is None:
        wedges_selector = lambda dist_arr: max_repeated_wedges(dist_arr, strict_mingen=False)
    n_repeats, wedges = wedges_selector(dist_arr)
    return _generate_solution_from_selection(n_loci, dist_arr, n_repeats, wedges)


def _generate_solution_from_selection(
    n_loci: int, dist_arr: List[int], n_repeats: int, wedges: List[bool]
):
    # which wedges are repeated
    n_pop = max(dist_arr) + 1

    n_cross = n_loci - n_repeats
    n_nodes = n_cross + n_pop
    tree_data = [[[0] * n_loci for _ in range(2)] for _ in range(n_nodes)]
    tree_type = ["Node"] * n_cross + ["Leaf"] * n_pop
    tree_left = [0] * n_nodes
    tree_right = [0] * n_nodes
    objective = n_cross

    # populate tree_data for initial population
    for j in range(n_loci):
        i = dist_arr[j]
        i_data = n_cross + i
        tree_data[i_data][0][j] = 1
        tree_data[i_data][1][j] = 1

    #####################
    #  graph structure  #
    #####################

    n_gen1_crossings = sum(wedges) - n_repeats

    q = []

    wedge_map = {}
    wedge_key = lambda l, r: l * n_pop + r if l <= r else r * n_pop + l

    i_selection = 0
    i_wedge = n_cross - n_gen1_crossings

    # Add each wedge and single to the queue
    while i_selection < n_loci - 1:
        if wedges[i_selection]:

            i_l = dist_arr[i_selection] + n_cross
            i_r = dist_arr[i_selection + 1] + n_cross

            if wedge_key(i_l, i_r) not in wedge_map:
                wedge_map[wedge_key(i_l, i_r)] = i_wedge
                tree_left[i_wedge] = i_l + ONE_INDEX
                tree_right[i_wedge] = i_r + ONE_INDEX
                tree_data[i_wedge] = [tree_data[i_l][0], tree_data[i_r][0]]
                i_wedge += 1

            s = i_selection
            e = i_selection + 2  # endpoints mirroring python endpoints
            g = (
                tree_data[i_l][0][: i_selection + 1]
                + tree_data[i_r][0][i_selection + 1 :]
            )
            assert all(g[j] in (1, True) for j in range(s, e))
            q.append(((s, e, g), wedge_map[wedge_key(i_l, i_r)]))
            i_selection += 2

        else:
            s = i_selection
            e = i_selection + 1
            g = tree_data[n_cross + dist_arr[i_selection]][0]
            assert all(g[j] in (1, True) for j in range(s, e))
            q.append(((s, e, g), n_cross + dist_arr[i_selection]))
            i_selection += 1

    if i_selection == n_loci - 1:
        s = i_selection
        e = i_selection + 1
        g = tree_data[n_cross + dist_arr[i_selection]][0]
        assert all(g[j] in (1, True) for j in range(s, e))
        q.append(((s, e, g), n_cross + dist_arr[i_selection]))
        i_selection += 1

    i_node = n_cross - 1 - n_gen1_crossings
    while i_node > 0:
        assert len(q) > 1
        assert len(q) == i_node + 1
        q_new = []
        while len(q) > 1:
            pack_r = q.pop()
            pack_l = q.pop()
            if pack_r[0][0] < pack_l[0][0]:
                pack_l, pack_r = pack_r, pack_l
            (s_l, e_l, g_l), i_l = pack_l
            (s_r, e_r, g_r), i_r = pack_r

            assert s_l < s_r <= e_l < e_r
            tree_left[i_node] = i_l + ONE_INDEX
            tree_right[i_node] = i_r + ONE_INDEX
            tree_data[i_node] = [g_l, g_r]

            g = g_l[:e_l] + g_r[e_l:]

            assert all(g[j] in (1, True) for j in range(s_l, e_r))
            q_new.append(((s_l, e_r, g), i_node))
            i_node -= 1
        q = q + list(reversed(q_new))

    tree_left[0] = tree_right[0] = 1 + ONE_INDEX
    tree_data[0] = [[1] * n_loci] * 2

    return BaseSolution(tree_data, tree_type, tree_left, tree_right, objective,)


def min_crossings_wedges_heuristic(dist_arr: List[int]) -> int:
    res, res_xs = max_repeated_wedges(dist_arr)
    return len(dist_arr) - res


def min_crossings_wedges_heuristic_w_solution(dist_arr: List[int]):
    pass


def min_crossings_wedges_heuristic_w_bounded_objective(dist_arr: List[int], bound: int):
    return min(min_crossings_wedges_heuristic(dist_arr), bound)
