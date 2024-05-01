from typing import List
import eugene.utils as eu
import eugene.solvers.base_min_crossings_minizinc as em
from eugene.solution import BaseSolution


def breeding_program_distribute(
    n_loci: int, dist_array: List[int], sub_solver=None
) -> BaseSolution:
    dist_array_simple = dist_array.copy()
    dist_array_extend = eu.sanitise_distribute_array(dist_array)

    n_pop_orig = max(dist_array_extend) + 1

    tree_data = em.instance_array_genotype_homo(dist_array, max_crossovers=1)[
        "genotypes"
    ]
    n_leaves = len(tree_data)
    tree_type = ["Leaf"] * n_leaves
    tree_left = [0] * n_leaves
    tree_right = [0] * n_leaves
    objective = 0

    parent_index = [i + 1 for i in range(n_pop_orig)]

    while len(dist_array_simple) > 1:
        dist_array_simple = eu.sanitise_distribute_array(dist_array_simple)
        ranges = eu.distribute_to_ranges(dist_array_simple)

        # find join
        full_join_exists = False
        for i in range(len(ranges)):
            _, e = ranges[i]
            for j in range(i + 1, len(ranges)):
                s, _ = ranges[j]
                if e + 1 == s:
                    full_join_exists = True
                    break
            if e + 1 == s:
                full_join_exists = True
                break

        if not full_join_exists:
            break

        # make child
        child = [
            [1 if x == i else 0 for x in dist_array_extend],
            [1 if x == j else 0 for x in dist_array_extend],
        ]
        tree_data.append(child)
        tree_type.append("Node")
        tree_left.append(parent_index[i])
        tree_right.append(parent_index[j])
        objective += 1

        # replace the j's with i's in the dist_array_extend
        for k, x in enumerate(dist_array_extend):
            if dist_array_extend[k] == j:
                dist_array_extend[k] = i
            elif dist_array_extend[k] > j:
                dist_array_extend[k] -= 1
        # replace the j's with i's in the dist_array_simple
        for k in range(ranges[j][0], ranges[j][1] + 1):
            if dist_array_simple[k] == j:
                dist_array_simple[k] = i
        dist_array_simple.pop(s)

        parent_index[i] = len(tree_data)
        parent_index.pop(j)

    # map loci in simplified to segments in extended
    segments = [None] * len(dist_array_simple)
    big_to_sub = [0] * len(dist_array_extend)
    i = e = 0
    while i < len(dist_array_simple) and e < len(dist_array_extend):
        assert dist_array_simple[i] == dist_array_extend[e]
        s = e
        while (
            e < len(dist_array_extend) and dist_array_extend[e] == dist_array_simple[i]
        ):
            big_to_sub[e] = i
            e += 1
        segments[i] = (s, e)
        i += 1

    # solve subproblem
    if sub_solver is None:
        sub_solver = em.breeding_program_distribute
    sub_solution = sub_solver(len(dist_array_simple), dist_array_simple)

    n_crossings = sub_solution.crossings()
    n_plants = sub_solution.n_plants

    # translate solution back to original problem
    sub_sol_rev = sub_solution.permute_base_solution(
        list(reversed(range(n_crossings))) + list(range(n_crossings, n_plants))
    )
    old_tree_data_len = len(tree_data)
    tree_data.extend(
        [
            [
                [chrom[big_to_sub[j]] for j in range(len(dist_array_extend))]
                for chrom in x
            ]
            for x in sub_sol_rev.tree_data[:n_crossings]
        ]
    )
    tree_type.extend(["Node"] * n_crossings)
    objective += n_crossings
    tree_left.extend(
        [
            # None
            # if sub_sol_rev.tree_left[i] > sub_sol_rev.crossings()
            # else
            old_tree_data_len + sub_sol_rev.tree_left[i]
            if sub_sol_rev.tree_left[i] > 0
            else 0
            for i in range(sub_sol_rev.crossings())
        ]
    )
    tree_right.extend(
        [
            # None
            # if sub_sol_rev.tree_right[i] > sub_sol_rev.crossings()
            # else
            old_tree_data_len + sub_sol_rev.tree_right[i]
            if sub_sol_rev.tree_right[i] > 0
            else 0
            for i in range(sub_sol_rev.crossings())
        ]
    )

    for i in range(sub_sol_rev.crossings()):
        if (
            tree_left[old_tree_data_len + i]
            > old_tree_data_len + sub_sol_rev.crossings()
        ):
            data_left = sub_sol_rev.tree_data[sub_sol_rev.tree_left[i] - 1]
            j = 0
            while j < len(dist_array_simple) and data_left[0][j] == 0:
                j += 1

            gx_in_sub_dist_array = dist_array_simple[j]

            p_idx = parent_index[gx_in_sub_dist_array]
            tree_left[old_tree_data_len + i] = p_idx

        if (
            tree_right[old_tree_data_len + i]
            > old_tree_data_len + sub_sol_rev.crossings()
        ):
            data_right = sub_sol_rev.tree_data[sub_sol_rev.tree_right[i] - 1]
            j = 0
            while j < len(dist_array_simple) and data_right[0][j] == 0:
                j += 1

            gx_in_sub_dist_array = dist_array_simple[j]

            p_idx = parent_index[gx_in_sub_dist_array]
            tree_right[old_tree_data_len + i] = p_idx

    for i in range(sub_sol_rev.crossings(), sub_sol_rev.n_plants):
        assert sub_sol_rev.tree_type[i] == "Leaf"
        assert sub_sol_rev.tree_data[i] == sub_solution.tree_data[i]
        assert sub_sol_rev.tree_data[i][0] == sub_sol_rev.tree_data[i][1]

        j = 0
        while j < len(dist_array_simple) and sub_sol_rev.tree_data[i][0][j] == 0:
            j += 1

        gx_in_sub_dist_array = dist_array_simple[j]
        for j in range(j + 1, len(dist_array_simple)):
            if (dist_array_simple[j] == gx_in_sub_dist_array) != bool(
                sub_sol_rev.tree_data[i][0][j]
            ):
                raise ValueError("solution Leaf and dist_array_simple are incongruent")

    return BaseSolution(
        tree_data=tree_data,
        tree_type=tree_type,
        tree_left=tree_left,
        tree_right=tree_right,
        objective=objective,
    ).permute_base_solution(
        list(reversed(range(len(tree_data) - objective, len(tree_data),)))
        + list(range(len(tree_data) - objective)),
    )


if __name__ == "__main__":
    print(breeding_program_distribute(2, [0, 1]))
    print(breeding_program_distribute(3, [0, 1, 0]))
    print(breeding_program_distribute(3, [0, 1, 2]))
    print(breeding_program_distribute(4, [0, 1, 0, 2]))
    print(breeding_program_distribute(6, [0, 1, 2, 3, 0, 2]))
    print(breeding_program_distribute(10, [0, 1, 0, 2, 0, 3, 4, 5, 6, 5]))
