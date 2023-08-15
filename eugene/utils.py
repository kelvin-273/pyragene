from random import randint, randrange, sample, choice
from typing import List, Optional, NewType, Callable, Tuple
from collections import namedtuple


def count_ones(x: int):
    """
    Returns the number of one bits in the integer
    """
    if x < 0:
        raise ValueError(
            "Gave negative int and don't know how to deal with this yet"
        )
    out = 0
    while x > 0:
        x, r = x >> 1, x & 1
        out += r
    return out


def random_distribute_array(n_loci: int) -> List[int]:
    out = [0] * n_loci
    x_max = 0
    x_pre = 0
    for i in range(1, n_loci):
        x = randint(0, x_max)
        if x >= x_pre:
            x += 1
        x_max = max(x_max, x)
        out[i] = x
        x_pre = x
    return out


def max_repeated_wedges(dist_arr: List[int]) -> int:
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
    # __import__('pprint').pprint([[int(x) for x in row] for row in adjacent])

    xs = [False] * n_diff

    def aux(i, prev_skip=False, obj_curr=0):
        if i >= n_diff:
            # list is fully allocated
            return obj_curr
        else:
            out = obj_curr

            # choose current and continue
            xs[i] = True
            repeated = any(adjacent[i][j] and xs[j] for j in range(i))
            # print(i, repeated, sep='\t')
            obj_new = obj_curr + repeated
            res = aux(i + 2, prev_skip=False, obj_curr=obj_new)
            out = max(out, res)

            # remove current from choice
            xs[i] = False
            # skip
            if not prev_skip:
                res = aux(i + 1, prev_skip=True, obj_curr=obj_curr)
                out = max(out, res)
            # print(f"return {out}")
            return out

    return aux(0)


def test_max_wedges():
    case = [0, 1, 2, 1, 2, 1, 2, 0, 2, 1]
    check = 3
    guess = max_repeated_wedges(case)
    assert guess == check, f"guess({guess}) != check({check})"


if __name__ == "__main__":
    test_max_wedges()
