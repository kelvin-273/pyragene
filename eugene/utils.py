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


def isolated_subproblems_ranges(xs: List[Tuple[int, int]]):
    assert len(xs) > 0
    assert all(s <= e for s, e in xs)
    xs = sorted(xs)
    s_curr, e_curr = xs[0]
    for i in range(1, len(xs)):
        s, e = xs[i]
        pass


def gen_distribute_instances(n_loci):
    state = [0] * n_loci

    def aux(i=1, max_val=0):
        if i == n_loci:
            yield state.copy()
        else:
            for val_i in range(max_val + 1):
                if val_i != state[i - 1]:
                    state[i] = val_i
                    yield from aux(i + 1, max_val)
            state[i] = max_val + 1
            yield from aux(i + 1, max_val + 1)

    yield from aux()


def random_distribute_instance(n_loci):
    state = [0] * n_loci
    max_val = 0
    for i in range(1, n_loci):
        next_val = randint(0, max_val)
        if next_val >= state[i-1]:
            next_val += 1
        state[i] = next_val
        max_val = max(max_val, next_val)
    return state


def distribute_instance_to_ranges(instance_array: List[int]) -> List[Tuple[int, int]]:
    """
    Given a distribute instance, returns an array of ranges
    [(s_0, e_0), ..., (s_{m-1}, e_{m-1})] where m is the number of genotypes in
    the initial population.
    """
    d = {}
    for i, x in enumerate(instance_array):
        if x in d:
            s, _ = d[x]
            d[x] = (s, i)
        else:
            d[x] = (i, i)
    return list(d.values())


def distribute_instance_to_index_lists(instance_array: List[int]):
    d = {}
    for i, x in enumerate(instance_array):
        if x in d:
            d[x].append(i)
        else:
            d[x] = [i]
    return list(d.values())


def print_lines(lines):
    assert all(s <= e for s, e in lines)
    n = max(e for s, e in lines) + 1
    lines = sorted(lines)
    for s, e in lines:
        print(" " * s + "-" * (e - s + 1) + " " * (n - e))


def print_lines_with_markings(instance_array):
    indexes = distribute_instance_to_index_lists(instance_array)
    n = max(gen[-1] for gen in indexes) + 1
    indexes.sort()
    for gen in indexes:
        s = gen[0]
        e = gen[-1]
        string = " " * s
        string += (
            "x"
            if len(gen) == 1
            else (
                "x"
                + "x".join(["-" * (b - a - 1) for a, b in zip(gen[:-1], gen[1:])])
                + "x"
            )
        )
        string += " " * (n - e - 1)
        print(string)


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


def main1():
    for ia in gen_distribute_instances(20):
        print(ia)
        # print_lines(distribute_instance_to_ranges(ia))
        print_lines_with_markings(ia)
        print()


def main2():
    while True:
        print_lines_with_markings(random_distribute_instance(80))
        print()
        print()


if __name__ == "__main__":
    main2()
