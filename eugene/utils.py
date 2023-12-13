from math import ceil
from random import randint
from typing import List, Tuple


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


def gen_distribute_instances_with_n_pop(n_loci, n_pop):
    state = [0] * n_loci

    def aux(i=1, max_val=0):
        if i == n_loci:
            yield state.copy()
        else:
            for val_i in range(max_val + 1):
                if val_i != state[i - 1]:
                    state[i] = val_i
                    yield from aux(i + 1, max_val)
            if max_val + 1 < n_pop:
                state[i] = max_val + 1
                yield from aux(i + 1, max_val + 1)

    yield from aux()


def random_distribute_instance(n_loci):
    state = [0] * n_loci
    max_val = 0
    for i in range(1, n_loci):
        next_val = randint(0, max_val)
        if next_val >= state[i - 1]:
            next_val += 1
        state[i] = next_val
        max_val = max(max_val, next_val)
    return state


def distribute_to_ranges(instance_array: List[int]) -> List[Tuple[int, int]]:
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


def distribute_to_index_lists(instance_array: List[int]):
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
    indexes = distribute_to_index_lists(instance_array)
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
                + "x".join(
                    ["-" * (b - a - 1) for a, b in zip(gen[:-1], gen[1:])]
                )
                + "x"
            )
        )
        string += " " * (n - e - 1)
        print(string)


def distribute_to_isolated_subproblems(
    instance: List[int],
) -> List[Tuple[int, int]]:
    n_pop = max(instance) + 1
    s = {}
    e = {}
    for i, x in enumerate(instance):
        e[x] = i
        if x not in s:
            s[x] = i

    ranges = sorted([(s[x], e[x]) for x in range(n_pop)])
    chosen = [True] * n_pop
    for i in range(n_pop - 1):
        s1, e1 = ranges[i]
        s2, e2 = ranges[i + 1]
        if e1 > s2:
            ranges[i + 1] = (s1, max(e2, e1))
            chosen[i] = False
    return [se for i, se in enumerate(ranges) if chosen[i]]


def distribute_lower_bound(instance_array):
    n_loci = len(instance_array)
    n_pop = len(set(instance_array))
    return ceil((n_loci + n_pop) / 2)


def verify_distribute_array(array):
    """
    Throws an assert error if array is invalid distribute array

    >>> sanitise_distribute_array([0, 1, 0, 2, 0])
    >>> sanitise_distribute_array([0, 1, 0, 3, 0])
    AssertionError
    """
    assert array[0] == 0
    x_pre = 0
    x_max = 1
    for x in array[1:]:
        assert x != x_pre
        assert x in range(0, x_max + 1)
        x_pre = x
        if x == x_max:
            x_max += 1
    return array


def sanitise_distribute_array(dist_array: List[int]) -> List[int]:
    """
    Translates any substring of a dist_array to a dist_array

    >>> sanitise_distribute_array([3, 5, 2, 6, 5, 3])
    [0, 1, 2, 3, 1, 0]
    """
    assert len(dist_array) > 0
    out = dist_array.copy()
    d = {dist_array[0]: 0}
    dx_max = 0
    for i in range(len(dist_array)):
        x = dist_array[i]
        if x not in d:
            dx_max += 1
            d[x] = dx_max
        out[i] = d[x]
    return out


def gen_covering_subsets(n_loci: int, segments: list):
    n = len(segments)
    segments.sort()
    selection = [False] * n

    def aux(i, current_end):
        if i == n or current_end == n_loci:
            yield [seg for seg, b in zip(segments, selection) if b]
        elif segments[i][0] > current_end:
            return
        else:
            yield from aux(i + 1, current_end)
            if segments[i][1] > current_end:
                selection[i] = True
                yield from aux(i + 1, max(current_end, segments[i][1]))
                selection[i] = False

    return aux(0, 0)


def main1():
    for ia in gen_distribute_instances(20):
        print(ia)
        # print_lines(distribute_to_ranges(ia))
        print_lines_with_markings(ia)
        print()


def main2():
    while True:
        print_lines_with_markings(random_distribute_instance(80))
        print()
        print()


def main3():
    for n_loci in range(1, 14):
        for case in gen_distribute_instances(n_loci):
            iso_probs = distribute_to_isolated_subproblems(case)
            if len(iso_probs) == 1:
                print(case, iso_probs, sep="\t")


if __name__ == "__main__":
    main3()
