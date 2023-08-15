from typing import List, Tuple


def dp_cost(x: List[Tuple[int]]):
    n_loci = len(x)
    out = [0] * n_loci
    i = 0
    while i < n_loci:
        # INV: i >= 0 or x[i-1] == (0, 0)
        while i < n_loci and x[i] == (0, 0):
            out[i] = 0
            i += 1
        if i == n_loci:
            break

        out[i] = 1
        while i < n_loci and x[i] == (1, 1):
            out[i] = 1
            i += 1
        if i < n_loci and x[i] != (0, 0):
            out[i] = 1

        while i < n_loci and x[i] != (0, 0):
            j = i + 1
            while j < n_loci and (x[j] == x[i] or x[j] == (1, 1)):
                out[j] = out[i]
                j += 1
            if j == n_loci or x[j] == (0, 0):
                i = j
                break

            out[j] = out[i] + 1
            i = j

    return [(z + 1) // 2 for z in out]


if __name__ == "__main__":
    print(dp_cost([(0, 0)] * 2 + [(1, 1)] * 3))
    print(dp_cost(list(zip(
        [0, 1, 1, 1, 0, 0, 1, 1],
        [0, 1, 1, 0, 1, 0, 0, 1],
    ))))
    x = list(zip(
        [0, 1, 0, 1, 0, 1, 0, 1],
        [0, 1, 1, 0, 1, 0, 0, 1],
    ))
    print(dp_cost(x))
    print(list(reversed(dp_cost(list(reversed(x))))))
