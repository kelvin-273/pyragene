def generate_instances(n_loci):
    key = [0] * n_loci

    def aux(i, c_max):
        if i == n_loci:
            yield key.copy()
        else:
            for c in range(c_max + 1):
                if key[i - 1] != c:
                    key[i] = c
                    yield from aux(i + 1, c_max)

            key[i] = c_max + 1
            yield from aux(i + 1, c_max + 1)

    yield from aux(1, 0)


if __name__ == "__main__":
    import sys
    n_loci = int(sys.argv[1])
    for x in generate_instances(n_loci):
        print(x)
