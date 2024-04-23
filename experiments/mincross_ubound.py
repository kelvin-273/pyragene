from eugene import utils as eu


def max_repeated_wedges_faster_1(dist_arr):
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
    __import__('pprint').pprint([[int(x) for x in row] for row in adjacent])

    xs = [False] * n_diff

    # get sorted list of groups
    group_counts = []
    checked = [False] * n_diff
    for i in range(n_diff):
        if not checked[i]:
            count_extra = sum(adjacent[i])
            group_counts.append((i, count_extra))
            checked[i] = True
            for j in range(i + 1, n_diff):
                if adjacent[i][j]:
                    checked[j] = True
    __import__('pprint').pprint(group_counts)
    group_counts.sort(key=lambda ic: ic[1])


if __name__ == "__main__":
    print("hello")
    for _ in range(10):
        case = eu.random_distribute_instance(20)
        check = eu.max_repeated_wedges(case)
        guess = max_repeated_wedges_faster_1(case)
        print(check)
