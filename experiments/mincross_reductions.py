"""
This file defines functions that encode reductions to MinCross.

The decision problem from MinCross is defined as follows:
    Given initial population P_0, number of loci n_l and integer k,
    Does there exist a crossing schedule exist from genotypes in P_0
    where the number of internal nodes is less than or equal to k
"""


def min_cross(n_loci: int, pop_0: list, k: int):
    pass


def reduction_istrue(x: bool):
    """
    Given boolean variables x and y,
    does there exist a variable z that is true if and only if x and y are true?

    This answers the equivalent question:
    Given formula Φ = { x } does there exist a formula τ that satisfies Φ?
    """
    pop_0 = [
        ("1", "1" if x else "0"),
    ]
    return min_cross(1, pop_0, 0)


def reduction_and_pair_findx(x: bool, y: bool):
    """
    Given boolean variables x and y,
    does there exist a variable z that is true if and only if x and y are true?
    """
    return reduction_and_many_findx([x, y])


def reduction_and_many_findx(xs: list):
    """
    Given a list of boolean variables xs,
    does there exist a variable z that is true if and only if
    all variables in xs are true?
    """
    n_loci = len(xs)
    pop_0 = [
        (
            "0" * (i) + ("1" if xs[i] else "0") + "0" * (n_loci - i - 1),
            "0" * (i) + ("1" if xs[i] else "0") + "0" * (n_loci - i - 1),
        )
        for i in range(n_loci)
    ]
    k = n_loci
    min_cross(n_loci, pop_0, k)


def reduction_not(x: bool):
    """
    Given boolean variable x,
    does there exist a variable z that is true if and only if x is false?

    This answers the equivalent question:
    Given formula Φ = { ¬x } does there exist a formula τ that satisfies Φ?
    """
    return reduction_istrue(not x)


def reduction_and_pair(x: bool, y: bool):
    """
    Given formula Φ = { x ∧ y },
    does there exist an assignment τ that satisfies Φ?

    This is trivially satisfiable.
    """
    pass


def reduction_and_many(xs: list):
    """
    Given a formula Φ of literals of xs
    does there exist an assignment τ that satisfies Φ?
    """
    n_loci = max(abs(x) for x in xs)
    xs = sorted(xs, key=lambda x: abs(x))
    pop_0 = []
    return min_cross(n_loci, pop_0)


def reduction_or_pair(x: bool, y: bool):
    """
    Given boolean variables x and y,
    does there exist a variable z that is true if and only if either of x or y are true?
    """
    pass


def reduction_cnf(formula: list):
    """
    Given a CNF formula Φ = [[x1, x3, ...], [x5, ...], ...]
    where xi ∈ {-n, ..., -1, 1, ..., n}, does there exist an assignment
    τ : {x1, ..., xn} -> {0, 1} that satisfies Φ|τ = 1?

    We reduce this to MinCross, we prove it's NP-Hard.
    """
    pass
