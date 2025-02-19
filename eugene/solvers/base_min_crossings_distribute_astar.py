import dataclasses
import eugene.solvers.base_min_crossings_minizinc as em
from typing import List, Optional, Tuple
import heapq

BRANCHING_CTX = em.MinizincContext.from_solver_and_model_file(
    "gecode", "./eugene/solvers/minizinc/distribute_recomb.mzn"
)


@dataclasses.dataclass
class Node:
    xs: List[int]
    parent_gametes: Optional[Tuple[int, int]]
    parent_node: Optional[object]
    n_gametes: int
    n_segments: int
    g: int
    f: int

    def __lt__(self, other):
        return self.n_gametes + self.n_segments

    @staticmethod
    def from_dist_array(dist_array):
        n_gametes = max(dist_array) + 1
        n_segments = len(dist_array)
        return Node(
            xs=dist_array,
            parent_gametes=None,
            parent_node=Node,
            n_gametes=n_gametes,
            n_segments=n_segments,
            g=0,
            f=(n_segments + n_gametes) / 2
        )

    def __str__(self):
        return f"{self.xs} -> ({self.n_gametes}, {self.n_segments})"


def astar(dist_array, ctx):
    open_list = []
    open_dict = {}
    closed_list = set()

    while len(open_list) > 0:
        node = heapq.heappop(open_list)
        if success(node):
            # TODO: Success case #
            return
        open_list.extend(branching(node, ctx=ctx, simplify_results=True))


def branching(state: Node, ctx) -> List[Node]:
    
    with ctx.instance.branch() as instance:
        instance["instance"] = state.xs
        instance["nLoci"] = state.n_segments
        instance["nPop"] = state.n_gametes
        return [
            Node(
                xs=_simplify_dist_array(res.xs, res.nSegments)
                if simplify_results
                else res.xs,
                parent_gametes=(res.gx, res.gy),
                parent_node=state,
                n_gametes=res.nGametes,
                n_segments=res.nSegments,
                g=state.g + 2,
                f=state.g + 2 + (res.nGametes + res.nSegments) / 2
            )
            # for res in (lambda x: (print(x), x)[1])(instance.solve(all_solutions=True))
            for res in instance.solve(all_solutions=True)
        ]


def success(state: Node):
    len(state.xs) > 0 and all(x == state.xs[0] for x in state.xs)


def _simplify_dist_array(xs, n_segments):
    xs = xs.copy() # get rid of this for more performance (and instability)
    mapping = [None] * max(xs) # assuming that xs is 1-indexed
    mapping[xs[0]] = 0
    x_max = 0
    i = j = 1
    while j < n_segments:
        if xs[j] != xs[j - 1]:
            if mapping[xs[j]] is not None:
                xs[i] = mapping[xs[j]]
            else:
                x_max += 1
                mapping[xs[j]] = x_max
                xs[i] = x_max
            i += 1
        j += 1
    return xs[:i]


def first_full_join_raw(xs):
    n_loci = len(xs)
    n_pop = max(xs) + 1 # assumes xs is 0-indexed

    start_points = [n_loci - 1] * n_pop
    end_points = [0] * n_pop
    for i, gx in enumerate(xs):
        end_points[gx] = i
    for i, gx in reversed(enumerate(xs)):
        start_points[gx] = i

    for (gx, ex) in enumerate(end_points):
        if ex < n_loci - 1 and start_points[xs[ex + 1]] == ex + 1:
            gy = xs[ex + 1]
            out = xs.copy()
            for i in range(start_points[gy], end_points[gy] + 1):
                if xs[i] == gy:
                    out[i] = gx
            return (out, gx, gy)
    return None
