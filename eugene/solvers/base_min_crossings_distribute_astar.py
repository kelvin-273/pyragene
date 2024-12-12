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


def branching(state: Node, ctx, simplify_results=True) -> List[Node]:
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


def _simplify_dist_array(dist_array, n_segments):
    out = [None] * n_segments
    n_loci = len(dist_array)
    i = j = 0
    while i < n_segments:
        x = dist_array[j]
        while j < n_loci and dist_array[j] == x:
            j += 1
        out[i] = x
        i += 1
    return out
