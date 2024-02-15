from dataclasses import dataclass
from typing import List


@dataclass
class BaseSolution:
    """Solution for BaseMinCrossing and BaseMinGenerations"""

    tree_data: List[List[List[int]]]
    tree_type: List[str]
    tree_left: List[int]
    tree_right: List[int]
    objective: int

    @staticmethod
    def from_dict(d: dict):
        """Returns a BaseSolution from a dictionary."""
        return BaseSolution(
            tree_data=d["treeData"],
            tree_type=d["treeType"],
            tree_left=d["treeLeft"],
            tree_right=d["treeRight"],
            objective=d["objective"],
        )

    def to_dict(self) -> dict:
        """Returns a dictionary from a BaseSolution."""
        return {
            "treeData": self.tree_data,
            "treeType": self.tree_type,
            "treeLeft": self.tree_left,
            "treeRight": self.tree_right,
            "objective": self.objective,
        }

    @property
    def n_loci(self):
        if len(self.tree_data) == 0:
            raise ValueError("tree_data is empty")
        return len(self.tree_data[0][0])

    @property
    def n_plants(self):
        return len(self.tree_data)

    def __str__(self):
        return "\n".join([
            f"{self.__class__.__name__}(",
            "\ttree_data=[",
            *["\t\t" + repr(genotype) for genotype in self.tree_data],
            "\t],",
            f"\ttree_type={repr(self.tree_type)},",
            f"\ttree_left={repr(self.tree_left)},",
            f"\ttree_right={repr(self.tree_right)},",
            f"\tobjective={repr(self.objective)}",
            ")"
        ])

    def generations(self):
        gen = [0] * self.n_plants

        def _aux(i):
            if self.tree_type[i] == "Null":
                raise ValueError("looking for the generation number of a Null")
            if self.tree_type[i] == "Leaf":
                return 0
            if gen[i] > 0:
                return gen[i]
            res_l = _aux(self.tree_left[i] - 1)     # Because one-indexed
            res_r = _aux(self.tree_right[i] - 1)
            res = 1 + max(res_l, res_r)
            gen[i] = res
            return res

        return _aux(0)

    def crossings(self):
        return sum(t == 'Node' for t in self.tree_type)
