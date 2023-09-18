from abc import ABC, abstractmethod


class Genotype(ABC):

    """Abstract class for Genotype"""

    @staticmethod
    @abstractmethod
    def from_gametes(gx, gy):
        pass


class Diploid(ABC, Genotype):

    """Abstract class for Diploid genotypes"""

    @abstractmethod
    def upper(self):
        pass

    def lower(self):
        pass


class Dominance(ABC):

    """Abstract class for dominance relations"""

    @staticmethod
    @abstractmethod
    def dom(x, y) -> bool:
        pass


class Unionable(ABC):

    """Abstract class for feasibility checks"""

    @staticmethod
    @abstractmethod
    def is_feasible(pop) -> bool:
        pass

    @staticmethod
    @abstractmethod
    def goal(n_loci):
        pass


class InitialPop(ABC):

    """Abstract class for initial population generators"""

    @staticmethod
    @abstractmethod
    def initial_pop_random(n_loci, n_elements):
        pass

    @staticmethod
    @abstractmethod
    def initial_pop_trait_introgression(n_loci, n_holes, n_elite, n_donor):
        pass

    @staticmethod
    @abstractmethod
    def initial_pop_singles_homo(n_loci):
        pass

    @staticmethod
    @abstractmethod
    def initial_pop_singles_hetero(n_loci):
        pass


class Crossover(ABC):

    """Abstract class for crossovers. These are classes that implement cross"""

    @staticmethod
    @abstractmethod
    def cross(n_loci, x: "Genotype"):
        pass
