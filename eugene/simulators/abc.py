from abc import ABC, abstractmethod, abstractclassmethod

class BreedingProgram(ABC):

    """Abstract base class for Breeding Program"""

    @abstractmethod
    def __init__(self, n_loci, plant_type):
        self._n_loci = n_loci
        self._plant_type = plant_type
        self._pruning = True
        self._pop_0 = None
        self._ideotype = None
        self._constraint_time = None

    def set_init_pop(self, pop_0):
        """Set the initial population unwrapped"""

        if not type(pop_0) is list:
            raise TypeError("pop_0 must be a list")

        if not all(isinstance(x, self._plant_type) for x in pop_0):
            raise TypeError("Element in pop_0 mismatches the plant type")

        self._pop_0 = pop_0
        self._pop_current = self._pop_0

    def set_ideotype(self, ideotype):
        self._ideotype = ideotype

    @abstractmethod
    def run():
        pass
