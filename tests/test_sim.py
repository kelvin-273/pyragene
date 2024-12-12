import unittest

from eugene.solvers import sim
from eugene.plant_models import plant
from eugene.solvers.sim import PopulationGenerators


class TestSimulation(unittest.TestCase):

    """Test case docstring."""

    def test_simulation_with_all_alleles(self):
        for n_loci in range(10):
            args = sim.Args(
                n_loci,
                0,
                10,
                0.5,
                sim.choose_parents_uniform,
                sim.choose_intermediate_target_uniform,
                PopulationGenerators.generate_initial_pop_trait_introgression,
                pruning=True,
            )
            for _ in range(10):
                res = sim.seq_breeding(args)
                assert res.success
                assert res.n_generations == 0
                assert res.n_plants_max == 0

    def test_single_plants_worst_case_bounds(self):
        for n_loci in range(1, 50):
            for n_remaining_loci in range(n_loci+1):
                args = sim.Args(
                    n_loci,
                    n_remaining_loci,
                    1,
                    0.5,
                    sim.choose_parents_uniform,
                    sim.choose_intermediate_target_uniform,
                    sim.PopulationGenerators.generate_initial_pop_single,
                    pruning=True
                )
                for _ in range(10):
                    res = sim.seq_breeding(args)
                    assert res.n_generations <= n_remaining_loci
        args.n_loci = 4
        args.n_remaining_loci = 1
        for _ in range(50):
            res = sim.seq_breeding(args)
            assert res.n_generations <= args.n_remaining_loci

    def test_generate_single_plant(self):
        from eugene.solvers.sim import PopulationGenerators
        n_loci = 8
        for n_remaining_loci in range(n_loci):
            args = sim.Args(
                n_loci=n_loci,
                n_remaining_loci=n_remaining_loci,
                n_initial_pop=1,
                gamma=0.5,
                choose_parents=sim.choose_parents_uniform,
                choose_intermediate_target=sim.choose_intermediate_target_uniform,
                generate_initial_pop=PopulationGenerators.generate_initial_pop_single,
                pruning=True,
            )
            for _ in range(1000):
                pop = PopulationGenerators.generate_initial_pop_single(args)
                assert len(pop) == 1
                x = pop[0]
                xu = x.chrom1
                xl = x.chrom2
                assert 2*args.n_loci - args.n_remaining_loci == bin(xu).count('1') + bin(xl).count('1'), f"{x} has the wrong number of alleles"

    def test_population_transition(self):
        pop = [plant.Plant(4, 13, 15)]
        x, y = sim.choose_parents_uniform(pop)
        assert x is y and y is pop[0]
        z = sim.choose_intermediate_target_uniform(x, y)

    def test_generate_initial_pop_trait_introgression(self):
        from eugene.solvers.sim import PopulationGenerators
        from eugene.utils import count_ones
        n_loci = 8
        for n_remaining_loci in range(n_loci // 2):
            args = sim.Args(
                n_loci=n_loci,
                n_remaining_loci=n_remaining_loci,
                n_initial_pop=2,
                gamma=0.5,
                choose_parents=sim.choose_parents_uniform,
                choose_intermediate_target=sim.choose_intermediate_target_uniform,
                generate_initial_pop=PopulationGenerators.generate_initial_pop_trait_introgression,
                pruning=True,
            )
            for i in range(100):
                pop = PopulationGenerators.generate_initial_pop_trait_introgression(args)
                elite = pop[0]
                donor = pop[1]
                assert count_ones(elite.chrom1 & elite.chrom2) \
                    == n_loci - n_remaining_loci, f"""
                {elite.chrom1}-{elite.chrom2}
                {n_loci - n_remaining_loci}
                """

    def test_generate_initial_pop_trait_introgression2(self):
        from eugene.solvers.sim import PopulationGenerators
        from eugene.utils import count_ones
        n_loci = 8
        for n_remaining_loci in range(n_loci // 2):
            args = sim.Args(
                n_loci=n_loci,
                n_remaining_loci=n_remaining_loci,
                n_initial_pop=2,
                gamma=0.5,
                choose_parents=sim.choose_parents_uniform,
                choose_intermediate_target=sim.choose_intermediate_target_uniform,
                generate_initial_pop=PopulationGenerators.generate_initial_pop_trait_introgression,
                pruning=True,
            )
            for i in range(100):
                pop = PopulationGenerators.generate_initial_pop_trait_introgression(args)
                elite = pop[0]
                donor = pop[1]
                assert (elite.chrom1 & elite.chrom2 & donor.chrom1 & donor.chrom2) == 0, f"\n\t{(elite, donor)}\n\t{( elite.format(), donor.format() )}"
