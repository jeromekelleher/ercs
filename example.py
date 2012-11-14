"""
Example using the ercs module
"""
import ercs
import sys
import math
import random
"""
import pickle
import numpy as np
import multiprocessing

# TEMP: for use on non-X server machine.
import matplotlib
matplotlib.use('Agg')

from matplotlib import ticker 
from matplotlib import pyplot
"""
#################################

def first_example(seed):
    sim = ercs.Simulator(20)
    sim.sample =  [None, (0, 0), (0, 5), (0, 10)]
    sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
    return sim.run(seed) 

def oriented_forest_example(seed):
    sim = ercs.Simulator(20)
    sim.sample = [None] + [(j, j) for j in range(10)]
    sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
    sim.max_time = 1e5
    pi, tau = sim.run(seed)
    return pi[0]

def mrca_example(seed, n):
    sim = ercs.Simulator(20)
    sim.sample = [None] + [(j, j) for j in range(n)]
    sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
    pi, tau = sim.run(seed)
    mc = ercs.MRCACalculator(pi[0])
    print("pi  = ", pi)
    print("tau = ", tau)
    for j in range(1, n + 1):
        for k in range(j + 1, n + 1):
            mrca = mc.get_mrca(j, k)
            t = tau[0][mrca]
            print("\tmrca({0}, {1}) = {2} @ {3:.2f}".format(j, k, mrca, t))

def two_locus_example(seed):
    mu = 1e-7
    sim = ercs.Simulator(40)
    sim.sample = [None] + [(10, 10), (20, 10)]
    sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
    sim.recombination_probabilities = [0.1]
    pi, tau = sim.run(seed)
    return math.exp(-2 * mu * tau[0][3]) * math.exp(-2 * mu * tau[1][3])

def out_of_memory_example():
    sim = ercs.Simulator(40)
    sim.sample = [None] + [(10, 10), (20, 10)]
    sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
    sim.max_lineages = 10 
    sim.recombination_probabilities = [0.1 for j in range(500)]
    pi, tau = sim.run(1)


#################################

class SingleLocusIdentitySimulator(ercs.Simulator):
    """
    Class that calculates identity in state for genes separated by a range 
    of distances.
    """
    def setup(self, num_points, max_distance, mutation_rate): 
        """
        Sets up the simulation so that we calculate identity at the specified 
        number of points, the maximum distance between points is 
        max_distance and mutation happens at the specified rate.
        """
        self.mutation_rate = mutation_rate
        self.distances = np.linspace(0, max_distance, num_points)
        self.sample = [None, (0, 0)] + [(0, x) for x in self.distances] 
    
    def set_max_time(self, accuracy_goal, num_replicates):
        """
        Sets the maximum amount of time to run the simulation based on having
        the absolute specified accuracy goal over the specified number of 
        replicates.
        """
        t = math.log(num_replicates * accuracy_goal) / (-2 * self.mutation_rate)
        self.max_time = t 

    def get_identity(self, seed):
        """
        Returns the probability of identity at all distance classes 
        in this replicate.
        """
        pi, tau = self.run(seed)
        mc = ercs.MRCACalculator(pi[0])
        n = len(self.distances)
        F = [0.0 for j in range(n)] 
        for j in range(n):
            mrca = mc.get_mrca(1, j + 2)
            if mrca != 0:
                F[j] = math.exp(-2 * self.mutation_rate * tau[0][mrca]) 
        return F


def subprocess_runner(t):
    sim, seed = t
    return sim.get_identity(seed)

def run_replicates(sim, filename, num_replicates, pool):
    args = [(sim, random.randint(0, sys.maxsize)) for j in range(num_replicates)]
    replicates = np.array(pool.map(subprocess_runner, args))
    mean_identity = np.mean(replicates, axis=0)
    mean_identity.tofile(filename)
    print(mean_identity)

def run_simulations(num_replicates):
    sim = SingleLocusIdentitySimulator(100)
    sim.setup(50, 20, 1e-6)
    sim.set_max_time(1e-8, num_replicates)
    small_events = ercs.DiscEventClass(rate=1.0, r=1, u=0.5)
    large_events = ercs.DiscEventClass(rate=0.1, r=10, u=0.05)
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())       
    sim.event_classes = [small_events]
    run_replicates(sim, "small.dat", num_replicates, pool)
    sim.event_classes = [large_events]
    run_replicates(sim, "large.dat", num_replicates, pool)
    sim.event_classes = [small_events, large_events]
    run_replicates(sim, "mixed.dat", num_replicates, pool)
    with open("simulator.dat", "wb") as f:
        pickle.dump(sim, f)


def generate_plot():
    small = np.fromfile("small.dat") 
    mixed = np.fromfile("mixed.dat") 
    large = np.fromfile("large.dat") 
    with open("simulator.dat", "rb") as f:
        sim = pickle.load(f)
    x = sim.distances 
    pyplot.plot(x, small, label="small")
    pyplot.plot(x, mixed, label="mixed")
    pyplot.plot(x, large, label="large")
    pyplot.yscale('log')
    pyplot.ylim(min(large), max(small))
    pyplot.gca().yaxis.set_minor_formatter(ticker.ScalarFormatter())
    pyplot.xlabel("x")
    pyplot.ylabel("F(x)")
    pyplot.legend(loc="upper right")
    pyplot.savefig("identity.png", dpi=72)


def main():
    #tmp()
    #multiprocessing_example()
    #nca_example()
    #print(first_example(3))
    #print(oriented_forest_example(5))
    #mrca_example(10292, 4)
    #print(two_locus_example(30))
    #
    out_of_memory_example()
    #run_simulations(10)
    #generate_plot()

if __name__ == "__main__":
    main()
