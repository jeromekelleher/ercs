"""
Example using the ercs module
"""
import ercs
import sys
import math
import random

import pickle
import numpy as np
import multiprocessing

from matplotlib import ticker 
from matplotlib import pyplot

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
    sim.max_lineage_memory = 1 
    sim.recombination_probabilities = [0.1 for j in range(500)]
    pi, tau = sim.run(1)


#################################

class SingleLocusIdentitySimulator(ercs.Simulator):
    """
    Class that extends the simulator class and calculates identity 
    in state at a set of distance classes for a replicates.
    """
    def setup(self, mutation_rate, spacing, max_distance):
        """
        Sets up the simulation so that we have a mutation rate mu, a
        distance of spacing between samples and the maximum distance 
        between samples is max_distance.
        """
        self.mutation_rate = mutation_rate
        self.num_distances = int(max_distance / spacing)
        s = [(0, j * spacing) for j in range(self.num_distances)]
        self.sample = [None, (0, 0)] + s 
    
    def set_max_time(self, accuracy_goal, num_replicates):
        """
        Sets the maximum amount of time to run the simulation based on having
        the absolute specified accuracy goal over the specified number of 
        replicates.
        """
        t = math.log(num_replicates * accuracy_goal) / (-2 * self.mutation_rate)
        self.max_time = t 

    def get_distances(self):
        """
        Returns the list of distances between sample[1] and all other elements.
        """
        return [self.sample[j + 2][1] for j in range(self.num_distances)]

    def get_identity(self, seed):
        """
        Returns the probability of identity at all distance classes 
        in this replicate.
        """
        pi, tau = self.run(seed)
        mc = ercs.MRCACalculator(pi[0])
        F = [0.0 for j in range(self.num_distances)]
        for j in range(self.num_distances):
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
    print("wrote ", filename)

def full_example():
    """
    A full example of using ercs to see the effect of mixed events 
    on the probability of identity in state against distance.
    """
    num_replicates = 100000 
    small_events = ercs.DiscEventClass(rate=1.0, r=1, u=0.5)
    large_events = ercs.DiscEventClass(rate=0.1, r=10, u=0.05)
    sim = SingleLocusIdentitySimulator(100)
    sim.setup(1e-6, 0.25, 20)
    sim.set_max_time(1e-6, num_replicates)
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())       
    sim.event_classes = [small_events]
    run_replicates(sim, "small.dat", num_replicates, pool)
    sim.event_classes = [large_events]
    run_replicates(sim, "large.dat", num_replicates, pool)
    sim.event_classes = [small_events, large_events]
    run_replicates(sim, "mixed.dat", num_replicates, pool)
    with open("simulator.dat", "wb") as f:
        pickle.dump(sim, f)
    
def label_form(x, pos): 
    return str(float(x))

def plot():
    small = np.fromfile("small.dat") 
    mixed = np.fromfile("mixed.dat") 
    large = np.fromfile("large.dat") 
    with open("simulator.dat", "rb") as f:
        sim = pickle.load(f)
    x = np.array(sim.get_distances())
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
    #out_of_memory_example()
    #full_example()
    plot()

if __name__ == "__main__":
    main()
