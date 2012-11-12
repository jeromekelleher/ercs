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


def parallel_run(seed):
    sim = ercs.Simulator(100)
    sim.sample =  [(1, 1), (2, 2)]
    sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
    pi, tau = sim.run(seed) 
    coal_time = tau[0][3]
    return coal_time 

def multiprocessing_example():
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())       
    num_replicates = 100
    
    
    seeds = [random.randint(0, sys.maxsize) for j in range(num_replicates)]
    coal_times = pool.map(parallel_run, seeds) 
    
    
    print("Mean coalescence time =", sum(coal_times) / num_replicates)


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


def tmp():
    sim = ercs.Simulator(10)
    sim.sample =  [(1, 1), (2, 2), (3, 3)]
    sim.recombination_probabilities = [0.1 for j in range(5)]
    sim.event_classes = [ercs.DiscEventClass(rate=1.0, u=0.5, r=1), 
        ercs.GaussianEventClass(rate=0.1, theta=1.5, alpha=0.6, u0=0.75)]
    sim.max_lineages = 4

    for i in range(5):
        pi, tau = sim.run(i)
        print(pi, tau)
    
class SingleLocusIdentitySimulator(ercs.Simulator):
    """
    Class that extends the simulator class and calculates identity 
    in state at a set of distance classes for a replicates.
    """
    def get_identity(self, seed):
        """
        Returns the probability of identity between all sample[0] 
        and sample[1], sample[0] and sample[2], ..., 
        sample[0] and sample[n].
        """
        pi, tau = self.run(seed)
        mc = ercs.MRCACalculator(pi[0])
        F = [0.0 for x in self.sample]
        for j in range(2, len(self.sample)):
            mrca = mc.get_mrca(1, j)
            F[j] = 0 if mrca == 0 else math.exp(-2 * self.mu * tau[0][mrca]) 
        return F

    def set_max_time(self, accuracy_goal, num_replicates):
        """
        Sets the maximum amount of time to run the simulation based on having
        the absolute specified accuracy goal over the specified number of 
        replicates.
        """
        t = math.log(num_replicates * accuracy_goal) / (-2 * self.mu)
        self.max_time = t 

    def get_distances(self):
        """
        Returns the list of distances between sample[1] and all other elements.
        """
        x = self.sample[1]
        d = [ercs.torus_distance(x, y, self.torus_diameter) for y in self.sample[1:]]
        return [-1] + d 

def subprocess_runner(t):
    sim, seed = t
    return sim.get_identity(seed)

def get_mean(replicates):
    """
    Returns the columnwise mean of the specified matrix.
    """
    a = np.array(replicates) 
    return np.mean(a, axis=0)

def run_replicates(sim, filename, num_replicates, pool):
    args = [(sim, random.randint(0, sys.maxsize)) for j in range(num_replicates)]
    mean_identity = get_mean(pool.map(subprocess_runner, args))
    small_result = zip(sim.get_distances(), mean_identity)[2:]
    with open(filename, "wb") as f:
        pickle.dump(small_result, f)
    print("wrote ", filename)

def full_example():
    """
    A full example of using ercs to see the effect of mixed events 
    on the probability of identity in state against distance.
    """
    num_distances = 100
    num_replicates = 1000
    small_events = ercs.DiscEventClass(rate=1.0, r=1, u=0.5)
    large_events = ercs.DiscEventClass(rate=0.1, r=10, u=0.05)
    sim = SingleLocusIdentitySimulator(100)
    sim.mu = 1e-6
    sim.set_max_time(1e-10, num_replicates)
    s = [(0, j * 20.0 / num_distances) for j in range(num_distances)]
    sim.sample = [None] + s 
    # Setup our worker pool and do some work 
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())       
    sim.event_classes = [small_events]
    run_replicates(sim, "small.dat", num_replicates, pool)
    sim.event_classes = [large_events]
    run_replicates(sim, "large.dat", num_replicates, pool)
    sim.event_classes = [small_events, large_events]
    run_replicates(sim, "mixed.dat", num_replicates, pool)


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
    full_example()
if __name__ == "__main__":
    main()
