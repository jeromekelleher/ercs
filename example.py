"""
Example using the ercs module
"""
import ercs
import sys
import math
import random

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
    
   


def main():
    #tmp()
    #multiprocessing_example()
    #nca_example()
    #print(first_example(3))
    #print(oriented_forest_example(5))
    #mrca_example(10292, 4)
    print(two_locus_example(30))
    #
    #out_of_memory_example()

if __name__ == "__main__":
    main()
