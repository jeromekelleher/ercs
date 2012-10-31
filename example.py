"""
Example using the ercs module
"""
import ercs
import sys
import math
import random

import multiprocessing


def parallel_simulate(seed):
    sim = ercs.Simulator(100)
    sim.sample =  [(1, 1), (2, 2)]
    sim.event_classes = [ercs.DiscEventClass(rate=1.0, u=0.5, r=1)]
    pi, tau = sim.simulate(seed) 
    coal_time = tau[0][3]
    return coal_time 

def multiprocessing_example():
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())       
    num_replicates = 100
    seeds = [random.randint(0, sys.maxsize) for j in range(num_replicates)]
    coal_times = pool.map(parallel_simulate, seeds) 
    print("Mean coalescence time =", sum(coal_times) / num_replicates)



def first_example(seed):
    sim = ercs.Simulator(20)
    sim.sample =  [(0, 0), (0, 5), (0, 10)]
    sim.event_classes = [ercs.DiscEventClass(rate=1.0, u=0.5, r=1)]
    return sim.simulate(seed) 

def oriented_forest_example(seed):
    L = 20
    sim = ercs.Simulator(L)
    sim.sample = [(j, j) for j in range(10)]
    sim.event_classes = [ercs.DiscEventClass(rate=1.0, u=0.5, r=1)]
    sim.max_time = 1e5
    pi, tau = sim.simulate(seed)
    return pi[0]

def mrca_example(seed):
    sim = ercs.Simulator(40)
    sim.sample =  [(0, j) for j in range(10)]
    sim.event_classes = [ercs.DiscEventClass(rate=1.0, u=0.5, r=1)]
    pi, tau = sim.simulate(seed)
    sv = ercs.MRCACalculator(pi[0])
    print("d", "\t", "coal_time")
    for j in range(2, 11):
        mrca = sv.get_mrca(1, j)
        coal_time = tau[0][mrca]
        print(j, "\t", coal_time)


def tmp():
    sim = ercs.Simulator(10)
    sim.sample =  [(1, 1), (2, 2), (3, 3)]
    sim.recombination_probabilities = [0.1 for j in range(5)]
    sim.event_classes = [ercs.DiscEventClass(rate=1.0, u=0.5, r=1), 
        ercs.GaussianEventClass(rate=0.1, theta=1.5, alpha=0.6, u0=0.75)]
    sim.max_lineages = 4


    for i in range(5):
        pi, tau = sim.simulate(i)
        print(pi, tau)
    
   


def main():
    #tmp()
    #multiprocessing_example()
    #nca_example()
    #print(first_example(3))
    #print(oriented_forest_example(5))
    mrca_example(2)

if __name__ == "__main__":
    main()
