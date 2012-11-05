Extinction/recolonisation coalescent simulator
==============================================

.. toctree::
   :maxdepth: 2

.. contents::


------------
Introduction
------------

This document provides a reference to the ``ercs`` Python module, which provides 
a straightforward interface to coalescent simulations of the 
extinction/recolonisation model. See [E08]_, [BEV10]_, [BKE10]_  and [BEV12]_

-----------
Examples
-----------

Simulating the coalescent for the extinction/recolonisation model using 
:mod:`ercs` follows a basic pattern:

1. Allocate an instance of :class:`ercs.Simulator`;

2. Set the properties of the desired simulation by setting values to the 
   attributes of this object;

3. Run the simulation for a given random seed using the :meth:`ercs.Simulator.run`
   method;

4. Analyse the resulting genealogies to obtain the information of interest.

In the following examples we look at the parameters of the simulation, the 
structure of the simulated genealogies (and how we can analyse them) and 
how we can use these tools to estimate values of interest.

**************************
Basic usage
**************************

To simulate the history of a set of individuals at a sample of 
locations on a 2D torus, we first allocate an instance of the 
:class:`ercs.Simulator` class. This class has a number of attributes
which can be set to describe the parameters of the desired simulation.
Most of these parameters have sensible defaults, but we must specify
at least three of these before we can run a simulation. Here is a 
simple example::

    import ercs
    
    def first_example(seed):
        sim = ercs.Simulator(20)
        sim.sample =  {1:(0, 0), 2:(0, 5), 3:(0, 10)}
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        return sim.run(seed) 

In this example we allocate a simulator on a torus of diameter 20,  
set up our sample and event classes, and the run the simulation
returning the resulting genealogy. The size of the torus is
rather arbitrary, as the real size of the habitat that 
we imagine our population evolving on is determined by the scale 
of events relative to the size of the torus. Thus, the size of the 
torus can be any value you wish, once the spatial extent of events 
is scaled appropriately. For the following examples we'll tend to 
use rather large events, as it's useful to have examples that 
run quickly. These are very unrealistic evolutionary scenarios. 

The initial locations of the lineages whose ancestry we wish to simulate
are specified using the :attr:`ercs.Simulator.sample` attibute. These
are 2-tuples describing locations in two dimensions on the torus.
Here, we simulate the history of three locations, 
``(0, 0)``, ``(0, 5)`` and ``(0, 10)``.

Before we can simulate the history of this sample, we must describe the 
model under which we imagine the population has evolved. This is done by 
allocating some objects that describe the type of events that we 
are interested in, and assigning these to the 
`ercs.Simulator.event_classes` attribute. In the example above, 
we state that all events in the simulation are from the Disc model, 
and they have radius ``r =  1`` and impact ``u = 0.5``.
There can be any number of event classes happening at different 
rates: see `Event Classes`_ for details.

After we have completed setting up the parameters of the simulation
we can then run the simulation by calling the 
:meth:`ercs.Simulator.run` method for a given random seed. 
This returns the simulated history of the sample.

**************************
Oriented trees and forests
**************************

The most part of :mod:`ercs` to understand is the way in which we encode 
genealogies. Running the example above, we get

>>> first_example(3)
([[-1, 4, 4, 5, 5, 0]], [[-1, 0.0, 0.0, 0.0, 30441.574004183603, 46750.11224375103]])

(Note there is nothing special about the seed 3 here---it is just a value which 
produced a neat example to discuss).
This output completely describes the ancestry of the sample, although it's not immediately
obvious how. In ``ercs`` we use *oriented trees* to represent the genealogy of a sample.
In an oriented tree, we are only interested in the parent-child relationships between nodes,
and don't care about the order of children at a node. Therefore, in an oriented tree
``pi``, the parent of node ``j`` is ``pi[j]``. If we map each node in the tree to a unique 
positive integer and adopt the convention that any node whose parent is the 
special "null node" ``0`` is a root, 
we can then represent an oriented tree very simply as a list of integers.

In our example above, we have a list of three locations as our sample, and so we map 
these to the integers ``1``, ``2`` and ``3`` (i.e., lineage ``1`` is sampled at 
location ``(0, 0)`` and so on). The :meth:`ercs.Simulator.run` 
method returns a tuple, ``(pi, tau)``;
``pi`` is a list of oriented forests (one for each locus) and ``tau`` is a list of 
node times (one for each locus). In the example, we are dealing with a single locus 
only, so ``pi`` is a list consisting of one list, ``[-1, 4, 4, 5, 5, 0]``, that 
encodes the following tree:

.. image::  ../images/oriented-tree.png
   :align: center 
   :alt: An oriented tree
   :width: 15cm

It may be easier to see this if we explicity map the nodes to their parents:

>>> pi, tau = first_example(3)
>>> [(node, pi[0][node]) for node in range(1, len(pi[0]))]
[(1, 4), (2, 4), (3, 5), (4, 5), (5, 0)]

.. note:: The zero'th element of an oriented forest and its associated node 
   time list is not used and is set to -1 by convention, 
   following Knuth (Algorithm O, section 7.2.1.6) [K11]_. 

The times labelled on the tree are derived from the node times list for this 
locus, ``tau[0]``. The node times list associated with an oriented tree 
records the time that the associated lineage entered the sample, looking 
backwards in time (hence, for each node in the sample the time is ``0.0``).


Oriented *forests* occur when there is more than one root in a list ``pi``, and 
so we have a set of disconnected trees. This can happen when we specify
the ``max_time`` attribute, potentially stopping the simulation before
the sample has completely coalesced. Consider the following example::
    
    def oriented_forest_example(seed):
        L = 20
        sim = ercs.Simulator(L)
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        sim.sample = [(j, j) for j in range(10)]
        sim.max_time = 1e5
        pi, tau = sim.run(seed)
        return pi[0]

Here we allocate a Simulator on a torus of diameter 20 as before and
use the usual event class. This time we allocate a sample of size 10,
arranged regularly along a line, and stipulate that the simulation should 
continue for no more then 10000 time units. As we're only 
interested in the structure of the genealogy this time, we just 
return the oriented forest at the first locus. Running this, we get 

>>> oriented_forest_example(5)
[-1, 0, 15, 0, 12, 12, 13, 11, 13, 11, 16, 16, 14, 14, 15, 0, 0]

This corresponds to the forest:

.. image::  ../images/oriented-forest.png
   :align: center 
   :alt: An oriented forest 
   :width: 15cm

In this forest there are *four* roots: 1, 3, 15 and 16.

.. note:: This forest is **not** a correct representation of the 
    node times; in any simultation, node ``n + 1`` cannot be more recent 
    than node ``n``.


****************************
Coalescence times and MRCAs 
****************************

The most important quantity in coalescent simulations is the 
*coalescence time* for a set of individuals, or the time back 
to their most recent common ancestor (MRCA). This is
straightforward to do in :mod:`ercs` using the 
:class:`ercs.MRCACalculator` class to find the most recent
common ancestor of two nodes and then looking up the node 
times list to find the time that this node entered the sample.

Suppose we wished to find the coalescence time for a set of  
lineages sampled at a regular set of distances. We can set our 
sample so that the lineages are located at the relevent distances,
but it's not clear how we can get::
    
    def mrca_example(seed):
        L = 40
        sim = ercs.Simulator(L)
        sim.sample = [None] + [(0, j) for j in range(1, 10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        pi, tau = sim.run(seed)
        sv = ercs.MRCACalculator(pi[0])
        for j in range(2, 10):
            mrca = sv.get_mrca(1, j)
            coal_time = tau[0][mrca]
            distance = ercs.torus_distance(sim.sample[1], sim.sample[j], L) 
            print(distance, "\t", coal_time)


Running this

>>> mrca_example(2)
1.0      293.516221072
2.0      1240.09792256
3.0      1505.42515816
4.0      247645.676128
5.0      247645.676128
6.0      263761.554043
7.0      263761.554043
8.0      263761.554043

****************************
Multiple event classes
****************************

Up to this point we have considered models in which a single class 
of event occurs at rate 1.0. We can simulate an arbitrary number 
of event classes happening at different rates, however; we simply
set the :attr:`ercs.Simulator.event_classes` attribute to a list 
consisting of :class:`ercs.EventClass` instances with the required rates and 
properties.

For example, Figure 2 of [BEV12]_  plots the probability of identity in
state for three event regimes: one in which we have small events only;
another in which we have large events only; and finally a regime
in which we have a mixture of the two. The following example 
returns the probability of identity for evenly spaced samples in a 
single replicate::
    
    def mixed_events_example(seed):
        L = 40
        sim = ercs.Simulator(L)
        sim.sample = [None] + [(0, j) for j in range(1, 10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        pi, tau = sim.run(seed)
        sv = ercs.MRCACalculator(pi[0])
        # TODO  Fill me in. 
        
****************************
Multiple loci
****************************

For meaningful multilocus simulations we must specifify the rate at which 
recombination occurs between loci. In the extinction/recolonisation model 
this is done by describing the probability that two adjacent loci ``l`` 
and ``l + 1`` descend from different parents at an event. Therefore, in 
a system of ``m`` loci, we need a list of ``m - 1`` recombination 
probabilities to describe the system in a general way.

The :attr:`ercs.Simulator.recombination_probabilities` attribute is then used 
to describe both the number of loci and the recombination rates between them.
By default, this attibute is set to the empty list. In the following example
we compute the joint probability of identity in state given a mutation 
rate ``mu`` at two loci in a single 
replicate::
    
    def two_locus_example(seed, mu):
        sim = ercs.Simulator(40)
        sim.sample = [None] + [(10, 10), (20, 10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        sim.recombination_probabilities = [0.1]
        pi, tau = sim.run(seed)
        t1 = tau[0][ercs.MRCACalculator(pi[0]).get_mrca(1, 2)]
        t2 = tau[1][ercs.MRCACalculator(pi[1]).get_mrca(1, 2)]
        return math.exp(-2 * mu * t1) * math.exp(-2 * mu * t2)

(Since our sample is of size two, the MRCA of nodes ``1`` and ``2`` must 
be ``3``, so we don't really need to allocate :class:`ercs.MRCACalculator`
objects in this instance. However, in general we need to allocate a different
:class:`ercs.MRCACalculator` for each locus.)

Running this example gives us:

>>> two_locus_example(30, 1e-7)
0.06931300943428219

.. note:: Loci are zero-indexed in the usual Python way, unlike individuals 
    in the sample.

The example above shows how we can simulate two loci. This can be generalised
to larger numbers of loci in a straightforward way. However, since the number 
of lineages can grow to be very large when we deal with large numbers of loci,
it is necessary to become more familiar with some more advanced properties 
of the simulation.

The first issue is to decide how much memory you are willing to dedicate 
to the task of tracking lineages. The C library allocates all its memory in 
advance and fails in a predictable way when it reaches the point in the
simulation where there are too many lineages to fit into this space. 
This is illustrated in the following example::

    def out_of_memory_example():
        sim = ercs.Simulator(40)
        sim.sample = [None] + [(10, 10), (20, 10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        sim.max_lineage_memory = 1
        sim.recombination_probabilities = [0.1 for j in range(500)]
        pi, tau = sim.run(1)

>>> out_of_memory_example()
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "example.py", line 77, in out_of_memory_example
    pi, tau = sim.run(1)
  File "ercs.py", line 116, in run
     self.max_time, 0)
_ercs.LibraryError: Out of lineage memory

The example fails because we try to simulate 500 locus individuals with 
only one MiB of memory. There are two things we can do to alleviate 
this problem:

1. Increase the amount of memory available to the library;

2. Stop the simulation before it runs out of memory using the 
   :attr:`ercs.Simulator.max_time` attribute.



**************************
Parallelism
**************************

The most common use of coalescent simulation is to estimate the distribution 
of some quantity by aggregating over many different replicates. This is 
done in ``ercs`` by running the ``run`` method with different random 
seeds, one for each replicate. Since each replicate is then completely 
independant, we can easily parallise the process. One possible way 
to this is using the :mod:`multiprocessing` module::
    
    import ercs
    import multiprocessing

    def parallel_run(seed):
        sim = ercs.Simulator(50)
        sim.sample = [(1, 1), (2, 2)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        pi, tau = sim.run(seed) 
        coal_time = tau[0][3]
        return coal_time 

    def parallel_example(num_replicates):
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())       
        coal_times = pool.map(parallel_run, range(1, num_replicates + 1)) 
        return sum(coal_times) / num_replicates

In this example we are working on a torus of diameter 100, so the simulations 
require a lot longer to run. On most modern systems we have many CPU cores, 
and so we use the multiprocessing module to distribute the work of many
replicates across these cores.

>>> parallel_example(100)
2968953.9276501946

This is the mean coalescence time among 100 replicates. The multiprocessing 
module runs ``parallel_run`` function for each of the seeds in 
a subprocess and collects the coalescence times into the list 
``coal_times``. We then take the mean of this list and return it.
The random seeds are simply the integers from 1 to 100. This is a perfectly
legitimate way to choose seeds for a single example, since the random 
sequences for adjacent seeds are not correlated. If, however, we are 
are doing lots of simulations with different parameter values or are
distributing our simulations over several machines, it would be better 
to spread our choice of seeds out more evenly across the possible space.
One way to do this is::
    
    import random
    seeds = [random.randint(1, 2**31 - 1) for j in range(100)]

There is no issue with using the Python random generator within your
code, as the ``ercs`` C library generates random numbers independantly 
of Python (using the ``mt19937`` random number generator from the 
`GNU Scientific Library 
<http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html>`_).


-----------------------------------
:mod:`ercs` -- Module reference
-----------------------------------

.. automodule:: ercs

**********************
Event Classes
**********************

The classes of event in a given simulation are specified by 
providing a list of Event Class instances 
in the :attr:`ercs.Simulator.event_classes` attribute. 
Two classes of event are currently supported:
the Disc model and the Gaussian model. See [E08]_, [BEV10]_, [BEV12]_
and several other articles for details of the Disc model, and see 
[BKE10]_ and [BEV12]_ for details of the Gaussian model.

.. autoclass:: ercs.EventClass

.. autoclass:: ercs.DiscEventClass

.. autoclass:: ercs.GaussianEventClass


***********************
:class:`ercs.Simulator`
***********************
.. autoclass:: ercs.Simulator

    .. attribute:: sample
       
        The location of lineages at the beginning of the simulation. This 
        must be a non-empty list of two-tuples describing locations 
        within the 2D space defined by the torus.

        **Default value:** None.

    .. attribute:: event_classes
        
        The event classes to simulate. This must be a list of 
        :class:`ercs.EventClass`
        instances. There must be at least one event class specified.
        
        **Default value:** None.

    .. attribute:: torus_diameter

        The diameter of the torus we are simulating on. This defines the 
        size of the space that lineages can move around in. 
        
        **Default value:** Specified at instantiation time.

    .. attribute:: num_parents 

       The number of parents in each event. For a single locus simulation 
       there must be at least one parent and for multi-locus simulations 
       at least two. 
        
       **Default value:** 1 if the simulation is single locus, otherwise 2.

    .. attribute:: recombination_probabilities
    	
       The list of inter-locus recombination probabilities; the length of 
       this list also determines the number of loci for each individual.
       At an event, the probability of locus ``j`` and ``j + 1`` descending 
       from different parents is ``recombination_probablities[j]``.
       The number of loci in the simulation is therefore 
       ``len(recombination_probablities) + 1``.
       
       **Default value:** The empty list [] (so, we have a single locus 
       simulation by default).

    .. attribute:: max_time 

       The maximum amount of time (in simulation units) that we simulate. If 
       this is set to `0.0` the simulation continues until all loci have 
       coalesced.
       
       **Default value:** 0.0

    .. attribute:: max_lineage_memory

       The maximum amount of memory used for tracking lineages in MiB
       (i.e., 2^20 bytes). If the number of lineages 
       we are tracking grows so much that we exceed this limit, 
       the simulation aborts and raises an :exc:`_ercs.LibraryError`.  
       This is an only an approximate limit on the total 
       amount of memory used by the simulation.
       
       **Default value:** 32 MiB  

    .. attribute:: kdtree_bucket_size 

       The number of locations in a leaf node of the kdtree; must be a power of 
       two, greater than 0. The ``kdtree_bucket_size``
       is an advanced parameter that may be useful in tuning preformance when very 
       large numbers of lineages are involved. Larger values will result in less 
       time and memory spent indexing the lineages, but more lineages will need to
       be tested to see if they are within the critical radius of the event. **Note:**
       changing this parameter affects the outcome of simulations! That is, if we 
       change the value of the bucket size, we cannot expect the outcome of two 
       simulations with the same random seed to be the same. The reason for this 
       is that, although we are guaranteed to end up with the same set of lineages
       in an event in any case, the *order* in which they die may be different,
       pushing the simulation onto a different stochastic trajectory.
       
       **Default value:** 1 

    .. attribute:: max_kdtree_insertions

       The maximum number of insertions into the kdtree before a rebuild, or 
       0 if the tree is not to be rebuilt. This parameter is useful for tuning the 
       performance of the simulation when 
       we have large numbers of loci, particularly if we begin with a relatively
       small sample. In this case, as the number of lineages increases over time
       and they spread outwards to cover more and more of the torus, we need to 
       rebuild the index periodically. If we begin with a large sample uniformly
       distributed around the space then this can safely be set to 0.
       
       **Default value:** 0 

    .. automethod:: run



**********************
Utilities
**********************


.. autofunction:: ercs.torus_distance

.. autoclass:: ercs.MRCACalculator 
   :members:


-----------------------------------
:mod:`_ercs` -- Module reference 
-----------------------------------

The :mod:`ercs` module delegates the majority of its work to the low-level :mod:`_ercs`
extension module, which is written in C. It is not recommended to call this 
module directly - the :mod:`ercs` module provides all of the functionality with a 
much more straightforward interface. In the interested of completeness,
however, the low-level module is documented here.

.. automodule:: _ercs

********************
Simulation interface
********************

.. autofunction:: _ercs.simulate 

------------
Bibliography
------------
.. [E08] A. Etheridge. Drift, draft and structure: some mathematical models of evolution,
   *Banach Center Publications* 80, pp 121--144, 2008.
.. [BEV10] N. H.  Barton,  A. M. Etheridge and A. V\ |eacute|\ ber.
    A new model for evolution in a spatial continuum, *Electronic Journal of 
    Probability* 15:7, 2010.
.. [BKE10] N. H.  Barton,  J. Kelleher and A. M. Etheridge.
    A new model for extinction and recolonisation in two dimensions: quantifying phylogeography, 
    *Evolution*, 64(9), pp 2701--2715, 2010.
.. [BEV12] N. H.  Barton,  A. M. Etheridge and A. V\ |eacute|\ ber. 
    Modelling Evolution in a Spatial continuum,
    *J. Stat. Mech.*, to appear, 2012.
.. [K11] D. E. Knuth, Combinatorial Algorithms, Part 1; Volume 4A of *The Art of Computer
    Programming*, 2011.

.. include:: <isolat1.txt>




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

