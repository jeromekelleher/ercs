Extinction/recolonisation coalescent simulator
==============================================

.. toctree::
   :maxdepth: 2

.. contents::


------------
Introduction
------------

The extinction/recolonisation model (or the spatial Lambda-Fleming-Viot process)
is a recently-introduced model of populations evolving in continuous space.
In this model, individuals occupy a fixed location and all movement and 
reproduction occur as a result of extinction/recolonisation *events*. Events 
fall randomly in space, independent of the location of extant individuals.
In each event, some of the nearby individuals  
die and are replaced by the descendants of a small
number of parents, drawn randomly from the population immediately 
before the event.  See [E08]_, [BEV10]_, [BKE10]_  and [BEV12]_ for 
discussions of the model, its background, and detailed mathematical 
results.

This is the documentation for the :mod:`ercs` Python module, 
which provides a straightforward interface to coalescent simulations 
for this model. The basic approach to simulating 
the history of a sample and interpreting the results is explained through 
a series of straightforward examples in the `Examples`_  section 
along with thorough API documentation for the :mod:`ercs` module.

-----------
Examples
-----------

Simulating the coalescent for the extinction/recolonisation model using 
:mod:`ercs` follows a basic pattern:

1. Allocate an instance of :class:`ercs.Simulator`;

2. Set the properties of the desired simulation by assigning values to the 
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
        sim.sample =  [None, (0, 0), (0, 5), (0, 10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        return sim.run(seed) 

In this example we allocate a simulator on a torus of diameter 20,  
set up our sample and event classes, run the simulation
and the resulting genealogy. The size of the torus is
rather arbitrary, as the real size of the habitat that 
we imagine our population evolving on is determined by the scale 
of events relative to the size of the torus. Thus, the size of the 
torus can be any value you wish, once the spatial extent of events 
is scaled appropriately.

The details of the model under which we imagine the population has evolved
are specified via one or more *event classes*. 
We allocate a list of objects 
that describe the type of events that we 
are interested in, and assign these to the 
:attr:`ercs.Simulator.event_classes` attribute. In the example above, 
we state that all events in the simulation are from the Disc model, 
and they have radius ``r =  1`` and impact ``u = 0.5``.
There can be any number of event classes happening at different 
rates: see `Event Classes`_ for details.


.. note:: In these examples we'll tend to 
    use rather large events, as it's useful to have examples that 
    run quickly. These represent highly unrealistic evolutionary scenarios. 

**************************
Samples
**************************

The initial locations of the lineages whose ancestry we wish to simulate
are specified using the :attr:`ercs.Simulator.sample` attribute, which is a list of
2-tuples describing locations on the torus.
There is a slightly awkward issue to deal with about how we use this list,
however.

In the simulation algorithm, locations in the sample are mapped to the 
positive integers, ``1`` to ``n``. These integers represent the nodes 
that the lineage occupies in the tree describing the history of the 
sample, and this relationship is most simply described using a list.
The awkwardness arises in the discrepancy between this convention of 
counting from ``1`` and Python's convention of making the index of 
the first element of a list ``0``.

To avoid inconsistencies between the input and output of the simulations,
we therefore adopt the convention that the zero'th element of the sample 
must be ``None``. So, in the sample above we simulate the history of 
three locations ``(0, 0)``, ``(0, 5)`` and ``(0, 10)``, which are
mapped to nodes ``1``, ``2`` and ``3`` respectively. Writing the
mapping out explicitly we get:

>>> s = [None, (0, 0), (0, 5), (0, 10)]
>>> {j: s[j] for j in range(1, 4)}
{1: (0, 0), 2: (0, 5), 3: (0, 10)}




**************************
Oriented trees and forests
**************************

After we have completed setting up the parameters of the simulation
we can then run the simulation by calling the 
:meth:`ercs.Simulator.run` method for a given random seed. 
Running the example above, we get

>>> first_example(3)
([[-1, 4, 4, 5, 5, 0]], [[-1, 0.0, 0.0, 0.0, 30441.574004183603, 46750.11224375103]])

.. note:: There is nothing special about the seed 3 here---it is just a value which 
    produced a neat example to discuss.

This output completely describes the ancestry of the sample, although it's not immediately
obvious how. In ``ercs`` we use *oriented trees* to represent the genealogy of a sample.
In an oriented tree, we are only interested in parent-child relationships
and don't care about the order of children at a node. Therefore, in an oriented tree
``pi``, the parent of node ``j`` is ``pi[j]``. If we map each node in the tree to a unique 
positive integer and adopt the convention that any node whose parent is the 
special "null node" ``0`` is a root, 
we can then represent an oriented tree very simply as a list of integers.

The :meth:`ercs.Simulator.run` method returns a tuple, ``(pi, tau)``;
``pi`` is a list of oriented forests (one for each locus) and ``tau`` is a list of 
node times (one for each locus). In the example, we are dealing with a single locus 
only, so ``pi`` is a list consisting of one list, ``[-1, 4, 4, 5, 5, 0]``, that 
encodes the following tree:

.. image::  ../images/oriented-tree.png
   :align: center 
   :alt: An oriented tree
   :width: 15cm

The leaf nodes ``1``, ``2`` and ``3`` correspond to the locations in our sample
as discussed above.
Mapping all nodes in the oriented tree to their parents explicitly, we get:

>>> pi, tau = first_example(3)
>>> {node: pi[0][node] for node in range(1, len(pi[0]))}
{1: 4, 2: 4, 3: 5, 4: 5, 5: 0}

.. note:: The zero'th element of an oriented forest and its associated node 
   time list is not used and is set to -1 by convention, 
   following Knuth (Algorithm O, section 7.2.1.6) [K11]_. 

The times labelled on the tree are derived from the node times list for this 
locus, ``tau[0]``. The node times list associated with an oriented tree 
records the time that the associated lineage entered the sample, looking 
backwards in time (hence, for each node in the sample the time is ``0.0``).


Oriented *forests* occur when there is more than one root in a list ``pi``, and 
so we have a set of disconnected trees. This can happen when we specify
the :attr:`ercs.Simulator.max_time` attribute, potentially stopping the simulation before
the sample has completely coalesced. Consider the following example::
   

    def oriented_forest_example(seed):
        sim = ercs.Simulator(20)
        sim.sample = [None] + [(j, j) for j in range(10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        sim.max_time = 1e5
        pi, tau = sim.run(seed)
        return pi[0]

Here we allocate a Simulator on a torus of diameter 20 as before and
use the usual event class. This time we allocate a sample of size 10,
arranged regularly along a line, and stipulate that the simulation should 
continue for no more then 10\ :sup:`5` time units. As we're only 
interested in the structure of the genealogy this time, we just 
return the oriented forest at the first locus. Running this, we get 

>>> oriented_forest_example(5)
[-1, 0, 15, 0, 12, 12, 13, 11, 13, 11, 16, 16, 14, 14, 15, 0, 0]

In this forest there are *four* roots: 1, 3, 15 and 16.

.. image::  ../images/oriented-forest.png
   :align: center 
   :alt: An oriented forest 
   :width: 15cm

.. note:: This forest is **not** a correct representation of the 
    node times; in any simulation, node ``n + 1`` cannot be more recent 
    than node ``n``.


****************************
Coalescence times and MRCAs 
****************************

In coalescent simulations we are usually interested in the  
coalescence time of two or more lineages in our sample,
which is the time back to their most recent
recent common ancestor (MRCA).
This is straightforward to do in :mod:`ercs` using the 
:class:`ercs.MRCACalculator` class to find the most recent
common ancestor of two nodes and then looking up the node 
times list to find the time that this node entered the sample.

In the following example, we print out the pairwise coalescence 
time for a sample of size 4::
    
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

>>> mrca_example(10292, 4) 
pi  =  [[-1, 5, 5, 6, 6, 7, 7, 0]]
tau =  [[-1, 0.0, 0.0, 0.0, 0.0, 1495.7013764597423, 51935.87882804881, 231859.81270041558]]
        mrca(1, 2) = 5 @ 1495.70
        mrca(1, 3) = 7 @ 231859.81
        mrca(1, 4) = 7 @ 231859.81
        mrca(2, 3) = 7 @ 231859.81
        mrca(2, 4) = 7 @ 231859.81
        mrca(3, 4) = 6 @ 51935.88

Here we can see, for example, that the MRCA of nodes ``1`` and ``2`` is ``5`` 
which can be read directly from ``pi``, as it is also their parent. Then, 
reading the time of node ``5`` from ``tau`` we see that their coalescence 
time is 1495.70. 

****************************
Multiple loci
****************************

In the extinction/recolonisation model 
recombination rates are specified by 
describing the probability that two adjacent loci ``l`` 
and ``l + 1`` descend from different parents at an event. Therefore, in 
a system of ``m`` loci, we need a list of ``m - 1`` recombination 
probabilities to describe the system in a general way.

The :attr:`ercs.Simulator.recombination_probabilities` attribute is used 
to describe both the number of loci and the recombination rates between them.
By default, this attribute is set to the empty list, specifying a single
locus system. In the following example
we compute the joint probability of identity in state given a mutation 
rate ``mu`` at two loci in a single 
replicate::
    
    def two_locus_example(seed):
        mu = 1e-7
        sim = ercs.Simulator(40)
        sim.sample = [None] + [(10, 10), (20, 10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        sim.recombination_probabilities = [0.1]
        pi, tau = sim.run(seed)
        return math.exp(-2 * mu * tau[0][3]) * math.exp(-2 * mu * tau[1][3])


Since our sample is of size two, the MRCA of nodes ``1`` and ``2`` must 
be ``3``. The probability of identity in state for two genes with coalescence
time ``t`` and mutation at rate ``mu`` 
under the infinite sites model is ``exp(-2 * mu * t)``. The joint probability 
of identity is then the product of the probability of identity at 
each locus.

Running this example gives us:

>>> two_locus_example(30)
0.06931300943428219

.. note:: Loci are zero-indexed in the usual Python way, unlike individuals 
    in the sample.

The example above shows how we can simulate two loci, and can be generalised
to larger numbers of loci in a straightforward way. However, since the number 
of lineages can grow to be very large when we deal with large numbers of loci,
it may be necessary to become more familiar with some more advanced properties 
of the simulation implementation.

The first issue is to decide how much memory you are willing to dedicate 
to the task of tracking lineages; this is done by specifiying the maximum 
number of lineages in the sample using the 
:attr:`ercs.Simulator.max_lineages` attribute.
When the number of lineages in the sample
exceeds this limit, the simulation fails with an exception, as 
illustrated in the following example::

    def out_of_memory_example():
        sim = ercs.Simulator(40)
        sim.sample = [None] + [(10, 10), (20, 10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        sim.max_lineages = 10
        sim.recombination_probabilities = [0.1 for j in range(499)]
        pi, tau = sim.run(1)

>>> out_of_memory_example()
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "example.py", line 77, in out_of_memory_example
    pi, tau = sim.run(1)
  File "ercs.py", line 116, in run
     self.max_time, 0)
_ercs.LibraryError: Out of lineage memory

The example fails because we try to simulate the
ancestry of two 500 locus individuals with a limit of only 10 extant individuals. 
If we wish to simulate the entire history of the sample, we must increase
the number maximum number of lineages. On the other hand, if we are only 
interested in the recent history of the sample, we can stop the simulation 
before the sample gets too large using the :attr:`ercs.Simulator.max_time` 
attribute.

Lineages are indexed using a kdtree to allow us to find 
the lineages that may potentially be hit by an event 
efficiently. The 
:attr:`ercs.Simulator.kdtree_bucket_size` 
:attr:`ercs.Simulator.max_kdtree_insertions` 
attributes provide a means of tuning this data structure for 
performance when very large numbers of lineages are involved.

****************************
DRAFT 
****************************
Everything below this point is in development.


****************************
A full example
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
 


The most common use of coalescent simulation is to estimate the distribution 
of some quantity by aggregating over many different replicates. This is 
done in ``ercs`` by running the ``run`` method with different random 
seeds, one for each replicate. Since each replicate is then completely 
independent, we can easily parallelise the process. One possible way 
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
code, as the ``ercs`` C library generates random numbers independently 
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
        must be a list of 2-tuples describing locations 
        within the 2D space defined by the torus. The zero'th element 
        of this list **must be** ``None``.

        **Default value:** ``None``.

    .. attribute:: event_classes
        
        The event classes to simulate. This must be a list of 
        :class:`ercs.EventClass`
        instances. There must be at least one event class specified.
        
        **Default value:** ``None``.

    .. attribute:: torus_diameter

        The diameter of the torus we are simulating on. This defines the 
        size of the space that lineages can move around in. 
        
        **Default value:** Specified at substantiation time.

    .. attribute:: num_parents 

       The number of parents in each event. For a single locus simulation 
       there must be at least one parent and for multi-locus simulations 
       at least two. 
        
       **Default value:** 1 if the simulation is single locus, otherwise 2.

    .. attribute:: recombination_probabilities
    	
       The list of inter-locus recombination probabilities; the length of 
       this list also determines the number of loci for each individual.
       At an event, the probability of locus ``j`` and ``j + 1`` descending 
       from different parents is ``recombination_probabilities[j]``.
       The number of loci in the simulation is therefore 
       ``len(recombination_probabilities) + 1``.
       
       **Default value:** The empty list [] (so, we have a single locus 
       simulation by default).

    .. attribute:: max_time 

       The maximum amount of time (in simulation units) that we simulate. If 
       this is set to `0.0` the simulation continues until all loci have 
       coalesced.
       
       **Default value:** 0.0

    .. attribute:: max_lineages

       The maximum number of extant lineages in the simulation. 
       If the number of lineages we are tracking exceeds this limit, 
       the simulation aborts and raises an :exc:`_ercs.LibraryError`.  
       
       **Default value:** 1000 

    .. attribute:: kdtree_bucket_size 

       The number of locations in a leaf node of the kdtree; must be a power of 
       two, greater than 0. The ``kdtree_bucket_size``
       is an advanced parameter that may be useful in tuning performance when very 
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

