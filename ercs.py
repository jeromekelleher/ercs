# 
# Copyright (C) 2012 Jerome Kelleher <jerome.kelleher@ed.ac.uk>
#
# This file is part of ercs.
# 
# ercs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# ercs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with ercs.  If not, see <http://www.gnu.org/licenses/>.
# 
"""
Simulate the coalescent in the extinction recolonisation model. 
"""

# Must be single quotes for parsing in setup.py
__version__ = '1.0.1'

import math
import _ercs

def torus_distance(x, y, L):
    """
    Returns the Euclidean distance between two points x and y on a 2D
    square torus with diameter L.
   
    :param x: first point
    :type x: two-tuple of numeric values
    :param y: second point
    :type y: two-tuple of numeric values
    :param L: torus diameter
    :rtype: floating point value
    """
    xabs = math.fabs(y[0] - x[0]);
    yabs = math.fabs(y[1] - x[1]);
    xd = min(xabs, L - xabs);
    yd = min(yabs, L - yabs);
    return math.sqrt(xd**2 + yd**2)


class Simulator(object):
    """
    Class representing a coalescent simulator for the extinction/recolonisation
    model on a torus of the specified diameter. 
    """
    def __init__(self, torus_diameter):
        self.torus_diameter = torus_diameter
        self.sample = None 
        self.event_classes = None 
        self.num_parents = None
        self.recombination_probabilities = None 
        self.kdtree_bucket_size = None
        self.max_kdtree_insertions = None
        self.max_lineages = None 
        self.max_time = None

    def __set_defaults(self):
        """
        Sets up the default values for instances that have not been 
        specified.
        """
        if self.recombination_probabilities == None:
            self.recombination_probabilities = []
        if self.max_time == None:
            self.max_time = 0
        if self.max_kdtree_insertions == None:
            self.max_kdtree_insertions = 0    
        if self.kdtree_bucket_size == None:
            self.kdtree_bucket_size = 1
        m = len(self.recombination_probabilities) + 1
        if self.max_lineages == None:
            self.max_lineages = 1000 
        if self.num_parents == None:
            self.num_parents = 1 if m == 1 else 2 

    def __convert_sample(self):
        """
        Coverts the sample of locations in the sample instance variable to 
        the format used by _ercs, a zero-indexed list. The zero'th element 
        of this list must be None
        """
        if self.sample[0] != None:
            raise ValueError("zeroth element of list samples must be None")
        sample = self.sample[1:] 
        return sample

    def run(self, random_seed):
        """
        Runs the coalescent simulation for the specified random seed, 
        and returns the simulated history, (pi, tau). The history consists
        of a list of oriented forests (one for each locus) and their
        corresponding node times (one for each locus).
        
        :param random_seed: the value to initialise the random number
            generator
        :type random_seed: integer.
        :return: the simulated history of the sample, (pi, tau)
        :rtype: a tuple ``(pi, tau)``; ``pi`` is a list of lists of integers, 
            and ``tau`` is a list of lists of doubles
        :raises: :exc:`_ercs.InputError` when the input is not correctly formed
        :raises: :exc:`_ercs.LibraryError` when the C library encounters an 
            error
        """
        self.__set_defaults()
        sample = self.__convert_sample()
        llec = [ec.get_low_level_representation() for ec in self.event_classes]
        pi, tau = _ercs.simulate(random_seed, self.torus_diameter, 
                self.num_parents, sample, llec, 
                self.recombination_probabilities, self.kdtree_bucket_size, 
                self.max_kdtree_insertions, self.max_lineages, 
                self.max_time, 0)
        # trim the output and set pi[0] = -1
        for j in range(len(pi)):
            pi[j][0] = -1
            tau[j][0] = -1
            last_node = max(len(self.sample), max(pi[j])) + 1
            pi[j] = pi[j][:last_node]
            tau[j] =tau[j][:last_node]
        return pi, tau


class EventClass(object):
    """
    Class representing the an Event Class in the extinction/recolonisation model. 
    Events of a particular class occur at a specific rate and have fixed parameters,
    the details of which depend on the specific model.
    """
    _TYPE = "type"
    _RATE = "rate"
    def __init__(self, rate=1.0):
        self._ll_representation = {self._RATE:rate}
        self.rate = rate
    
    def get_low_level_representation(self):
        """
        Returns the low-level dictionary representation of this event class. 
        """
        return self._ll_representation
    
    def __str__(self):
        return str(self._ll_representation)

class DiscEventClass(EventClass):
    """
    Class representing events from the Disc model, in which all individuals within
    distance *r* of the centre of an event have probability *u* of dying in the 
    event and parents are thrown down uniformly within this disc. 
    """
    _R = "r"
    _U = "u"
    def __init__(self, r, u, rate=1.0):
        super(DiscEventClass, self).__init__(rate)
        self.r = r
        self.u = u
        self._ll_representation[self._TYPE] = _ercs.DISC_EVENT_CLASS
        self._ll_representation[self._R] = r
        self._ll_representation[self._U] = u
        
class GaussianEventClass(EventClass):
    """
    Class representing events from the Gaussian model, in which an individual at  
    distance *d* of the centre of the event have probability 
    :math:`u_0\\exp(-d^2/(2\\theta^2))` of dying in the event. Parents are thrown down 
    around the centre of the event according to a 2D Gaussian with variance 
    :math:`\\theta^2\\alpha^2`.  
    """
    _THETA = "theta"
    _ALPHA = "alpha"
    _U0 = "u0"
    def __init__(self, theta, alpha, u0, rate=1.0):
        super(GaussianEventClass, self).__init__(rate)
        self.theta = theta
        self.alpha = alpha
        self.u = u0
        self._ll_representation[self._TYPE] = _ercs.GAUSSIAN_EVENT_CLASS
        self._ll_representation[self._THETA] = theta
        self._ll_representation[self._ALPHA] = alpha
        self._ll_representation[self._U0] = u0
        



class MRCACalculator(object):
    """
    Class to that allows us to compute the nearest common ancestor of arbitrary
    nodes in an oriented forest.
    
    This is an implementation of Schieber and Vishkin's nearest common ancestor 
    algorithm from TAOCP volume 4A, pg.164-167 [K11]_. Preprocesses the 
    input tree into a sideways heap in O(n) time and processes queries for the 
    nearest common ancestor between an arbitary pair of nodes in O(1) time.
    
    :param oriented_forest: the input oriented forest
    :type oriented_forest: list of integers
    """
    LAMBDA = 0

    def __init__(self, oriented_forest):
        self.__preprocess(oriented_forest)

    def __preprocess(self, oriented_forest):
        """
        Preprocess the oriented forest, so that we can answer mrca queries 
        in constant time.
        """
        n = len(oriented_forest)
        child = [self.LAMBDA for i in range(n)]
        parent = [self.LAMBDA for i in range(n)]
        sib = [self.LAMBDA for i in range(n)]
        self.__lambda = [0 for i in range(n)]
        self.__pi = [0 for i in range(n)]
        self.__tau = [0 for i in range(n)]
        self.__beta = [0 for i in range(n)]
        self.__alpha = [0 for i in range(n)]
        for u in range(n):
            v = oriented_forest[u]
            sib[u] = child[v]
            child[v] = u
            parent[u] = v
        p = child[self.LAMBDA]
        n = 0
        self.__lambda[0] = -1 
        while p != self.LAMBDA:
            notDone = True
            while notDone:
                n += 1
                self.__pi[p] = n
                self.__tau[n] = self.LAMBDA
                self.__lambda[n] = 1 + self.__lambda[n >> 1]
                if child[p] != self.LAMBDA: 
                    p = child[p]
                else:
                    notDone = False
            self.__beta[p] = n
            notDone = True
            while notDone: 
                self.__tau[self.__beta[p]] = parent[p]
                if sib[p] != self.LAMBDA:
                    p = sib[p]
                    notDone = False
                else:
                    p = parent[p]
                    if p != self.LAMBDA:
                        h = self.__lambda[n & -self.__pi[p]]
                        self.__beta[p] = ((n >> h) | 1) << h
                    else:
                        notDone = False 
        # Begin the second traversal
        self.__lambda[0] = self.__lambda[n]
        self.__pi[self.LAMBDA] = 0
        self.__beta[self.LAMBDA] = 0
        self.__alpha[self.LAMBDA] = 0
        p = child[self.LAMBDA]
        while p != self.LAMBDA:
            notDone = True
            while notDone:
                a = self.__alpha[parent[p]] | (self.__beta[p] & -self.__beta[p])
                self.__alpha[p] = a
                if child[p] != self.LAMBDA: 
                    p = child[p]
                else:
                    notDone = False
            notDone = True
            while notDone: 
                if sib[p] != self.LAMBDA:
                    p = sib[p]
                    notDone = False
                else:
                    p = parent[p]
                    notDone = p != self.LAMBDA

    def get_mrca(self, x, y):
        """
        Returns the most recent common ancestor of the nodes x and y, 
        or 0 if the nodes belong to different trees.

        :param x: the first node
        :type x: positive integer
        :param y: the second node
        :type y: positive integer
        :return: the MRCA of nodes x and y
        :type: non-negative integer
        """
        if self.__beta[x] <= self.__beta[y]:
            h = self.__lambda[self.__beta[y] & -self.__beta[x]]
        else:
            h = self.__lambda[self.__beta[x] & -self.__beta[y]]
        k = self.__alpha[x] & self.__alpha[y] & -(1 << h)
        h = self.__lambda[k & -k]
        j = ((self.__beta[x] >> h) | 1) << h
        if j == self.__beta[x]:
            xhat = x
        else:
            l = self.__lambda[self.__alpha[x] & ((1 << h) - 1)]
            xhat = self.__tau[((self.__beta[x] >> l) | 1) << l]
        if j == self.__beta[y]:
            yhat = y
        else:
            l = self.__lambda[self.__alpha[y] & ((1 << h) - 1)]
            yhat = self.__tau[((self.__beta[y] >> l) | 1) << l]
        if self.__pi[xhat] <= self.__pi[yhat]:
            z = xhat
        else:
            z = yhat
        return z    
    



