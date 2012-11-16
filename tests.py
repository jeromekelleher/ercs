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
Unit tests checking the API.
"""
import math
import random
import unittest

import ercs
import _ercs

def get_ll_disc_event(u=1.0, r=1.0, rate=1.0):
    """
    Convenience function to return a low-level disc event.
    """
    return {"type":_ercs.DISC_EVENT_CLASS, "r":r, "u":u, "rate":rate}

def get_ll_gaussian_event(u0=1.0, theta=1.0, alpha=1.0, rate=1.0):
    """
    Convenience function to return a low-level disc event.
    """
    return {"type":_ercs.GAUSSIAN_EVENT_CLASS, "theta":theta, "u0":u0, 
            "rate":rate, "alpha":alpha}


def random_point(L):
    """
    Returns a point from a 2D torus of diameter L.
    """
    return (random.uniform(0, L), random.uniform(0, L))



class TestLowLevelSimulate(unittest.TestCase):
    """
    Superclass of tests for the low-level api.
    """

    def setUp(self):
        """
        Set up the default values for the simulator.
        """
        self._random_seed = 1
        self._torus_diameter = 10
        self._num_parents = 1
        self._sample = [(0, 0)]
        self._event_classes = [get_ll_disc_event()]
        self._recombination_probabilities = []
        self._kdtree_bucket_size = 1
        self._max_kdtree_insertions = 0
        self._max_lineages = 1000
        self._max_time = 1.0

    def simulate(self):
        """
        Runs the low-level simulate function, returning the results.
        """
        #print("running :", self._event_classes)
        return _ercs.simulate(self._random_seed, self._torus_diameter, 
                self._num_parents, self._sample, self._event_classes, 
                self._recombination_probabilities, self._kdtree_bucket_size, 
                self._max_kdtree_insertions, self._max_lineages, 
                self._max_time, 0)
    
    def test_default_arguments(self):
        """
        Verify that the default arguments don't raise an error.
        """
        self.simulate()

class TestSampleError(TestLowLevelSimulate):
    """
    Tests the various ways in which we can specified erroneous samples.
    """

    def test_bad_types(self):
        self._sample = {} 
        self.assertRaises(TypeError, self.simulate)
        self._sample =  None 
        self.assertRaises(TypeError, self.simulate)
    
    def test_empty_list(self):
        self._sample = []
        self.assertRaises(_ercs.InputError, self.simulate)
        
    def test_bad_values(self):
        errors = [
            [1], # wrong dimensions
            [(1, 2, 3)], # dimension > 2
            [(1, "12")], # non numeric
            [(-1, -1)], # negative
            [(self._torus_diameter + 1, 1)], # outside the torus
        ]
        for error in errors:
            self._sample = error
            self.assertRaises(_ercs.InputError, self.simulate)
        # Mix these errors into a legal sample.
        sample = [(j, j) for j in range(self._torus_diameter - 1)]
        # verify this does not throw an error
        self._sample = sample
        self.simulate()
        for error in errors:
            self._sample = sample + error
            random.shuffle(self._sample)
            self.assertRaises(_ercs.InputError, self.simulate)
        

class TestRecombinationError(TestLowLevelSimulate):
    """
    Test the list of recombination probabilities to see if bad arguments 
    are caught correctly.
    """
    def test_bad_types(self):
        self._recombination_probabilities = {} 
        self.assertRaises(TypeError, self.simulate)
        self._recombination_probabilities = None 
        self.assertRaises(TypeError, self.simulate)
    
    def test_bad_values(self):
        errors = ["0.1", {}, [], None, -1, 1000, 1.01] 
        for error in errors:
            self._recombination_probabilities = [error]
            self.assertRaises(_ercs.InputError, self.simulate)
    
        good_values = [0.1, 0.2, 0.3]
        self._recombination_probabilities = good_values
        self.simulate()
        for error in errors:
            self._recombination_probabilities = good_values + [error]
            random.shuffle(self._recombination_probabilities)
            self.assertRaises(_ercs.InputError, self.simulate)
            
class TestEventClassError(TestLowLevelSimulate):
    """
    Test the list of event_classes to see if bad arguments 
    are caught correctly.
    """
    def test_empty_list(self):
        self._event_classes = []
        self.assertRaises(_ercs.InputError, self.simulate)

    def test_bad_types(self):
        self._event_classes = {} 
        self.assertRaises(TypeError, self.simulate)
        self._event_classes = None 
        self.assertRaises(TypeError, self.simulate)
    
    def test_bad_values(self):
        errors = [
            {},
            {"rate":1.0},
            {"type":500},
            {"type":_ercs.DISC_EVENT_CLASS, "u":1.0},
            get_ll_disc_event(u="xx"),
            get_ll_gaussian_event(u0="test")
        ]
        
        for error in errors:
            self._event_classes = [error]
            self.assertRaises(_ercs.InputError, self.simulate)
        good_values = [get_ll_disc_event(), get_ll_gaussian_event(rate=10),
            get_ll_gaussian_event(u0=0.1), get_ll_disc_event(r=0.1)] 
        self._event_classes = good_values
        self.simulate()
        for error in errors:
            self._event_classes = good_values + [error]
            random.shuffle(self._event_classes)
            self.assertRaises(_ercs.InputError, self.simulate)
           
class TestBadArguments(TestLowLevelSimulate):
    """
    Tests to see if bad arguments for the lesser used parameters to 
    simulate throw an error as required.
    """
    def test_kdtree_bucket_size(self):
        self._kdtree_bucket_size = 0
        self.assertRaises(_ercs.InputError, self.simulate)
        self._kdtree_bucket_size = 3
        self.assertRaises(_ercs.InputError, self.simulate)
        self._kdtree_bucket_size = -1 
        self.assertRaises(_ercs.InputError, self.simulate)

    def test_max_time(self):
        self._max_time = -1
        self.assertRaises(_ercs.InputError, self.simulate)
    
    def test_num_parents(self):
        self._num_parents = 0
        self.assertRaises(_ercs.InputError, self.simulate)
        self._num_parents = -1 
        self.assertRaises(_ercs.InputError, self.simulate)

    def test_torus_diameter(self):
        self._torus_diameter = 0
        self.assertRaises(_ercs.InputError, self.simulate)
        self._torus_diameter = -1 
        self.assertRaises(_ercs.InputError, self.simulate)

    def test_max_lineages(self):
        self._max_lineages = 0
        self.assertRaises(_ercs.InputError, self.simulate)
        self._max_lineages = -1 
        self.assertRaises(_ercs.InputError, self.simulate)
        self._max_lineages = 1 
        self.assertRaises(_ercs.InputError, self.simulate)

    def test_max_kdtree_insertions(self):
        self._max_kdtree_insertions = -1
        self.assertRaises(_ercs.InputError, self.simulate)


class TestOutput(TestLowLevelSimulate):
    """
    Tests the output of simulate to see if it has the right basic 
    properties. To do this, we first setup the test by 
    running the simulation for a range of sample sizes and 
    numbers of loci, and then run some tests on the output.
    """
    def setUp(self):
        TestLowLevelSimulate.setUp(self)
        L = self._torus_diameter
        self._event_classes = [get_ll_disc_event(u=1.0, r=L / 4)]
        self._max_time = 1000
        self._results = {}
        for n in range(2, 10):
            self._sample = [random_point(L) for j in range(n)]
            for m in range(4):
                self._recombination_probabilities = [0.0 for j in range(m)]
                pi, tau = self.simulate() 
                self._results[(n, m)] = pi, tau 
    

    def test_times(self):
        """
        All times should have the following properties: 
        1) tau[1, n] = 0.0
        2) for everything above n, it should be greater than 
            any previous values, if it has coalesced. 
        3) All times should be less than max_time.
        """
        for (n, m), (pi, tau) in self._results.items():
            for l in range(m):
                # Zero'th element doesn't mean anything and is zero.
                self.assertEqual(pi[l][0], 0)
                self.assertEqual(tau[l][0], 0.0) 
                for j in range(1, n + 1):
                    # times for the sample must be 0.0
                    self.assertEqual(tau[l][j], 0.0) 
                max_t = 0.0
                for j in range(n + 1, 2 * n):   
                    t = tau[l][j] 
                    if t != 0.0:
                        self.assertTrue(t >= max_t)
                        max_t = t
                    self.assertTrue(self._max_time >= max_t)


def oriented_forests(n):
    """
    Implementation of Algorithm O from TAOCP section 7.2.1.6. 
    Generates all canonical n-node oriented forests.
    """
    p = [k - 1 for k in range(0, n + 1)]
    k = 1
    while k != 0:
        yield p
        if p[n] > 0:
            p[n] = p[p[n]]
            yield p
        k = n
        while k > 0 and p[k] == 0:
            k -= 1
        if k != 0:
            j = p[k]
            d = k - j
            notDone = True
            while notDone:
                if p[k - d] == p[j]:
                    p[k] = p[j]
                else:
                    p[k] = p[k - d] + d
                if k == n:
                    notDone = False
                else:
                    k += 1


def get_mrca(pi, x, y):
    """
    Returns the most recent common ancestor of nodes x and y in the 
    oriented forest pi.
    """
    x_parents = [x]
    j = x
    while j != 0:
        j = pi[j]
        x_parents.append(j)
    y_parents = {y:None}
    j = y
    while j != 0:
        j = pi[j]
        y_parents[j] = None
    # We have the complete list of parents for x and y back to root.
    mrca = 0
    j = 0
    while x_parents[j] not in y_parents:
        j += 1
    mrca = x_parents[j]
    return mrca

class TestMRCACalculator(unittest.TestCase):
    """
    Class to test the Schieber-Vishkin algorithm. 
    """
        
    def test_all_oriented_forsts(self):
        """
        Runs through all possible oriented forests and checks all possible 
        node pairs using an inferior algorithm.
        """
        for n in range(2, 9):
            for pi in oriented_forests(n):
                sv = ercs.MRCACalculator(pi)
                for j in range(1, n + 1):
                    for k in range(1, j + 1):
                        mrca = get_mrca(pi, j, k)
                        self.assertEqual(mrca, sv.get_mrca(j, k))
        
    def test_simulated_oriented_forests(self):
        """  
        Tests some known oriented forests and MRCA values.
        """
        L = 20
        sim = ercs.Simulator(L)
        sim.sample = [None] + [(j, k) for j in range(L) for k in range(L)]
        sim.event_classes = [ercs.DiscEventClass(rate=1.0, u=1.0, r=1)]
        pi, tau = sim.run(1)
        sv = ercs.MRCACalculator(pi[0])
        n = len(sim.sample)
        for j in range(1, n):
            for k in range(1, j + 1):
                mrca = get_mrca(pi[0], j, k)
                self.assertEqual(mrca, sv.get_mrca(j, k))
        


def test_torus_distance(p1, p2, R):
    """
    Returns the distance between the two specified locations on a 
    square torus of side R.

    This is almost identical to the given definition, but has been 
    in use for a very long time, so can provide a reliable point 
    of reference in case any bugs creep in.
    """
    xabs = math.fabs(p2[0] - p1[0])
    yabs = math.fabs(p2[1] - p1[1])
    xd = min(xabs, R - xabs)
    yd = min(yabs, R - yabs)
    return math.sqrt(xd * xd + yd * yd)


class TestTorusDistance(unittest.TestCase):
    """
    Tests the torus distance utility.
    """
    def test_random_locations(self):
        L = random.uniform(1, 50)
        random_point = lambda: (random.uniform(0, L), random.uniform(0, L))
        for j in range(100):
            x = random_point()
            y = random_point()
            d1 = ercs.torus_distance(x, y, L)
            d2 = test_torus_distance(x, y, L)
            self.assertEqual(d1, d2)
    

if __name__ == '__main__':
    
    unittest.main()
    """
    # check for memory leaks
    mod = __import__("__main__")
    suite = unittest.TestLoader().loadTestsFromModule(mod)
    while True:
        unittest.TextTestRunner(verbosity=0).run(suite)
    """

