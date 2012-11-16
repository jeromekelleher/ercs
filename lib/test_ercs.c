/*
** Copyright (C) 2012 Jerome Kelleher <jerome.kelleher@ed.ac.uk>
**  
** This file is part of ercs.
** 
** ercs is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
** 
** ercs is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with ercs.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <float.h>
#include <limits.h>

#include "ercs.h"


static ercs_t *
alloc_ercs(double L, int num_points, int num_loci, int num_parents) 
{
    int j;
    ercs_t *sim= xcalloc(1, sizeof(ercs_t));
    sim->num_loci = num_loci;
    sim->num_parents = num_parents;
    sim->num_event_classes = 1;
    sim->kdtree_bucket_size = 1;
    sim->max_lineages = 1000; 
    sim->random_seed = 1; 
    sim->torus_diameter = L;
    sim->sample_size = num_points;
    sim->max_time = DBL_MAX;
    sim->max_kdtree_insertions = INT_MAX;
    sim->sample = xmalloc(2 * num_points * sizeof(double));
    sim->event_classes = xmalloc(sizeof(event_class_t));
    sim->recombination_probabilities = xmalloc(num_loci * sizeof(double));
    for (j = 0; j < num_points; j++) {
        sim->sample[2 * j] = j * L / num_points;
        sim->sample[2 * j + 1] = j * L / num_points;
    }
    for (j = 0; j < num_loci - 1; j++) {
        sim->recombination_probabilities[j] = 0.75; 
    }
    sim->recombination_probabilities[num_loci - 1] = 0.0;
    alloc_disc_event_class(sim->event_classes, 1.0, L / 10.0, 0.95); 
    return sim;
}

static void
free_ercs(ercs_t *sim)
{
    free(sim->event_classes);
    free(sim->sample);
    free(sim->recombination_probabilities);
    ercs_free(sim);
    free(sim);
}


/*
 * Compares the two specified simulators to see if their histories 
 * are equal.
 */
static int 
compare_histories(ercs_t *sim1, ercs_t *sim2)
{
    int ret = 0;
    int l, j;
    int n = sim1->sample_size; 
    for (l = 0; l < sim1->num_loci; l++) {
        for (j = 0; j < 2 * n; j++) {
            if (sim1->pi[l][j] != sim2->pi[l][j]) {
                ret = -576;
                goto out;
            }
            if (sim1->tau[l][j] != sim2->tau[l][j]) {
                ret = -577;
                goto out;
            }
        }
    }   

out:
    return ret;
}

/*
 * Checks the history of the specified simulator to see if it is 
 * consistent.
 */
static int 
check_history(ercs_t *sim)
{
    int ret = 0; 
    int l, j, k, sample_mrca;
    int n = sim->sample_size; 
    int *pi;
    double *tau;
    double sample_ct, max_ct;
    max_ct = -1.0;
    for (l = 0; l < sim->num_loci; l++) {
        pi = sim->pi[l];
        tau = sim->tau[l];
        sample_mrca = 2 * n - 1;
        while (tau[sample_mrca] == 0.0) {
            sample_mrca--;
        }
        if (sample_mrca <= n) {
            ret = -600;
            goto out;
        }
        sample_ct = tau[sample_mrca];
        if (sample_ct > max_ct) {
            max_ct = sample_ct;
        }
        for (j = 1; j <=  n; j++) {
            if (tau[j] != 0.0) {
                ret = -610;
                goto out;
            }
            /* trace the path of every lineage in the sample to root */
            k = j;
            while (pi[k] != 0) {
                k = pi[k];
            }
            if (k != sample_mrca) {
                ret = -611;
                goto out;
            }
        }
    }
    if (max_ct != sim->time) {
        ret = -650;
        goto out;

    }

out:
    return ret;
}


static int 
test_ercs_size(double L, int num_points, int m, int nu)
{
    int ret = 0;
    int not_done;
    int num_events;
    ercs_t *sim1 = alloc_ercs(L, num_points, m, nu);
    ercs_t *sim2 = alloc_ercs(L, num_points, m, nu);
    ret = ercs_initialise(sim1); 
    ERCS_ERROR_CHECK(ret, out);
    ret = ercs_initialise(sim2); 
    ERCS_ERROR_CHECK(ret, out);
    ret = ercs_simulate(sim1, INT_MAX); 
    ERCS_ERROR_CHECK(ret, out);
    not_done = 1;
    num_events = 0;
    while (not_done) {
        ret = ercs_simulate(sim2, num_events);
        ERCS_ERROR_CHECK(ret, out);
        not_done = ret == ERCS_SIM_NOT_DONE;
        num_events++;
    }
    ret = compare_histories(sim1, sim2);
    ERCS_ERROR_CHECK(ret, out);
    ret = check_history(sim1);
    free_ercs(sim1);
    free_ercs(sim2);
out:
    return ret;
}


static int 
test_ercs(void)
{
    int ret = 0;
    size_t j, k, n, m, nu;
    double torus_sizes[] = {1.9, 0.2, 50.3, 123.22, 1111.1};
    int num_points[] = {2, 11, 50, 212, 333, 515};
    printf("\trunning ercs tests");
    for (j = 0; j < sizeof(torus_sizes) / sizeof(double); j++) {
        for (k = 0; k < sizeof(num_points) / sizeof(int); k++) {
            for (m = 1; m < 4; m++) {
                for (nu = 1; nu < 5; nu++) {
                    printf(".");
                    fflush(stdout);
                    n = num_points[k]; 
                    ret = test_ercs_size(torus_sizes[j], n, m, nu);
                    ERCS_ERROR_CHECK(ret, out);
                }
            }
        }
    } 
out:
    return ret;

}
    
int 
main(void)
{
    int result;
    result = test_ercs();
    if (result < 0) {
        printf("Error in ercs tests: %d\n", result);
        exit(result);
    }
    return EXIT_SUCCESS;
}
