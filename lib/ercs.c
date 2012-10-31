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

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include "ercs.h"
#include "torus.h"


const char *
ercs_error_str(int err) 
{
    const char *ret = "Unknown error";
    switch (-err) {
        case ILLEGAL_ARGUMENT:
            ret = "Illegal argument";
            break;
        case OUT_OF_LINEAGES:
            ret = "Out of lineage memory";
            break;
        case OUT_OF_ANCESTRIES:
            ret = "Out of ancestry memory";
            break;
        /*  kdtree errors */
        case OUT_OF_KDT_INTERNAL_NODES:
            ret = "Out of KDT_INTERNAL_NODES";
            break;
        case OUT_OF_KDT_EXTERNAL_NODES:
            ret = "Out of KDT_INTERNAL_NODES";
            break;
        case OUT_OF_LIST_NODES:
            ret = "Out of LIST_NODES";
            break;
        case NOT_POWER_OF_TWO:
            ret = "Bucket size not a power of two";
            break;
        case OUT_OF_LIST_SPACE:
            ret = "Out of LIST_SPACE";
            break;
        case ITERATOR_ERROR:
            ret = "Iterator error";
            break;
    }

    return ret;
}


/*
 * Uses the specified random variate u to return the index of the specified 
 * list with the appropriate probability.
 */
int 
probability_list_select(double *probabilities, const unsigned int n, 
        const double u) 
{
    int ret = 0;
    int x;
    unsigned int i;
    double lhs, rhs;
    assert(n > 0);
    if (n == 1) {
        ret = 0;
    } else {
        x = -1;
        i = 0;
        lhs = 0.0; 
        rhs = probabilities[0];
        while (x == -1) {
            if (lhs <= u && u < rhs) {
                x = (int) i;
            } else {
                i++;
                lhs = rhs;
                rhs += probabilities[i];
            }
        }
        ret = x; 
    }
    return ret;
}

/*
 * Inserts a random point_t from the torus of side R into the specified 
 * double pointer.
 */
static inline void 
random_point(double *p, const double R, gsl_rng *generator) 
{
    p[0] = (double) gsl_ran_flat(generator, 0.0, R);
    p[1] = (double) gsl_ran_flat(generator, 0.0, R);
}

/* Gaussian model */

static double
gaussian_death_probability(event_class_t *event, double d2)
{
    double theta = event->state.gaussian_state.theta;
    double u0 = event->state.gaussian_state.u0;
    return u0 * exp(-d2 / (2 * gsl_pow_2(theta)));
}

static void
gaussian_parent_location(event_class_t *event, double *z, gsl_rng *rng, 
        double *p)
{
    double theta = event->state.gaussian_state.theta;
    double alpha = event->state.gaussian_state.alpha;
    p[0] = z[0] + gsl_ran_gaussian(rng, theta * alpha);  
    p[1] = z[1] + gsl_ran_gaussian(rng, theta * alpha);  
}

static void
gaussian_print_state(event_class_t *event)
{   
    printf("Gauss: theta = %f; alpha = %f u0 = %f", 
            event->state.gaussian_state.theta,
            event->state.gaussian_state.alpha,
            event->state.gaussian_state.u0);
}

static int
gaussian_sanity_check(event_class_t *event, ercs_t *sim) 
{
    int ret = 0;
    double u0 = event->state.gaussian_state.u0;
    double theta = event->state.gaussian_state.theta;
    double alpha = event->state.gaussian_state.alpha;
    if (u0 < 0.0 || u0 > 1.0) {
        ret = ILLEGAL_ARGUMENT;
        goto out;
    } 
    if (theta < 0.0 || theta >= sim->torus_edge / 2.0) {
        ret = ILLEGAL_ARGUMENT;
        goto out;
    } 
    if (alpha < 0.0) {
        ret = ILLEGAL_ARGUMENT;
        goto out;
    }
out: 
    return ret;
}   

void
alloc_gaussian_event_class(event_class_t *event, double rate, double theta,
        double alpha, double u0) 
{
    event->state.gaussian_state.theta = theta;
    event->state.gaussian_state.alpha = alpha;
    event->state.gaussian_state.u0 = u0;
    event->rate = rate;
    event->radius = 3 * theta;
    event->sanity_check = gaussian_sanity_check;
    event->print_state = gaussian_print_state;
    event->death_probability = gaussian_death_probability;
    event->parent_location = gaussian_parent_location;
}

/* Disc events */

static double
disc_death_probability(event_class_t *event, double d2)
{
    double r2 = gsl_pow_2(event->state.disc_state.r); 
    return d2 < r2 ? event->state.disc_state.u : 0.0; 
}

/*
 * Updates the specified point p such that it contains a uniformly distributed
 * within a circle of radius r, centered at z.
 * This uses the canonical method for generating a point uniformly distributed within
 * the unit circle (from TAOCP), then scales and translates the result. 
 */
static void
disc_parent_location(event_class_t *event, double *z, gsl_rng *rng, double *p)
{
    double x, y;
    double s = 1.1;
    double r = event->state.disc_state.r; 
    while (s >= 1) {   
        x = 2 * gsl_rng_uniform(rng) - 1;
        y = 2 * gsl_rng_uniform(rng) - 1;
        s = (x * x) + (y * y);
    }
    p[0] = z[0] + (x * r);
    p[1] = z[1] + (y * r);
}

static void
disc_print_state(event_class_t *event)
{   
    printf("Disc: r = %f; u = %f", event->state.disc_state.r, 
            event->state.disc_state.u);
}

static int
disc_sanity_check(event_class_t *event, ercs_t *sim) 
{
    int ret = 0;
    double u = event->state.disc_state.u;
    double r = event->state.disc_state.r;
    if (u < 0.0 || u > 1.0) {
        ret = ILLEGAL_ARGUMENT;
        goto out;
    } 
    if (r < 0.0 || r >= sim->torus_edge / 2.0) {
        ret = ILLEGAL_ARGUMENT;
        goto out;
    } 
out: 
    return ret;
}   

void
alloc_disc_event_class(event_class_t *event, double rate, double r, double u)
{
    event->state.disc_state.u = u;
    event->state.disc_state.r = r;
    event->rate = rate;
    event->radius = r;
    event->print_state = disc_print_state;
    event->sanity_check = disc_sanity_check;
    event->death_probability = disc_death_probability;
    event->parent_location = disc_parent_location;
}


/*
 * The linear ancestry algorithm.
 */
static void
aa_linear_print_state(ercs_t *self)
{
}

static void
aa_linear_print_ancestry(ercs_t *self, void *ap)
{
    unsigned int j;
    int *a = (int *) ap;
    for (j = 0; j < self->num_loci; j++) {
        printf("%3d ", a[j]); 
    }
    printf("\n");
}


static void *
aa_linear_initialise(ercs_t *self)
{
    aa_linear_t *aas = xmalloc(sizeof(aa_linear_t));
    unsigned int j;
    unsigned int m = self->num_loci;
    unsigned int nu = self->num_parents;
    /* Allocate the ancestry */
    aas->ancestry_memory = xmalloc(m * self->max_lineages * sizeof(int));
    aas->ancestry_heap = xmalloc(self->max_lineages * sizeof(int *));
    for (j = 0; j < self->max_lineages; j++) {
        aas->ancestry_heap[j] = &aas->ancestry_memory[j * m]; 
    }
    aas->ancestry_heap_top = (int) self->max_lineages - 1;
    aas->coalesced_loci = xmalloc(nu * m * sizeof(unsigned int));
    return (void *) aas;
}

static void
aa_linear_free(ercs_t *self)
{
    aa_linear_t *aas = self->aa_state;
    free(aas->ancestry_memory);
    free(aas->ancestry_heap);
    free(aas->coalesced_loci);
    free(aas);
}
static void
aa_linear_free_ancestry(aa_linear_t *aas, void *a)
{
    aas->ancestry_heap_top++;
    aas->ancestry_heap[aas->ancestry_heap_top] = a;
}

void *
aa_linear_alloc_ancestry(aa_linear_t *aas)
{
    int *ancestry = aas->ancestry_heap[aas->ancestry_heap_top];
    aas->ancestry_heap_top--;
    if (aas->ancestry_heap_top < 0) {
        ancestry = NULL; 
    }
    return ancestry;
}

void *
aa_linear_get_initial_ancestry(ercs_t *self, int value)
{
    unsigned int l;
    aa_linear_t *aas = (aa_linear_t *) self->aa_state; 
    int *a = aa_linear_alloc_ancestry(aas);
    if (a != NULL) {
        for (l = 0; l < self->num_loci; l++) {
            a[l] = value; 
        }
    }
    return a;
}


int 
aa_linear_coalesce(ercs_t *self, lineage_t **children, unsigned int num_children, 
        double t, lineage_t **parents, unsigned int *s_p)
{
    int ret = 0;
    aa_linear_t *aas = (aa_linear_t *) self->aa_state; 
    unsigned int s = 0;
    unsigned int m = self->num_loci;
    unsigned int nu = self->num_parents;
    double * rho = self->recombination_probabilities;
    unsigned int j, k, l, v;
    unsigned int *coalesced_loci = aas->coalesced_loci; 
    int *a, *p, g, h;
    for (k = 0; k < nu; k++) {
        a = aa_linear_alloc_ancestry(aas);
        if (a == NULL) {
            ret = -OUT_OF_LINEAGES;
            goto out;
        }
        for (l = 0; l < m; l++) {
            a[l] = 0;
            coalesced_loci[k * m + l] = 0;
        }
        parents[k]->ancestry = a; 
    }
    //printf("%d children:\n", num_children);
    for (j = 0; j < num_children; j++) {
        //printf("\t");
        //aa_linear_print_ancestry(self, children[j]->ancestry);
        a = (int *) children[j]->ancestry;
        k = gsl_rng_uniform_int(self->rng, nu);
        for (l = 0; l < m; l++) { 
            p = (int *) parents[k]->ancestry;
            if (a[l] != 0) {
                g = a[l];
                if (p[l] == 0) {
                    p[l] = g;
                } else {
                    /* coalesce */
                    s++;
                    v = k * m + l;
                    if (coalesced_loci[v] == 0) {
                        h = self->eta[l];
                        self->eta[l]++;
                        self->tau[l][h] = t;
                        self->pi[l][p[l]] = h;
                        coalesced_loci[v] = h;
                    } else {
                        h = coalesced_loci[v];
                    }
                    self->pi[l][g] = h;
                    p[l] = h;
                }
            }

            if (gsl_rng_uniform(self->rng) < rho[l]) {
                // choose a different parent;
                // FIXME
                k++;
                k = k % nu;
            }
        }
        /* we're finished with this child's ancestry now */
        aa_linear_free_ancestry(aas, a);
    }
    for (k = 0; k < nu; k++) {
        a = (int *) parents[k]->ancestry; 
        j = 0;
        //printf("parent %d: ", k);
        for (l = 0; l < m; l++) {
            j += a[l];
            //printf("%d ", a[l]);
        }
        //printf("\n");
        if (j == 0) {
            aa_linear_free_ancestry(aas, a);
            parents[k]->ancestry = NULL; 
        }
    }
    *s_p = s;
    
    //printf("s = %d\n", s);
    //ercs_print_state(self);
out: 
    return ret;
}
void 
ercs_print_state(ercs_t *self)
{
    unsigned int j, k;
    int *pi;
    double *tau;
    double z[] = {0.0, 0.0};
    double L = self->torus_edge;
    lineage_t *lin;
    kri_t *iter = xmalloc(sizeof(kri_t));
    printf("num_parents = %u\n", self->num_parents);
    printf("sample_size = %u\n", self->sample_size);
    printf("max_lineages = %u\n", self->max_lineages);
    printf("kdtree_bucket_size = %u\n", self->kdtree_bucket_size);
    printf("max_kdtree_insertions = %u\n", self->max_kdtree_insertions);
    printf("num_loci = %u\n", self->num_loci);
    printf("random_seed = %ld\n", self->random_seed);
    printf("torus_edge = %f\n", self->torus_edge);
    printf("max_time = %G\n", self->max_time);
    printf("recombination_probabilities = [");
    for (j = 0; j < self->num_loci - 1; j++) {
        printf("%f, ", self->recombination_probabilities[j]);
    }
    printf("]\n");
    printf("events = \n");
    for (j = 0; j < self->num_event_classes; j++) {
        printf("\t%f:", self->event_classes[j].rate);
        self->event_classes[j].print_state(&self->event_classes[j]);
        printf("\n");
    }   
    
    printf("sample = [");
    for (j = 0; j < self->sample_size; j++) {
        printf(" (%f, %f),  ", self->sample[j * 2], self->sample[j * 2 + 1]);
    }
    printf("]\n");
    
    kdtree_get_torus_region_iterator(self->kdtree, z, L, L, iter);
    while ((lin = (lineage_t*) kri_next(iter)) != NULL) {
        printf("(%6.2f, %6.2f);\t::\t", lin->location[0], lin->location[1]);
        self->aa_print_ancestry(self, lin->ancestry);
    }
    
    self->aa_print_state(self); 
    printf("history = \n");
    for (j = 0; j < self->num_loci; j++) {
        pi = self->pi[j];
        tau = self->tau[j];
        printf("\t pi = ");
        for (k = 0; k < 2 * self->sample_size; k++) {
            printf("%3d ", pi[k]); 
        }
        printf(" tau = ");
        for (k = 0; k < 2 * self->sample_size; k++) {
            printf("%f ", tau[k]); 
        }
        printf("\n");
    }
    free(iter);
}

static void
ercs_free_lineage(ercs_t *self, lineage_t *l)
{
    self->lineage_heap_top++;
    assert(self->lineage_heap_top < (int) self->max_lineages);
    self->lineage_heap[self->lineage_heap_top] = l;
}

static lineage_t *
ercs_alloc_lineage(ercs_t *self)
{
    lineage_t *lineage= self->lineage_heap[self->lineage_heap_top];
    self->lineage_heap_top--;
    if (self->lineage_heap_top < 1) {
        lineage = NULL;
    }
    return lineage;
}



int
ercs_initialise(ercs_t *self)
{
    int ret = 0;
    double *x;
    unsigned int j;
    unsigned int m = self->num_loci;
    unsigned int n = self->sample_size;
    lineage_t **sample = xmalloc(n * sizeof(lineage_t *));
    const gsl_rng_type *rng_type = gsl_rng_mt19937;
    self->rng = gsl_rng_alloc(rng_type);
    gsl_rng_set(self->rng, self->random_seed);                                               
    self->lineage_mem = xmalloc(self->max_lineages * sizeof(lineage_t));
    self->lineage_heap = xmalloc(self->max_lineages * sizeof(lineage_t *));
    for (j = 0; j < self->max_lineages; j++) {
        self->lineage_heap[j] = self->lineage_mem + j;
    }    
    self->lineage_heap_top = self->max_lineages - 1;
    self->pi = xmalloc(m * sizeof(int *));
    self->tau = xmalloc(m * sizeof(double *));
    self->eta = xmalloc(m * sizeof(int));
    for (j = 0; j < m; j++) {
        self->pi[j] = xcalloc(2 * n, sizeof(int));
        self->tau[j] = xcalloc(2 * n, sizeof(double));
        self->eta[j] = (int) n + 1;
    }
    self->kdtree = xmalloc(sizeof(kdtree_t));
    self->aa_initialise = aa_linear_initialise;
    self->aa_get_initial_ancestry = aa_linear_get_initial_ancestry;
    self->aa_print_state = aa_linear_print_state;
    self->aa_print_ancestry = aa_linear_print_ancestry;
    self->aa_free = aa_linear_free;
    self->aa_coalesce = aa_linear_coalesce;
    self->aa_state = self->aa_initialise(self);
    /* All mallocing is done - we can now do things that might fail */ 
    ret = kdtree_init(self->kdtree, self->max_lineages, self->kdtree_bucket_size, 
            self->random_seed); 
    ERROR_CHECK(ret, out);
    for (j = 0; j < n; j++) {
        x =  self->sample + 2 * j;
        sample[j] = ercs_alloc_lineage(self);
        if (sample[j] == NULL) {
            ret = -OUT_OF_LINEAGES;
            goto out;
        }
        sample[j]->location[0] = x[0];
        sample[j]->location[1] = x[1];
        sample[j]->ancestry = self->aa_get_initial_ancestry(self, (int) j + 1);
        if (sample[j]->ancestry == NULL) {
            ret = -OUT_OF_ANCESTRIES;
            goto out;
        }

    }
    ret = kdtree_build(self->kdtree, NULL, (point_t **) sample, n);
    ERROR_CHECK(ret, out);
out:
    free(sample);
    return ret;
}

void
ercs_free(ercs_t *self)
{
    unsigned int j;
    for (j = 0; j < self->num_loci; j++) {
        free(self->pi[j]);
        free(self->tau[j]);
    }
    free(self->pi);
    free(self->tau); 
    free(self->eta);
    
    self->aa_free(self);
    free(self->lineage_heap);
    free(self->lineage_mem);
    kdtree_free(self->kdtree);
    free(self->kdtree);
    gsl_rng_free(self->rng);
}

int 
ercs_sanity_check(ercs_t *self) 
{
    int ret = 0;
    unsigned int j, k;
    double v;
    event_class_t *e;
    if (self->torus_edge <= 0.0) {
        ret = -ILLEGAL_ARGUMENT;
        goto out;
    }
    if (self->sample_size == 0) {
        ret = -ILLEGAL_ARGUMENT;
        goto out;
    }
    for (j = 0; j < self->sample_size; j++) {
        for (k = 0; k < 2; k++) {
            v = self->sample[2 * j + k];
            if (v < 0.0 || v > self->torus_edge) {
                ret = -ILLEGAL_ARGUMENT;
                goto out;
            }
        }
    }
    for (j = 0; j < self->num_loci - 1; j++) {
        v = self->recombination_probabilities[j];
        if (v < 0.0 || v > 1.0) {
            ret = -ILLEGAL_ARGUMENT;
            goto out;
        }
    }
    if (self->num_parents == 0 || self->num_parents > 2) { /* FIXME */
        ret = -ILLEGAL_ARGUMENT;
        goto out;
    }
    if (self->num_event_classes == 0) {
        ret = -ILLEGAL_ARGUMENT;
        goto out;
    }
    for (j = 0; j < self->num_event_classes; j++) {
        e = &self->event_classes[j];
        ret = e->sanity_check(e, self);
        ERROR_CHECK(ret, out);
    }
    if (self->max_lineages == 0) {
        ret = -ILLEGAL_ARGUMENT;
        goto out;
    }
    if (self->kdtree_bucket_size == 0) {
        ret = -ILLEGAL_ARGUMENT;
        goto out;
    }
    if (self->max_kdtree_insertions == 0) {
        ret = -ILLEGAL_ARGUMENT;
        goto out;
    }
    if (self->max_time < 0.0) {
        ret = -ILLEGAL_ARGUMENT;
        goto out;
    }
out:
    return ret;
}


int 
ercs_simulate(ercs_t *self)
{
    int ret = 0;
    unsigned int j, k;
    unsigned int n = self->sample_size;
    unsigned int m = self->num_loci;
    unsigned int num_children, num_inserted;
    double *event_probabilities = xmalloc(
            self->num_event_classes * sizeof(double));
    double L = self->torus_edge;
    double z[2] = {0.0, 0.0};
    double d2, u;
    unsigned int kappa = n * m;
    unsigned int kdtree_insertions = 0; 
    double t = 0;
    gsl_rng *rng = self->rng;
    lineage_t *lin;
    kri_t *iter = xmalloc(sizeof(kri_t));
    kdtree_t *space = self->kdtree; 
    point_t **inserted = xmalloc(self->num_parents * sizeof(point_t *));
    lineage_t **parents = xmalloc(self->num_parents * sizeof(lineage_t *));
    lineage_t **children = xmalloc(self->max_lineages * sizeof(lineage_t *));
    event_class_t *event;
    unsigned int num_lineages = n;
    double total_rate;
    total_rate = 0.0;
    for (j = 0; j < self->num_event_classes; j++) {
        total_rate += self->event_classes[j].rate;
    }
    for (j = 0; j < self->num_event_classes; j++) {
        event_probabilities[j] = self->event_classes[j].rate / total_rate;
    }
    ret = ercs_sanity_check(self);
    ERROR_CHECK(ret, out);
    while (kappa > m && t < self->max_time) {
        ret = probability_list_select(event_probabilities, 
                self->num_event_classes, gsl_rng_uniform(self->rng));
        ERROR_CHECK(ret, out);
        event = &self->event_classes[ret];
        t += gsl_ran_exponential(self->rng, 1.0 / total_rate);
        random_point(z, L, rng); 
        ret = kdtree_get_torus_region_iterator(space, z, event->radius, L, 
                iter);
        ERROR_CHECK(ret, out);
        num_children = 0;
        while ((lin = (lineage_t*) kri_next(iter)) != NULL) {
            d2 = torus_squared_distance(z, lin->location, L);
            if (d2 < gsl_pow_2(event->radius)) {
                u = event->death_probability(event, d2);
                if (gsl_rng_uniform(rng) < u) {
                    ret = kri_delete(iter);
                    ERROR_CHECK(ret, out);
                    num_lineages--;
                    children[num_children] = lin;
                    num_children++;
                }
            }
        }
        if (num_children > 0) {
            for (k = 0; k < self->num_parents; k++) {
                parents[k] = ercs_alloc_lineage(self);
                if (parents[k] == NULL) {
                    ret = -OUT_OF_LINEAGES; 
                    goto out;
                }
            }
            ret = self->aa_coalesce(self, children, num_children, t, parents, 
                    &k); 
            ERROR_CHECK(ret, out);
            kappa -= k;
            num_inserted = 0;
            for (j = 0; j < self->num_parents; j++) {
                lin = parents[j];
                assert(lin != NULL);
                if (lin->ancestry != NULL) {
                    event->parent_location(event, z, rng, lin->location);
                    torus_wrap(lin->location, L);
                    inserted[num_inserted] = (point_t *) lin;
                    num_inserted++;
                } else {
                    ercs_free_lineage(self, parents[j]);
                }
            }
            ret = kdtree_insert_points(space, inserted, num_inserted);
            ERROR_CHECK(ret, out);
            kdtree_insertions += num_inserted;
            num_lineages += num_inserted; 
            for (j = 0; j < num_children; j++) {
                ercs_free_lineage(self, children[j]);
            }
        }
        if (kdtree_insertions > self->max_kdtree_insertions) {
            kdtree_insertions = 0;
            j = kdtree_copy_points(self->kdtree, (point_t **) children);
            assert(j == num_lineages);
            kdtree_clear(self->kdtree);
            ret = kdtree_build(self->kdtree, NULL, (point_t **) children, 
                    num_lineages);
            ERROR_CHECK(ret, out);
        }
    } 

out:
    free(iter);
    free(children);
    free(inserted);
    free(parents);
    free(event_probabilities);
    return ret;
}



