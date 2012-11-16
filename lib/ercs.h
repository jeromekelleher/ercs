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

#include "kdtree.h"
#include "util.h"

#include <gsl/gsl_rng.h>

#define ERCS_SIM_DONE 0
#define ERCS_SIM_NOT_DONE 1


/* Errors */
#define OUT_OF_LINEAGES     201
#define OUT_OF_ANCESTRIES   202
#define ILLEGAL_ARGUMENT    203


typedef struct {
    double location[2];
    void *ancestry; 
} lineage_t;

/* ancestry algorithms */
typedef struct {
    int **ancestry_heap;
    int ancestry_heap_top;
    int *ancestry_memory;
    int *coalesced_loci;
} aa_linear_t;

/* Event classes */

typedef struct {
    double r;
    double u;
} disc_event_state_t;

typedef struct {
    double theta;
    double alpha;
    double u0;
} gaussian_event_state_t;


struct ercs_t_t;

typedef struct event_class_t_t {
    double rate;
    double radius;
    union {
        disc_event_state_t disc_state;
        gaussian_event_state_t gaussian_state;
    } state;
    int (*sanity_check)(struct event_class_t_t *, struct ercs_t_t *);
    void (*print_state)(struct event_class_t_t *);
    double (*death_probability)(struct event_class_t_t *, double);
    void (*parent_location)(struct event_class_t_t *, double *, gsl_rng*, 
            double *);
} event_class_t;


typedef struct ercs_t_t {
    /* Input */
    int num_loci;
    int num_parents;
    int sample_size;
    int num_event_classes;
    int kdtree_bucket_size;
    int max_kdtree_insertions;
    int max_lineages;
    long random_seed; 
    double torus_diameter;
    double max_time;
    double *sample;
    event_class_t *event_classes;
    double *recombination_probabilities;
    /* state */
    gsl_rng *rng;
    kdtree_t *kdtree;
    kri_t *kdtree_iterator;
    point_t **insertion_buffer;
    lineage_t **parents_buffer;
    lineage_t **children_buffer;
    int kdtree_insertions;
    double *event_probabilities;
    double total_event_rate; 
    double time;
    /* algorithm state */ 
    int kappa;
    int *eta;
    int **pi;
    double **tau;
    /* lineage memory management */
    lineage_t *lineage_mem;
    lineage_t **lineage_heap;
    int lineage_heap_top;
    /* ancestry algorithm */ 
    void * (*aa_initialise)(struct ercs_t_t *);
    void * (*aa_get_initial_ancestry)(struct ercs_t_t *, int);
    void (*aa_print_state)(struct ercs_t_t *); 
    void (*aa_print_ancestry)(struct ercs_t_t *, void *);
    void (*aa_free)(struct ercs_t_t *); 
    int (*aa_coalesce)(struct ercs_t_t *, lineage_t**, int, double, 
            lineage_t**, int *); 
    void *aa_state;
} ercs_t;

void alloc_gaussian_event_class(event_class_t *event, double, double, double,
        double);
void alloc_disc_event_class(event_class_t *, double, double, double);

int ercs_initialise(ercs_t *);
int ercs_simulate(ercs_t *, unsigned int);
void ercs_print_state(ercs_t *);
void ercs_free(ercs_t *);
const char *ercs_error_str(int);


