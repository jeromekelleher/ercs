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

#ifndef KDTREE_H
#define KDTREE_H

#define OUT_OF_KDT_INTERNAL_NODES 101 
#define OUT_OF_KDT_EXTERNAL_NODES 102
#define OUT_OF_LIST_NODES 103 
#define NOT_POWER_OF_TWO 104
#define OUT_OF_LIST_SPACE 105
#define ITERATOR_ERROR 106


/* The point_t type is the interface to the kdtree_t. Any input points to 
 * the kdtree_t must be castable to the point_t type and have location as its 
 * first member
 */
typedef struct {
    double location[2];
} point_t;

typedef struct list_node_t {
    point_t *point;
    struct list_node_t *next;
} list_node;

typedef struct list_set_t {
    int num_keys;
    int max_keys;
    void **list;
    int total_memory;
} list_set;

typedef struct kd_internal_node_t {
    char flags;
    int cut_dimension;
    double cut_value;
    void *left_child;
    void *right_child;
} kd_internal_node;

typedef struct kd_external_node_t {
    char flags;
    int num_points;
    list_node *head;
} kd_external_node;

typedef struct insertion_stack_entry_t {
    void *node;
    int left;
    int right;
} insertion_stack_entry;

typedef struct {
    void *root;
    int max_points;
    int bucket_max;
    int depth;
    int total_memory;
    
    void **search_stack;
    insertion_stack_entry *insertion_stack;
    list_set *node_set;

    /* Heaps for allocated memory */
    int current_list_node;
    int max_list_node;
    list_node **list_node_heap;
    void *list_node_mem;

    int current_kd_external_node;
    int max_kd_external_node;
    kd_external_node **kd_external_node_heap;
    void *kd_external_node_mem;
    
    int current_kd_internal_node;
    int max_kd_internal_node;
    kd_internal_node **kd_internal_node_heap;
    void *kd_internal_node_mem;
} kdtree_t;


/* kdtree_t region iterator */
typedef struct {
    list_set *node_set;
    kdtree_t *kdtree;
    int current_tree_node_index;
    list_node *previous_list_node;
    list_node *current_list_node;
    list_node *deleted_list_node; /* hack */
    int delete_possible;
} kri_t;


int kdtree_init(kdtree_t *tree, int max_points, int bucket_max, 
        int random_seed);
void kdtree_free(kdtree_t *tree); 
void kdtree_clear(kdtree_t *tree); 
int kdtree_build(kdtree_t *tree, kd_internal_node *root, point_t **points, 
        const int num_points); 
int kdtree_insert_points(kdtree_t *tree, point_t **points, int num_points);
int kdtree_insert_point(kdtree_t *tree, point_t *ind);
int kdtree_copy_points(kdtree_t *tree, point_t **points);
int kdtree_get_torus_region_iterator(kdtree_t *tree,  const double *p, const double r, 
        const double R, kri_t *iterator);
point_t* kri_next(kri_t *iterator);
int kri_delete(kri_t *iterator);

void kdtree_print_heaps(kdtree_t *self); 
void kdtree_print(kdtree_t *tree);


#endif /*KDTREE_H*/
