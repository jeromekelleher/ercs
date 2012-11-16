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

#include "kdtree.h"
#include "util.h"

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>

#define LEFT_CHILD 1 
#define RIGHT_CHILD 2
#define BUCKET 4 

#define IS_BUCKET(node) (*((char *) (node)) & BUCKET)

#define SWAP(a, i, j) {\
    tmp = a[i];  \
    a[i] = a[j]; \
    a[j] = tmp;  \
}

#define TORUS_ADD(v, r, R) fmod(((v) + (r)),(R))
#define TORUS_SUBTRACT(v, r, R) fmod(((v) - (r) + (R)), R)
#define DIMENSIONS 2
#define X_DIMENSION 0
#define Y_DIMENSION 1

/*
 * Returns the log base 2 of the specified integer n >= 1. Finds the 
 * position of the right-most 1 in the binary representation. 
 */
static int 
binary_logarithm(const int n)
{
    int x = n;
    int ret = -1;
    while (x != 0) {
        ret++;
        x >>= 1;
    }
    return ret;
}


/*
 * Returns the dimension with the maximum spread accross the set of points from 
 * left to right inclusive.
 */
static inline int 
get_max_spread(point_t **points, int left, int right) 
{
    int i;
    double min_x = DBL_MAX; 
    double min_y = DBL_MAX; 
    double max_x = DBL_MIN; 
    double max_y = DBL_MIN;
    double *p; 
    for (i = left; i <= right; i++) {
        p = points[i]->location;
        min_x = GSL_MIN(p[X_DIMENSION], min_x);
        min_y = GSL_MIN(p[Y_DIMENSION], min_y);
        max_x = GSL_MAX(p[X_DIMENSION], max_x);
        max_y = GSL_MAX(p[Y_DIMENSION], max_y);
    }
    return max_x - min_x > max_y - min_y? X_DIMENSION: Y_DIMENSION;
}

/*
 * Permutes the input array of point_t pointers a such that a[k] contains
 * a point_t that is not greater in the d'th dimension than any point_t to its
 * left, and is similarly not greater than the points to its right.
 */
static inline void 
distribute(point_t **a, int i, int j, int k, 
        int d) 
{
    int left = i;
    int right = j;
    int p = 0;
    int pivot_index;
    int r;

    double pivot_value;
    point_t *tmp; /* required for the swap macro */
    while (k != p) {
        pivot_index = left == right ? right
                :(int) rand() % (right - left) + left; 
        SWAP(a, left, pivot_index); 
        pivot_value = a[left]->location[d];
        p = left;   
        for (r = left + 1; r <= right; ++r) {
            if (a[r]->location[d] < pivot_value) {
                ++p;
                SWAP(a, p, r);
            }       
        }
        SWAP(a, left, p);
        if (k < p) {
            right = p - 1;
        } else if (k > p) {
            left = p + 1;
        }
    }
}

/*
 * Initialises the specified list_set of keys.
 */
static int 
list_set_init(list_set *set, int max_keys) 
{
    int ret = 0;
    set->list = xmalloc(max_keys * sizeof(void *));
    set->num_keys = 0u;
    set->max_keys = max_keys;
    set->total_memory = (int) (sizeof(list_set) 
            + (max_keys * sizeof(void *)));
    return ret;
}

/*
 * Frees the resources held by the specified list_set. 
 */
static int 
list_set_free(list_set *set) 
{
    int ret = 0;
    free(set->list);
    return ret;
}

/*
 * Inserts the specified key into the set. 
 */
static inline int 
list_set_insert(list_set *set, void *key) 
{
    int ret = 0;
    int i;
    int j;
    void **a = set->list;
    int max_keys = (int) set->max_keys;
    int num_keys = (int) set->num_keys;
    if (num_keys == max_keys) {
        ret = -OUT_OF_LIST_SPACE;
        goto out;
    }
    j = 0;
    while (j < num_keys && a[j] < key) {
        j++;
    }
    if (j == num_keys) {
        a[j] = key;
        set->num_keys++;
    } else if (a[j] != key) {
        i = num_keys;
        while (i > j) {
            a[i] = a[i - 1];
            i--;
        }
        a[j] = key;
        set->num_keys++;
    }
out:
    return ret;
}

/*
 * Removes all keys from the specified list set.
 */ 
static inline int 
list_set_clear(list_set *set)
{
    set->num_keys = 0;
    return 0;
}

/*
 * Allocates and builds a new kdtree_t object on the specified maximum number of points 
 * into the specified kdtree_t object. 
 */
int 
kdtree_init(kdtree_t *self, int max_points, int bucket_max, 
        int random_seed)
{
    int num_buckets;
    int ret = 0;
    int i;
    kd_external_node *kdenp;
    kd_internal_node *kdinp;
    list_node *lnp;
    srand(random_seed);
    if (bucket_max == 0) {
        ret = -NOT_POWER_OF_TWO;
        num_buckets = 1; 
    } else {
        num_buckets = (int) (max_points < bucket_max? 1 :
                (max_points / bucket_max) * 2);    
    }
    self->node_set = NULL;
    self->search_stack = NULL;
    self->insertion_stack = NULL;
    if ((bucket_max & (bucket_max - 1)) != 0) {
        ret = -NOT_POWER_OF_TWO;
        /* continue through so that there is valid memory allocated */
        num_buckets = 1; 
    }
    self->max_points = max_points;
    self->bucket_max = bucket_max;
    self->depth = 0;
    self->root = NULL;
    self->total_memory = sizeof(kdtree_t);
    
    self->max_kd_external_node = (int) num_buckets;
    self->current_kd_external_node = num_buckets - 1;
    self->max_kd_internal_node = (int) num_buckets;
    self->current_kd_internal_node = num_buckets - 1;
    self->max_list_node = max_points + 2; 
    self->current_list_node = (int) self->max_list_node - 1; 
    
    self->node_set = xmalloc(sizeof(list_set));
    list_set_init(self->node_set, self->max_kd_external_node + 1);
    self->total_memory += self->node_set->total_memory;
    self->kd_external_node_heap = xmalloc(
            self->max_kd_external_node * sizeof(void *));
    self->kd_external_node_mem = xmalloc(
            self->max_kd_external_node * sizeof(kd_external_node));
    kdenp = self->kd_external_node_mem;
    for (i = 0; i < self->max_kd_external_node; i++) {
        self->kd_external_node_heap[i] = kdenp;
        kdenp++;
    } 
    self->kd_internal_node_heap = xmalloc(
            self->max_kd_internal_node * sizeof(void *));
    self->kd_internal_node_mem = xmalloc(
            self->max_kd_internal_node * sizeof(kd_internal_node));
    kdinp = self->kd_internal_node_mem;
    for (i = 0; i < self->max_kd_internal_node; i++) {
        self->kd_internal_node_heap[i] = kdinp;
        kdinp++;
    } 
    self->list_node_heap = xmalloc(self->max_list_node * sizeof(void *));
    self->list_node_mem = xmalloc(self->max_list_node * sizeof(list_node));
    lnp = self->list_node_mem;
    for (i = 0; i < self->max_list_node; i++) {
        self->list_node_heap[i] = lnp;
        lnp++;
    }
    self->total_memory += (int) (
            self->max_kd_external_node * sizeof(void *)
            + self->max_kd_external_node * sizeof(kd_external_node)
            + self->max_kd_internal_node * sizeof(void *)
            + self->max_kd_internal_node * sizeof(kd_internal_node)
            + self->max_list_node * sizeof(void *)
            + self->max_list_node * sizeof(list_node));
    return ret;
}

/* 
 * Clear out this kdtree_t, ready for build to be called again. 
 */
void 
kdtree_clear(kdtree_t *self)
{
    int i;
    kd_external_node *kdenp;
    kd_internal_node *kdinp;
    list_node *lnp;
    self->current_kd_external_node = self->max_kd_external_node - 1; 
    self->current_kd_internal_node = self->max_kd_internal_node - 1;
    self->current_list_node = (int) self->max_list_node - 1; 
    /* we must also rebuild the heaps, as the pointers can be anywhere. */
    kdenp = self->kd_external_node_mem;
    for (i = 0; i < self->max_kd_external_node; i++) {
        self->kd_external_node_heap[i] = kdenp;
        kdenp++;
    } 
    kdinp = self->kd_internal_node_mem;
    for (i = 0; i < self->max_kd_internal_node; i++) {
        self->kd_internal_node_heap[i] = kdinp;
        kdinp++;
    } 
    lnp = self->list_node_mem;
    for (i = 0; i < self->max_list_node; i++) {
        self->list_node_heap[i] = lnp;
        lnp++;
    }
    /* these are alloced in build, so we must free here to avoid a leak */
    free(self->search_stack);
    free(self->insertion_stack);
}

/*
 * Debug method to print the current state of the heaps.
 */
void 
kdtree_print_heaps(kdtree_t *self) 
{
    int i;
    printf("LIST NODE HEAP\n");
    for (i = (int) self->max_list_node - 1; i >= 0; i--) {
        if (i == self->current_list_node) {
            printf(" == ");
        }
        printf("\t%p\n", (void *) self->list_node_heap[i]);    
    }
    printf("INTERNAL NODE HEAP\n");
    for (i = (int) self->max_kd_internal_node - 1; i >= 0; i--) {
        if (i == self->current_kd_internal_node) {
            printf(" == ");
        }
        printf("\t%p\n", (void *) self->kd_internal_node_heap[i]);    
    }
    printf("EXTERNAL NODE HEAP\n");
    for (i = (int) self->max_kd_external_node - 1; i >= 0; i--) {
        if (i == self->current_kd_external_node) {
            printf(" == ");
        }
        printf("\t%p\n", (void *) self->kd_external_node_heap[i]);    
    }
}


/*
 * Allocates a new list node from the heap.
 */
static inline list_node* 
kdtree_alloc_list_node(kdtree_t *self)
{
    list_node *node = NULL;
    if (self->current_list_node >= 0) { 
        node = self->list_node_heap[self->current_list_node];    
        self->current_list_node--;
        node->next = NULL;
    }  
    return node;
}

static inline void 
kdtree_free_list_node(kdtree_t *self, list_node *node)
{
    self->current_list_node++;
    self->list_node_heap[self->current_list_node] = node;
}

/*
 * Allocates a new kdtree_t external node from the heap.
 */
static inline kd_external_node* 
kdtree_alloc_external_node(kdtree_t *self)
{
    kd_external_node *node = NULL;
    if (self->current_kd_external_node >= 0) { 
        node = self->kd_external_node_heap[self->current_kd_external_node];    
        self->current_kd_external_node--;
        node->flags = BUCKET;
        node->head = NULL;
    }    
    return node;
}

static inline void 
kdtree_free_kd_external_node(kdtree_t *self, kd_external_node *node)
{
    self->current_kd_external_node++;
    self->kd_external_node_heap[self->current_kd_external_node] = node;
}

/*
 * Allocates a new kdtree_t internal node from the heap.
 */
static inline kd_internal_node* 
kdtree_alloc_internal_node(kdtree_t *self)
{
    kd_internal_node* node = NULL;
    if (self->current_kd_internal_node >= 0) { 
        node = self->kd_internal_node_heap[self->current_kd_internal_node];    
        self->current_kd_internal_node--;
        node->flags = 0;
        node->left_child = NULL;
        node->right_child = NULL;
    }    
    return node;
}

static inline void 
kdtree_free_kd_internal_node(kdtree_t *self, kd_internal_node *node)
{
    self->current_kd_internal_node++;
    self->kd_internal_node_heap[self->current_kd_internal_node] = node;
}

/* 
 * Frees all resources used by the specified self.
 */
void 
kdtree_free(kdtree_t *self) 
{
    if (self->search_stack != NULL) {
        free(self->search_stack);
    }
    if (self->insertion_stack != NULL) {
        free(self->insertion_stack);
    }
    list_set_free(self->node_set);
    free(self->node_set);
    free(self->kd_external_node_heap);
    free(self->kd_internal_node_heap);
    free(self->list_node_heap);
    free(self->kd_external_node_mem);
    free(self->kd_internal_node_mem);
    free(self->list_node_mem);
}

/*
 * Prints out the specified kdtree_t rooted at the specified node and recursively
 * prints out the children.
 */
static void 
kdtree_print_node(kdtree_t *self, void *node, int depth, int label) {
    
    kd_internal_node *itnode;
    kd_external_node *etnode;
    list_node *next_lnode;
    int i;
    int n;
    for (i = 0; i < depth; i++) {
        printf("  ");
    }
    printf("(%d)", label);
    if (IS_BUCKET(node)) {
        etnode = node;
        printf("bucket (%d) [", etnode->num_points);
        next_lnode = etnode->head;
        n = 0;
        while (next_lnode != NULL) {
            n++;
            printf("(%p), ", next_lnode->point); 
            next_lnode = next_lnode->next;
        }
        printf("]\n");
        if (n != etnode->num_points) {
            printf("ERROR!!!! %d != %d \n", n, etnode->num_points);
        }
    } else {
        itnode = node;
        printf("internal v=%.2f d=%d\n", itnode->cut_value, itnode->cut_dimension);
        kdtree_print_node(self, itnode->left_child, depth + 1, 2 * label);
        kdtree_print_node(self, itnode->right_child, depth + 1, 2 * label + 1);
    }
}

void 
kdtree_print(kdtree_t *self) 
{
    kdtree_print_node(self, self->root, 0, 0);
}


/*
 * Builds a new kdtree rooted at the specified node (null if a new tree is 
 * required).
 */
int 
kdtree_build(kdtree_t *self, kd_internal_node *root, point_t **points, 
        const int num_points)
{
    int ret = 0;
    int s = 0;
    int left, right;
    int left_child, right_child;
    int i;
    int d, m;
    double v;
    void *node;
        
    kd_internal_node *parent;
    kd_internal_node *itnode;
    kd_external_node *etnode;
    list_node *next_lnode;
    list_node *prev_lnode;

    typedef struct stack_entry_t {
        kd_internal_node *parent;
        int left;
        int right;
        int flags;
    } stack_entry;
    int stack_size = GSL_MAX(128, binary_logarithm(num_points) + 1);
    stack_entry *stack = xmalloc(stack_size * sizeof(stack_entry));
    stack_entry *se;
    se = &stack[0];
    se->parent = root;
    se->left = 0;
    se->right = num_points - 1;
    se->flags = 0;
    while (s >= 0) {
        assert(s < (int) stack_size);
        se = &stack[s];                    
        self->depth = (int) s > self->depth? 
                (int) s: self->depth;
        s--;
        parent = se->parent;
        left = se->left;
        right = se->right;
        left_child = se->flags & LEFT_CHILD; 
        right_child = se->flags & RIGHT_CHILD; 
        if (right - left <= self->bucket_max) {
            etnode = kdtree_alloc_external_node(self); 
            if (etnode == NULL) {
                ret = -OUT_OF_KDT_EXTERNAL_NODES;
                goto cleanup;
            }
            i = left;
            prev_lnode = NULL;
            next_lnode = kdtree_alloc_list_node(self);
            if (next_lnode == NULL) {
                ret = -OUT_OF_LIST_NODES;
                goto cleanup;
            }
            etnode->head = next_lnode;
            while (i <= right) {
                next_lnode->point = points[i];
                prev_lnode = next_lnode;
                next_lnode = kdtree_alloc_list_node(self);
                if (next_lnode == NULL) {
                    ret = -OUT_OF_LIST_NODES;
                    goto cleanup;
                }
                prev_lnode->next = next_lnode;
                i++;
            }
            if (prev_lnode != NULL) {
                prev_lnode->next = NULL;
            }
            kdtree_free_list_node(self, next_lnode);
            etnode->num_points = right - left + 1;
            if (etnode->num_points == 0) {
                etnode->head = NULL;
            }
            node = etnode;
        } else {
            itnode = kdtree_alloc_internal_node(self);
            if (itnode == NULL) {
                ret = -OUT_OF_KDT_INTERNAL_NODES;
                goto cleanup;
            }
            d = get_max_spread(points, left, right);
            m = left + (right - left) / 2;
            distribute(points, left, right, m, d);
            v = points[m]->location[d];
            itnode->cut_value = v;
            itnode->cut_dimension = d;
            
            s++;
            se = &stack[s];
            se->parent = itnode;
            se->left = left;
            se->right = m - 1;
            se->flags = LEFT_CHILD;

            s++;
            se = &stack[s];
            se->parent = itnode;
            se->left = m;
            se->right = right;
            se->flags = RIGHT_CHILD;
            node = itnode;
            
        }
        if (parent == NULL) {
            self->root = node; 
        } else {
            if (left_child) {
                parent->left_child = node;
            } else if (right_child) {
                parent->right_child = node;
            }
        }
    }
    self->depth++;
    self->search_stack = xmalloc(self->depth * sizeof(void *));
    self->insertion_stack = xmalloc(self->depth 
            * sizeof(insertion_stack_entry));
cleanup:
    free(stack);
    return ret;
}

/*
 * Finds all bucket nodes in the specified kdtree_t that are within the specified bounds,
 * and adds these nodes to the self's node_set.
 */
static int 
kdtree_get_buckets_in_region(kdtree_t *self, double *min_bounds, double *max_bounds) 
{
    void *node;
    kd_internal_node *itnode;
    int d;
    double v;
    void **stack = self->search_stack;
    int s = 0;
    int ret = 0; 
    stack[0] = self->root;
    while (s >= 0) {
        node = stack[s];                    
        s--;
        if (IS_BUCKET(node)) {
            ret = list_set_insert(self->node_set, node);
            ERCS_ERROR_CHECK(ret, out);
        } else {
            itnode = node;
            d = itnode->cut_dimension;
            v = itnode->cut_value;
            if (v < min_bounds[d]) {
                s++;
                stack[s] = itnode->right_child;
            } else if(min_bounds[d] <= v && v <= max_bounds[d]) {
                s++;
                stack[s] = itnode->right_child;
                s++;
                stack[s] = itnode->left_child;
            } else {
                s++;
                stack[s] = itnode->left_child;
            }
        }
    }
out:   
    return ret;
}

/*
 * Helper function: set the bounds appropriately.
 */
static inline void 
set_bounds(double *min_bounds, double *max_bounds, double x_min, 
        double y_min, double x_max, double y_max) 
{
        min_bounds[X_DIMENSION] = x_min;
        min_bounds[Y_DIMENSION] = y_min;
        max_bounds[X_DIMENSION] = x_max;
        max_bounds[Y_DIMENSION] = y_max;
}

/*
 * Adds all external nodes that are within the specified radius of the specified point_t 
 * on a square torus with edge R to the node set of the kdtree.
 */
static int 
kdtree_get_buckets_in_torus_radius(kdtree_t *self, const double *p, const double r, 
        const double torus_edge)
{ 
    int ret = 0;
    double R = torus_edge;
    double x = p[X_DIMENSION];
    double y = p[Y_DIMENSION];
    int x_middle_full = (x + r < R) && (x - r > 0);
    int y_middle_full = (y + r < R) && (y - r > 0);
    double min_bounds[DIMENSIONS];
    double max_bounds[DIMENSIONS];
    if (r >= R / 2) {
        /* just look at the whole thing - it's not worth finding the holes. */
        set_bounds(min_bounds, max_bounds, 0.0, 0.0, R, R);
        ret = kdtree_get_buckets_in_region(self, min_bounds, max_bounds);
        ERCS_ERROR_CHECK(ret, out);
    } else if (x_middle_full && y_middle_full) {
        set_bounds(min_bounds, max_bounds, x - r, y - r, x + r, y + r);
        ret = kdtree_get_buckets_in_region(self, min_bounds, max_bounds);
        ERCS_ERROR_CHECK(ret, out);
    } else if (x_middle_full && !y_middle_full) {
        set_bounds(min_bounds, max_bounds, x - r, 0.0, x + r, TORUS_ADD(y, r, R));
        ret = kdtree_get_buckets_in_region(self, min_bounds, max_bounds);
        ERCS_ERROR_CHECK(ret, out);
        set_bounds(min_bounds, max_bounds, x - r, TORUS_SUBTRACT(y, r, R), x + r, R); 
        ret = kdtree_get_buckets_in_region(self, min_bounds, max_bounds);
        ERCS_ERROR_CHECK(ret, out);
    } else if (y_middle_full && !x_middle_full) {
        set_bounds(min_bounds, max_bounds, 0.0, y - r, TORUS_ADD(x, r, R), y + r);
        ret = kdtree_get_buckets_in_region(self, min_bounds, max_bounds);
        ERCS_ERROR_CHECK(ret, out);
        set_bounds(min_bounds, max_bounds, TORUS_SUBTRACT(x, r, R), y - r, R, y + r); 
        ret = kdtree_get_buckets_in_region(self, min_bounds, max_bounds);
        ERCS_ERROR_CHECK(ret, out);
    } else {
        /* Define our four regions to search, defined clockwise from origin */
        set_bounds(min_bounds, max_bounds, 0.0, 0.0, TORUS_ADD(x, r, R), TORUS_ADD(y, r, R)); 
        ret = kdtree_get_buckets_in_region(self, min_bounds, max_bounds);
        ERCS_ERROR_CHECK(ret, out);
        set_bounds(min_bounds, max_bounds, 0.0, TORUS_SUBTRACT(y, r, R), TORUS_ADD(x, r, R), R); 
        ret = kdtree_get_buckets_in_region(self, min_bounds, max_bounds);
        ERCS_ERROR_CHECK(ret, out);
        set_bounds(min_bounds, max_bounds, TORUS_SUBTRACT(x, r, R), TORUS_SUBTRACT(y, r, R), R, R); 
        ret = kdtree_get_buckets_in_region(self, min_bounds, max_bounds);
        ERCS_ERROR_CHECK(ret, out);
        set_bounds(min_bounds, max_bounds, TORUS_SUBTRACT(x, r, R), 0.0, R, TORUS_ADD(y, r, R)); 
        ret = kdtree_get_buckets_in_region(self, min_bounds, max_bounds);
        ERCS_ERROR_CHECK(ret, out);
    }
out:    
    return ret;
}

/*
 * Inserts all points between left and right inclusive into the specified 
 * kd_external_node;
 */
static inline int 
kdtree_external_node_insert_points(kdtree_t *self, kd_external_node *tree_node, 
        point_t **points, int left, int right)
{
    int ret = 0;
    list_node *last = tree_node->head;
    list_node *next = kdtree_alloc_list_node(self);
    list_node *prev = NULL;
    int i = left;
    tree_node->head = next;
    if (next == NULL) {
        ret = -OUT_OF_LIST_NODES;
        goto out;
    }    
    while (i <= right) {
        next->point = points[i];
        prev = next;
        next = kdtree_alloc_list_node(self);
        if (next == NULL) {
            ret = -OUT_OF_LIST_NODES;
            goto out;
        }    
        prev->next = next;
        i++;
    }
    if (prev != NULL) {
        prev->next = last;
    } else {
        tree_node->head = last;
    }
    kdtree_free_list_node(self, next);
    tree_node->num_points += right - left + 1; 
    
out:
    return ret;
}

/*
 * Inserts the specified point_t into the kdtree_t.
 */
int 
kdtree_insert_point(kdtree_t *self, point_t *ind) 
{
    point_t *points[] = {ind};
    return kdtree_insert_points(self, points, 1);
}

/*
 * Inserts the specified set of points into the specified kdtree_t.
 */
int 
kdtree_insert_points(kdtree_t *self, point_t **points, 
        int num_points)
{
    int ret = 0;
    insertion_stack_entry *se;
    int d;
    double v;
    int pivot, left, right;
    int s = 0;
    int j;
    point_t *tmp;
    void *node;
    kd_internal_node *itnode;
    kd_external_node *etnode;
    insertion_stack_entry *stack = self->insertion_stack; 
    if (num_points == 0) {
        goto out;
    }
    /* traverse the entire set of points */
    se = &stack[s];                    
    se->node = self->root; 
    se->left = 0;
    se->right= num_points - 1;
    while (s >= 0) {
        se = &stack[s];                    
        s--;
        node = se->node;
        left = se->left;
        right = se->right;
        if (IS_BUCKET(node)) {
            etnode = (kd_external_node *) node; 
            ret = kdtree_external_node_insert_points(self, etnode, points, 
                    left, right);
            ERCS_ERROR_CHECK(ret, out);
       } else {
            itnode = (kd_internal_node *) node; 
            d = itnode->cut_dimension;
            v = itnode->cut_value;
            /* Permutes the points such that for all left <= j < pivot
             * we have a[j][d] <= v and pivot <= j <= right we have a[j][d] > v.
             */
            pivot = left;
            for (j = left; j <= right; j++) {
                if (points[j]->location[d] <= v) {
                    SWAP(points, pivot, j);
                    pivot += 1;
                }
            }
            if (pivot > left) {
                s++;
                se = &stack[s];                    
                se->node = itnode->left_child;
                se->left = left;
                se->right = pivot - 1;
            }
            if (pivot <= right) {
                s++;
                se = &stack[s];                    
                se->node = itnode->right_child;
                se->left = pivot;
                se->right= right;
            }    
        }
    }
out:
    return ret;
}


int 
kdtree_copy_points(kdtree_t *self, point_t **points) 
{
    int n = 0;
    int ret;
    double min_bounds[2] = {0.0, 0.0};
    double max_bounds[2] = {DBL_MAX, DBL_MAX};
    void **a;
    int i;
    kd_external_node *tnode;
    list_node *next_lnode;
    ret = kdtree_get_buckets_in_region(self, min_bounds, max_bounds); 
    ERCS_ERROR_CHECK(ret, out);
    a = self->node_set->list;
    for (i = 0; i < self->node_set->num_keys; i++) {
        tnode = (kd_external_node *) a[i];
        next_lnode = (list_node *) tnode->head;
        while (next_lnode != NULL) {
            points[n] = next_lnode->point;
            next_lnode = next_lnode->next;
            n++;
        }
    }    
    ret = n;
out:
    list_set_clear(self->node_set);
    return ret;
}
/*
 * Increments the iter->current_tree_node_index until a non-empty node
 * is found, and sets the current_list_node to the head of this tree node's
 * list. The current_list_node is set to NULL if no non empty nodes 
 * are left.
 */
static inline void
kri_goto_next_tree_node(kri_t *iter)
{
    int i = iter->current_tree_node_index;
    int n = (int) iter->node_set->num_keys - 1;
    kd_external_node *tnode;
    list_node *next = NULL;
    while (next == NULL && i < n) {
        i++;
        tnode = (kd_external_node *) iter->node_set->list[i];
        next = (list_node *) tnode->head;
    }
    iter->current_tree_node_index = i;
    iter->previous_list_node = NULL;
    iter->current_list_node = next;
}

int 
kdtree_get_torus_region_iterator(kdtree_t *self, const double *z, const double r, 
        const double R, kri_t *iter)
{
    int ret = 0;
    double p[] = {0.0, 0.0};
    p[0] = z[0];
    p[1] = z[1];
    list_set_clear(self->node_set);
    ret = kdtree_get_buckets_in_torus_radius(self, p, r, R);
    ERCS_ERROR_CHECK(ret, out);
    iter->delete_possible = 0;
    iter->kdtree = self;
    iter->node_set = self->node_set;
    iter->previous_list_node = NULL;
    iter->current_list_node = NULL;
    iter->current_tree_node_index = -1;
    iter->deleted_list_node = NULL;
out:
    return ret;
}

point_t* 
kri_next(kri_t *iter) 
{
    point_t *ret = NULL;
    if (iter->deleted_list_node == NULL) {
        iter->previous_list_node = iter->current_list_node;
    } 
    if (iter->current_list_node == NULL) {
        kri_goto_next_tree_node(iter);
    } else { 
        iter->current_list_node = iter->current_list_node->next;
        if (iter->current_list_node == NULL) {
            kri_goto_next_tree_node(iter);
        }
    }
    if (iter->current_list_node != NULL) {
        ret = iter->current_list_node->point;
    }
    if (iter->deleted_list_node != NULL) {
        kdtree_free_list_node(iter->kdtree, iter->deleted_list_node);
        iter->deleted_list_node->next = NULL;
        iter->deleted_list_node = NULL;
    } 
    iter->delete_possible = ret != NULL;
    return ret;
}
/*
 * Deletes the point_t last returned by the specified iterator. Each call
 * to delete must be preceded by at lease one call to next.
 */
int 
kri_delete(kri_t *iter) {
    int ret = 0;
    kd_external_node *tnode;
    list_node *deleted, *prev;
    int i = iter->current_tree_node_index; 
    if (!iter->delete_possible) {
        ret = -ITERATOR_ERROR;
        goto out;
    }
    iter->delete_possible = 0;
    tnode = (kd_external_node *) iter->node_set->list[i];
    prev = iter->previous_list_node;
    deleted = iter->current_list_node;
    if (prev == NULL) {
        tnode->head = deleted->next;
    } else {
        prev->next = deleted->next;
    }
    tnode->num_points--;
    iter->deleted_list_node = deleted;
out:
    return ret;
}


