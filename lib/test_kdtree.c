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

#include "util.h"
#include "torus.h"
#include "kdtree.h"

#include <stdio.h>
#include <stdint.h>
#include <gsl/gsl_math.h>

/*
 * Fills in the specified point with values in the range 0 .. R
 */
static void 
random_point_rand(double *p, double R) 
{
    p[0] = rand() / (((double) RAND_MAX + 1) / R); 
    p[1] = rand() / (((double) RAND_MAX + 1) / R); 
}

/*
 * Returns a random integer from l to h inclusive.
 */
static int 
random_int(int l, int h) {
    return l + (int) (rand() / (((double) RAND_MAX + 1) / h));
}


/*
 * Comparison function on pointers for library functions
 */
static int 
pointer_compare(const void *ip1, const void *ip2) 
{
    const intptr_t *a = ((const intptr_t *) ip1);
    const intptr_t *b = ((const intptr_t *) ip2);
    return (*a > *b) - (*a < *b);
}

static int
compare_point_lists(kdtree_t *tree, point_t **f1, point_t **f2, 
        int num_found1, int num_found2, int error_code)
{
    int error = 0;
    int ret = 0;
    int i;
    if (num_found2 != num_found1) {
        printf("ERROR %d != %d\n", num_found2, num_found1);
        error = 1;
    } else {
        qsort(f1, (size_t) num_found1, sizeof(point_t *), &pointer_compare); 
        qsort(f2, (size_t) num_found1, sizeof(point_t *), &pointer_compare); 
        for (i = 0; i < num_found2; i++) {
            if ((! error) && f1[i] != f2[i]) {
                printf("ERROR!!\n");
                error = 1;
            } 
        }
    }
    if (error) { 
        for (i = 0; i < GSL_MAX_INT(num_found2, num_found1); i++) {
            printf("f1[%d] = %p \tf2[%d] = %p \t%d\n", 
                    i, f1[i], i, f2[i], f1[i] == f2[i]);
        }
        kdtree_print(tree);
        ret = error_code;
    }   
    return ret;
}
/*
 * Tests to see if the region iterator will correctly return all points in the 
 * tree.
 */
static int 
kdtree_test_region_iterator(kdtree_t *tree, int num_points, double R) 
{
    int ret = 0;
    int i, j, u;
    point_t **f1 = xcalloc((size_t) num_points, sizeof(point_t *));
    point_t **f2 = xcalloc((size_t) num_points, sizeof(point_t *));
    kri_t *iter = xmalloc(sizeof(kri_t)); 
    point_t *ind;
    double p[2] = {0.0, 0.0};
    for (i = num_points; i > 0; i--) {
        ret = kdtree_copy_points(tree, f1);
        ERCS_ERROR_CHECK(ret, out);
        ret = kdtree_get_torus_region_iterator(tree, p, R, R, iter);
        ERCS_ERROR_CHECK(ret, out);
        j = 0;
        while ((ind = kri_next(iter)) != NULL) {
            f2[j] = ind;
            j++;
        }
        ret = compare_point_lists(tree, f1, f2, i, j, -124);
        ERCS_ERROR_CHECK(ret, out);
        /* delete a random element of the list */
        u = (int) random_int(0, (int) i - 1); 
        ret = kdtree_get_torus_region_iterator(tree, p, R, R, iter);
        ERCS_ERROR_CHECK(ret, out);
        j = 0;
        while ((ind = kri_next(iter)) != NULL && j < u) {
            j++;
        }
        ret = kri_delete(iter);
        ERCS_ERROR_CHECK(ret, out);
        while ((ind = kri_next(iter)) != NULL) {
            j++;
        }
        if (j != i - 1) {
            printf("Error in iterator: not completed correctly: %d %d\n", j, i);
            ret = -125;
            goto out;
        }

    }
out:
    free(f1);
    free(f2);
    free(iter);
    return ret;
}
/*
 * Tests the specified kdtree to see if all points are correctly found by doing an 
 * exhaustive search.
 */
static int 
kdtree_test_torus_search(kdtree_t *tree, double *p, double r, int R, 
        int num_points, point_t **points) 
{

    int ret = 0;
    int i;
    point_t **f1 = xcalloc((size_t) tree->max_points, sizeof(point_t *));
    point_t **f2 = xcalloc((size_t) tree->max_points, sizeof(point_t *));
    point_t *ind;
    int num_found1 = 0;
    int num_found2 = 0;
    double r2 = r * r;
    kri_t *iter = xmalloc(sizeof(kri_t)); 
    for (i = 0; i < num_points; i++) {
        ind = points[i];
        if (torus_squared_distance(p, ind->location, (double) R) <= r2) {
            f1[num_found1] = ind;
            num_found1++;
        } 
    }
    ret = kdtree_get_torus_region_iterator(tree, p, r, (double) R, iter);
    ERCS_ERROR_CHECK(ret, out);
    i = 0;
    while ((ind = kri_next(iter)) != NULL) {
        if (torus_squared_distance(p, ind->location, (double) R) <= r2) {
            f2[i] = ind;
            i++;
        }
    }
    num_found2 = i;
    ret = compare_point_lists(tree, f1, f2, num_found1, num_found2, -123);
out:
    free(f1);
    free(f2);
    free(iter);
    return ret;

}

static int 
kdtree_test_search(kdtree_t *tree, int R, int num_points, 
        point_t **points) 
{
    int ret = 0;
    int i;
    double r; 
    const int num_radii= 10;
    double p[2] = {0.0, 0.0};
    for (i = 0; i < num_radii; i++) {
        r = 0.1 + ((R / 2) / ((double) num_radii)) * i;
        random_point_rand(p, (double) R);
        ret = kdtree_test_torus_search(tree, p, r, R, num_points, 
                points); 
        ERCS_ERROR_CHECK(ret, out);
    } 
    random_point_rand(p, (double) R);
    ret = kdtree_test_torus_search(tree, p, (R / 2.0) + 1, R, num_points,
            points); 
    ERCS_ERROR_CHECK(ret, out);
    ret = kdtree_test_torus_search(tree, p, (double) R, R, num_points,
            points); 
    ERCS_ERROR_CHECK(ret, out);
out:            
    return ret;
}

static int 
test_kdtree_size(int R, int num_points, 
        int bucket_size)
{
    int ret = 0;
    int i;
    kdtree_t *tree = xmalloc(sizeof(kdtree_t));
    point_t **points = xmalloc(num_points * sizeof(point_t *));
    ret = kdtree_init(tree, num_points, bucket_size, 0); 
    ERCS_ERROR_CHECK(ret, out);
    for (i = 0; i < num_points; i++) {
        points[i] = xmalloc(sizeof(point_t)); 
        random_point_rand(points[i]->location, (double) R);
    }   
    ret = kdtree_build(tree, tree->root, points, num_points); 
    ERCS_ERROR_CHECK(ret, out);
    /* OK, now we can do stuff */
    ret = kdtree_test_search(tree, R, num_points, points);
    ERCS_ERROR_CHECK(ret, out);
    ret = kdtree_test_region_iterator(tree, num_points, (double) R);
    ERCS_ERROR_CHECK(ret, out);
out:
    kdtree_free(tree);
    for (i = 0; i < num_points; i++) {
        free(points[i]);
    }
    free(points);    
    free(tree);
    return ret;
}

static int 
test_kdtree(void)
{
    int ret = 0;
    size_t j, k, L, n;
    int torus_sizes[] = {1, 2, 5, 123, 1111};
    int num_points[] = {1, 2, 11, 50, 212, 333, 515};
    printf("\trunning kdtree tests");
    for (j = 0; j < sizeof(torus_sizes) / sizeof(int); j++) {
        L = torus_sizes[j];
        for (k = 0; k < sizeof(num_points) / sizeof(int); k++) {
            printf(".");
            fflush(stdout);
            n = num_points[k]; 
            ret = test_kdtree_size(L, n, 64);
            ERCS_ERROR_CHECK(ret, out);
            ret = test_kdtree_size(L, n, 8);
            ERCS_ERROR_CHECK(ret, out);
            ret = test_kdtree_size(L, n, 4);
            ERCS_ERROR_CHECK(ret, out);
            ret = test_kdtree_size(L, n, 1);
            ERCS_ERROR_CHECK(ret, out);
        }
    }  
    printf("\n");

out:
    return ret;

}

int 
main(void)
{
    int result;
    result = test_kdtree();
    if (result < 0) {
        printf("Error in kdtree tests: %d\n", result);
        exit(result);
    }
    return EXIT_SUCCESS;

}
