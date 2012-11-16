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
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>

/* Euclidean squared distance */
static double 
esd(double *x, double *y)
{
    return gsl_pow_2(y[0] - x[0]) + gsl_pow_2(y[1] - x[1]);
}

/* Alternative version of the squared torus distance */
static double
alternative_tsd(double *a, double *b, double R)
{
    double x[] = {0.0, 0.0};
    double origin[] = {0.0, 0.0};
    /* translate a to the origin, and reflect so that x is positive */
    x[0] = fabs(b[0] - a[0]);
    x[1] = fabs(b[1] - a[1]);
    /* wrap around so that x is as near to the origin as possible */
    x[0] = GSL_MIN_DBL(x[0], R - x[0]);
    x[1] = GSL_MIN_DBL(x[1], R - x[1]); 
    return esd(origin, x); 
}

/*
 * Heaviside Lambda function
 */
static double
heaviside_lambda(double x)
{
    double x_abs = fabs(x);
    return x_abs >= 1 ? 0 : 1.0 - x_abs;
}

/* 
 * Squared distance between two points x and y on a circle of length L.
 */
static double 
squared_circle_distance(double x, double y, double L) 
{
    double z = L * heaviside_lambda(2 * fabs(x - y) / L - 1) / 2;
    return z * z;
}



/*
 * Squared distance on a torus of side L using the Heaviside Lambda function.
 */
static double
heaviside_tsd(double *x, double *y, double L)
{
    return squared_circle_distance(x[0], y[0], L) 
            + squared_circle_distance(x[1], y[1], L);
}


static int 
test_torus_distance_known(void)
{
    int ret = 0;
    double p1[2] = {0.0, 0.0};
    double p2[2] = {10.0, 10.0};
    int num_tests = 13;
    double R[13] = {10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
            20.0, 21.0, 22.0};
    double result[13] = {0.0, 2.0, 8.0, 18.0, 32.0, 50.0, 72.0, 98.0, 128.0, 
            162.0, 200.0, 200.0, 200.0};
    double d21, d22; 
    int i;

    for (i = 0; i < num_tests; i++) {
        d21 = torus_squared_distance(p1, p2, R[i]);
        if (gsl_fcmp(d21, result[i], 1E-8) != 0) {
            ret = -601;
            goto out;
        }
        d22 = alternative_tsd(p1, p2, R[i]);
        if (gsl_fcmp(d22, result[i], 1E-8) != 0) {
            printf("%f:%f %f\n", R[i], d22, result[i]);
            ret = -603;
            goto out;
        }
        d22 = heaviside_tsd(p1, p2, R[i]);
        if (gsl_fcmp(d22, result[i], 1E-8) != 0) {
            printf("%f:%f %f\n", R[i], d22, result[i]);
            ret = -604;
            goto out;
        }


    }     
out:
    return ret;
    
}

static int 
test_grid(double *p1, double R, int grid_size)
{
    int ret = 0;
    int i, j;
    double d21, d22;
    double p2[2] = {0.0, 0.0};
    for (i = 0; i < grid_size; i++) {
        p2[0] = i * R / grid_size;
        for (j = 0; j < grid_size; j++) {
            p2[1] = j * R / grid_size;
            d21 = torus_squared_distance(p1, p2, R);
            d22 = alternative_tsd(p1, p2, R);
            if (gsl_fcmp(d21, d22, 1E-8) != 0) {
                ret = -605;
                goto out;
            }
            d22 = heaviside_tsd(p1, p2, R);
            if (gsl_fcmp(d21, d22, 1E-8) != 0) {
                ret = -606;
                goto out;
            }
        }
    }
out:
    return ret;
}

static int 
test_torus_distance_systematic(void)
{
    int ret = 0;
    double p1[2] = {0.0, 0.0};
    double R = 8;
    int grid_size = 32;
    test_grid(p1, R, grid_size);    
    p1[0] = 4;
    p1[1] = 4;
    test_grid(p1, R, grid_size);    
    p1[0] = 7;
    p1[1] = 7;
    test_grid(p1, R, grid_size);    
    p1[0] = 7.99;
    p1[1] = 7.99;
    test_grid(p1, R, grid_size);    
    p1[0] = 0.11;
    p1[1] = 7.99;
    test_grid(p1, R, grid_size);    
    for (R = 2.0; R < 16.0; R += 1.0) {
        p1[0] = 0.0;
        p1[1] = 0.0;
        test_grid(p1, R, grid_size);    
    }
    return ret;
    
}


int 
main(void)
{
    int result;
    printf("\trunning torus tests");
    result = test_torus_distance_known();
    if (result < 0) {
        printf("Error in known torus_distance tests: %d\n", result);
        exit(result);
    }
    printf(".");
    fflush(stdout);
    result = test_torus_distance_systematic();
    if (result < 0) {
        printf("Error in known torus_distance tests: %d\n", result);
        exit(result);
    }
    printf(".");
    printf("\n"); 
    return EXIT_SUCCESS;

}
