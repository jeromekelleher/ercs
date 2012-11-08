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

#include "ercs.h"
#include <stdio.h>


static ercs_t *
alloc_ercs(double L, unsigned int num_points) 
{
    ercs_t *sim = xmalloc(sizeof(ercs_t));
    sim->torus_diameter = L;
    sim->sample_size = num_points;
    return sim;
}

static int 
test_ercs_size(double L, unsigned int num_points)
{
    int ret = 0;
    ercs_t *sim = alloc_ercs(L, num_points);
    if (sim == NULL) {
        ret = -111;
        goto out;
    }
out:
    return ret;
}




static int 
test_ercs(void)
{
    int ret = 0;
    unsigned int j, k, n;
    double torus_sizes[] = {1.9, 0.2, 50.3, 123.22, 1111.1};
    unsigned int num_points[] = {1, 2, 11, 50, 212, 333, 515};
    printf("\trunning ercs tests");
    for (j = 0; j < sizeof(torus_sizes) / sizeof(double); j++) {
        for (k = 0; k < sizeof(num_points) / sizeof(unsigned int); k++) {
            printf(".");
            fflush(stdout);
            n = num_points[k]; 
            ret = test_ercs_size(torus_sizes[j], n);
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
    result = test_ercs();
    if (result < 0) {
        printf("Error in ercs tests: %d\n", result);
        exit(result);
    }
    return EXIT_SUCCESS;
}
