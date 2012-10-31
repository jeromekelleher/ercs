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

/*
 * Collection of miscellaneous functions shared throughout source.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "torus.h"

#include <assert.h>
#include <math.h>
#include <gsl/gsl_math.h>

/*
 * Returns the square of the distance between the two specified points on a 
 * square torus of side R.
 */
double 
torus_squared_distance(const double *p1, const double *p2, double R) 
{
    double xabs = fabs(p2[0] - p1[0]);
    double yabs = fabs(p2[1] - p1[1]);
    double xd = GSL_MIN_DBL(xabs, R - xabs);
    double yd = GSL_MIN_DBL(yabs, R - yabs);
    return xd * xd + yd * yd;
}

/*
 * Wraps the specified point around such that it lies on a torus of side R.
 */
void 
torus_wrap(double *p, double R)
{
    p[0] = fmod(p[0] + R, R);
    p[1] = fmod(p[1] + R, R); 
    assert(0 <= p[0] && p[0] < R && 0 <= p[1] && p[1] < R);

}

