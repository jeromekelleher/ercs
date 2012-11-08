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

#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdlib.h>


#ifdef DEBUG_ERRORS
#define ERCS_ERROR_CHECK(r, label) \
    printf("ERROR_CHECK: %d: %d: %s\n", r, __LINE__, __FILE__); \
    if (r < 0) { goto label;}
#else
#define ERCS_ERROR_CHECK(r, label) \
    if (r < 0) { goto label;}
#endif

void * xmalloc(size_t size);
void * xcalloc(size_t count, size_t eltsize);
void * xrealloc(void *ptr, size_t size);


#endif /*__UTIL_H__*/
