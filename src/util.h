/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#ifndef __UTIL_H_
#define __UTIL_H_
#define HLINE                                                                            \
    "----------------------------------------------------------------------------\n"

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif
#ifndef ABS
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#endif

#define P(i, j)   p[(j) * (imax + 2) + (i)]
#define F(i, j)   f[(j) * (imax + 2) + (i)]
#define G(i, j)   g[(j) * (imax + 2) + (i)]
#define U(i, j)   u[(j) * (imax + 2) + (i)]
#define V(i, j)   v[(j) * (imax + 2) + (i)]
#define RHS(i, j) rhs[(j) * (imax + 2) + (i)]

#endif // __UTIL_H_
