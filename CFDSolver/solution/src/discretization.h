/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#ifndef __DISCRETIZATION_H_
#define __DISCRETIZATION_H_
#include "grid.h"
#include "parameter.h"


typedef struct {
    /* geometry and grid information */
    Grid grid;
    /* arrays */
    double *p, *rhs;
    double *f, *g;
    double *u, *v;
    /* parameters */
    double rho;
    double re, tau, gamma;
    double gx, gy;
    /* time stepping */
    double dt, te;
    double dtBound;
    char* problem;
    int bcLeft, bcRight, bcBottom, bcTop;
} Discretization;

extern void initDiscretization(Discretization*, Parameter*);
extern void setObjectBoundaryCondition(Discretization*);
extern void print(Discretization*, double*);
extern void printGrid(Discretization*, int*);
#endif
