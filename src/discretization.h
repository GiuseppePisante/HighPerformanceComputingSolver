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

enum BC { NOSLIP = 1, SLIP, OUTFLOW, PERIODIC };

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
extern void computeRHS(Discretization*);
extern void normalizePressure(Discretization*);
extern void computeTimestep(Discretization*);
extern void setBoundaryConditions(Discretization*);
extern void setSpecialBoundaryCondition(Discretization*);
extern void setObjectBoundaryCondition(Discretization*);
extern void computeFG(Discretization*);
extern void adaptUV(Discretization*);
extern void writeResult(Discretization*);
extern void print(Discretization*, double*);
extern void printGrid(Discretization*, int*);
#endif
