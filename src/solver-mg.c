/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>

#include "allocate.h"
#include "solver.h"
#include "util.h"

#define FINEST_LEVEL   0
#define COARSEST_LEVEL (s->levels - 1)
#define S(i, j)        s[(j) * (imax + 2) + (i)]
#define E(i, j)        e[(j) * (imax + 2) + (i)]
#define R(i, j)        r[(j) * (imax + 2) + (i)]
#define OLD(i, j)      old[(j) * (imax + 2) + (i)]

static void restrictMG(Solver* s, int level, int imax, int jmax)
{
    double* r   = s->r[level + 1];
    double* old = s->r[level];

    for (int j = 1; j < jmax + 1; j++) {
        for (int i = 1; i < imax + 1; i++) {
            R(i, j) = (OLD(2 * i - 1, 2 * j - 1) + OLD(2 * i, 2 * j - 1) * 2 +
                          OLD(2 * i + 1, 2 * j - 1) + OLD(2 * i - 1, 2 * j) * 2 +
                          OLD(2 * i, 2 * j) * 4 + OLD(2 * i + 1, 2 * j) * 2 +
                          OLD(2 * i - 1, 2 * j + 1) + OLD(2 * i, 2 * j + 1) * 2 +
                          OLD(2 * i + 1, 2 * j + 1)) /
                      16.0;
        }
    }
}

static void prolongate(Solver* s, int level, int imax, int jmax)
{
    double* old = s->r[level + 1];
    double* e   = s->r[level];

    for (int j = 2; j < jmax + 1; j += 2) {
        for (int i = 2; i < imax + 1; i += 2) {
            E(i, j) = OLD(i / 2, j / 2);
        }
    }
}

static void correct(Solver* s, double* p, int level, int imax, int jmax)
{
    double* e = s->e[level];

    for (int j = 1; j < jmax + 1; ++j) {
        for (int i = 1; i < imax + 1; ++i) {
            P(i, j) += E(i, j);
        }
    }
}

static void setBoundaryCondition(double* p, int imax, int jmax)
{
    for (int i = 1; i < imax + 1; i++) {
        P(i, 0)        = P(i, 1);
        P(i, jmax + 1) = P(i, jmax);
    }

    for (int j = 1; j < jmax + 1; j++) {
        P(0, j)        = P(1, j);
        P(imax + 1, j) = P(imax, j);
    }
}

static void smooth(Solver* s, double* p, double* rhs, int level, int imax, int jmax)
{
    double dx2    = s->grid->dx * s->grid->dx;
    double dy2    = s->grid->dy * s->grid->dy;
    double idx2   = 1.0 / dx2;
    double idy2   = 1.0 / dy2;
    double factor = s->omega * 0.5 * (dx2 * dy2) / (dx2 + dy2);
    double* r     = s->r[level];
    double res    = 1.0;
    int pass, jsw, isw;

    jsw = 1;

    for (pass = 0; pass < 2; pass++) {
        isw = jsw;

        for (int j = 1; j < jmax + 1; j++) {
            for (int i = isw; i < imax + 1; i += 2) {

                P(i, j) -= factor * (RHS(i, j) -
                          ((P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) * idx2 +
                              (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) * idy2));

            }
            isw = 3 - isw;
        }
        jsw = 3 - jsw;
    }

}

static double calculateResidual(Solver* s, double* p, double* rhs, int level, int imax, int jmax)
{
    double dx2    = s->grid->dx * s->grid->dx;
    double dy2    = s->grid->dy * s->grid->dy;
    double idx2   = 1.0 / dx2;
    double idy2   = 1.0 / dy2;
    double factor = s->omega * 0.5 * (dx2 * dy2) / (dx2 + dy2);
    double* r     = s->r[level];
    double res    = 1.0;
    int pass, jsw, isw;

    jsw = 1;

    for (pass = 0; pass < 2; pass++) {
        isw = jsw;

        for (int j = 1; j < jmax + 1; j++) {
            for (int i = isw; i < imax + 1; i += 2) {

                R(i, j) = RHS(i, j) -
                          ((P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) * idx2 +
                              (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) * idy2);

                res += (R(i, j) * R(i, j));
            }
            isw = 3 - isw;
        }
        jsw = 3 - jsw;
    }

    res = res / (double)(imax * jmax);
    return res;
}

static double multiGrid(Solver* s, double* p, double* rhs, int level, int imax, int jmax)
{
    double res = 0.0;

    // coarsest level
    if (level == COARSEST_LEVEL) {
        for (int i = 0; i < 5; i++) {
            smooth(s, p, rhs, level, imax, jmax);
        }
        return res;
    }

    // pre-smoothing
    for (int i = 0; i < s->presmooth; i++) {
        smooth(s, p, rhs, level, imax, jmax);
        if (level == FINEST_LEVEL) setBoundaryCondition(p, imax, jmax);
    }

    res = calculateResidual(s, p, rhs, level, imax, jmax);

    // restrict
    restrictMG(s, level, imax, jmax);

    // MGSolver on residual and error.
    multiGrid(s, s->e[level + 1], s->r[level + 1], level + 1, imax / 2, jmax / 2);

    // prolongate
    prolongate(s, level, imax, jmax);

    // correct p on finer level using residual
    correct(s, p, level, imax, jmax);
    if (level == FINEST_LEVEL) setBoundaryCondition(p, imax, jmax);

    // post-smoothing
    for (int i = 0; i < s->postsmooth; i++) {
        smooth(s, p, rhs, level, imax, jmax);
        if (level == FINEST_LEVEL) setBoundaryCondition(p, imax, jmax);
    }

    return res;
}

void initSolver(Solver* s, Discretization* d, Parameter* p)
{
    s->eps     = p->eps;
    s->omega   = p->omg;
    s->itermax = p->itermax;
    s->levels  = p->levels;
    s->grid    = &d->grid;
    s->presmooth = p->presmooth;
    s->postsmooth = p->postsmooth;

    int imax   = s->grid->imax;
    int jmax   = s->grid->jmax;
    int levels = s->levels;
    printf("Using Multigrid solver with %d levels\n", levels);

    s->r = malloc(levels * sizeof(double*));
    s->e = malloc(levels * sizeof(double*));

    size_t size = (imax + 2) * (jmax + 2) * sizeof(double);

    for (int j = 0; j < levels; j++) {
        s->r[j] = allocate(64, size);
        s->e[j] = allocate(64, size);

        for (int i = 0; i < (imax + 2) * (jmax + 2); i++) {
            s->r[j][i] = 0.0;
            s->e[j][i] = 0.0;
        }
    }
}

double solve(Solver* s, double* p, double* rhs)
{
    double res = multiGrid(s, p, rhs, 0, s->grid->imax, s->grid->jmax);

#ifdef VERBOSE
    printf("Residuum: %.6f\n", res);
#endif

return res;
}
