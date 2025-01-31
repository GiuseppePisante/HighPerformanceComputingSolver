/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include "grid.h"
#include "solver.h"
#include "util.h"

void initSolver(Solver* s, Discretization* d, Parameter* p)
{
    s->grid    = &d->grid;
    s->itermax = p->itermax;
    s->eps     = p->eps;
    s->omega   = p->omg;

    Grid* g  = s->grid;
    int imax = s->grid->imax;
    int jmax = s->grid->jmax;
}

double solve(Solver* s, double* p, double* rhs)
{
    int imax      = s->grid->imax;
    int jmax      = s->grid->jmax;
    double eps    = s->eps;
    int itermax   = s->itermax;
    double dx2    = s->grid->dx * s->grid->dx;
    double dy2    = s->grid->dy * s->grid->dy;
    double idx2   = 1.0 / dx2;
    double idy2   = 1.0 / dy2;
    double factor = s->omega * 0.5 * (dx2 * dy2) / (dx2 + dy2);
    double epssq  = eps * eps;
    int it        = 0;
    Grid* g       = s->grid;
    double res    = 1.0;

    while ((res >= epssq) && (it < itermax)) {
        res = 0.0;

        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                double r = RHS(i, j) -
                           ((P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) * idx2 +
                               (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) * idy2);

                P(i, j) -= (factor * r);
                res += (r * r);
            }
        }

        for (int i = 1; i < imax + 1; i++) {
            P(i, 0)        = P(i, 1);
            P(i, jmax + 1) = P(i, jmax);
        }

        for (int j = 1; j < jmax + 1; j++) {
            P(0, j)        = P(1, j);
            P(imax + 1, j) = P(imax, j);
        }

        res = res / (double)(imax * jmax);
#ifdef DEBUG
        printf("%d Residuum: %e\n", it, res);
#endif
        it++;
    }

#ifdef VERBOSE
    printf("Solver took %d iterations to reach %f\n", it, sqrt(res));
#endif

    return res;
}
