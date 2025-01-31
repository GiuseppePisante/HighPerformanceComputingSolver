/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocate.h"
#include "discretization.h"
#include "grid.h"
#include "parameter.h"
#include "util.h"

#define S(i, j) s[(j) * (imax + 2) + (i)]

static double distance(double i, double j, double iCenter, double jCenter)
{
    return sqrt(pow(iCenter - i, 2) + pow(jCenter - j, 2) * 1.0);
}

void print(Discretization* d, double* grid)
{
    int imax = d->grid.imax;
    int jmax = d->grid.jmax;

    for (int j = 0; j < jmax + 2; j++) {
        printf("%02d: ", j);
        for (int i = 0; i < imax + 2; i++) {
            printf("%3.2f  ", grid[j * (imax + 2) + i]);
        }
        printf("\n");
    }
    fflush(stdout);
}

void printGrid(Discretization* d, int* grid)
{
    int imax = d->grid.imax;
    int jmax = d->grid.jmax;

    for (int j = 0; j < jmax + 2; j++) {
        printf("%02d: ", j);
        for (int i = 0; i < imax + 2; i++) {
            printf("%2d  ", grid[j * (imax + 2) + i]);
        }
        printf("\n");
    }
    fflush(stdout);
}

static void printConfig(Discretization* d)
{
    printf("Parameters for #%s#\n", d->problem);
    printf("Boundary conditions Left:%d Right:%d Bottom:%d Top:%d\n",
        d->bcLeft,
        d->bcRight,
        d->bcBottom,
        d->bcTop);
    printf("\tReynolds number: %.2f\n", d->re);
    printf("\tGx Gy: %.2f %.2f\n", d->gx, d->gy);
    printf("Geometry data:\n");
    printf("\tDomain box size (x, y): %.2f, %.2f\n", d->grid.xlength, d->grid.ylength);
    printf("\tCells (x, y): %d, %d\n", d->grid.imax, d->grid.jmax);
    printf("Timestep parameters:\n");
    printf("\tDefault stepsize: %.2f, Final time %.2f\n", d->dt, d->te);
    printf("\tdt bound: %.6f\n", d->dtBound);
    printf("\tTau factor: %.2f\n", d->tau);
    printf("Iterative d parameters:\n");
    printf("\tgamma factor: %f\n", d->gamma);
}

void initDiscretization(Discretization* d, Parameter* p)
{
    d->problem      = p->name;
    d->bcLeft       = p->bcLeft;
    d->bcRight      = p->bcRight;
    d->bcBottom     = p->bcBottom;
    d->bcTop        = p->bcTop;
    d->grid.imax    = p->imax;
    d->grid.jmax    = p->jmax;
    d->grid.xlength = p->xlength;
    d->grid.ylength = p->ylength;
    d->grid.dx      = p->xlength / p->imax;
    d->grid.dy      = p->ylength / p->jmax;
    d->re           = p->re;
    d->gx           = p->gx;
    d->gy           = p->gy;
    d->dt           = p->dt;
    d->te           = p->te;
    d->tau          = p->tau;
    d->gamma        = p->gamma;

    int imax    = d->grid.imax;
    int jmax    = d->grid.jmax;
    size_t size = (imax + 2) * (jmax + 2) * sizeof(double);
    d->u        = allocate(64, size);
    d->v        = allocate(64, size);
    d->grid.s   = allocate(64, size);
    d->p        = allocate(64, size);
    d->rhs      = allocate(64, size);
    d->f        = allocate(64, size);
    d->g        = allocate(64, size);

    for (int i = 0; i < (imax + 2) * (jmax + 2); i++) {

        d->u[i]      = p->u_init;
        d->v[i]      = p->v_init;
        d->p[i]      = p->p_init;
        d->rhs[i]    = 0.0;
        d->f[i]      = 0.0;
        d->g[i]      = 0.0;
        d->grid.s[i] = FLUID;
    }

    double dx        = d->grid.dx;
    double dy        = d->grid.dy;
    double invSqrSum = 1.0 / (dx * dx) + 1.0 / (dy * dy);
    d->dtBound       = 0.5 * d->re * 1.0 / invSqrSum;

    double xCenter = 0, yCenter = 0, radius = 0;
    double x1 = 0, x2 = 0, y1 = 0, y2 = 0;

    int* s = d->grid.s;

    switch (p->shape) {
    case NOSHAPE:
        break;
    case RECT:
        x1 = p->xCenter - p->xRectLength / 2;
        x2 = p->xCenter + p->xRectLength / 2;
        y1 = p->yCenter - p->yRectLength / 2;
        y2 = p->yCenter + p->yRectLength / 2;

        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                if ((x1 <= (i * dx)) && ((i * dx) <= x2) && (y1 <= (j * dy)) &&
                    ((j * dy) <= y2)) {
                    S(i, j) = OBSTACLE;
                }
            }
        }
        break;
    case CIRCLE:
        xCenter = p->xCenter;
        yCenter = p->yCenter;
        radius  = p->circleRadius;

        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                if (distance((i * dx), (j * dy), xCenter, yCenter) <= radius) {
                    S(i, j) = OBSTACLE;
                }
            }
        }
        break;
    }

    if (p->shape != NOSHAPE) {
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {

                if (S(i, j - 1) == FLUID && S(i, j + 1) == OBSTACLE &&
                    S(i, j) == OBSTACLE)
                    S(i, j) = BOTTOM; // TOP
                if (S(i - 1, j) == FLUID && S(i + 1, j) == OBSTACLE &&
                    S(i, j) == OBSTACLE)
                    S(i, j) = LEFT; // LEFT
                if (S(i + 1, j) == FLUID && S(i - 1, j) == OBSTACLE &&
                    S(i, j) == OBSTACLE)
                    S(i, j) = RIGHT; // RIGHT
                if (S(i, j + 1) == FLUID && S(i, j - 1) == OBSTACLE &&
                    S(i, j) == OBSTACLE)
                    S(i, j) = TOP; // BOTTOM
                if (S(i - 1, j - 1) == FLUID && S(i, j - 1) == FLUID &&
                    S(i - 1, j) == FLUID && S(i + 1, j + 1) == OBSTACLE &&
                    (S(i, j) == OBSTACLE || S(i, j) == LEFT || S(i, j) == BOTTOM))
                    S(i, j) = BOTTOMLEFT; // TOPLEFT
                if (S(i + 1, j - 1) == FLUID && S(i, j - 1) == FLUID &&
                    S(i + 1, j) == FLUID && S(i - 1, j + 1) == OBSTACLE &&
                    (S(i, j) == OBSTACLE || S(i, j) == RIGHT || S(i, j) == BOTTOM))
                    S(i, j) = BOTTOMRIGHT; // TOPRIGHT
                if (S(i - 1, j + 1) == FLUID && S(i - 1, j) == FLUID &&
                    S(i, j + 1) == FLUID && S(i + 1, j - 1) == OBSTACLE &&
                    (S(i, j) == OBSTACLE || S(i, j) == LEFT || S(i, j) == TOP))
                    S(i, j) = TOPLEFT; // BOTTOMLEFT
                if (S(i + 1, j + 1) == FLUID && S(i + 1, j) == FLUID &&
                    S(i, j + 1) == FLUID && S(i - 1, j - 1) == OBSTACLE &&
                    (S(i, j) == OBSTACLE || S(i, j) == RIGHT || S(i, j) == TOP))
                    S(i, j) = TOPRIGHT; // BOTTOMRIGHT
            }
        }
    }

#ifdef VERBOSE
    printConfig(solver);
#endif
}

static double maxElement(Discretization* d, double* m)
{
    int size      = (d->grid.imax + 2) * (d->grid.jmax + 2);
    double maxval = DBL_MIN;

    for (int i = 0; i < size; i++) {
        maxval = MAX(maxval, fabs(m[i]));
    }

    return maxval;
}

void computeRHS(Discretization* d)
{
    int imax    = d->grid.imax;
    int jmax    = d->grid.jmax;
    double idx  = 1.0 / d->grid.dx;
    double idy  = 1.0 / d->grid.dy;
    double idt  = 1.0 / d->dt;
    double* rhs = d->rhs;
    double* f   = d->f;
    double* g   = d->g;
    int* s      = d->grid.s;

    for (int j = 1; j < jmax + 1; j++) {
        for (int i = 1; i < imax + 1; i++) {
            RHS(i, j) = idt *
                        ((F(i, j) - F(i - 1, j)) * idx + (G(i, j) - G(i, j - 1)) * idy);
        }
    }
}

void normalizePressure(Discretization* d)
{
    int size    = (d->grid.imax + 2) * (d->grid.jmax + 2);
    double* p   = d->p;
    double avgP = 0.0;

    for (int i = 0; i < size; i++) {
        avgP += p[i];
    }
    avgP /= size;

    for (int i = 0; i < size; i++) {
        p[i] = p[i] - avgP;
    }
}

void computeTimestep(Discretization* d)
{
    double dt   = d->dtBound;
    double dx   = d->grid.dx;
    double dy   = d->grid.dy;
    double umax = maxElement(d, d->u);
    double vmax = maxElement(d, d->v);

    if (umax > 0) {
        dt = (dt > dx / umax) ? dx / umax : dt;
    }
    if (vmax > 0) {
        dt = (dt > dy / vmax) ? dy / vmax : dt;
    }

    d->dt = dt * d->tau;
}

void setBoundaryConditions(Discretization* d)
{
    int imax  = d->grid.imax;
    int jmax  = d->grid.jmax;
    double* u = d->u;
    double* v = d->v;

    // Left boundary
    switch (d->bcLeft) {
    case NOSLIP:
        for (int j = 1; j < jmax + 1; j++) {
            U(0, j) = 0.0;
            V(0, j) = -V(1, j);
        }
        break;
    case SLIP:
        for (int j = 1; j < jmax + 1; j++) {
            U(0, j) = 0.0;
            V(0, j) = V(1, j);
        }
        break;
    case OUTFLOW:
        for (int j = 1; j < jmax + 1; j++) {
            U(0, j) = U(1, j);
            V(0, j) = V(1, j);
        }
        break;
    case PERIODIC:
        break;
    }

    // Right boundary
    switch (d->bcRight) {
    case NOSLIP:
        for (int j = 1; j < jmax + 1; j++) {
            U(imax, j)     = 0.0;
            V(imax + 1, j) = -V(imax, j);
        }
        break;
    case SLIP:
        for (int j = 1; j < jmax + 1; j++) {
            U(imax, j)     = 0.0;
            V(imax + 1, j) = V(imax, j);
        }
        break;
    case OUTFLOW:
        for (int j = 1; j < jmax + 1; j++) {
            U(imax, j)     = U(imax - 1, j);
            V(imax + 1, j) = V(imax, j);
        }
        break;
    case PERIODIC:
        break;
    }

    // Bottom boundary
    switch (d->bcBottom) {
    case NOSLIP:
        for (int i = 1; i < imax + 1; i++) {
            V(i, 0) = 0.0;
            U(i, 0) = -U(i, 1);
        }
        break;
    case SLIP:
        for (int i = 1; i < imax + 1; i++) {
            V(i, 0) = 0.0;
            U(i, 0) = U(i, 1);
        }
        break;
    case OUTFLOW:
        for (int i = 1; i < imax + 1; i++) {
            U(i, 0) = U(i, 1);
            V(i, 0) = V(i, 1);
        }
        break;
    case PERIODIC:
        break;
    }

    // Top boundary
    switch (d->bcTop) {
    case NOSLIP:
        for (int i = 1; i < imax + 1; i++) {
            V(i, jmax)     = 0.0;
            U(i, jmax + 1) = -U(i, jmax);
        }
        break;
    case SLIP:
        for (int i = 1; i < imax + 1; i++) {
            V(i, jmax)     = 0.0;
            U(i, jmax + 1) = U(i, jmax);
        }
        break;
    case OUTFLOW:
        for (int i = 1; i < imax + 1; i++) {
            U(i, jmax + 1) = U(i, jmax);
            V(i, jmax)     = V(i, jmax - 1);
        }
        break;
    case PERIODIC:
        break;
    }
}

void setSpecialBoundaryCondition(Discretization* d)
{
    int imax   = d->grid.imax;
    int jmax   = d->grid.jmax;
    double mDy = d->grid.dy;
    double* u  = d->u;
    int* s     = d->grid.s;

    if (strcmp(d->problem, "dcavity") == 0) {
        for (int i = 1; i < imax; i++) {
            U(i, jmax + 1) = 2.0 - U(i, jmax);
        }
    } else if (strcmp(d->problem, "canal") == 0) {
        double ylength = d->grid.ylength;
        double y;

        for (int j = 1; j < jmax + 1; j++) {
            y       = mDy * (j - 0.5);
            U(0, j) = y * (ylength - y) * 4.0 / (ylength * ylength);
        }
    } else if (strcmp(d->problem, "backstep") == 0) {
        for (int j = 1; j < jmax + 1; j++) {
            if (S(0, j) == FLUID) U(0, j) = 1.0;
        }
    } else if (strcmp(d->problem, "karman") == 0) {
        for (int j = 1; j < jmax + 1; j++) {
            U(0, j) = 1.0;
        }
    }
}

void setObjectBoundaryCondition(Discretization* d)
{
    int imax  = d->grid.imax;
    int jmax  = d->grid.jmax;
    double* u = d->u;
    double* v = d->v;
    int* s    = d->grid.s;

    for (int j = 1; j < jmax + 1; j++) {
        for (int i = 1; i < imax + 1; i++) {
            switch (S(i, j)) {
            case TOP:
                U(i, j)     = -U(i, j + 1);
                U(i - 1, j) = -U(i - 1, j + 1);
                V(i, j)     = 0.0;
                break;
            case BOTTOM:
                U(i, j)     = -U(i, j - 1);
                U(i - 1, j) = -U(i - 1, j - 1);
                V(i, j)     = 0.0;
                break;
            case LEFT:
                U(i - 1, j) = 0.0;
                V(i, j)     = -V(i - 1, j);
                V(i, j - 1) = -V(i - 1, j - 1);
                break;
            case RIGHT:
                U(i, j)     = 0.0;
                V(i, j)     = -V(i + 1, j);
                V(i, j - 1) = -V(i + 1, j - 1);
                break;
            case TOPLEFT:
                U(i, j)     = -U(i, j + 1);
                U(i - 1, j) = 0.0;
                V(i, j)     = 0.0;
                V(i, j - 1) = -V(i - 1, j - 1);
                break;
            case TOPRIGHT:
                U(i, j)     = 0.0;
                U(i - 1, j) = -U(i - 1, j + 1);
                V(i, j)     = 0.0;
                V(i, j - 1) = -V(i + 1, j - 1);
                break;
            case BOTTOMLEFT:
                U(i, j)     = -U(i, j - 1);
                U(i - 1, j) = 0.0;
                V(i, j)     = -V(i - 1, j);
                V(i, j - 1) = 0.0;
                break;
            case BOTTOMRIGHT:
                U(i, j)     = 0.0;
                U(i - 1, j) = -U(i - 1, j - 1);
                V(i, j)     = -V(i, j + 1);
                V(i, j - 1) = 0.0;
                break;
            }
        }
    }
}

void computeFG(Discretization* d)
{
    double* u        = d->u;
    double* v        = d->v;
    double* f        = d->f;
    double* g        = d->g;
    int* s           = d->grid.s;
    int imax         = d->grid.imax;
    int jmax         = d->grid.jmax;
    double gx        = d->gx;
    double gy        = d->gy;
    double gamma     = d->gamma;
    double dt        = d->dt;
    double inverseRe = 1.0 / d->re;
    double inverseDx = 1.0 / d->grid.dx;
    double inverseDy = 1.0 / d->grid.dy;
    double du2dx, dv2dy, duvdx, duvdy;
    double du2dx2, du2dy2, dv2dx2, dv2dy2;

    for (int j = 1; j < jmax + 1; j++) {
        for (int i = 1; i < imax + 1; i++) {
            if (S(i, j) == FLUID) {
                du2dx = inverseDx * 0.25 *
                            ((U(i, j) + U(i + 1, j)) * (U(i, j) + U(i + 1, j)) -
                                (U(i, j) + U(i - 1, j)) * (U(i, j) + U(i - 1, j))) +
                        gamma * inverseDx * 0.25 *
                            (fabs(U(i, j) + U(i + 1, j)) * (U(i, j) - U(i + 1, j)) +
                                fabs(U(i, j) + U(i - 1, j)) * (U(i, j) - U(i - 1, j)));

                duvdy = inverseDy * 0.25 *
                            ((V(i, j) + V(i + 1, j)) * (U(i, j) + U(i, j + 1)) -
                                (V(i, j - 1) + V(i + 1, j - 1)) *
                                    (U(i, j) + U(i, j - 1))) +
                        gamma * inverseDy * 0.25 *
                            (fabs(V(i, j) + V(i + 1, j)) * (U(i, j) - U(i, j + 1)) +
                                fabs(V(i, j - 1) + V(i + 1, j - 1)) *
                                    (U(i, j) - U(i, j - 1)));

                du2dx2 = inverseDx * inverseDx *
                         (U(i + 1, j) - 2.0 * U(i, j) + U(i - 1, j));
                du2dy2 = inverseDy * inverseDy *
                         (U(i, j + 1) - 2.0 * U(i, j) + U(i, j - 1));
                F(i, j) = U(i, j) +
                          dt * (inverseRe * (du2dx2 + du2dy2) - du2dx - duvdy + gx);

                duvdx = inverseDx * 0.25 *
                            ((U(i, j) + U(i, j + 1)) * (V(i, j) + V(i + 1, j)) -
                                (U(i - 1, j) + U(i - 1, j + 1)) *
                                    (V(i, j) + V(i - 1, j))) +
                        gamma * inverseDx * 0.25 *
                            (fabs(U(i, j) + U(i, j + 1)) * (V(i, j) - V(i + 1, j)) +
                                fabs(U(i - 1, j) + U(i - 1, j + 1)) *
                                    (V(i, j) - V(i - 1, j)));

                dv2dy = inverseDy * 0.25 *
                            ((V(i, j) + V(i, j + 1)) * (V(i, j) + V(i, j + 1)) -
                                (V(i, j) + V(i, j - 1)) * (V(i, j) + V(i, j - 1))) +
                        gamma * inverseDy * 0.25 *
                            (fabs(V(i, j) + V(i, j + 1)) * (V(i, j) - V(i, j + 1)) +
                                fabs(V(i, j) + V(i, j - 1)) * (V(i, j) - V(i, j - 1)));

                dv2dx2 = inverseDx * inverseDx *
                         (V(i + 1, j) - 2.0 * V(i, j) + V(i - 1, j));
                dv2dy2 = inverseDy * inverseDy *
                         (V(i, j + 1) - 2.0 * V(i, j) + V(i, j - 1));
                G(i, j) = V(i, j) +
                          dt * (inverseRe * (dv2dx2 + dv2dy2) - duvdx - dv2dy + gy);
            } else {
                switch (S(i, j)) {
                case TOP:
                    G(i, j) = V(i, j);
                    break;
                case BOTTOM:
                    G(i, j - 1) = V(i, j - 1);
                    break;
                case LEFT:
                    F(i - 1, j) = U(i - 1, j);
                    break;
                case RIGHT:
                    F(i, j) = U(i, j);
                    break;
                case TOPLEFT:
                    F(i - 1, j) = U(i - 1, j);
                    G(i, j)     = V(i, j);
                    break;
                case TOPRIGHT:
                    F(i, j) = U(i, j);
                    G(i, j) = V(i, j);
                    break;
                case BOTTOMLEFT:
                    F(i - 1, j) = U(i - 1, j);
                    G(i, j - 1) = V(i, j - 1);
                    break;
                case BOTTOMRIGHT:
                    F(i, j)     = U(i, j);
                    G(i, j - 1) = V(i, j - 1);
                    break;
                }
            }
        }
    }

    /* ---------------------- boundary of F --------------------------- */
    for (int j = 1; j < jmax + 1; j++) {
        F(0, j)    = U(0, j);
        F(imax, j) = U(imax, j);
    }

    /* ---------------------- boundary of G --------------------------- */
    for (int i = 1; i < imax + 1; i++) {
        G(i, 0)    = V(i, 0);
        G(i, jmax) = V(i, jmax);
    }
}

void adaptUV(Discretization* d)
{
    int imax       = d->grid.imax;
    int jmax       = d->grid.jmax;
    double* p      = d->p;
    double* u      = d->u;
    double* v      = d->v;
    double* f      = d->f;
    double* g      = d->g;
    double factorX = d->dt / d->grid.dx;
    double factorY = d->dt / d->grid.dy;

    for (int j = 1; j < jmax + 1; j++) {
        for (int i = 1; i < imax + 1; i++) {
            U(i, j) = F(i, j) - (P(i + 1, j) - P(i, j)) * factorX;
            V(i, j) = G(i, j) - (P(i, j + 1) - P(i, j)) * factorY;
        }
    }
}

void writeResult(Discretization* d)
{
    int imax  = d->grid.imax;
    int jmax  = d->grid.jmax;
    double dx = d->grid.dx;
    double dy = d->grid.dy;
    double* p = d->p;
    double* u = d->u;
    double* v = d->v;
    double x = 0.0, y = 0.0;

    FILE* fp;
    fp = fopen("pressure.dat", "w");

    if (fp == NULL) {
        printf("Error!\n");
        exit(EXIT_FAILURE);
    }

    for (int j = 1; j < jmax + 1; j++) {
        y = (double)(j - 0.5) * dy;
        for (int i = 1; i < imax + 1; i++) {
            x = (double)(i - 0.5) * dx;
            fprintf(fp, "%.2f %.2f %f\n", x, y, P(i, j));
        }
        fprintf(fp, "\n");
    }

    fclose(fp);

    fp = fopen("velocity.dat", "w");

    if (fp == NULL) {
        printf("Error!\n");
        exit(EXIT_FAILURE);
    }

    for (int j = 1; j < jmax + 1; j++) {
        y = dy * (j - 0.5);
        for (int i = 1; i < imax + 1; i++) {
            x           = dx * (i - 0.5);
            double velU = (U(i, j) + U(i - 1, j)) / 2.0;
            double velV = (V(i, j) + V(i, j - 1)) / 2.0;
            double len  = sqrt((velU * velU) + (velV * velV));
            fprintf(fp, "%.2f %.2f %f %f %f\n", x, y, velU, velV, len);
        }
    }

    fclose(fp);
}
