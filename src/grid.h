/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#ifndef __GRID_H_
#define __GRID_H_

#define S(i, j) s[(j) * (imax + 2) + (i)]

enum OBJECTBOUNDARY {
    FLUID = 0,
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
    TOPLEFT,
    BOTTOMLEFT,
    TOPRIGHT,
    BOTTOMRIGHT,
    OBSTACLE
};

enum SHAPE { NOSHAPE = 0, RECT, CIRCLE };

typedef struct {
    double dx, dy;
    int imax, jmax;
    double xlength, ylength;
    int* s;
} Grid;

static inline int gridIsFluid(Grid* g, int i, int j)
{
    return g->s[j * (g->imax + 2) + i] == FLUID;
}
#endif // __GRID_H_
