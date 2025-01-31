/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#ifndef __PARAMETER_H_
#define __PARAMETER_H_

typedef struct {
    double xlength, ylength;
    int imax, jmax;
    int itermax, levels, presmooth, postsmooth;
    double eps, omg, rho;
    double re, tau, gamma;
    double te, dt;
    double gx, gy;
    char* name;
    int bcLeft, bcRight, bcBottom, bcTop;
    double u_init, v_init, p_init;

    int numberOfParticles;
    double startTime, injectTimePeriod, writeTimePeriod;

    double x1, y1, x2, y2;

    int shape;
    double xCenter, yCenter, xRectLength, yRectLength, circleRadius;

} Parameter;

void initParameter(Parameter*);
void readParameter(Parameter*, const char*);
void printParameter(Parameter*);
#endif
