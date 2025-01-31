/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "discretization.h"
#include "parameter.h"
#include "particletracing.h"
#include "progress.h"
#include "solver.h"
#include "timing.h"

static FILE* initResidualWriter()
{
    FILE* fp;
    fp = fopen("residual.dat", "w");

    if (fp == NULL) {
        printf("Error!\n");
        exit(EXIT_FAILURE);
    }

    return fp;

}

static void writeResidual(FILE* fp, double ts, double res)
{
    fprintf(fp, "%f, %f\n", ts, res);
}

int main(int argc, char** argv)
{
    double timeStart, timeStop;
    Parameter p;
    Discretization d;
    Solver s;
    ParticleTracer particletracer;

    initParameter(&p);
    FILE* fp;
    fp = initResidualWriter();

    if (argc != 2) {
        printf("Usage: %s <configFile>\n", argv[0]);
        exit(EXIT_SUCCESS);
    }

    readParameter(&p, argv[1]);
    printParameter(&p);
    initDiscretization(&d, &p);
    initSolver(&s, &d, &p);
    initParticleTracer(&particletracer, &d.grid, &p);
    printParticleTracerParameters(&particletracer);

#ifndef VERBOSE
    initProgress(d.te);
#endif

    double tau = d.tau;
    double te  = d.te;
    double t   = 0.0;
    int nt     = 0;
    double res = 0.0;

    timeStart  = getTimeStamp();

    while (t <= te) {
        if (tau > 0.0) computeTimestep(&d);
        setBoundaryConditions(&d);
        setSpecialBoundaryCondition(&d);
        setObjectBoundaryCondition(&d);

        computeFG(&d);
        computeRHS(&d);
        if (nt % 100 == 0) normalizePressure(&d);
        res = solve(&s, d.p, d.rhs);
        adaptUV(&d);
        trace(&particletracer, d.u, d.v, d.dt, t);

        writeResidual(fp, t, res);

        t += d.dt;
        nt++;

#ifdef VERBOSE
        printf("TIME %f , TIMESTEP %f\n", t, solver.dt);
#else
        printProgress(t);
#endif
    }

    timeStop = getTimeStamp();

    fclose(fp);
    stopProgress();
    freeParticles(&particletracer);
    printf("Solution took %.2fs\n", timeStop - timeStart);
    writeResult(&d);
    return EXIT_SUCCESS;
}
