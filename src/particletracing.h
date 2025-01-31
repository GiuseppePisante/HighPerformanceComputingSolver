/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#ifndef __PARTICLETRACING_H_
#define __PARTICLETRACING_H_
#include <stdbool.h>

#include "grid.h"
#include "parameter.h"

typedef enum COORD { X = 0, Y, NCOORD } COORD;

typedef struct {
    double x, y;
    bool flag;
} Particle;

typedef struct {
    int numParticlesInLine, removedParticles, totalParticles;
    double startTime, injectTimePeriod, writeTimePeriod;
    double lastInjectTime, lastUpdateTime, lastWriteTime;

    int numAllocatedParticles;

    Particle* linSpaceLine;
    Particle* particlePool;

    int pointer;
    double x1, y1, x2, y2;
    Grid* grid;
} ParticleTracer;

extern void initParticleTracer(ParticleTracer*, Grid*, Parameter*);
extern void freeParticles(ParticleTracer*);
extern void writeParticles(ParticleTracer*);
extern void printParticleTracerParameters(ParticleTracer*);
extern void trace(ParticleTracer*, double*, double*, double, double);
#endif
