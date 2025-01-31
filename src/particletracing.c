/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grid.h"
#include "vtkWriter.h"

#define U(i, j) u[(j) * (imax + 2) + (i)]
#define V(i, j) v[(j) * (imax + 2) + (i)]
#define S(i, j) s[(j) * (imax + 2) + (i)]

void printParticles(ParticleTracer* p)
{
    for (int i = 0; i < p->totalParticles; ++i) {
        printf("Particle position X : %.2f, Y : %.2f, flag : %d\n",
            p->particlePool[i].x,
            p->particlePool[i].y,
            p->particlePool[i].flag);
    }
}

static void injectParticles(ParticleTracer* p)
{
    if (p->totalParticles + p->numParticlesInLine > p->numAllocatedParticles) {
        return;
    }
    for (int i = 0; i < p->numParticlesInLine; ++i) {
        p->particlePool[p->pointer].x    = p->linSpaceLine[i].x;
        p->particlePool[p->pointer].y    = p->linSpaceLine[i].y;
        p->particlePool[p->pointer].flag = true;
        p->pointer++;
        p->totalParticles++;
    }
}

static void advanceParticles(
    ParticleTracer* p, double* restrict u, double* restrict v, double dt)
{
    int imax       = p->grid->imax;
    int jmax       = p->grid->jmax;
    double dx      = p->grid->dx;
    double dy      = p->grid->dy;
    double xlength = p->grid->xlength;
    double ylength = p->grid->ylength;

    for (int i = 0; i < p->totalParticles; ++i) {
        if (p->particlePool[i].flag == true) {
            double x = p->particlePool[i].x;
            double y = p->particlePool[i].y;

            int iCoord = (int)(x / dx) + 1;
            int jCoord = (int)((y + 0.5 * dy) / dy) + 1;

            double x1 = (double)(iCoord - 1) * dx;
            double y1 = ((double)(jCoord - 1) - 0.5) * dy;
            double x2 = (double)iCoord * dx;
            double y2 = ((double)jCoord - 0.5) * dy;

            double intU = (1.0 / (dx * dy)) *
                          ((x2 - x) * (y2 - y) * U(iCoord - 1, jCoord - 1) +
                              (x - x1) * (y2 - y) * U(iCoord, jCoord - 1) +
                              (x2 - x) * (y - y1) * U(iCoord - 1, jCoord) +
                              (x - x1) * (y - y1) * U(iCoord, jCoord));

            double newX          = x + dt * intU;
            p->particlePool[i].x = newX;

            iCoord = (int)((x + 0.5 * dx) / dx) + 1;
            jCoord = (int)(y / dy) + 1;

            x1 = ((double)(iCoord - 1) - 0.5) * dx;
            y1 = (double)(jCoord - 1) * dy;
            x2 = ((double)iCoord - 0.5) * dx;
            y2 = (double)jCoord * dy;

            double intV = (1.0 / (dx * dy)) *
                          ((x2 - x) * (y2 - y) * V(iCoord - 1, jCoord - 1) +
                              (x - x1) * (y2 - y) * V(iCoord, jCoord - 1) +
                              (x2 - x) * (y - y1) * V(iCoord - 1, jCoord) +
                              (x - x1) * (y - y1) * V(iCoord, jCoord));

            double newY          = y + dt * intV;
            p->particlePool[i].y = newY;

            if (((newX < 0.0) || (newX > xlength) || (newY < 0.0) || (newY > ylength))) {
                p->particlePool[i].flag = false;
                p->removedParticles++;
            }

            int newI = newX / dx, newJ = newY / dy;

            if (!gridIsFluid(p->grid, newI, newJ)) {
                p->particlePool[i].flag = false;
                p->removedParticles++;
                // printf("Forbidden movement of particle into obstacle!\n");
            }
        }
    }
}

static void compress(ParticleTracer* p)
{
    Particle* memPool = p->particlePool;
    Particle tempPool[p->totalParticles];
    int totalParticles = 0;

    // printf("Performing compression ...");

    for (int i = 0; i < p->totalParticles; i++) {
        if (memPool[i].flag == 1) {
            tempPool[totalParticles].x    = memPool[i].x;
            tempPool[totalParticles].y    = memPool[i].y;
            tempPool[totalParticles].flag = memPool[i].flag;
            totalParticles++;
        }
    }

    // printf(" remove %d particles\n", p->totalParticles - totalParticles);
    p->totalParticles   = totalParticles;
    p->removedParticles = 0;
    p->pointer          = totalParticles + 1;
    memcpy(p->particlePool, tempPool, totalParticles * sizeof(Particle));
}

void writeParticles(ParticleTracer* p)
{
    static int ts = 0;
    compress(p);
    
    VtkOptions opts = { .particletracer = p };

    char filename[50];
    // snprintf(filename, 50, "vtk_files/particles%d.vtk", ts);
    // vtkOpen(&opts, filename, ts);
    // vtkParticle(&opts, "particle");
    // vtkClose(&opts);

    FILE* fp;
    Particle* particlePool = p->particlePool;

    snprintf(filename, 50, "vis_files/particles_%d.dat", ts);
    fp = fopen(filename, "w");

    if (fp == NULL) {
        printf("Error!\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < p->totalParticles; ++i) {
        double x = particlePool[i].x;
        double y = particlePool[i].y;
        fprintf(fp, "%f %f\n", x, y);
    }
    fclose(fp);

    ++ts;
}

void initParticleTracer(ParticleTracer* pt, Grid* g, Parameter* p)
{
    pt->numParticlesInLine = p->numberOfParticles;
    pt->startTime          = p->startTime;
    pt->injectTimePeriod   = p->injectTimePeriod;
    pt->writeTimePeriod    = p->writeTimePeriod;
    pt->grid               = g;

    pt->x1 = p->x1;
    pt->y1 = p->y1;
    pt->x2 = p->x2;
    pt->y2 = p->y2;

    pt->lastInjectTime = p->startTime;
    pt->lastWriteTime  = p->startTime;

    pt->pointer          = 0;
    pt->removedParticles = 0;
    pt->totalParticles   = 0;

    if (p->te > p->startTime) {
        pt->numAllocatedParticles = ((p->te - p->startTime) / p->injectTimePeriod) *
                                    p->numberOfParticles;
        pt->numAllocatedParticles += (2 * p->numberOfParticles);

        pt->particlePool = malloc(sizeof(Particle) * pt->numAllocatedParticles);
        pt->linSpaceLine = malloc(sizeof(Particle) * pt->numParticlesInLine);

        for (int i = 0; i < pt->numParticlesInLine; ++i) {
            double spacing           = (double)i / (double)(pt->numParticlesInLine - 1);
            pt->linSpaceLine[i].x    = spacing * pt->x1 + (1.0 - spacing) * pt->x2;
            pt->linSpaceLine[i].y    = spacing * pt->y1 + (1.0 - spacing) * pt->y2;
            pt->linSpaceLine[i].flag = true;
        }
    } else {
        pt->particlePool = NULL;
        pt->linSpaceLine = NULL;
    }
}

void printParticleTracerParameters(ParticleTracer* p)
{
    printf("Particle Tracing data:\n");
    printf("\tNumber of particles : %d being injected for every period of %.2f\n",
        p->numParticlesInLine,
        p->injectTimePeriod);
    printf("\tstartTime : %.2f\n", p->startTime);
    printf("\t(Line along which the particles are to be injected) \n\tx1 : %.2f, y1 : "
           "%.2f, x2 : %.2f, y2 : %.2f\n",
        p->x1,
        p->y1,
        p->x2,
        p->y2);
    printf("\tPointer : %d, TotalParticles : %d\n", p->pointer, p->totalParticles);
}

void trace(ParticleTracer* p, double* u, double* v, double dt, double time)
{
    if (time >= p->startTime) {
        if ((time - p->lastInjectTime) >= p->injectTimePeriod) {
            injectParticles(p);
            p->lastInjectTime = time;
        }

        if ((time - p->lastWriteTime) >= p->writeTimePeriod) {
            writeParticles(p);
            p->lastWriteTime = time;
        }

        advanceParticles(p, u, v, dt);

        if (p->removedParticles > (p->totalParticles * 0.2)) {
            compress(p);
        }
    }
}

void freeParticles(ParticleTracer* p)
{
    if (p->particlePool != NULL) {
        free(p->particlePool);
        free(p->linSpaceLine);
    }
}
