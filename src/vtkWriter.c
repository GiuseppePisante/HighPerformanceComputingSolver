/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vtkWriter.h"

static float floatSwap(float f)
{
    union {
        float f;
        char b[4];
    } dat1, dat2;

    dat1.f    = f;
    dat2.b[0] = dat1.b[3];
    dat2.b[1] = dat1.b[2];
    dat2.b[2] = dat1.b[1];
    dat2.b[3] = dat1.b[0];
    return dat2.f;
}

static void writeHeader(VtkOptions* o, int ts)
{
    fprintf(o->fh, "# vtk DataFile Version 3.0\n");
    fprintf(o->fh, "PAMPI cfd solver particle tracing file\n");
    if (o->fmt == ASCII) {
        fprintf(o->fh, "ASCII\n");
    } else if (o->fmt == BINARY) {
        fprintf(o->fh, "BINARY\n");
    }

    fprintf(o->fh, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(o->fh, "FIELD FieldData 2\n");
    fprintf(o->fh, "TIME 1 1 double\n");
    fprintf(o->fh, "%d\n", ts);
    fprintf(o->fh, "CYCLE 1 1 int\n");
    fprintf(o->fh, "1\n");
}

void vtkOpen(VtkOptions* o, char* problem, int ts)
{
    o->fh = fopen(problem, "w");

    if (o->fh == NULL) {
        printf("vtkWriter not initialize! Call vtkOpen first!\n");
        exit(EXIT_FAILURE);
    }

    writeHeader(o, ts);

    printf("Writing VTK output for %s\n", problem);
}

void vtkParticle(VtkOptions* o, char* name)
{
    Particle* particlePool = o->particletracer->particlePool;

    int imax = o->particletracer->grid->imax;
    int jmax = o->particletracer->grid->jmax;

    if (o->fh == NULL) {
        printf("vtkWriter not initialize! Call vtkOpen first!\n");
        exit(EXIT_FAILURE);
    }

    fprintf(o->fh, "POINTS %d float\n", o->particletracer->totalParticles);

    for (int i = 0; i < o->particletracer->totalParticles; ++i) {
        double x = particlePool[i].x;
        double y = particlePool[i].y;
        fprintf(o->fh, "%.2f %.2f 0.0\n", x, y);
    }

    fprintf(o->fh,
        "CELLS %d %d\n",
        o->particletracer->totalParticles,
        2 * o->particletracer->totalParticles);

    for (int i = 0; i < o->particletracer->totalParticles; ++i) {
        fprintf(o->fh, "1 %d\n", i);
    }

    fprintf(o->fh, "CELL_TYPES %d\n", o->particletracer->totalParticles);

    for (int i = 0; i < o->particletracer->totalParticles; ++i) {
        fprintf(o->fh, "1\n");
    }

    /*
    for (int k = 0; k < kmax; k++) {
        for (int j = 0; j < jmax; j++) {
            for (int i = 0; i < imax; i++) {
                if (o->fmt == ASCII) {
                    fprintf(o->fh,
                        "%f %f %f\n",
                        G(vec.u, i, j, k),
                        G(vec.v, i, j, k),
                        G(vec.w, i, j, k));
                } else if (o->fmt == BINARY) {
                    fwrite((float[3]) { floatSwap(G(vec.u, i, j, k)),
                               floatSwap(G(vec.v, i, j, k)),
                               floatSwap(G(vec.w, i, j, k)) },
                        sizeof(float),
                        3,
                        o->fh);
                }
            }
        }
    }
    if (o->fmt == BINARY) fprintf(o->fh, "\n");

    */
}

void vtkClose(VtkOptions* o)
{
    fclose(o->fh);
    o->fh = NULL;
}
