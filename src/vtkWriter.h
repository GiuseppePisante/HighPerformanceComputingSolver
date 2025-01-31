/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#ifndef __VTKWRITER_H_
#define __VTKWRITER_H_
#include <stdio.h>

#include "particletracing.h"
#include "solver.h"

typedef enum VtkFormat { ASCII = 0, BINARY } VtkFormat;

typedef struct VtkOptions {
    VtkFormat fmt;
    FILE* fh;
    ParticleTracer* particletracer;
} VtkOptions;

extern void vtkOpen(VtkOptions* opts, char* filename, int ts);
extern void vtkParticle(VtkOptions* opts, char* name);
extern void vtkClose(VtkOptions* opts);
#endif // __VTKWRITER_H_
