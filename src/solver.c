/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocate.h"
#include "parameter.h"
#include "solver.h"
#include "util.h"

#define P(i, j) p[(j) * (imaxLocal + 2) + (i)]
#define F(i, j) f[(j) * (imaxLocal + 2) + (i)]
#define G(i, j) g[(j) * (imaxLocal + 2) + (i)]
#define U(i, j) u[(j) * (imaxLocal + 2) + (i)]
#define V(i, j) v[(j) * (imaxLocal + 2) + (i)]
#define RHS(i, j) rhs[(j) * (imaxLocal + 2) + (i)]
#define Pall(i, j) Pall[(j) * (imax + 2) + (i)]
#define Uall(i, j) Uall[(j) * (imax + 2) + (i)]
#define Vall(i, j) Vall[(j) * (imax + 2) + (i)]
#define R(i, j) r[(j) * (imaxLocal + 2) + (i)]
#define OLD(i, j) old[(j) * (imaxLocal + 2) + (i)]
#define E(i, j) e[(j) * (imaxLocal + 2) + (i)]

static double multiGrid(Solver *s, double *p, double *rhs, int level, int imax, int jmax);
static double calculateResidual(Solver *s, double *p, double *rhs, int level, int imaxLocal, int jmaxLocal);
static void correct(Solver *s, double *p, int level, int imaxLocal, int jmaxLocal);
static double maxElement(Solver *solver, double *m);
static void restrictMG(Solver *s, int level, int imaxLocal, int jmaxLocal);
static void prolongate(Solver *s, int level, int imaxLocal, int jmaxLocal);
static void smooth(Solver *s, double *p, double *rhs, int level, int imaxLocal, int jmaxLocal);

// ! AGGIUNTA
#define FINEST_LEVEL 0
#define COARSEST_LEVEL (s->levels - 1)
//

static int sizeOfRank(int rank, int size, int N)
{
    return N / size + ((N % size > rank) ? 1 : 0);
}

static void printExchange(Solver *solver, double *grid)
{
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;
    double *p = grid;
    for (int j = 0; j < jmaxLocal + 2; j++)
    {
        for (int i = 0; i < imaxLocal + 2; i++)
        {
            P(i, j) = solver->rank;
        }
    }

    exchange(solver, p);

    for (int i = 0; i < solver->size; i++)
    {
        if (i == solver->rank)
        {
            printf("### RANK %d "
                   "#######################################################\n",
                   solver->rank);
            for (int j = 0; j < jmaxLocal + 2; j++)
            {
                printf("%02d: ", j);
                for (int i = 0; i < imaxLocal + 2; i++)
                {
                    printf("%12.2f  ", grid[j * (imaxLocal + 2) + i]);
                }
                printf("\n");
            }
            fflush(stdout);
        }
        MPI_Barrier(solver->comm);
    }
}

static void printShift(Solver *solver)
{
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;
    double *f = solver->f;
    double *g = solver->g;

    for (int j = 0; j < jmaxLocal + 2; j++)
    {
        for (int i = 0; i < imaxLocal + 2; i++)
        {
            F(i, j) = solver->rank;
            G(i, j) = solver->rank;
        }
    }

    shift(solver);

    printf("Vector F\n");
    for (int i = 0; i < solver->size; i++)
    {
        if (i == solver->rank)
        {
            printf("### RANK %d "
                   "#######################################################\n",
                   solver->rank);
            for (int j = 0; j < jmaxLocal + 2; j++)
            {
                printf("%02d: ", j);
                for (int i = 0; i < imaxLocal + 2; i++)
                {
                    printf("%12.2f  ", F(i, j));
                }
                printf("\n");
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    printf("Vector G\n");
    for (int i = 0; i < solver->size; i++)
    {
        if (i == solver->rank)
        {
            printf("### RANK %d "
                   "#######################################################\n",
                   solver->rank);
            for (int j = 0; j < jmaxLocal + 2; j++)
            {
                printf("%02d: ", j);
                for (int i = 0; i < imaxLocal + 2; i++)
                {
                    printf("%12.2f  ", G(i, j));
                }
                printf("\n");
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// subroutines local to this module
static int sum(int *sizes, int init, int offset, int coord)
{
    int sum = 0;

    for (int i = init - offset; coord > 0; i -= offset, --coord)
    {
        sum += sizes[i];
    }

    return sum;
}

void exchange(Solver *solver, double *grid)
{
    /* TODO ghost cell exchange with respective neighbors */
    int counts[NDIRS] = {1, 1, 1, 1};
    MPI_Neighbor_alltoallw(grid,
                           counts,
                           solver->sdispls,
                           solver->bufferTypes,
                           grid,
                           counts,
                           solver->rdispls,
                           solver->bufferTypes,
                           solver->comm);
}

void shift(Solver *solver)
{
    /* TODO shift (send) f ghost layer from left to right */
    /* TODO shift (send) g ghost layer from bottom to top */
    double *f = solver->f;
    double *g = solver->g;

    MPI_Request requests[4] = {MPI_REQUEST_NULL,
                               MPI_REQUEST_NULL,
                               MPI_REQUEST_NULL,
                               MPI_REQUEST_NULL};

    /* shift G */
    /* receive ghost cells from bottom neighbor */
    double *buf = g + 1;
    MPI_Irecv(buf,
              1,
              solver->bufferTypes[BOTTOM],
              solver->neighbours[BOTTOM],
              0,
              solver->comm,
              &requests[0]);

    /* send ghost cells to top neighbor */
    /* Use information from solver->neighbours */
    buf = g + (solver->jmaxLocal) * (solver->imaxLocal + 2) + 1;
    MPI_Isend(buf, 1, solver->bufferTypes[TOP], solver->neighbours[TOP], 0, solver->comm, &requests[1]);

    /* shift F */
    /* receive ghost cells from left neighbor */
    buf = f + (solver->imaxLocal + 2);
    MPI_Irecv(buf,
              1,
              solver->bufferTypes[LEFT],
              solver->neighbours[LEFT],
              1,
              solver->comm,
              &requests[2]);

    /* send ghost cells to right neighbor */
    buf = f + (solver->imaxLocal + 2) + (solver->imaxLocal);
    MPI_Isend(buf,
              1,
              solver->bufferTypes[RIGHT],
              solver->neighbours[RIGHT],
              1,
              solver->comm,
              &requests[3]);

    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
}

static int isBoundary(Solver *solver, int direction)
{
    switch (direction)
    {
    case LEFT:
        return solver->coords[IDIM] == 0;
        break;
    case RIGHT:
        return solver->coords[IDIM] == (solver->dims[IDIM] - 1);
        break;
    case BOTTOM:
        return solver->coords[JDIM] == 0;
        break;
    case TOP:
        return solver->coords[JDIM] == (solver->dims[JDIM] - 1);
        break;
    }
    return 1;
}

static void assembleResult(Solver *solver, double *src, double *dst)
{
    int imax = solver->imax;
    int jmax = solver->jmax;

    MPI_Request *requests;
    int numRequests = 1;

    if (solver->rank == 0)
    {
        numRequests = solver->size + 1;
    }
    else
    {
        numRequests = 1;
    }

    requests = (MPI_Request *)malloc(numRequests * sizeof(MPI_Request));

    /* all ranks send their bulk array, including the external boundary layer */
    MPI_Datatype bulkType;
    int oldSizes[NDIMS] = {solver->jmaxLocal + 2, solver->imaxLocal + 2};
    int newSizes[NDIMS] = {solver->jmaxLocal, solver->imaxLocal};
    int starts[NDIMS] = {1, 1};

    if (isBoundary(solver, LEFT))
    {
        newSizes[CIDIM] += 1;
        starts[CIDIM] = 0;
    }
    if (isBoundary(solver, RIGHT))
    {
        newSizes[CIDIM] += 1;
    }
    if (isBoundary(solver, BOTTOM))
    {
        newSizes[CJDIM] += 1;
        starts[CJDIM] = 0;
    }
    if (isBoundary(solver, TOP))
    {
        newSizes[CJDIM] += 1;
    }

    MPI_Type_create_subarray(NDIMS,
                             oldSizes,
                             newSizes,
                             starts,
                             MPI_ORDER_C,
                             MPI_DOUBLE,
                             &bulkType);
    MPI_Type_commit(&bulkType);
    MPI_Isend(src, 1, bulkType, 0, 0, solver->comm, &requests[0]);

    int newSizesI[solver->size];
    int newSizesJ[solver->size];
    MPI_Gather(&newSizes[CIDIM], 1, MPI_INT, newSizesI, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&newSizes[CJDIM], 1, MPI_INT, newSizesJ, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* rank 0 assembles the subdomains */
    if (solver->rank == 0)
    {
        for (int i = 0; i < solver->size; i++)
        {
            MPI_Datatype domainType;
            int oldSizes[NDIMS] = {jmax + 2, imax + 2};
            int newSizes[NDIMS] = {newSizesJ[i], newSizesI[i]};
            int coords[NDIMS];
            MPI_Cart_coords(solver->comm, i, NDIMS, coords);
            int starts[NDIMS] = {sum(newSizesJ, i, 1, coords[JDIM]),
                                 sum(newSizesI, i, solver->dims[JDIM], coords[IDIM])};
            /*printf(
                "Rank: %d, Coords(i,j): %d %d, Size(i,j): %d %d, Target Size(i,j): %d %d "
                "Starts(i,j): %d %d\n",
                i,
                coords[IDIM],
                coords[JDIM],
                oldSizes[CIDIM],
                oldSizes[CJDIM],
                newSizes[CIDIM],
                newSizes[CJDIM],
                starts[CIDIM],
                starts[CJDIM]);*/

            MPI_Type_create_subarray(NDIMS,
                                     oldSizes,
                                     newSizes,
                                     starts,
                                     MPI_ORDER_C,
                                     MPI_DOUBLE,
                                     &domainType);
            MPI_Type_commit(&domainType);

            MPI_Irecv(dst, 1, domainType, i, 0, solver->comm, &requests[i + 1]);
            MPI_Type_free(&domainType);
        }
    }

    MPI_Waitall(numRequests, requests, MPI_STATUSES_IGNORE);
}

void collectResult(Solver *solver)
{
    size_t bytesize = (solver->imax + 2) * (solver->jmax + 2) * sizeof(double);

    double *Pall = allocate(64, bytesize);
    double *Uall = allocate(64, bytesize);
    double *Vall = allocate(64, bytesize);
    /* collect P */
    assembleResult(solver, solver->p, Pall);

    /* collect U */
    assembleResult(solver, solver->u, Uall);

    /* collect V */
    assembleResult(solver, solver->v, Vall);

    if (solver->rank == 0)
    {
        writeResult(solver, Pall, Uall, Vall);
    }

    free(Pall);
    free(Uall);
    free(Vall);
}

static void printConfig(Solver *solver)
{
    if (solver->rank == 0)
    {
        printf("Parameters for #%s#\n", solver->problem);
        printf("Boundary conditions Left:%d Right:%d Bottom:%d Top:%d\n",
               solver->bcLeft,
               solver->bcRight,
               solver->bcBottom,
               solver->bcTop);
        printf("\tReynolds number: %.2f\n", solver->re);
        printf("\tGx Gy: %.2f %.2f\n", solver->gx, solver->gy);
        printf("Geometry data:\n");
        printf("\tDomain box size (x, y): %.2f, %.2f\n",
               solver->xlength,
               solver->ylength);
        printf("\tCells (x, y): %d, %d\n", solver->imax, solver->jmax);
        printf("Timestep parameters:\n");
        printf("\tDefault stepsize: %.2f, Final time %.2f\n", solver->dt, solver->te);
        printf("\tdt bound: %.6f\n", solver->dtBound);
        printf("\tTau factor: %.2f\n", solver->tau);
        printf("Iterative solver parameters:\n");
        printf("\tMax iterations: %d\n", solver->itermax);
        printf("\tepsilon (stopping tolerance) : %f\n", solver->eps);
        printf("\tgamma factor: %f\n", solver->gamma);
        printf("\tomega (SOR relaxation): %f\n", solver->omega);
        printf("Communication parameters:\n");
    }
    for (int i = 0; i < solver->size; i++)
    {
        if (i == solver->rank)
        {
            printf("\tRank %d of %d\n", solver->rank, solver->size);
            printf("\tNeighbours (bottom, top, left, right): %d %d, %d, %d\n",
                   solver->neighbours[BOTTOM],
                   solver->neighbours[TOP],
                   solver->neighbours[LEFT],
                   solver->neighbours[RIGHT]);
            printf("\tIs boundary:\n");
            printf("\t\tLEFT: %d\n", isBoundary(solver, LEFT));
            printf("\t\tRIGHT: %d\n", isBoundary(solver, RIGHT));
            printf("\t\tBOTTOM: %d\n", isBoundary(solver, BOTTOM));
            printf("\t\tTOP: %d\n", isBoundary(solver, TOP));
            printf("\tCoordinates (i,j) %d %d\n", solver->coords[IDIM], solver->coords[JDIM]);
            printf("\tDims (i,j) %d %d\n", solver->dims[IDIM], solver->dims[JDIM]);
            printf("\tLocal domain size (i,j) %dx%d\n", solver->imaxLocal, solver->jmaxLocal);
            fflush(stdout);
        }
        MPI_Barrier(solver->comm);
    }
}

void initSolver(Solver *solver, Parameter *params) // !Aggiunta
{
    // Comment it out
    MPI_Comm_rank(MPI_COMM_WORLD, &(solver->rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(solver->size));
    solver->problem = params->name;
    solver->bcLeft = params->bcLeft;
    solver->bcRight = params->bcRight;
    solver->bcBottom = params->bcBottom;
    solver->bcTop = params->bcTop;
    solver->imax = params->imax;
    solver->jmax = params->jmax;
    solver->xlength = params->xlength;
    solver->ylength = params->ylength;
    solver->dx = params->xlength / params->imax;
    solver->dy = params->ylength / params->jmax;
    solver->eps = params->eps;
    solver->omega = params->omg;
    solver->itermax = params->itermax;
    solver->re = params->re;
    solver->gx = params->gx;
    solver->gy = params->gy;
    solver->dt = params->dt;
    solver->te = params->te;
    solver->tau = params->tau;
    solver->gamma = params->gamma;

    // ! AGGIUNTA:
    int imax = solver->imax;
    int jmax = solver->jmax;

    solver->levels = 3;         // Number of levels in the multigrid
    solver->presmooth = 5;   // Number of pre-smoothing steps
    solver->postsmooth = 5; // Number of post-smoothing steps
    //

    /* Use virtual topologies to get dimension and coordinates of the rank */
    int dims[NDIMS] = {0, 0};
    int periods[NDIMS] = {0, 0};

    /* Retrieves best possible #ranks in X and Y dimensions from the given communicator size. */
    /* Say for 6 ranks, best possible #ranks in in X and Y dimension would be 3 and 2         */
    /* Store in dims array                                                                    */
    /* https://rookiehpc.org/mpi/docs/mpi_dims_create/index.html                              */
    MPI_Dims_create(solver->size, NDIMS, dims);
    memcpy(&solver->dims, &dims, NDIMS * sizeof(int));

    /* Create a new communicator based on the dimensions             */
    /* https://rookiehpc.org/mpi/docs/mpi_cart_create/index.html     */
    MPI_Cart_create(MPI_COMM_WORLD, NDIMS, dims, periods, 0, &solver->comm);

    /* Use newly created communicator solver->comm                   */
    /* For left and right neighbors                                  */
    /* Store in solver->neighbors[LEFT] and solver->neighbors[RIGHT] */
    /* https://rookiehpc.org/mpi/docs/mpi_cart_shift/index.html      */
    MPI_Cart_shift(solver->comm, IDIM, 1, &solver->neighbours[LEFT], &solver->neighbours[RIGHT]);

    /* Use newly created communicator solver->comm                   */
    /* For bottom and top neighbors                                  */
    /* Store in solver->neighbors[BOTTOM] and solver->neighbors[TOP] */
    MPI_Cart_shift(solver->comm, JDIM, 1, &solver->neighbours[BOTTOM], &solver->neighbours[TOP]);

    /* Retrieve the cartesion topology of the rank in solver->coords */
    /* https://rookiehpc.org/mpi/docs/mpi_cart_get/index.html        */
    MPI_Cart_get(solver->comm, NDIMS, solver->dims, periods, solver->coords);

    solver->imaxLocal = sizeOfRank(solver->coords[IDIM], dims[IDIM], imax);
    solver->jmaxLocal = sizeOfRank(solver->coords[JDIM], dims[JDIM], jmax);

    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;
    size_t bytesize = (imaxLocal + 2) * (jmaxLocal + 2) * sizeof(double);

    MPI_Datatype jBufferType;
    MPI_Type_contiguous(imaxLocal, MPI_DOUBLE, &jBufferType);
    MPI_Type_commit(&jBufferType);

    MPI_Datatype iBufferType;
    MPI_Type_vector(jmaxLocal, 1, imaxLocal + 2, MPI_DOUBLE, &iBufferType);
    MPI_Type_commit(&iBufferType);

    solver->bufferTypes[LEFT] = iBufferType;
    solver->bufferTypes[RIGHT] = iBufferType;
    solver->bufferTypes[BOTTOM] = jBufferType;
    solver->bufferTypes[TOP] = jBufferType;

    size_t dblsize = sizeof(double);
    solver->sdispls[LEFT] = ((imaxLocal + 2) + 1) * dblsize;
    solver->sdispls[RIGHT] = ((imaxLocal + 2) + imaxLocal) * dblsize;
    solver->sdispls[BOTTOM] = ((imaxLocal + 2) + 1) * dblsize;
    solver->sdispls[TOP] = (jmaxLocal * (imaxLocal + 2) + 1) * dblsize;

    solver->rdispls[LEFT] = (imaxLocal + 2) * dblsize;
    solver->rdispls[RIGHT] = ((imaxLocal + 2) + (imaxLocal + 1)) * dblsize;
    solver->rdispls[BOTTOM] = 1 * dblsize;
    solver->rdispls[TOP] = ((jmaxLocal + 1) * (imaxLocal + 2) + 1) * dblsize;

    solver->u = allocate(64, bytesize);
    solver->v = allocate(64, bytesize);
    solver->p = allocate(64, bytesize);
    solver->rhs = allocate(64, bytesize);
    solver->f = allocate(64, bytesize);
    solver->g = allocate(64, bytesize);

    // ! AGGIUNTA:
    size_t size = (imaxLocal + 2) * (jmaxLocal + 2) * sizeof(double);
    //

    /* change to imaxLocal & jmaxLocal */
    for (int i = 0; i < solver->size; i++)
    {
        solver->u[i] = params->u_init;
        solver->v[i] = params->v_init;
        solver->p[i] = params->p_init;
        solver->rhs[i] = 0.0;
        solver->f[i] = 0.0;
        solver->g[i] = 0.0;
    }

    // ! AGGIUNTA:
    solver->r = malloc(solver->levels * sizeof(double *)); // residual vector
    solver->e = malloc(solver->levels * sizeof(double *)); // error vector
                                                           //

    // ! AGGIUNTA:
    for (int j = 0; j < solver->levels; j++)
    {
        solver->r[j] = allocate(64, size);
        solver->e[j] = allocate(64, size);

        for (int i = 0; i < (imaxLocal + 2) * (jmaxLocal + 2); i++)
        {
            solver->r[j][i] = 0.0;
            solver->e[j][i] = 0.0;
        }
    }
    //

    double dx = solver->dx;
    double dy = solver->dy;
    double invSqrSum = 1.0 / (dx * dx) + 1.0 / (dy * dy);
    solver->dtBound = 0.5 * solver->re * 1.0 / invSqrSum;
#ifdef VERBOSE
    printConfig(solver);
#endif

}

void computeRHS(Solver *solver)
{
    int imax = solver->imax;
    int jmax = solver->jmax;
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;
    double idx = 1.0 / solver->dx;
    double idy = 1.0 / solver->dy;
    double idt = 1.0 / solver->dt;
    double *rhs = solver->rhs;
    double *f = solver->f;
    double *g = solver->g;

    /* Shift function is similar to exchange. Just that you need to exchange values for F and G array */
    /* A simple Send Recv function for F and G array can be implemented. */
    /* You can look at the formulation below to guess which values will be required. */
    /* For F, you just have to send the values from left neighbor to right neighbor's ghost cells */
    /* For G, you just have to send the values from bottom neighbor to top neighbor's ghost cells */

    shift(solver);

    /* change to imaxLocal & jmaxLocal */
    for (int j = 1; j < jmaxLocal + 1; j++)
    {
        for (int i = 1; i < imaxLocal + 1; i++)
        {
            RHS(i, j) = ((F(i, j) - F(i - 1, j)) * idx + (G(i, j) - G(i, j - 1)) * idy) *
                        idt;
        }
    }
}

void solve(Solver *solver)
{
    int imax = solver->imax;
    int jmax = solver->jmax;
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;
    double eps = solver->eps;
    int itermax = solver->itermax;
    double dx2 = solver->dx * solver->dx;
    double dy2 = solver->dy * solver->dy;
    double idx2 = 1.0 / dx2;
    double idy2 = 1.0 / dy2;
    double factor = solver->omega * 0.5 * (dx2 * dy2) / (dx2 + dy2);
    double *p = solver->p;
    double *rhs = solver->rhs;
    double epssq = eps * eps;
    int it = 0;
    double res = 1.0;

    while ((res >= epssq) && (it < itermax))
    {
        res = 0.0;

        /* Exchange the ghost cells with respective neighbors */
        /* Remember this is a 2D domain decomposition.        */
        exchange(solver, p);

        // ! Commentato perché qua lo risolve in maniera diretta, noi lo vogliamo con il multigrid
        /* change to imaxLocal & jmaxLocal                    */
        /* for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {

                double r = RHS(i, j) -
                           ((P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) * idx2 +
                               (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) * idy2);

                P(i, j) -= (factor * r);
                res += (r * r);
            }
        } */

        // ! AGGIUNTA:
        double res = multiGrid(solver, p, rhs, 0, solver->imaxLocal, solver->jmaxLocal);
        //

        /* Extend boundary condition handling from Poisson assignment */
        if (isBoundary(solver, BOTTOM))
        { // set bottom bc
            for (int i = 1; i < imaxLocal + 1; i++)
            {
                P(i, 0) = P(i, 1);
            }
        }

        if (isBoundary(solver, TOP))
        { // set top bc
            for (int i = 1; i < imaxLocal + 1; i++)
            {
                P(i, jmaxLocal + 1) = P(i, jmaxLocal);
            }
        }

        if (isBoundary(solver, LEFT))
        { // set left bc
            for (int j = 1; j < jmaxLocal + 1; j++)
            {
                P(0, j) = P(1, j);
            }
        }

        if (isBoundary(solver, RIGHT))
        { // set right bc
            for (int j = 1; j < jmaxLocal + 1; j++)
            {
                P(imaxLocal + 1, j) = P(imaxLocal, j);
            }
        }

        /* add global reduction for the residual */
        MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        res = res / (double)(imax * jmax);
#ifdef DEBUG
        printf("%d Residuum: %e\n", it, res);
#endif
        it++;
    }

#ifdef VERBOSE
    if (solver->rank == 0)
        printf("Solver took %d iterations to reach %f\n", it, sqrt(res));
#endif
}

static double multiGrid(Solver *s, double *p, double *rhs, int level, int imax, int jmax)
{
    double res = 0.0;

    // coarsest level
    if (level == COARSEST_LEVEL)
    {
        for (int i = 0; i < 5; i++)
        {
            smooth(s, p, rhs, level, imax, jmax);
        }
        return res;
    }

    // pre-smoothing
    for (int i = 0; i < s->presmooth; i++)
    {
        smooth(s, p, rhs, level, imax, jmax);
        if (level == FINEST_LEVEL)
            setBoundaryConditions(s);
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
    if (level == FINEST_LEVEL)
        setBoundaryConditions(s);
    // post-smoothing
    for (int i = 0; i < s->postsmooth; i++)
    {
        smooth(s, p, rhs, level, imax, jmax);
        if (level == FINEST_LEVEL)
            setBoundaryConditions(s);
    }

    return res;
}
//! Aggiunta
static void smooth(Solver *s, double *p, double *rhs, int level, int imaxLocal, int jmaxLocal)
{
    double dx2 = s->dx * s->dx;
    double dy2 = s->dy * s->dy;
    double idx2 = 1.0 / dx2;
    double idy2 = 1.0 / dy2;
    double factor = s->omega * 0.5 * (dx2 * dy2) / (dx2 + dy2);
    double *r = s->r[level];
    double res = 1.0;
    int pass, jsw, isw;

    jsw = 1;

    for (pass = 0; pass < 2; pass++)
    {
        isw = jsw;

        for (int j = 1; j < jmaxLocal + 1; j++)
        {
            for (int i = isw; i < imaxLocal + 1; i += 2)
            {

                P(i, j) -= factor * (RHS(i, j) -
                                     ((P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) * idx2 +
                                      (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) * idy2));
            }
            isw = 3 - isw;
        }
        jsw = 3 - jsw;
    }
}
//

// ! Aggiunta
static double calculateResidual(Solver *s, double *p, double *rhs, int level, int imaxLocal, int jmaxLocal)
{
    double dx2 = s->dx * s->dx;
    double dy2 = s->dy * s->dy;
    double idx2 = 1.0 / dx2;
    double idy2 = 1.0 / dy2;
    double factor = s->omega * 0.5 * (dx2 * dy2) / (dx2 + dy2);
    double *r = s->r[level];
    double res = 1.0;
    int pass, jsw, isw;

    jsw = 1;

    for (pass = 0; pass < 2; pass++)
    {
        isw = jsw;

        for (int j = 1; j < jmaxLocal + 1; j++)
        {
            for (int i = isw; i < imaxLocal + 1; i += 2)
            {

                R(i, j) = RHS(i, j) -
                          ((P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) * idx2 +
                           (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) * idy2);

                res += (R(i, j) * R(i, j));
            }
            isw = 3 - isw;
        }
        jsw = 3 - jsw;
    }

    res = res / (double)(imaxLocal * jmaxLocal);
    return res;
}
//

// ! Aggiunta
static void restrictMG(Solver *s, int level, int imaxLocal, int jmaxLocal)
{
    double *r = s->r[level + 1];
    double *old = s->r[level];

    for (int j = 1; j < jmaxLocal + 1; j++)
    {
        for (int i = 1; i < imaxLocal + 1; i++)
        {
            R(i, j) = (OLD(2 * i - 1, 2 * j - 1) + OLD(2 * i, 2 * j - 1) * 2 +
                       OLD(2 * i + 1, 2 * j - 1) + OLD(2 * i - 1, 2 * j) * 2 +
                       OLD(2 * i, 2 * j) * 4 + OLD(2 * i + 1, 2 * j) * 2 +
                       OLD(2 * i - 1, 2 * j + 1) + OLD(2 * i, 2 * j + 1) * 2 +
                       OLD(2 * i + 1, 2 * j + 1)) /
                      16.0;
        }
    }
}
//

// ! Aggiunta:
static void prolongate(Solver *s, int level, int imaxLocal, int jmaxLocal)
{
    double *old = s->r[level + 1];
    double *e = s->r[level];

    for (int j = 2; j < jmaxLocal + 1; j += 2)
    {
        for (int i = 2; i < imaxLocal + 1; i += 2)
        {
            E(i, j) = OLD(i / 2, j / 2);
        }
    }
}
//

// !Aggiunta:

static void correct(Solver *s, double *p, int level, int imaxLocal, int jmaxLocal)
{
    double *e = s->e[level];

    for (int j = 1; j < jmaxLocal + 1; ++j)
    {
        for (int i = 1; i < imaxLocal + 1; ++i)
        {
            P(i, j) += E(i, j);
        }
    }
}
//

static double maxElement(Solver *solver, double *m)
{
    int size = (solver->imaxLocal + 2) * (solver->jmaxLocal + 2);
    double maxval = DBL_MIN;

    for (int i = 0; i < size; i++)
    {
        maxval = MAX(maxval, fabs(m[i]));
    }

    /* TODO add global MPI reduction for maxval */
    MPI_Allreduce(MPI_IN_PLACE, &maxval, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    return maxval;
}

void normalizePressure(Solver *solver)
{
    /* adapt to local size imaxLocal & jmaxLocal */
    int size = (solver->imaxLocal + 2) * (solver->jmaxLocal + 2);
    double *p = solver->p;
    double avgP = 0.0;

    for (int i = 0; i < size; i++)
    {
        avgP += p[i];
    }
    avgP /= size;

    /* TODO add global MPI reduction for avgP */
    MPI_Allreduce(MPI_IN_PLACE, &avgP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    avgP /= (solver->imax + 2) * (solver->jmax + 2);

    for (int i = 0; i < size; i++)
    {
        p[i] = p[i] - avgP;
    }
}

void computeTimestep(Solver *solver)
{
    double dt = solver->dtBound;
    double dx = solver->dx;
    double dy = solver->dy;
    double umax = maxElement(solver, solver->u);
    double vmax = maxElement(solver, solver->v);

    if (umax > 0)
    {
        dt = (dt > dx / umax) ? dx / umax : dt;
    }
    if (vmax > 0)
    {
        dt = (dt > dy / vmax) ? dy / vmax : dt;
    }

    solver->dt = dt * solver->tau;
}

// TODO adapt to local sizes
void setBoundaryConditions(Solver *solver)
{
    int imax = solver->imax;
    int jmax = solver->jmax;
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;
    double *u = solver->u;
    double *v = solver->v;

    if (isBoundary(solver, TOP))
    {
        switch (solver->bcTop)
        {
        case NOSLIP:
            for (int i = 1; i < imaxLocal + 1; i++)
            {
                V(i, jmaxLocal) = 0.0;
                U(i, jmaxLocal + 1) = -U(i, jmaxLocal);
            }
            break;
        case SLIP:
            for (int i = 1; i < imaxLocal + 1; i++)
            {
                V(i, jmaxLocal) = 0.0;
                U(i, jmaxLocal + 1) = U(i, jmaxLocal);
            }
            break;
        case OUTFLOW:
            for (int i = 1; i < imaxLocal + 1; i++)
            {
                U(i, jmaxLocal + 1) = U(i, jmaxLocal);
                V(i, jmaxLocal) = V(i, jmaxLocal - 1);
            }
            break;
        case PERIODIC:
            break;
        }
    }

    if (isBoundary(solver, BOTTOM))
    {
        switch (solver->bcBottom)
        {
        case NOSLIP:
            for (int i = 1; i < imaxLocal + 1; i++)
            {
                V(i, 0) = 0.0;
                U(i, 0) = -U(i, 1);
            }
            break;
        case SLIP:
            for (int i = 1; i < imaxLocal + 1; i++)
            {
                V(i, 0) = 0.0;
                U(i, 0) = U(i, 1);
            }
            break;
        case OUTFLOW:
            for (int i = 1; i < imaxLocal + 1; i++)
            {
                U(i, 0) = U(i, 1);
                V(i, 0) = V(i, 1);
            }
            break;
        case PERIODIC:
            break;
        }
    }

    if (isBoundary(solver, RIGHT))
    {
        switch (solver->bcRight)
        {
        case NOSLIP:
            for (int j = 1; j < jmaxLocal + 1; j++)
            {
                U(imaxLocal, j) = 0.0;
                V(imaxLocal + 1, j) = -V(imaxLocal, j);
            }
            break;
        case SLIP:
            for (int j = 1; j < jmaxLocal + 1; j++)
            {
                U(imaxLocal, j) = 0.0;
                V(imaxLocal + 1, j) = V(imaxLocal, j);
            }
            break;
        case OUTFLOW:
            for (int j = 1; j < jmaxLocal + 1; j++)
            {
                U(imaxLocal, j) = U(imaxLocal - 1, j);
                V(imaxLocal + 1, j) = V(imaxLocal, j);
            }
            break;
        case PERIODIC:
            break;
        }
    }

    if (isBoundary(solver, LEFT))
    {
        switch (solver->bcLeft)
        {
        case NOSLIP:
            for (int j = 1; j < jmaxLocal + 1; j++)
            {
                U(0, j) = 0.0;
                V(0, j) = -V(1, j);
            }
            break;
        case SLIP:
            for (int j = 1; j < jmaxLocal + 1; j++)
            {
                U(0, j) = 0.0;
                V(0, j) = V(1, j);
            }
            break;
        case OUTFLOW:
            for (int j = 1; j < jmaxLocal + 1; j++)
            {
                U(0, j) = U(1, j);
                V(0, j) = V(1, j);
            }
            break;
        case PERIODIC:
            break;
        }
    }
}

// TODO adapt to local sizes
void setSpecialBoundaryCondition(Solver *solver)
{
    int imax = solver->imax;
    int jmax = solver->jmax;
    double mDy = solver->dy;
    double *u = solver->u;
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;

    if (strcmp(solver->problem, "dcavity") == 0)
    {
        if (isBoundary(solver, TOP))
        {
            for (int i = 1; i < imaxLocal + 1; i++)
            {
                U(i, jmaxLocal + 1) = 2.0 - U(i, jmaxLocal);
            }
        }
    }
    else if (strcmp(solver->problem, "canal") == 0)
    {
        if (isBoundary(solver, LEFT))
        {
            double ylength = solver->ylength;
            double dy = solver->dy;
            int rest = solver->jmax % solver->dims[JDIM];
            int yc = solver->rank * (solver->jmax / solver->dims[JDIM]) +
                     MIN(rest, solver->rank);
            double ys = dy * (yc + 0.5);
            double y;

            // printf("RANK %d yc: %d ys: %f\n", s->comm.rank, yc, ys);

            for (int j = 1; j < jmaxLocal + 1; j++)
            {
                y = ys + dy * (j - 0.5);
                U(0, j) = y * (ylength - y) * 4.0 / (ylength * ylength);
                U(1, j) = U(0, j);
            }
        }
    }
}

// TODO adapt to local sizes
void computeFG(Solver *solver)
{
    double *u = solver->u;
    double *v = solver->v;
    double *f = solver->f;
    double *g = solver->g;
    int imax = solver->imax;
    int jmax = solver->jmax;
    double gx = solver->gx;
    double gy = solver->gy;
    double gamma = solver->gamma;
    double dt = solver->dt;
    double inverseRe = 1.0 / solver->re;
    double inverseDx = 1.0 / solver->dx;
    double inverseDy = 1.0 / solver->dy;
    double du2dx, dv2dy, duvdx, duvdy;
    double du2dx2, du2dy2, dv2dx2, dv2dy2;
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;

    exchange(solver, u);
    exchange(solver, v);

    for (int j = 1; j < jmaxLocal + 1; j++)
    {
        for (int i = 1; i < imaxLocal + 1; i++)
        {
            du2dx = inverseDx * 0.25 *
                        ((U(i, j) + U(i + 1, j)) * (U(i, j) + U(i + 1, j)) -
                         (U(i, j) + U(i - 1, j)) * (U(i, j) + U(i - 1, j))) +
                    gamma * inverseDx * 0.25 *
                        (fabs(U(i, j) + U(i + 1, j)) * (U(i, j) - U(i + 1, j)) +
                         fabs(U(i, j) + U(i - 1, j)) * (U(i, j) - U(i - 1, j)));

            duvdy = inverseDy * 0.25 *
                        ((V(i, j) + V(i + 1, j)) * (U(i, j) + U(i, j + 1)) -
                         (V(i, j - 1) + V(i + 1, j - 1)) * (U(i, j) + U(i, j - 1))) +
                    gamma * inverseDy * 0.25 *
                        (fabs(V(i, j) + V(i + 1, j)) * (U(i, j) - U(i, j + 1)) +
                         fabs(V(i, j - 1) + V(i + 1, j - 1)) *
                             (U(i, j) - U(i, j - 1)));

            du2dx2 = inverseDx * inverseDx * (U(i + 1, j) - 2.0 * U(i, j) + U(i - 1, j));
            du2dy2 = inverseDy * inverseDy * (U(i, j + 1) - 2.0 * U(i, j) + U(i, j - 1));
            F(i, j) = U(i, j) + dt * (inverseRe * (du2dx2 + du2dy2) - du2dx - duvdy + gx);

            duvdx = inverseDx * 0.25 *
                        ((U(i, j) + U(i, j + 1)) * (V(i, j) + V(i + 1, j)) -
                         (U(i - 1, j) + U(i - 1, j + 1)) * (V(i, j) + V(i - 1, j))) +
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

            dv2dx2 = inverseDx * inverseDx * (V(i + 1, j) - 2.0 * V(i, j) + V(i - 1, j));
            dv2dy2 = inverseDy * inverseDy * (V(i, j + 1) - 2.0 * V(i, j) + V(i, j - 1));
            G(i, j) = V(i, j) + dt * (inverseRe * (dv2dx2 + dv2dy2) - duvdx - dv2dy + gy);
        }
    }

    /* ----------------------------- boundary of F --------------------------- */
    if (isBoundary(solver, LEFT))
    {
        for (int j = 1; j < jmaxLocal + 1; j++)
        {
            F(0, j) = U(0, j);
        }
    }

    if (isBoundary(solver, RIGHT))
    {
        for (int j = 1; j < jmaxLocal + 1; j++)
        {
            F(imaxLocal, j) = U(imaxLocal, j);
        }
    }

    /* ----------------------------- boundary of G --------------------------- */
    if (isBoundary(solver, BOTTOM))
    {
        for (int i = 1; i < imaxLocal + 1; i++)
        {
            G(i, 0) = V(i, 0);
        }
    }

    if (isBoundary(solver, TOP))
    {
        for (int i = 1; i < imaxLocal + 1; i++)
        {
            G(i, jmaxLocal) = V(i, jmaxLocal);
        }
    }
}

// TODO adapt to local sizes
void adaptUV(Solver *solver)
{
    int imax = solver->imax;
    int jmax = solver->jmax;
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;
    double *p = solver->p;
    double *u = solver->u;
    double *v = solver->v;
    double *f = solver->f;
    double *g = solver->g;
    double factorX = solver->dt / solver->dx;
    double factorY = solver->dt / solver->dy;

    for (int j = 1; j < jmaxLocal + 1; j++)
    {
        for (int i = 1; i < imaxLocal + 1; i++)
        {
            U(i, j) = F(i, j) - (P(i + 1, j) - P(i, j)) * factorX;
            V(i, j) = G(i, j) - (P(i, j + 1) - P(i, j)) * factorY;
        }
    }
}

void writeResult(Solver *solver, double *Pall, double *Uall, double *Vall)
{

    int imax = solver->imax;
    int jmax = solver->jmax;
    double dx = solver->dx;
    double dy = solver->dy;
    double x = 0.0, y = 0.0;

    FILE *fp;
    fp = fopen("pressure.dat", "w");

    if (fp == NULL)
    {
        printf("Error!\n");
        exit(EXIT_FAILURE);
    }

    for (int j = 1; j < jmax + 1; j++)
    {
        y = (double)(j - 0.5) * dy;
        for (int i = 1; i < imax + 1; i++)
        {
            x = (double)(i - 0.5) * dx;
            fprintf(fp, "%.2f %.2f %f\n", x, y, Pall(i, j));
        }
        fprintf(fp, "\n");
    }

    fclose(fp);

    fp = fopen("velocity.dat", "w");

    if (fp == NULL)
    {
        printf("Error!\n");
        exit(EXIT_FAILURE);
    }

    for (int j = 1; j < jmax + 1; j++)
    {
        y = dy * (j - 0.5);
        for (int i = 1; i < imax + 1; i++)
        {
            x = dx * (i - 0.5);
            double velU = (Uall(i, j) + Uall(i - 1, j)) / 2.0;
            double velV = (Vall(i, j) + Vall(i, j - 1)) / 2.0;
            double len = sqrt((velU * velU) + (velV * velV));
            fprintf(fp, "%.2f %.2f %f %f %f\n", x, y, velU, velV, len);
        }
    }

    fclose(fp);
}
