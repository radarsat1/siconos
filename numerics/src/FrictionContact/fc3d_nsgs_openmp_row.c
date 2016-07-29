/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "fc3d_onecontact_nonsmooth_Newton_solvers.h"
#include "fc3d_Path.h"
#include "fc3d_NCPGlockerFixedPoint.h"
#include "fc3d_projection.h"
#include "fc3d_unitary_enumerative.h"
#include "fc3d_compute_error.h"
#include "NCP_Solvers.h"
#include "SiconosBlas.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <alloca.h>
#include "op3x3.h"
#pragma GCC diagnostic ignored "-Wmissing-prototypes"

#include <SiconosConfig.h>
#if defined(WITH_OPENMP) && defined(_OPENMP)
#define USE_OPENMP 1
#include <omp.h>
#endif

void snPrintf(int level, SolverOptions* opts, const char *fmt, ...);

void fc3d_nsgs_computeqLocal_rowpar(FrictionContactProblem * problem, FrictionContactProblem * localproblem, double *reaction, int contact)
{

  double *qLocal = localproblem->q;
  int n = 3 * problem->numberOfContacts;


  int in = 3 * contact, it = in + 1, is = it + 1;

  /* qLocal computation*/
  qLocal[0] = problem->q[in];
  qLocal[1] =  problem->q[it];
  qLocal[2] =  problem->q[is];

  if (problem->M->storageType == 0)
  {
    double * MM = problem->M->matrix0;
    int incx = n, incy = 1;
    /* reaction current block set to zero, to exclude current contact block */
    double rin = reaction[in] ;
    double rit = reaction[it] ;
    double ris = reaction[is] ;
    reaction[in] = 0.0;
    reaction[it] = 0.0;
    reaction[is] = 0.0;
    qLocal[0] += cblas_ddot(n , &MM[in] , incx , reaction , incy);
    qLocal[1] += cblas_ddot(n , &MM[it] , incx , reaction , incy);
    qLocal[2] += cblas_ddot(n , &MM[is] , incx , reaction , incy);
    reaction[in] = rin;
    reaction[it] = rit;
    reaction[is] = ris;
  }
  else if (problem->M->storageType == 1)
  {
    /* qLocal += rowMB * reaction
     * with rowMB the row of blocks of MGlobal which corresponds
     * to the current contact
     */
    rowProdNoDiagSBM3x3par(n, 3, contact, problem->M->matrix1, reaction, qLocal);
  }
}

void fc3d_nsgs_openmp_row(FrictionContactProblem* problem, double *reaction,
                          double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  unsigned int nc = problem->numberOfContacts;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];
  double normq = cblas_dnrm2(nc*3 , problem->q , 1);

  if (*info == 0)
    return;

  if (options->numberOfInternalSolvers < 1)
  {
    numericsError("fc3d_nsgs_redblack_openmp", "The NSGS method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >= 1");
  }
  assert(options->internalSolvers);

  SolverOptions * localsolver_options = options->internalSolvers;

  SolverPtr local_solver = NULL;
  UpdatePtr update_localproblem = NULL;
  FreeSolverNSGSPtr freeSolver = NULL;
  ComputeErrorPtr computeError = NULL;
  FrictionContactProblem *localproblem = malloc(sizeof(FrictionContactProblem));

  if (verbose > 0) printf("----------------------------------- number of threads %i\n", omp_get_max_threads()  );
  if (verbose > 0) printf("----------------------------------- number of contacts %i\n", nc );

  /* Connect local solver and local problem*/
  localproblem->numberOfContacts = 1;
  localproblem->dimension = 3;
  localproblem->q = (double*)malloc(3 * sizeof(double));
  localproblem->mu = (double*)malloc(sizeof(double));

  if (problem->M->storageType == NM_DENSE || problem->M->storageType == NM_SPARSE)
  {
    localproblem->M = createNumericsMatrixFromData(NM_DENSE, 3, 3,
                                                   malloc(9 * sizeof(double)));
  }
  else /* NM_SPARSE_BLOCK */
  {
    localproblem->M = createNumericsMatrixFromData(NM_DENSE, 3, 3, NULL);
  }

  fc3d_nsgs_initialize_local_solver(&local_solver, &update_localproblem,
                                    (FreeSolverNSGSPtr *)&freeSolver, &computeError,
                                    problem, localproblem,
                                    options, localsolver_options);

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  double error_delta_reaction=0.0;
  double error_nat=0.0;
  unsigned int *scontacts = NULL;

  printf("threads = %d\n", omp_get_max_threads());

  double t0 = omp_get_wtime();

  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    printf("Iteration %d.. ", iter);
    error_delta_reaction=0.0;

    /* Loop through the contact points */
    for ( unsigned int contact = 0 ; contact < nc ; contact+=1)
    {
      if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
      (*update_localproblem)(contact, problem, localproblem,
                             reaction, localsolver_options);

      localsolver_options->iparam[4] = contact;

      (*local_solver)(localproblem, &(reaction[3 * contact]), localsolver_options);
    }
    /* printf("done.\n"); */

    /* **** Criterium convergence **** */
    if (iparam[8] >0)
    {
      if (iter % iparam[8] == 0) {
        /* printf("Computing error.."); */
        /* fflush(stdout); */
        (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);
        /* printf("done.\n"); */
        if (error > tolerance && iparam[9] == 1)
          iparam[8] *= 2;
      }
    }
    else {
      /* printf("Computing error.."); */
      /* fflush(stdout); */
      (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);
      /* printf("done.\n"); */
    }
    printf("error = %g", error);
    double t = omp_get_wtime();
    printf("    time per iter = %g", (t - t0)/iter);
    printf("      \r");
    fflush(stdout);

    if (error < tolerance)
    {
      hasNotConverged = 0;
      if (verbose > 0)
        printf("----------------------------------- FC3D - NSGS - Iteration %i Residual = %14.7e < %7.3e\n", iter, error, options->dparam[0]);
    }
    else
    {
      if (verbose > 0)
        printf("----------------------------------- FC3D - NSGS - Iteration %i Residual = %14.7e > %7.3e\n", iter, error, options->dparam[0]);
    }

    *info = hasNotConverged;

    if (options->callback)
    {
      options->callback->collectStatsIteration(options->callback->env, 3 * nc,
                                               reaction, velocity,
                                               error, NULL);
    }
  }

  dparam[0] = tolerance;
  dparam[1] = error;
  iparam[7] = iter;
  printf("\n");

  /***** Free memory *****/
  (*freeSolver)(problem,localproblem,localsolver_options);
  if (problem->M->storageType == NM_DENSE && localproblem->M->matrix0)
  {
    free(localproblem->M->matrix0);
  }
  localproblem->M->matrix0 = NULL;
  freeFrictionContactProblem(localproblem);
  solver_options_delete(localsolver_options);
  free(localsolver_options);

  if (scontacts) /* shuffle */
  {
    free(scontacts);
  }
}
