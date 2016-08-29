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
#include <float.h>

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void fake_compute_error_nsgs(FrictionContactProblem* problem, double *reaction, double *velocity, double tolerance, SolverOptions  *options,  double* error)
{
  int n = 3 * problem->numberOfContacts;
  *error = 0.;
  int i, m;
  m = 5 * n / 3;
  double err = INFINITY;
  for (i = 0 ; i < m ; ++i)
  {
    *error += Compute_NCP_error1(i, err);
  }
}
void fc3d_nsgs_update(int contact, FrictionContactProblem* problem, FrictionContactProblem* localproblem, double * reaction, SolverOptions* options)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */
  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  fc3d_nsgs_fillMLocal(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
     excluding the block corresponding to the current contact. ****/
  fc3d_nsgs_computeqLocal(problem, localproblem, reaction, contact);

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];


}

void fc3d_nsgs_initialize_local_solver(SolverPtr* solve, UpdatePtr* update,
                                       FreeSolverNSGSPtr* freeSolver,
                                       ComputeErrorPtr* computeError,
                                       FrictionContactProblem* problem,
                                       FrictionContactProblem* localproblem,
                                       SolverOptions * options,
                                       SolverOptions * localsolver_options)
{

  /** Connect to local solver */
  switch (localsolver_options->solverId)
  {
    /* Projection */
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization:
  {
    *solve = &fc3d_projectionWithDiagonalization_solve;
    *update = &fc3d_projectionWithDiagonalization_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_projection_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_projection_initialize(problem, localproblem);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone:
  {
    *solve = &fc3d_projectionOnCone_solve;
    *update = &fc3d_projection_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_projection_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_projection_initialize(problem, localproblem);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration:
  {
    *solve = &fc3d_projectionOnConeWithLocalIteration_solve;
    *update = &fc3d_projection_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_projectionOnConeWithLocalIteration_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_projectionOnConeWithLocalIteration_initialize(problem, localproblem,localsolver_options );
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization:
  {
    *solve = &fc3d_projectionOnCone_solve;
    *update = &fc3d_projection_update_with_regularization;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_projection_with_regularization_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_projection_initialize_with_regularization(problem, localproblem);
    break;
  }
  /* Newton solver (Alart-Curnier) */
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_AC:
  {
    *solve = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
    *update = &fc3d_onecontact_nonsmooth_Newton_AC_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_onecontact_nonsmooth_Newton_solvers_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP:
  {
    *solve = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
    *update = &fc3d_onecontact_nonsmooth_Newton_AC_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_onecontact_nonsmooth_Newton_solvers_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP_P:
  {
    *solve = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
    *update = &fc3d_onecontact_nonsmooth_Newton_AC_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_onecontact_nonsmooth_Newton_solvers_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
    break;
  }  /* Newton solver (Glocker-Fischer-Burmeister)*/
  case SICONOS_FRICTION_3D_NCPGlockerFBNewton:
  {
    *solve = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
    *update = &NCPGlocker_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_onecontact_nonsmooth_Newton_solvers_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    // *computeError = &fake_compute_error;
    fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
    break;
  }
  /* Path solver (Glocker Formulation) */
  case SICONOS_FRICTION_3D_NCPGlockerFBPATH:
  {
    *solve = &fc3d_Path_solve;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_Path_free;
    *update = &NCPGlocker_update;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    // *computeError = &fake_compute_error;
    fc3d_Path_initialize(problem, localproblem, localsolver_options);
    break;
  }

  /* Fixed Point solver (Glocker Formulation) */
  case SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint:
  {
    *solve = &fc3d_FixedP_solve;
    *update = &NCPGlocker_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_FixedP_free;
    /* *computeError = &fake_compute_error_nsgs; */
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_FixedP_initialize(problem, localproblem, localsolver_options);
    break;
  }
  /*iparam[4] > 10 are reserved for Tresca resolution */
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder:
  {
    *solve = &fc3d_projectionOnCylinder_solve;
    *update = &fc3d_projectionOnCylinder_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_projectionOnCylinder_free;
    *computeError = (ComputeErrorPtr)&fc3d_Tresca_compute_error;
    fc3d_projectionOnCylinder_initialize(problem, localproblem, options );
    break;
  }
    /*iparam[4] > 10 are reserved for Tresca resolution */
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration:
  {
    *solve = &fc3d_projectionOnCylinderWithLocalIteration_solve;
    *update = &fc3d_projectionOnCylinder_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_projectionOnCylinderWithLocalIteration_free;
    *computeError = (ComputeErrorPtr)&fc3d_Tresca_compute_error;
    fc3d_projectionOnCylinderWithLocalIteration_initialize(problem, localproblem, options, localsolver_options );
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC:
  {
    *solve = &fc3d_unitary_enumerative_solve;
    *update = &fc3d_nsgs_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_unitary_enumerative_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_unitary_enumerative_initialize(localproblem);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU:
  {
    *solve = &fc3d_unitary_enumerative_solve;
    *update = &fc3d_nsgs_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_unitary_enumerative_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_unitary_enumerative_initialize(localproblem);
    break;
  }
  default:
  {
    fprintf(stderr, "Numerics, fc3d_nsgs failed. Unknown internal solver : %s.\n", solver_options_id_to_name(localsolver_options->solverId));
    exit(EXIT_FAILURE);
  }
  }
}
void fc3d_nsgs_computeqLocal(FrictionContactProblem * problem, FrictionContactProblem * localproblem, double *reaction, int contact)
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
    SBM_row_prod_no_diag_3x3(n, 3, contact, problem->M->matrix1, reaction, qLocal);
  }



}

void fc3d_nsgs_fillMLocal(FrictionContactProblem * problem, FrictionContactProblem * localproblem, int contact)
{

  NumericsMatrix * MGlobal = problem->M;

  int n = 3 * problem->numberOfContacts;


  int storageType = MGlobal->storageType;
  if (storageType == 0)
    // Dense storage
  {
    int in = 3 * contact, it = in + 1, is = it + 1;
    int inc = n * in;
    double * MM = MGlobal->matrix0;
    double * MLocal =  localproblem->M->matrix0;

    /* The part of MM which corresponds to the current block is copied into MLocal */
    MLocal[0] = MM[inc + in];
    MLocal[1] = MM[inc + it];
    MLocal[2] = MM[inc + is];
    inc += n;
    MLocal[3] = MM[inc + in];
    MLocal[4] = MM[inc + it];
    MLocal[5] = MM[inc + is];
    inc += n;
    MLocal[6] = MM[inc + in];
    MLocal[7] = MM[inc + it];
    MLocal[8] = MM[inc + is];
  }
  else if (storageType == 1)
  {
    int diagPos = getDiagonalBlockPos(MGlobal->matrix1, contact);
    localproblem->M->matrix0 = MGlobal->matrix1->block[diagPos];

  }
  else
    numericsError("fc3d_projection -", "unknown storage type for matrix M");
}


/* swap two indices */
void uint_swap (unsigned int *a, unsigned int *b)
{
    unsigned int temp = *a;
    *a = *b;
    *b = temp;
}

/* shuffle an unsigned array */
void uint_shuffle (unsigned int *a, unsigned int n) {

  for (unsigned int i = 0; i < n - 1; i++)
  {
    uint_swap  (&a[i], &a[i + rand()%(n - i)]);
  }
}

FrictionContactProblem* allocLocalProblem(FrictionContactProblem* problem)
{
  /* Connect local solver and local problem*/
  FrictionContactProblem* localproblem =
    (FrictionContactProblem*)malloc(sizeof(FrictionContactProblem));
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
  return localproblem;
}

static
void freeLocalProblem(FrictionContactProblem* localproblem,
                      FrictionContactProblem* problem)
{
  if (problem->M->storageType == NM_DENSE && localproblem->M->matrix0)
  {
    free(localproblem->M->matrix0);
  }
  localproblem->M->matrix0 = NULL;
  freeFrictionContactProblem(localproblem);
}

static
unsigned int* allocShuffledContacts(FrictionContactProblem *problem,
                                    SolverOptions *options)
{
  unsigned int *scontacts = 0;
  unsigned int nc = problem->numberOfContacts;
  if (options->iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE||
      options->iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP)
  {
    if (options->iparam[6] >0)
    {
      srand((unsigned int)options->iparam[6]);
    }
    else
      srand(1);
    scontacts = (unsigned int *) malloc(nc * sizeof(unsigned int));
    for (unsigned int i = 0; i < nc ; ++i)
    {
      scontacts[i] = i;
    }
    uint_shuffle(scontacts, nc);
  }
  return scontacts;
}

static
void dump_reaction_history(char *fn)
{
  printf("Dumping reaction history..\n");
  FILE *f = 0;
  f = fopen(fn, "w");
  if (!f) {
    printf("file problem\n");
    exit(1);
  }
  for (int i=0; i<reaction_history_pos; i++)
  {
    fprintf(f, "%g,%g,%g\n", reaction_history[i][0],
            reaction_history[i][1],
            reaction_history[i][2]);
  }
  fclose(f);
}

static
void dump_problem_volume_and_quit(char *fn,
                                  FrictionContactProblem *localproblem,
                                  SolverOptions *localsolver_options)
{
  printf("Dumping problem volume..\n");
  FILE *f = 0;
  f = fopen(fn, "w");
  if (!f) {
    printf("file problem\n");
    exit(1);
  }
  for (int i=0; i<100; i++) {
    for (int j=0; j<100; j++) {
      for (int k=0; k<100; k++) {
        double range = 1000.0;
        double r0 = (i/100.0 * 2.0 - 1.0) * range;
        double r1 = (j/100.0 * 2.0 - 1.0) * range;
        double r2 = (k/100.0 * 2.0 - 1.0) * range;
        double R[3] = {r0, r1, r2};
        double val =
          fc3d_onecontact_nonsmooth_Newton_solvers_solve_damped_oneitereval(
            localproblem, R, localsolver_options->iparam,
            localsolver_options->dparam);
        //fprintf(f, "%g,%g,%g,%g\n", r0, r1, r2, val);
        fprintf(f, "%g\n", log10(val));
      }
    }
  }
  fclose(f);
}

static
int solveLocalReaction(UpdatePtr update_localproblem, SolverPtr local_solver,
                       unsigned int contact, FrictionContactProblem *problem,
                       FrictionContactProblem *localproblem, double *reaction,
                       SolverOptions *localsolver_options, double localreaction[3])
{
  (*update_localproblem)(contact, problem, localproblem,
                         reaction, localsolver_options);

  localsolver_options->iparam[4] = contact;

  localreaction[0] = reaction[contact*3 + 0];
  localreaction[1] = reaction[contact*3 + 1];
  localreaction[2] = reaction[contact*3 + 2];

  return (*local_solver)(localproblem, localreaction, localsolver_options);
}

static
void performRelaxation(double localreaction[3], double *oldreaction, double omega)
{
  localreaction[0] = omega*localreaction[0]+(1.0-omega)*oldreaction[0];
  localreaction[1] = omega*localreaction[1]+(1.0-omega)*oldreaction[1];
  localreaction[2] = omega*localreaction[2]+(1.0-omega)*oldreaction[2];
}

static
void accumulateLightErrorSum(double *light_error_sum, double localreaction[3],
                             double *oldreaction)
{
  *light_error_sum += ( pow(oldreaction[0] - localreaction[0], 2) +
                        pow(oldreaction[1] - localreaction[1], 2) +
                        pow(oldreaction[2] - localreaction[2], 2) );
}

static
void acceptLocalReactionFiltered(FrictionContactProblem *localproblem,
                                 SolverOptions *localsolver_options,
                                 unsigned int contact, unsigned int iter,
                                 double *reaction, double localreaction[3])
{
  if (isnan(localsolver_options->dparam[1])
      || isinf(localsolver_options->dparam[1])
      || localsolver_options->dparam[1] > 1.0)
  {
    DEBUG_EXPR(frictionContact_display(localproblem));
    DEBUG_PRINTF("Discard local reaction for contact %i at iteration %i "
                 "with local_error = %e\n",
                 contact, iter, localsolver_options->dparam[1]);
  }
  else
    memcpy(&reaction[contact*3], localreaction, sizeof(double)*3);
}

static
void acceptLocalReactionProjected(int *rc, FrictionContactProblem *problem,
                                  FrictionContactProblem *localproblem,
                                  SolverPtr local_solver,
                                  SolverOptions *localsolver_options,
                                  unsigned int contact, unsigned int iter,
                                  double *reaction, double localreaction[3])
{
  int nan1 = isnan(localsolver_options->dparam[1]) || isinf(localsolver_options->dparam[1]);
  if (nan1 || localsolver_options->dparam[1] > 1.0)
  {
    /*
    if (isnan(localsolver_options->dparam[1])) {
      printf("Problem resulted in NaN\n");
      dump_problem_volume_and_quit("badproblem.csv", localproblem, localsolver_options);
      exit(1);
    }
    else {
      dump_problem_volume_and_quit("goodproblem.csv", localproblem, localsolver_options);
    }
    */

    if (isnan(localsolver_options->dparam[1])) {
      printf("Problem resulted in NaN\n");
      dump_reaction_history("badreaction.csv");
      exit(1);
    }
    else {
      dump_reaction_history("goodreaction.csv");
    }

    /* DEBUG_EXPR( */
    /*   frictionContact_display(localproblem); */
    #define DEBUG_PRINTF printf

      // In case of bad solution, re-solve with a projection-on-cone solver
      DEBUG_PRINTF("Discard local reaction for contact %i at iteration %i with local_error = %e\n", contact, iter, localsolver_options->dparam[1]);
      memcpy(localreaction, &reaction[contact*3], sizeof(double)*3);
      fc3d_projectionOnConeWithLocalIteration_initialize(problem, localproblem, localsolver_options );
      fc3d_projectionOnConeWithLocalIteration_solve(localproblem, localreaction, localsolver_options);
      int nan2 = isnan(localsolver_options->dparam[1]) || isinf(localsolver_options->dparam[1]);
      if (nan2) {
        DEBUG_PRINTF("No hope for contact %d, setting to zero.\n", contact);
        reaction[3*contact+0] = 0;
        reaction[3*contact+1] = 0;
        reaction[3*contact+2] = 0;
        return;
      }

      // Save the result in case next step fails
      double localreaction_proj[3];
      double error_proj = localsolver_options->dparam[1];
      memcpy(localreaction_proj, localreaction, sizeof(double)*3);

      // Complete it further with the original solver
      DEBUG_PRINTF("Try local fc3d_projectionOnConeWithLocalIteration_solve with local_error = %e\n", localsolver_options->dparam[1]);
      (*local_solver)(localproblem, localreaction , localsolver_options);
      int nan3 = isnan(localsolver_options->dparam[1]) || isinf(localsolver_options->dparam[1]);
      DEBUG_PRINTF("Try local another newton solve with local_error = %e\n", localsolver_options->dparam[1]);

      // If we produced a bad value doing this, keep the projectionOnCone solution
      // (Note: this has been observed, particularly when mu=1.0)
      if (nan3)
      {
        DEBUG_PRINTF("Keep the projectionOnCone local solution = %e\n", error_proj);
        memcpy(&reaction[contact*3], localreaction_proj, sizeof(double)*3);
      }
      else
      {
        if (nan1 || localsolver_options->dparam[1] <= localsolver_options->dparam[0])
        {
          DEBUG_PRINTF("Keep the new local solution = %e\n", localsolver_options->dparam[1]);
          memcpy(&reaction[contact*3], localreaction, sizeof(double)*3);
          //getchar();
          *rc = 0;
        }
        else
        {
          DEBUG_PRINTF("Keep the previous local solution = %e\n", localsolver_options->dparam[1]);
        }
      }
    /* ); */
      #undef DEBUG_PRINTF
  }
  else
    memcpy(&reaction[contact*3], localreaction, sizeof(double)*3);
}

static
void acceptLocalReactionUnconditionally(unsigned int contact,
                                        double *reaction, double localreaction[3])
{
  memcpy(&reaction[contact*3], localreaction, sizeof(double)*3);
}

static
double calculateLightError(double light_error_sum, unsigned int nc, double *reaction)
{
  double error = sqrt(light_error_sum);
  double norm_r = cblas_dnrm2(nc*3 , reaction , 1);
  if (fabs(norm_r) > DBL_EPSILON)
    error /= norm_r;
  return error;
}

static
double calculateFullErrorAdaptiveInterval(FrictionContactProblem *problem,
                                          ComputeErrorPtr computeError,
                                          SolverOptions *options, int iter,
                                          double *reaction, double *velocity,
                                          double tolerance, double normq)
{
  double error;
  if (options->iparam[8] >0)
  {
    if (iter % options->iparam[8] == 0) {
      (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);
      if (error > tolerance
          && options->iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE)
        options->iparam[8] *= 2;
    }
    if (verbose > 0)
      printf("----------------------------------- FC3D - NSGS - Iteration %i "
             "options->iparam[8] = %i, options->iparam[1] = % i \n",
             iter, options->iparam[8], options->iparam[1]);
  }
  else
    (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);

  return error;
}

static
int determineConvergence(double error, double tolerance, int iter,
                         SolverOptions *options)
{
  int hasNotConverged = 1;
  if (error < tolerance)
  {
    hasNotConverged = 0;
    if (verbose > 0)
      printf("----------------------------------- FC3D - NSGS - Iteration %i "
             "Residual = %14.7e < %7.3e\n", iter, error, options->dparam[0]);
  }
  else
  {
    if (verbose > 0)
      printf("----------------------------------- FC3D - NSGS - Iteration %i "
             "Residual = %14.7e > %7.3e\n", iter, error, options->dparam[0]);
  }
  return hasNotConverged;
}

static
double calculateFullErrorFinal(FrictionContactProblem *problem, SolverOptions *options,
                               ComputeErrorPtr computeError,
                               double *reaction, double *velocity, double tolerance,
                               double normq, double light_error)
{
  double absolute_error;
  (*computeError)(problem, reaction , velocity, tolerance,
                  options, normq, &absolute_error);

  if (verbose > 0 && absolute_error > light_error)
  {
    printf("--------------------------- FC3D - NSGS - Warning absolute "
           "Residual = %14.7e is larger than incremental error = %14.7e\n",
           absolute_error, light_error);
  }
  return absolute_error;
}

static
void statsIterationCallback(FrictionContactProblem *problem,
                            SolverOptions *options,
                            double *reaction, double *velocity, double error)
{
  if (options->callback)
  {
    options->callback->collectStatsIteration(options->callback->env,
                                             problem->numberOfContacts * 3,
                                             reaction, velocity,
                                             error, NULL);
  }
}

void fc3d_nsgs(FrictionContactProblem* problem, double *reaction,
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
  double omega = dparam[8];

  SolverOptions * localsolver_options = options->internalSolvers;
  SolverPtr local_solver = NULL;
  UpdatePtr update_localproblem = NULL;
  FreeSolverNSGSPtr freeSolver = NULL;
  ComputeErrorPtr computeError = NULL;

  FrictionContactProblem* localproblem;
  double localreaction[3];

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  unsigned int contact; /* Number of the current row of blocks in M */
  unsigned int *scontacts = NULL;

  if (*info == 0)
    return;

  if (options->numberOfInternalSolvers < 1)
  {
    numericsError("fc3d_nsgs",
                  "The NSGS method needs options for the internal solvers, "
                  "options[0].numberOfInternalSolvers should be >= 1");
  }
  assert(options->internalSolvers);

  /*****  Initialize various solver options *****/
  localproblem = allocLocalProblem(problem);

  fc3d_nsgs_initialize_local_solver(&local_solver, &update_localproblem,
                             (FreeSolverNSGSPtr *)&freeSolver, &computeError,
                             problem, localproblem, options, localsolver_options);

  scontacts = allocShuffledContacts(problem, options);

  /*****  Check solver options *****/
  if (! (iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE
         || iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE
         || iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP))
  {
    numericsError(
      "fc3d_nsgs", "iparam[5] must be equal to "
      "SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE (0), "
      "SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE (1) or "
      "SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP (2)");
    return;
  }

  if (! (iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL
         || iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL
         || iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT
         || iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE))
  {
    numericsError(
      "fc3d_nsgs", "iparam[1] must be equal to "
      "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL (0), "
      "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL (1), "
      "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT (2) or "
      "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE (3)");
    return;
  }

  /*****  NSGS Iterations *****/

#define BIGREACTION(x) ((x)>1e10 || (x)<-1e10)
#define BADREACTION(x) (isnan(x) || isinf(x))

  /* A special case for the most common options (should correspond
   * with mechanics_io.py **/
  if (iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE
      && iparam[4] == SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE
      && iparam[14] == SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE
      && iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
  {
    int *badi = calloc(nc, sizeof(int));
    double *badr = calloc(nc*4, sizeof(double));
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      double light_error_sum = 0.0;
      for (unsigned int i = 0 ; i < nc ; ++i)
      {
        contact = i;
        double localreaction[3] = {
          reaction[i*3+0], reaction[i*3+1], reaction[i*3+2] };
        if (BADREACTION(localreaction[0]) ||
            BADREACTION(localreaction[1]) ||
            BADREACTION(localreaction[2])) {
          printf("bad reaction found on _input_! for contact %d\n", i);
          printf(": (%g, %g, %g) <- (%g, %g, %g)\n", localreaction[0],
                 localreaction[1], localreaction[2],
                 reaction[i*3+0], reaction[i*3+1], reaction[i*3+2]);
        }
        int rc =
        solveLocalReaction(update_localproblem, local_solver, contact,
                           problem, localproblem, reaction, localsolver_options,
                           localreaction);

        accumulateLightErrorSum(&light_error_sum, localreaction, &reaction[contact*3]);

        #if 0
          acceptLocalReactionUnconditionally(contact, reaction, localreaction);
        #else
        #if 0
        acceptLocalReactionFiltered(localproblem, localsolver_options,
                                    contact, iter, reaction, localreaction);
        #else
        // Experimental
        acceptLocalReactionProjected(&rc, problem, localproblem, local_solver, localsolver_options,
                                     contact, iter, reaction, localreaction);
        #endif
        #endif
        #if 1
        if (rc==1
	    || BADREACTION(reaction[i*3])
            || BADREACTION(reaction[i*3+1])
            || BADREACTION(reaction[i*3+2])
	    || BIGREACTION(reaction[i*3])
            || BIGREACTION(reaction[i*3+1])
            || BIGREACTION(reaction[i*3+2]))
        {
          badi[i] = iter;
          badr[i*4+0] = reaction[i*3+0];
          badr[i*4+1] = reaction[i*3+1];
          badr[i*4+2] = reaction[i*3+2];
          badr[i*4+3] = dparam[1];
          printf("bad value (%g, %g, %g)  (rc=%d) found for contact %d\n",
                 reaction[i*3],
                 reaction[i*3+1],
                 reaction[i*3+2], rc, contact);
        }
        #endif

        #if 0
        if (BADREACTION(reaction[i*3])
            || BADREACTION(reaction[i*3+1])
            || BADREACTION(reaction[i*3+2]))
        {
          printf("bad value (%g, %g, %g) found, resetting to %g, %g, %g\n",
                 reaction[i*3],
                 reaction[i*3+1],
                 reaction[i*3+2],
                 localreaction[0],
                 localreaction[1],
                 localreaction[2]);
          reaction[i*3]=localreaction[0];
          reaction[i*3+1]=localreaction[1];
          reaction[i*3+2]=localreaction[2];

          if (BADREACTION(reaction[i*3])
              || BADREACTION(reaction[i*3+1])
              || BADREACTION(reaction[i*3+2]))
          { // bad values on input! without this, normq>0 error
            reaction[i*3] = 0;
            reaction[i*3+1] = 0;
            reaction[i*3+2] = 0;
            printf("still bad, now set to %g, %g, %g\n",
                   reaction[i*3],
                   reaction[i*3+1],
                   reaction[i*3+2]);
          }
          //abort();
        }
        #endif
      }

      error = calculateLightError(light_error_sum, nc, reaction);

      hasNotConverged = determineConvergence(error, tolerance, iter, options);

      statsIterationCallback(problem, options, reaction, velocity, error);
    }
    for (int i=0; i<nc; i++) {
      if (badi[i] > 0) {
	printf("Bad reaction for contact %d on interation %d\n", i, badi[i]);
	printf(": %g, %g, %g (solver dparam[1]=%g)\n",
	       badr[i*4], badr[i*4+1], badr[i*4+2], badr[i*4+3]);
      }
    }
    free(badi);
    free(badr);
  }

  /* All other cases, we put all the ifs inline.. otherwise, too many
   * variations to have dedicated loops, but add more if there are
   * common cases to avoid checking booleans on every iteration. **/
  else
  {
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      double light_error_sum = 0.0;
      for (unsigned int i = 0 ; i < nc ; ++i)
      {
        if (iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE
            || iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP)
        {
          if (iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP)
            uint_shuffle(scontacts, nc);
          contact = scontacts[i];
        }
        else
          contact = i;

        solveLocalReaction(update_localproblem, local_solver, contact,
                           problem, localproblem, reaction, localsolver_options,
                           localreaction);

        if (iparam[4] == SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE)
          performRelaxation(localreaction, &reaction[contact*3], omega);

        accumulateLightErrorSum(&light_error_sum, localreaction, &reaction[contact*3]);

        if (iparam[14] == SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE)
          acceptLocalReactionFiltered(localproblem, localsolver_options,
                                      contact, iter, reaction, localreaction);
        else
          acceptLocalReactionUnconditionally(contact, reaction, localreaction);
      }
      if (iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT ||
          iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
        error = calculateLightError(light_error_sum, nc, reaction);
      else if (iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL)
        error = calculateFullErrorAdaptiveInterval(problem, computeError, options,
                                                   iter, reaction, velocity,
                                                   tolerance, normq);

      hasNotConverged = determineConvergence(error, tolerance, iter, options);

      statsIterationCallback(problem, options, reaction, velocity, error);
    }
  }

  *info = hasNotConverged;

  /* Full criterium */
  if (iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
    error = calculateFullErrorFinal(problem, options, computeError, reaction, velocity,
                                    tolerance, normq, error /* = light error*/);

  /** return parameter values */
  dparam[0] = tolerance;
  dparam[1] = error;
  iparam[7] = iter;

  /** Free memory **/
  (*freeSolver)(problem,localproblem,localsolver_options);
  freeLocalProblem(localproblem, problem);
  if (scontacts) free(scontacts);
}

int fc3d_nsgs_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSGS Solver\n");
  }

  /*  strcpy(options->solverName,"NSGS");*/
  options->solverId = SICONOS_FRICTION_3D_NSGS;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 20;
  options->dSize = 20;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-4;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  fc3d_onecontact_nonsmooth_Newtow_setDefaultSolverOptions(options->internalSolvers);

  return 0;
}
